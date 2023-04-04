#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import os
import pandas as pd

# own modules
from ExonSurfer.ensembl import ensembl
from ExonSurfer.readFiles import readGBK, readFasta, commonFunctions
from ExonSurfer.specificity import blast, annotate, offtargets_file_input, g_offtargets
from ExonSurfer.resources import resources
from ExonSurfer.primerDesign import designPrimers
from ExonSurfer.primerDesign import penalizePrimers, designConfig


def CreatePrimers(file, species = "homo_sapiens_masked",
                  release = 108, design_dict = designConfig.design_dict, 
                  path_out = ".", save_files = True, e_value = 0.8,
                  i_cutoff = 70, max_sep = 700, opt_prod_size = 200, NPRIMERS = 100, 
                  d_option = 1):
    """
    This function is the main function of the pipeline. It takes a gene name and
    a transcript name and returns a list of primers.
    Args:
        gene [in] (str)        Gene name
        transcripts [in] (str) Transcript name or ALL
        species [in] (str)     Species for the data (human, mus_musculus or 
                                                      rattus_novegicus)
        release [in] (str)     Ensembl release number. Change with caution
        design_dict [in] (d)   Dict with primer3 parameters
        path_out [in] (str)    Path to output direct
        d_option [in] (int)    1 for primers only on the junctions or "ALL"
        logresult [out] (str)  False if everything OK, str if design error
        blast_df [out] (df)    Dataframe with blast results
        final_df [out] (df)    Dataframe with primer design results
    """ 
    ###########################################################################
    #                 STEP 1. RETRIEVE INFORMATION & CHOOSE TARGET            #
    ########################################################################### 
    print("Extracting ensemble info")
    data = ensembl.create_ensembl_data(release, 
                                       species.replace("_masked", ""))
    
    # GenBank file format
    if file[-3:] == ".gb": 
        gb_record = readGBK.read_genbank_file(file)
        target = readGBK.extract_cdna(gb_record)
        junction = readGBK.extract_junctions(gb_record)

    # Fasta file format 
    else: 
        header, target = readFasta.read_fasta_file(file)
        junction = readFasta.extract_junctions(header)
        
    jdict = commonFunctions.get_junc_dict(junction)
    elen = commonFunctions.get_elen(junction, target)
    
    ###########################################################################
    #                          STEP 2:  DESIGN PRIMERS                        #
    ###########################################################################    
    # Get sequence and junction index
    print("Starting primer design")
    
    # Create dataframe for design
    cols = ("option", "junction", "junction_description", "forward", "reverse", 
            "amplicon_size", "forward_tm", "reverse_tm", "forward_gc", "reverse_gc", 
            "amplicon_tm", "pair_penalty")
    df = pd.DataFrame(columns = cols)
    
    if len(junction) == 0: # only one exon
        design_dict["PRIMER_NUM_RETURN"] = NPRIMERS  
        # Design primers
        c2 = designPrimers.call_primer3(target, int(len(target)/2), design_dict, 
                                        d_option, enum = 1)
        
        item = ["EXON1", "one exon"]            

        df = designPrimers.report_one_exon_design(c2, elen, item, df)
        
    else: # Normal design (more than 1 exon)

        for item in jdict: 
            # Decide number of primers to design
            if d_option == 1: 
                num_primers = int(NPRIMERS / len(junction))
            else: 
                num_primers = int(NPRIMERS / (len(junction)*2))
            
            # avoid 0 primer design
            if num_primers == 0: 
                num_primers = 1
            design_dict["PRIMER_NUM_RETURN"] = num_primers
            
            # Design primers
            c1, c2 = designPrimers.call_primer3(target, jdict[item], design_dict, 
                                                d_option)

            df = designPrimers.report_design(c1, c2, elen, item, "", df)
    
    df["pair_num"] = ["Pair{}".format(x) for x in range(0, df.shape[0])]
    df = df.set_index('pair_num')
    
    if df.shape[0] > 1: # we designed primers
    
        #######################################################################
        #                       STEP 4:  BLAST AGAINST CDNA                   #
        #######################################################################          
        # Define three output files, with the same name as the gene and transcript
        if file[-3:] == ".gb": 
            t_string = gb_record.name
        else: 
            t_string = header.split(" ")[0]

        DESIGN_OUT = os.path.join(path_out, "{}_design.txt".format(t_string))
        BLAST_OUT = os.path.join(path_out, "{}_blast.txt".format(t_string))
        GBLAST_OUT = os.path.join(path_out, "{}_gblast.txt".format( t_string))
   
        # remove in case remaining from previous design
        for file in (DESIGN_OUT, BLAST_OUT): 
            if os.path.exists(file):
                os.remove(file)
    
        # Write blast input
        FASTA_F = resources.create_temp_fasta()
        resources.fillin_temp_fasta(df, FASTA_F)
        
        # Call blast
        blast_df = blast.run_blast_list(FASTA_F, BLAST_OUT, 
                                        resources.BLAST_DB(species), 
                                        species, i_cutoff, e_value)
        
        # Delete fasta file
        if os.path.exists(FASTA_F): os.remove(FASTA_F)
        
        # Filter blast results
        blast_df = blast.pre_filter_blast(blast_df, df, e_value, i_cutoff, False)
        blast_df, df = blast.filter_big_blast(blast_df, df)
        
        # Check blast results positions
        df = offtargets_file_input.check_specificity(blast_df, df, max_sep)
        
        #######################################################################
        #                        STEP 5:  FILTER DATAFRAME                    #
        #######################################################################   
        final_df = penalizePrimers.penalize_final_output(df, "transcripts", data)
        
        #######################################################################
        #                         STEP 6:  GENOMIC BLAST                      #
        #######################################################################           
        # Write genomic blast input
        resources.fillin_temp_fasta(final_df, FASTA_F)
        
        # Call genomic blast
        blast_df = blast.run_blast_list(FASTA_F, GBLAST_OUT, 
                                        resources.BLAST_GENOMIC_DB(species), 
                                        species, i_cutoff, e_value, 
                                        tomerge = False, genomic = True)    
        # Filter big blast if needed
        blast_df = blast.pre_filter_blast(blast_df, df, e_value, i_cutoff, False)
        blast_df, final_df = blast.filter_big_blast(blast_df, final_df)
        
        # Delete fasta file
        if os.path.exists(FASTA_F): os.remove(FASTA_F)
        
        # Check blast results positions
        final_df = g_offtargets.check_genomic_specificity(blast_df, final_df, max_sep)
        
        # Filter according to genomic specificity
        final_df = penalizePrimers.genomic_filter(final_df)

        #######################################################################
        #                      STEP 7: PREPARE FINAL RESULT                   #
        #######################################################################          
        # Make penalty score
        final_df = penalizePrimers.make_penalty_score(final_df)     
        
        # Annotate transcripts detected
        final_df["not_detected"] = ""
        
        # Annotate if there are off_targets or not
        final_df = annotate.show_off_targets(final_df, "transcripts")
        
        # Remove files if necessary
        if save_files == True:
            final_df.to_csv(DESIGN_OUT, sep = "\t")
        else:
          if os.path.exists(BLAST_OUT):
            os.remove(BLAST_OUT)
        
        logresult = False # not empty
    
    else: 
        # extract c2 reason
        if "PRIMER_PAIR_EXPLAIN" in c2: 
            msg = c2["PRIMER_PAIR_EXPLAIN"].split(",")
        else: 
            msg = c1["PRIMER_PAIR_EXPLAIN"].split(",")
        msg_f = [x for x in msg if "ok" not in x and "considered" not in x]
        
        if len(msg_f) > 1: 
            msg_f = [x for x in msg_f if "product size" not in x][0]
        else: 
            msg_f = msg_f[0]
        
        final_df = None
        blast_df = None
        logresult = "Main reason for design failure: {}".format(msg_f)
        
    return blast_df, final_df, logresult
    
    
    
if __name__ == "__main__":
    CreatePrimers()
