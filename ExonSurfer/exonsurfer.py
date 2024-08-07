#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import os
import pandas as pd

# own modules
from ExonSurfer.ensembl import ensembl, extractCanonical
from ExonSurfer.specificity import blast, offtargets, g_offtargets, annotate
from ExonSurfer.resources import resources
from ExonSurfer.primerDesign import chooseTarget, construct_cdna, designPrimers
from ExonSurfer.primerDesign import penalizePrimers, designConfig
from ExonSurfer.dimers import dimers

def CreatePrimers(gene, transcripts = "ALL", species = "homo_sapiens_masked",
                  release = 108, design_dict = designConfig.design_dict, 
                  path_out = ".", save_files = True, e_value = 0.8,
                  i_cutoff = 70, max_sep = 700, opt_prod_size = 200, NPRIMERS = 200, 
                  d_option = 1, min_3_overlap = 5, min_5_overlap = 6):
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
    #                        STEP 1. RETRIEVE INFORMATION                     #
    ########################################################################### 
    # Construct transcripts dictionary
    print("Extracting ensemble info")
    data = ensembl.create_ensembl_data(release, 
                                       species.replace("_masked", ""))
    try:
      gene_obj = ensembl.get_gene_by_symbol(gene, data)
    except:
      gene_obj = ensembl.get_gene_by_id(gene, data)

        
    d = ensembl.get_transcripts_dict(gene_obj, exclude_noncoding = False)
    
    cdna_d = ensembl.build_cdna_dict(data, gene_obj, resources.MASKED_SEQS(species))
    
    # Reverse transcripts list (it goes in reverse to importance)
    if transcripts != "ALL": 
        transcripts = sorted(transcripts)
    
    # If ALL transcripts are targeted and human species, get canonical
    if "homo_sapiens" in species: 
        canonical_t = extractCanonical.extract_canonical(gene_obj)
    else: 
        canonical_t = []
    
    # Put canonical transcript in the most important position
    if transcripts != "ALL" and canonical_t != [] and canonical_t in transcripts: 
        transcripts.remove(canonical_t)
        transcripts = [canonical_t] + transcripts
    
    # Set optimum junction length  
    opt_junc_len = design_dict["PRIMER_PRODUCT_SIZE_RANGE"][0][1]
    
    ###########################################################################
    #                           STEP 2: CHOOSE TARGET                         #
    ###########################################################################      
    junction = [0] # initialize for loop
    
    # Get best exonic junction
    print("Getting exon junction")
    junctions_d = chooseTarget.format_junctions(d, 
                                                transcripts, 
                                                opt_junc_len, 
                                                data, 
                                                canonical_t)
    
    if transcripts == "ALL" and canonical_t != []: 
        junction = chooseTarget.choose_canonical(junctions_d, canonical_t)

    else: 
        junction = chooseTarget.choose_target(d, 
                                            junctions_d, 
                                            transcripts, 
                                            canonical_t)
        
    print("Exon junctions: {}".format(junction))
    
    # Create dataframe for design
    cols = ("option", "junction", "junction_description", "forward", "reverse", 
            "amplicon_size", "forward_tm", "reverse_tm", "forward_gc", "reverse_gc", 
            "amplicon_tm", "pair_penalty")
    df = pd.DataFrame(columns = cols)
    
    # Check if one exon design
    if len(junctions_d) == 1 and "_" not in junctions_d[list(junctions_d.keys())[0]]: 
        enum = 1
        i = 0
        WINDOW = 100
        target_len = 100 # random value to initialize
    else: 
        enum, i, WINDOW = False, False, False # random values
    
    while (len(junction) > 0 and df.empty) or (enum == 1 and i*WINDOW < target_len): 
        
        #######################################################################
        #                          STEP 3:  DESIGN PRIMERS                    #
        #######################################################################      
        if enum == 1: # only one exon

            target, index, elen = construct_cdna.construct_one_exon_cdna(resources.MASKED_SEQS(species), 
                                                                         gene_obj, 
                                                                         data, 
                                                                         transcripts, i, WINDOW)
            
            # Define PRIMER_NUM_RETURN depending on target length
            num_primers = int(NPRIMERS / (len(target)/WINDOW)) * 2
            design_dict["PRIMER_NUM_RETURN"] = num_primers
            
            if len(target)-50 < index: 
                index = int(len(target) / 2)
            # Design primers
            c2 = designPrimers.call_primer3(target, index, design_dict, min_3_overlap, 
                                            min_5_overlap, d_option, enum = 1)
            if transcripts == "ALL": 
                item = [ensembl.get_transcript_from_gene(gene_obj)[0].exons[0].exon_id, 
                        "one exon"]            
            else: 
                item = [data.transcript_by_id(transcripts[0]).exons[0].exon_id, 
                        "one exon"]
            df = designPrimers.report_one_exon_design(c2, elen, item, df)
            
            target_len = len(target) # keep looping?
            i += 1 # keep looping?
            
            
        else: # Normal design (more than 1 exon)
            to_design = []
            for item in junction: 
                to_design += construct_cdna.construct_target_cdna(resources.MASKED_SEQS(species), 
                                                                  gene_obj,
                                                                  data, 
                                                                  transcripts, 
                                                                  item)
            for tupla in to_design: 
                # Decide number of primers to design
                if d_option == 1: 
                    num_primers = int(NPRIMERS / len(to_design))
                else: 
                    num_primers = int(NPRIMERS / (len(to_design)*2))
                # avoid 0 primer design
                if num_primers == 0: 
                    num_primers = 1
                design_dict["PRIMER_NUM_RETURN"] = num_primers
                
                # Design primers
                c1, c2 = designPrimers.call_primer3(tupla[0], tupla[1], design_dict, 
                                                    min_3_overlap, min_5_overlap, d_option)
                df = designPrimers.report_design(c1, c2, tupla[2], tupla[3], 
                                                 tupla[4], df)
                
                if df.empty: # no designed primers
                    # remove the chosen junction from the dict and try again
                    del junctions_d[junction[0][0]]
                    junction = chooseTarget.choose_target(d, 
                                      junctions_d, 
                                      transcripts, 
                                      canonical_t)
    
    if df.shape[0] > 1: # we designed primers
        
        # Calculate dimers
        df["dimers"] = df.apply(lambda row: dimers.get_max_comp(row["forward"],row["reverse"]),axis=1)
        # Filter dimers
        df = df[df['dimers'] < 5]   
        #Annotate pair number
        df["pair_num"] = ["Pair{}".format(x) for x in range(0, df.shape[0])]
    
        df = df.set_index('pair_num')
    
        #######################################################################
        #                       STEP 4:  BLAST AGAINST CDNA                   #
        #######################################################################          
        # Define three output files, with the same name as the gene and transcript
        if transcripts != "ALL": 
            t_string = "_".join(transcripts[:3])
        else: 
            t_string = transcripts
        
        DESIGN_OUT = os.path.join(path_out, 
                                  "{}_{}_design.txt".format(gene, t_string))
        BLAST_OUT = os.path.join(path_out, 
                                 "{}_{}_blast.txt".format(gene, t_string))
        GBLAST_OUT = os.path.join(path_out, 
                                  "{}_{}_gblast.txt".format(gene, t_string))
   
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
        blast_df, df = blast.filter_big_blast(gene, blast_df, df)
        
        # Check blast results positions
        df = offtargets.check_specificity(blast_df, df, gene_obj.name, transcripts, max_sep)
        
        #######################################################################
        #                        STEP 5:  FILTER DATAFRAME                    #
        #######################################################################   
        final_df = penalizePrimers.penalize_final_output(df, transcripts, data)
        
        #######################################################################
        #                         STEP 6:  GENOMIC BLAST                      #
        #######################################################################           
        # Write genomic blast input
        resources.fillin_temp_fasta(final_df, FASTA_F)
        
        # Call genomic blast
        gblast_df = blast.run_blast_list(FASTA_F, GBLAST_OUT, 
                                         resources.BLAST_GENOMIC_DB(species), 
                                         species, i_cutoff, e_value, tomerge = False, 
                                         genomic = True)    
        # Filter big blast if needed
        gblast_df, final_df = blast.filter_big_gblast(gblast_df, final_df)
        
        # Delete fasta file
        if os.path.exists(FASTA_F): os.remove(FASTA_F)
        
        # Check blast results positions
        final_df = g_offtargets.check_genomic_specificity(gblast_df, final_df, max_sep)
        
        # Filter according to genomic specificity
        final_df = penalizePrimers.genomic_filter(final_df)

        #######################################################################
        #                      STEP 7: PREPARE FINAL RESULT                   #
        #######################################################################          
        # Make penalty score
        final_df = penalizePrimers.make_penalty_score(final_df)     
        
        # Annotate transcripts detected
        final_df = annotate.annotate_notdetected(final_df, cdna_d, gene_obj)
        
        # Annotate if there are off_targets or not
        final_df = annotate.show_off_targets(final_df, transcripts)
        
        # Remove files if necessary
        if save_files:
            final_df.to_csv(DESIGN_OUT, sep = "\t")
            blast_df.to_csv(BLAST_OUT, sep = "\t")
            gblast_df.to_csv(GBLAST_OUT, sep = "\t")
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
        gblast_df = None
        logresult = "Main reason for design failure: {}".format(msg_f)
        
    return blast_df, gblast_df, final_df, logresult
    
    
    
if __name__ == "__main__":
    CreatePrimers()
