#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import os
import pandas as pd

# own modules
from ExonSurfer.ensembl import ensembl
from ExonSurfer.specificity import blast, annotate_detected
from ExonSurfer.resources import resources
from ExonSurfer.primerDesign import chooseTarget, construct_cdna, designPrimers
from ExonSurfer.primerDesign import penalizePrimers, designConfig


def CreatePrimers(gene, transcripts = "ALL", species = "homo_sapiens_masked",
                  release = 108, design_dict = designConfig.design_dict, 
                  path_out = ".", save_files = True, e_value = 0.8,
                  i_cutoff = 70, max_sep = 700, opt_prod_size = 200, NPRIMERS = 100):
    """
    This function is the main function of the pipeline. It takes a gene name and
    a transcript name and returns a list of primers.
    Args:
        gene [in] (str):        Gene name
        transcripts [in] (str): Transcript name or ALL
        species [in] (str):     Species for the data (human, mus_musculus or 
                                                      rattus_novegicus)
        release [in] (str):     Ensembl release number. Change with caution
        design_dict [in] (d):   Dict with primer3 parameters
        path_out [in] (str):    Path to output direct
        BLAST_OUT [out] (df):   Dataframe with blast results
        DESIGN_OUT [out] (df):  Ddataframe with primer design results
    """ 
    ###########################################################################
    #                        STEP 1. RETRIEVE INFORMATION                     #
    ########################################################################### 
    # Construct transcripts dictionary
    print("Extracting ensemble info")
    data = ensembl.create_ensembl_data(release, 
                                       species.replace("_masked", ""))
    gene_obj = ensembl.get_gene_by_symbol(gene, data)
        
    d = ensembl.get_transcripts_dict(gene_obj, exclude_noncoding = False)
    
    cdna_d = ensembl.build_cdna_dict(data, gene_obj, resources.MASKED_SEQS(species))
    
    # If ALL transcripts are targeted and human species, get canonical
    if "homo_sapiens" in species and transcripts == "ALL": 
        cfile = open(resources.CANONICAL(), "r")
        lines = cfile.read().split("\n")
        canonical_t = [  # list comprehension here
            l.split("\t")[1] 
            for l 
            in lines 
            if gene_obj.gene_id in l and "Ensembl Canonical" in l
            ]
    else: 
        canonical_t = []
        
    ###########################################################################
    #                           STEP 2: CHOOSE TARGET                         #
    ###########################################################################         
    # Get best exonic junction
    print("Getting exon junction")
    junctions_d = chooseTarget.format_junctions(d, 
                                                transcripts, 
                                                opt_prod_size, 
                                                data)
    
    junction = chooseTarget.choose_target(d, junctions_d, transcripts, canonical_t)
    print("Exon junctions: {}".format(junction))
    
    # Get sequence and junction index
    print("Starting primer design")
    
    # Create dataframe for design
    cols = ("option", "junction", "junction_description", "forward", "reverse", 
            "amplicon_size", "forward_tm", "reverse_tm", "forward_gc", "reverse_gc", 
            "amplicon_tm", "pair_penalty")
    df = pd.DataFrame(columns = cols)
    
    ###########################################################################
    #                          STEP 3:  DESIGN PRIMERS                        #
    ###########################################################################      
    if len(junction) == 0: # only one exon
        design_dict["PRIMER_NUM_RETURN"] = NPRIMERS
        target, index, elen = construct_cdna.construct_one_exon_cdna(resources.MASKED_SEQS(species), 
                                                                     gene_obj, 
                                                                     data, 
                                                                     transcripts)        
        # Design primers
        c2 = designPrimers.call_primer3(target, index, design_dict, enum = 1)
        if transcripts == "ALL": 
            item = [ensembl.get_transcript_from_gene(gene_obj)[0].exons[0].exon_id, 
                    "one exon"]            
        else: 
            item = [data.transcript_by_id(transcripts).exons[0].exon_id, 
                    "one exon"]
        df = designPrimers.report_one_exon_design(c2, elen, item, df)
        
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
            num_primers = int(NPRIMERS / (len(to_design)*2))
            design_dict["PRIMER_NUM_RETURN"] = num_primers
            
            # Design primers
            c1, c2 = designPrimers.call_primer3(tupla[0], tupla[1], design_dict)
    
            df = designPrimers.report_design(c1, c2, tupla[2], tupla[3], 
                                             tupla[4], df)
    
    df["pair_num"] = ["Pair{}".format(x) for x in range(0, df.shape[0])]
    df = df.set_index('pair_num')
    
    if df.shape[0] > 1: # we designed primers
    
        #######################################################################
        #                       STEP 4:  BLAST AGAINST CDNA                   #
        #######################################################################          
        # Define three output files, with the same name as the gene and transcript
        if transcripts != "ALL": 
            t_string = "-".join(transcripts[:3])
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
                                        species)
        
        # Delete fasta file
        if os.path.exists(FASTA_F): os.remove(FASTA_F)
        
        # Filter blast results
        blast_df = blast.pre_filter_blast(blast_df, transcripts, gene, df, 
                                          e_value, i_cutoff, False)
        blast_df, df = blast.filter_big_blast(blast_df, df)
        
        # Check blast results positions
        df = blast.check_specificity(blast_df, df, gene, transcripts, max_sep)
        
        #######################################################################
        #                        STEP 5:  FILTER DATAFRAME                    #
        #######################################################################   
        final_df = penalizePrimers.penalize_final_output(df, transcripts, data, 
                                                         gene_obj)
        
        #######################################################################
        #                         STEP 6:  GENOMIC BLAST                      #
        #######################################################################           
        # Write genomic blast input
        resources.fillin_temp_fasta(final_df, FASTA_F)
        
        # Call genomic blast
        blast_df = blast.run_blast_list(FASTA_F, GBLAST_OUT, 
                                        resources.BLAST_GENOMIC_DB(species), 
                                        species, tomerge = False)    
        # Filter big blast if needed
        blast_df = blast.pre_filter_blast(blast_df, transcripts, gene, df, 
                                          e_value, i_cutoff, False)
        blast_df, final_df = blast.filter_big_blast(blast_df, final_df)
        
        # Delete fasta file
        if os.path.exists(FASTA_F): os.remove(FASTA_F)
        
        # Check blast results positions
        final_df = blast.check_genomic_specificity(blast_df, final_df, max_sep)

        #######################################################################
        #                      STEP 7: PREPARE FINAL RESULT                   #
        #######################################################################          
        # Make penalty score
        final_df = penalizePrimers.make_penalty_score(final_df)     
        final_df = penalizePrimers.genomic_filter(final_df)
        
        # Annotate transcripts detected
        final_df = annotate_detected.annotate_notdetected(final_df, cdna_d, gene_obj)
        
        # Annotate if there are off_targets or not
        final_df = blast.show_off_targets(final_df, transcripts)
        
        if save_files == True:
            final_df.to_csv(DESIGN_OUT, sep = "\t")
        else:
          if os.path.exists(BLAST_OUT):
            os.remove(BLAST_OUT)
        
        logresult = False # not empty
    
    else: 
        # extract c2 reason
        msg = c2["PRIMER_PAIR_EXPLAIN"].split(",")
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
