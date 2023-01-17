#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import os
import pandas as pd

# own modules
from ExonSurfer.ensembl import ensembl
from ExonSurfer.blast import blast
from ExonSurfer.resources import resources
from ExonSurfer.primerDesign import designPrimers, chooseTarget


def CreatePrimers(gene, transcripts = "ALL", path_out = ".", release = 108):
    """
    This function is the main function of the pipeline. It takes a gene name and
    a transcript name and returns a list of primers.
        gene [in] (str): gene name
        transcripts [in] (str): transcript name
        path_out [in] (str): path to output directory
        BLAST_OUT [out] (df): dataframe with blast results
        DESIGN_OUT [out] (df): dataframe with primer design results
    """ 
        
    # Construct transcripts dictionary
    print("Extracting ensemble info")
    data = ensembl.create_ensembl_data(release = 108)
    gene_obj = ensembl.get_gene_by_symbol(gene, data)
    d = ensembl.get_transcripts_dict(gene_obj, exclude_noncoding = True)
    
    # Get best exonic junction
    print("Getting exon junction")
    junctions_d = chooseTarget.format_junctions(d)
    junction = chooseTarget.choose_target(d, junctions_d, transcripts)
    print("Exon junctions: {}".format(junction))
    
    # Get sequence and junction index
    print("Starting primer design")
    # Define three output files, with the same name as the gene and transcript
    DESIGN_OUT = os.path.join(path_out, "{}_{}_design.txt".format(gene, transcripts))
    BLAST_OUT = os.path.join(path_out, "{}_{}_blast.txt".format(gene, transcripts))
    
    # remove in case remaining from previous design
    for file in (DESIGN_OUT, BLAST_OUT): 
        if os.path.exists(file):
            os.remove(file)

    print("Design output: {}".format(DESIGN_OUT))
    
    # number of primers to design for each junction and option
    num_primers = int(100 / (len(junction)*2))
    
    for item in junction: 
        target, index = ensembl.construct_target_cdna(resources.MASKED_SEQS(), 
                                                      gene_obj,
                                                      data, 
                                                      transcripts, 
                                                      item)
        # Design primers
        c1, c2 = designPrimers.call_primer3(target, index, num_primers)
        designPrimers.report_design(c1, c2, item, DESIGN_OUT)

    # Add primer pair identifier
    cols = ("option", "junction", "junction_description", "forward", "reverse", 
            "amplicon_size", "forward_tm", "reverse_tm", "forward_gc", "reverse_gc", 
            "amplicon_tm", "pair_penalty")
    df = pd.read_csv(DESIGN_OUT, sep = "\t", names = cols, header = None)
    df["pair_num"] = ["Pair{}".format(x) for x in range(0, df.shape[0])]
    df = df.set_index('pair_num')
    df.to_csv(DESIGN_OUT, sep = "\t")
    
    # Write blast input
    FASTA_F = resources.create_temp_fasta()
    resources.fillin_temp_fasta(DESIGN_OUT, FASTA_F)
    
    # Call blast
    blast_df = blast.run_blast_list(FASTA_F, BLAST_OUT, resources.BLAST_DB())
    
    # Delete fasta file
    os.remove(FASTA_F)
    
    # Filter blast results
    blast_df = blast.pre_filter_blast(blast_df, transcripts, gene, df)
    
    # Check blast results positions
    blast.check_specificity(blast_df, df, gene)
    
    final_df = designPrimers.penalize_final_output(df, transcripts)
    final_df.to_csv(DESIGN_OUT, sep = "\t")
    
    
if __name__ == "__main__":
    CreatePrimers()