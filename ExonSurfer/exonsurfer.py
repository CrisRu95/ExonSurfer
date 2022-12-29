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


def CreatePrimers(gene, transcripts = "ALL", path_out = "."):
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
    gene_obj = ensembl.get_gene_by_symbol(gene)
    d = ensembl.get_transcripts_dict(gene_obj, exclude_coding = False)
    
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

    print("Design output: {}".format(DESIGN_OUT))
    
    for item in junction: 
        
        target, index = ensembl.constructu_target_cdna(resources.MASKED_SEQS(), 
                                                       gene_obj, 
                                                       transcripts, 
                                                       item)
        # Design primers
        c1, c2 = designPrimers.call_primer3(target, index)
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
    blast_df = blast.run_blast_list(FASTA_F, BLAST_OUT, resources.FILT_BLAST_DB())
    
    # Delete fasta file
    os.remove(FASTA_F)
    
    # Filter blast results
    blast_df = blast.pre_filter_blast(blast_df, transcripts, df)
    
    # Check blast results positions
    blast.check_specificity(blast_df, df, gene)


if __name__ == "__main__":
    CreatePrimers()