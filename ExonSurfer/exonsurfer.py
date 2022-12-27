#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import os

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
    DESIGN_OUT = os.path.join(path_out, 
                              gene + "_" + transcripts + "_design.txt")
    BLAST_OUT = os.path.join(path_out, 
                             gene + "_" + transcripts + "_blast.txt")
    print("Design output: {}".format(DESIGN_OUT))
    for item in junction: 
        FASTA_F = resources.create_temp_fasta()
        target, index = ensembl.constructu_target_cdna(resources.MASKED_SEQS(), 
                                                       gene_obj, 
                                                       transcripts, 
                                                       item)
        
        # Design primers
        print(target)
        c1, c2 = designPrimers.call_primer3(target, index)
        dfPrimer = designPrimers.report_design(c1, c2, item, DESIGN_OUT)
        
        # Write primers to fasta file 
        designPrimers.write_blast_fasta(c1, c2, FASTA_F)
        
        # Call blast
        dfBlast = blast.run_blast_list(FASTA_F, 
                             BLAST_OUT, 
                             "blastn", 
                             resources.BLAST_DB())
        
        # Delete fasta file
        os.remove(FASTA_F)

if __name__ == "__main__":
    CreatePrimers()