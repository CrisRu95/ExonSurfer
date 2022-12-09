#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import argparse 

# own modules
from ensembl import ensembl
from blast import blast
from resources import resources
from primerDesign import designPrimers, chooseTarget


def arguments_parser(): 
    """
    This function parses arguments
    """
    parser = argparse.ArgumentParser(description = "Design primers in the specified dir")
    parser.add_argument("--gene", "-g", dest = "gene", 
                        help = "Name of the gene to target (ENSEMBL ID)")
    parser.add_argument("--transcript", "-t", dest = "transcript", 
                        help = "Ensembl transcript identifier or ALL")
    args = parser.parse_args()
    
    return args


def pipe(): 
    
    args = arguments_parser()
    
    # Construct transcripts dictionary
    print("Extracting ensemble info")
    gene_obj = ensembl.get_gene_by_symbol(args.gene)
    d = ensembl.get_transcripts_dict(gene_obj, exclude_coding = False)
    
    # Get best exonic junction
    print("Getting exon junction")
    junctions_d = chooseTarget.format_junctions(d)
    junction = chooseTarget.choose_target(d, junctions_d, args.transcript)
    print("Exon junctions: {}".format(junction))
    
    # Get sequence and junction index
    print("Starting primer design")
    for item in junction: 
        
        target, index = ensembl.constructu_target_cdna(resources.MASKED_SEQS, 
                                                       gene_obj, 
                                                       args.transcript, 
                                                       item)
        
        # Design primers
        c1, c2 = designPrimers.call_primer3(target, index)
        designPrimers.report_design(c1, c2, item, resources.DESIGN_OUT)
        # Write primers to fasta file 
        designPrimers.write_blast_fasta(c1, c2, resources.FASTA_F)
        
        # Call blast
        blast.run_blast_list(resources.FASTA_F, 
                             resources.BLAST_OUT, 
                             resources.BLAST_PATH, 
                             resources.BLAST_DB)

if __name__ == "__main__":
    pipe()