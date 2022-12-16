#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import argparse 

# own modules
from ExonSurfer.exonsurfer import CreatePrimers


def arguments_parser(): 
    """
    This function parses arguments
    """
    parser = argparse.ArgumentParser(description = "Design primers in the specified dir")
    parser.add_argument("--gene", "-g", dest = "gene", 
                        help = "Name of the gene to target (ENSEMBL ID)")
    parser.add_argument("--transcript", "-t", dest = "transcript",default="ALL", 
                        help = "Ensembl transcript identifier or ALL")
    parser.add_argument("--out", "-o", dest = "path_out",default=".", 
                        help = "Out directory")                   
    args = parser.parse_args()
    
    return args


def main(): 
    
    args = arguments_parser()
    
    # Construct transcripts dictionary
    gene, transcripts, path_out = args.gene, args.transcript, args.path_out
    CreatePrimers(gene, transcripts, path_out)

if __name__ == "__main__":
    main()