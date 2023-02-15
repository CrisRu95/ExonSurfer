#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import argparse 

# own modules
from ExonSurfer.exonsurfer import CreatePrimers
from ExonSurfer.primerDesign import designConfig

def arguments_parser(): 
    """
    This function parses arguments
    """
    parser = argparse.ArgumentParser(description = "Design primers in the specified dir")
    parser.add_argument("--species", "-sp", dest = "species",
                       help = "Species of the gene", default = "homo_sapiens")
    
    parser.add_argument("--gene", "-g", dest = "gene", 
                        help = "Name of the gene to target (ENSEMBL ID)")
    
    parser.add_argument("--transcript", "-t", dest = "transcript", default="ALL", 
                        help = "Ensembl transcript identifier or ALL")
    
    parser.add_argument("--out", "-o", dest = "path_out", default = ".", 
                        help = "Out directory")

    parser.add_argument("--release", "-r", dest = "release", default = 108,
                        help = "Ensembl release")
    
    parser.add_argument("--design_dict", "-d", dest = "design_dict", 
                        help = "Dictionary with design parameters", default = None)
    
    parser.add_argument("--files", "-f", dest = "save_files", default = True, 
                        help = "True if to save files")

    parser.add_argument("--e_value", "-ev", dest = "e_value", default = 0.8, 
                        help = "E value cutoff for the blast search")

    parser.add_argument("--i_cutoff", "-i_c", dest = "i_cutoff", default = 70, 
                        help = "Identity cutoff for the blast search")

    parser.add_argument("--max_sep", "-ms", dest = "max_sep", default = 1000, 
                        help = "Maximum sep between alignments to consider amplicon")
             
    args = parser.parse_args()
    
    return args


def main(): 
    
    args = arguments_parser()
    
    if args.design_dict is None:
        args.design_dict = designConfig.design_dict
    # Construct transcripts dictionary
    CreatePrimers(gene = args.gene, 
                  transcripts = args.transcript, 
                  species = args.species, 
                  path_out = args.path_out,
                  release = args.release,
                  design_dict = args.design_dict,
                  save_files = args.save_files,
                  e_value = args.e_value,
                  i_cutoff = args.i_cutoff,
                  max_sep = args.max_sep
                  )

if __name__ == "__main__":
    main()
