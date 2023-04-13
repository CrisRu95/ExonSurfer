#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import argparse 

# own modules
from ExonSurfer.exonsurfer import CreatePrimers
from ExonSurfer.primerDesign import designConfig
from ExonSurfer.resources import resources

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
    
    parser.add_argument("--download", "-db", dest = "download", default = False, 
                        help = "Download all DB") 

    parser.add_argument("--numprimers", "-np", dest = "NPRIMERS", default = 200, 
                        help = "Number of primers") 

    parser.add_argument("--opt_prod_size", "-ops", dest = "opt_prod_size", default = 200, 
                        help = "Optimal product size") 

    parser.add_argument("--des_option", "-opt", dest = "d_option", default = 1, 
                        help = "1 if only primers on junctions") 

    parser.add_argument("--min_3_overlap", "-o3", dest = "min_3_overlap", default = 5, 
                        help = "Minimum overlap for the 3' end of the primers") 

    parser.add_argument("--min_5_overlap", "-o5", dest = "min_5_overlap", default = 6, 
                        help = "Minimum overlap for te 5' end of the primer") 
             
    args = parser.parse_args()
    
    return args


def main(): 
    
    args = arguments_parser()
    
    if args.design_dict is None:
        args.design_dict = designConfig.design_dict
        
    if args.download:
        resources.download_all_db()
    # Construct transcripts dictionary
    try:
        CreatePrimers(gene = args.gene, 
                      transcripts = args.transcript, 
                      species = args.species, 
                      path_out = args.path_out,
                      release = int(args.release),
                      design_dict = args.design_dict,
                      save_files = bool(args.save_files),
                      e_value = float(args.e_value),
                      i_cutoff = int(args.i_cutoff),
                      max_sep = int(args.max_sep), 
                      NPRIMERS = int(args.NPRIMERS), 
                      opt_prod_size = int(args.opt_prod_size), 
                      d_option = int(args.d_option), 
                      min_3_overlap = int(args.min_3_overlap), 
                      min_5_overlap = int(args.min_5_overlap),                       
                      )
    except Exception as error:
        print("[!] Error: ", error)
        pass

if __name__ == "__main__":
    main()
