#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 13:10:49 2022

@author: ecrisru
"""

# imported modules
import primer3

# my constants
from .designConfig import design_dict as DESIGN_DICT

# Constants
FLANK = 400
CASE2_PRIMERS = 100

###############################################################################
#                 primerDesign module FUNCTION DEFINITION SITE                #
###############################################################################

    
def call_primer3(target_seq, junction_i): 
    """
    This function calls primers3 for 2 design options: (1) primers are placed 
    ON the exon junction and (2) primers are placed FLANKING the exon junction. 
    Args: 
        target_seq [in] (str)      Full cDNA sequence of the transcript
        junction_i [in] (int)      Index of the exonic junction on the target_seq
        case1_primers [out] (dict) Dictionary of primers designed with option 1
        case2_primers [out] (dict) Dictionary of primers designed with option 2
    """
    
    # Case 1: place forward primer or right primer ON junction
    target_dict = {
            'SEQUENCE_ID': 'InternalID',
            'SEQUENCE_TEMPLATE': target_seq,
            'SEQUENCE_OVERLAP_JUNCTION_LIST': junction_i,
            'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION': 4,
            }

    # PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION=4]
    case1_primers = primer3.bindings.designPrimers(target_dict, DESIGN_DICT)
    
    # Case 2: place each primer on one exon
    target_dict = {
            'SEQUENCE_ID': 'InternalID',
            'SEQUENCE_TEMPLATE': target_seq,
            'SEQUENCE_TARGET': [junction_i, 1],
            }
    case2_primers = primer3.bindings.designPrimers(target_dict, DESIGN_DICT)
    
    return case1_primers, case2_primers
                
###############################################################################
                   
def report_design(c1, c2, exon_junction, ofile):    
    """
    This function writes a table with the designed primers
    Args: 
        c1 [in] (dict)             Dict of primers designed with option 1
        c2 [in] (dict)             Dict of primers designed with option 2
        exon_junction [in] (tuple) Tuple with two strings: (1) "Ensemble exon ID
                                   - Ensembl exon ID" (2) Design specification
        ofile [in] (str)           Path to the output file
    """

    out_open = open(ofile, "a+")
    for pdict in (c1, c2): 
        
        for n in range(pdict["PRIMER_PAIR_NUM_RETURNED"]): 
            
            patt = "_{}_".format(n)
            opt = "1" if pdict == c1 else "2" # design option
            
            l = (opt, 
                 exon_junction[0], 
                 exon_junction[1], 
                 pdict["PRIMER_LEFT{}SEQUENCE".format(patt)].upper(), 
                 pdict["PRIMER_RIGHT{}SEQUENCE".format(patt)].upper(), 
                 pdict["PRIMER_PAIR{}PRODUCT_SIZE".format(patt)], 
                 pdict["PRIMER_LEFT{}TM".format(patt)], 
                 pdict["PRIMER_LEFT{}TM".format(patt)], 
                 pdict["PRIMER_LEFT{}GC_PERCENT".format(patt)], 
                 pdict["PRIMER_LEFT{}GC_PERCENT".format(patt)], 
                 pdict["PRIMER_PAIR{}PRODUCT_TM".format(patt)], 
                 pdict["PRIMER_PAIR{}PENALTY".format(patt)]
                 )
            
            string = "\t".join([str(x) for x in l]) + "\n"
            
            out_open.write(string)
            
    out_open.close()

    
###############################################################################

def penalize_final_output(df, transcripts): 
    """
    This function penalizes the last data design DF (with blast information 
    appended) and returns a list of the best primer pairs for the task
    """
    if transcripts != "ALL": 
        # prioritize option 1, no other transcripts and no other genes
        if any(df.loc[(df['other_transcripts'] == "") & (df['other_genes'] == "") & (df["option"] == "1")]): 
            final_df = df.loc[(df['other_transcripts'] == "") & (df['other_genes'] == "") & (df["option"] == "1")]
        # any option, no other transcripts nor genes
        elif any(df.loc[(df['other_transcripts'] == "") & (df['other_genes'] == "")]): 
            final_df = df.loc[(df['other_transcripts'] == "") & (df['other_genes'] == "")]
        # try with no other transcripts
        elif any(df.loc[(df['other_transcripts'] == "")]): 
            final_df = df.loc[(df['other_transcripts'] == "")] 
        # try with no other genes at least
        elif any(df.loc[(df['other_genes'] == "")]): 
            final_df = df.loc[(df['other_genes'] == "")]  
        else: 
            final_df = df # whatever
    
    else: # do not need to prioritize option 1
        if any(df.loc[(df['other_transcripts'] == "") & (df['other_genes'] == "")]): 
            final_df = df.loc[(df['other_transcripts'] == "") & (df['other_genes'] == "")]
        elif any(df.loc[(df['other_transcripts'] == "")]): # try without other transcripts
            final_df = df.loc[(df['other_transcripts'] == "")] 
        elif any(df.loc[(df['other_genes'] == "")]): # try without other genes
            final_df = df.loc[(df['other_genes'] == "")]  
        else: 
            final_df = df # whatever
            
    return final_df
    
    