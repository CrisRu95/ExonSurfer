#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 13:10:49 2022

@author: ecrisru
"""

# imported modules
import primer3
import pandas as pd

# Constants
FLANK = 400
CASE2_PRIMERS = 100

###############################################################################
#                 primerDesign module FUNCTION DEFINITION SITE                #
###############################################################################

    
def call_primer3(target_seq, junction_i, design_dict, enum = 2): 
    """
    This function calls primers3 for 2 design options: (1) primers are placed 
    ON the exon junction and (2) primers are placed FLANKING the exon junction. 
    Args: 
        target_seq [in] (str)      Full cDNA sequence of the transcript
        junction_i [in] (int)      Index of the exonic junction on the target_seq
        design_dict [in] (str)     Dict of arguments for primer3
        case1_primers [out] (dict) Dictionary of primers designed with option 1
        case2_primers [out] (dict) Dictionary of primers designed with option 2
    """
    if enum == 1: # Only one exon per transcript (no juntion)
        target_dict = {
                'SEQUENCE_ID': 'InternalID',
                'SEQUENCE_TEMPLATE': target_seq,
                'SEQUENCE_TARGET': [junction_i, 1],
                }
        case2_primers = primer3.bindings.designPrimers(target_dict, design_dict)
        
        return case2_primers
        
    else:
        # Case 1: place forward primer or right primer ON junction
        target_dict = {
                'SEQUENCE_ID': 'InternalID',
                'SEQUENCE_TEMPLATE': target_seq,
                'SEQUENCE_OVERLAP_JUNCTION_LIST': junction_i,
                'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION': 4,
                }
    
        case1_primers = primer3.bindings.designPrimers(target_dict, design_dict)
        
        # Case 2: place each primer on one exon
        target_dict = {
                'SEQUENCE_ID': 'InternalID',
                'SEQUENCE_TEMPLATE': target_seq,
                'SEQUENCE_TARGET': [junction_i, 1],
                }
        case2_primers = primer3.bindings.designPrimers(target_dict, design_dict)
    
        return case1_primers, case2_primers
                
###############################################################################
                   
def report_design(c1, c2, exon_junction, df):    
    """
    This function writes a table with the designed primers
    Args: 
        c1 [in] (dict)             Dict of primers designed with option 1
        c2 [in] (dict)             Dict of primers designed with option 2
        exon_junction [in] (tuple) Tuple with two strings: (1) "Ensemble exon ID
                                   - Ensembl exon ID" (2) Design specification
        df [in|out] (ps.df)        Pandas df with the designed primers info
    """
    for pdict in (c1, c2): 
        
        for n in range(pdict["PRIMER_PAIR_NUM_RETURNED"]): 
            
            patt = "_{}_".format(n)            
            row = {"option": 1 if pdict == c1 else 2, # design option, 
                   "junction": exon_junction[0], 
                   "junction_description": exon_junction[1], 
                   "forward": pdict["PRIMER_LEFT{}SEQUENCE".format(patt)].upper(), 
                   "reverse": pdict["PRIMER_RIGHT{}SEQUENCE".format(patt)].upper(), 
                   "amplicon_size": pdict["PRIMER_PAIR{}PRODUCT_SIZE".format(patt)], 
                   "forward_tm": pdict["PRIMER_LEFT{}TM".format(patt)], 
                   "reverse_tm": pdict["PRIMER_LEFT{}TM".format(patt)], 
                   "forward_gc": pdict["PRIMER_LEFT{}GC_PERCENT".format(patt)],
                   "reverse_gc": pdict["PRIMER_LEFT{}GC_PERCENT".format(patt)], 
                   "amplicon_tm": pdict["PRIMER_PAIR{}PRODUCT_TM".format(patt)], 
                   "pair_penalty": pdict["PRIMER_PAIR{}PENALTY".format(patt)]}
            
            df = pd.concat([df, pd.DataFrame([row])], ignore_index = True)
            
    return df

###############################################################################

def report_one_exon_design(c2, exon_junction, df):    
    """
    This function writes a table with the designed primers, only to use for one
    exon designs. 
    Args: 
        c2 [in] (dict)             Dict of primers
        exon_junction [in] (tuple) Tuple with two strings: (1) "Ensemble exon ID
                                   - Ensembl exon ID" (2) Design specification
        df [in|out] (ps.df)        Pandas df with the designed primers info
    """
    

    for n in range(c2["PRIMER_PAIR_NUM_RETURNED"]): 
        
        patt = "_{}_".format(n)
        
        row = {"option": 2, 
               "junction": exon_junction[0], 
               "junction_description": exon_junction[1], 
               "forward": c2["PRIMER_LEFT{}SEQUENCE".format(patt)].upper(), 
               "reverse": c2["PRIMER_RIGHT{}SEQUENCE".format(patt)].upper(), 
               "amplicon_size": c2["PRIMER_PAIR{}PRODUCT_SIZE".format(patt)], 
               "forward_tm": c2["PRIMER_LEFT{}TM".format(patt)], 
               "reverse_tm": c2["PRIMER_LEFT{}TM".format(patt)], 
               "forward_gc": c2["PRIMER_LEFT{}GC_PERCENT".format(patt)],
               "reverse_gc": c2["PRIMER_LEFT{}GC_PERCENT".format(patt)], 
               "amplicon_tm": c2["PRIMER_PAIR{}PRODUCT_TM".format(patt)], 
               "pair_penalty": c2["PRIMER_PAIR{}PENALTY".format(patt)]}
        
        df = pd.concat([df, pd.DataFrame([row])], ignore_index = True)
    
    return df