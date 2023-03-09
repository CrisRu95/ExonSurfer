#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 11:22:50 2023

@author: q0791lp
"""
# own modules
from ExonSurfer.ensembl import ensembl

###############################################################################
#                annotate detected FUNCTION DEFINITION SITE                   #
###############################################################################

def return_not_detected_list(transcript_l, detected_str): 
    """
    This function returns the primers that a primer pair does NOT amplify perfectly. 
    Args: 
        transcript_l [in] (l) List of transcript ids
        detected_str [in] (str)   Detected transcripts
        not_detected [out] (str)  Not detected transcripts
    """
    not_detected = []
    for item in transcript_l: 
        if item not in detected_str: 
            not_detected.append(item)
    
    return ";".join(not_detected)
    
###############################################################################
#                      annotate_detected MAIN FUNCTION                        #
###############################################################################
    
def annotate_notdetected(final_df, cdna_d, gene_obj): 
    """
    This function annotates the detected and not detected transcripts in the 
    dataframe. Perfect matches to consider detected. 
    Args: 
        final_df [in] (pd.df) Design dataframe
        cdna_d [in] (d)       Dict with trans ids as keys and seqs as values
        gene_obj [in] (Gene obj) Pyensembl gene object
    """
    tlist = ensembl.get_transcript_from_gene(gene_obj, only_id = True)
    
    
    final_df["not_detected"] = final_df.apply(lambda row: return_not_detected_list(tlist, 
                                                                                   row["detected"]), 
                                              axis = 1)    
    return final_df