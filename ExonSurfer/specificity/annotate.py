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

###############################################################################

def show_ot_for_pair(transcripts, other_genes, other_transcripts, genomic_amp): 
    """
    This function says if a primer pair has off-targets to return or not. 
    Args: 
        transcripts [in] (str|l)      List of transcripts or ALL
        other_genes [in] (str)        List of other genes or empty
        other_genes [in] (str)        List of other transcripts or empty
        genomic_amp [in] (str)        List of possible genomic amplification
        to_return [out] (int)         1 if there are off_targets, 0 if not
    """
    to_return = 0 #  default not show 
    
    if transcripts == "ALL": 
        if other_genes != "" or genomic_amp != "": 
            to_return = 1
    else: 
        if other_genes != "" or other_transcripts != "" or genomic_amp != "":
            to_return = 1
            
    return to_return

###############################################################################

def show_off_targets(final_df, transcripts): 
    """
    This function says if a primer pair has off-targets to return or not. 
    Args: 
        final_df [in] (pd.df)  Final pandas dataframe with the design
        transcripts [in] (str|l)      List of transcripts or ALL
    """
    final_df["off_targets"] = final_df.apply(lambda row: show_ot_for_pair(transcripts, 
                                                                          row["other_genes"], 
                                                                          row["other_transcripts"], 
                                                                          row["genomic_amp"]), 
                                             axis = 1)
    
    return final_df