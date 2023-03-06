#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 11:22:50 2023

@author: q0791lp
"""
# imported module
import re

# own modules
from ExonSurfer.resources import resources
from ExonSurfer.ensembl import ensembl

###############################################################################
#                annotate detected FUNCTION DEFINITION SITE                   #
###############################################################################

def search_primers(seq, forward, reverse): 
    """
    This function searches primers in a sequence. 
    Args: 
        seq [in] (str)      Transcript CDNA
        forward [in] (str) Forward primer sequence (all caps)
        reverse [in] (str) Reverse primer sequence (all caps)
        found [out] (bool)  True if both primers in sequence
    """
    # make reverse complement
    forward_rc = resources.reverse_complement(forward)
    reverse_rc = resources.reverse_complement(reverse)
    
    # forward found
    f_found = True if re.search(forward, seq, re.I) or re.search(forward_rc, seq, re.I) else False

    # forward found
    r_found = True if re.search(reverse, seq, re.I) or re.search(reverse_rc, seq, re.I) else False
    
    found = all((f_found, r_found))
    
    return found
        
###############################################################################

def return_detected_list(cdna_d, transcript_l, forward, reverse, 
                         other_transcripts_str): 
    """
    This function searches the transcripts that a primer pair amplifies perfectly
    (allowing 0 mismatches) + the identifiers found by BLAST. 
    Args: 
        cdna_d [in] (d)       Dict with trans ids as keys and seqs as values
        transcript_l [in] (l)     List of transcript ids
        forward [in] (str)        Forward primer sequence (all caps)
        reverse [in] (str)        Reverse primer sequence (all caps)
        other_transcripts_str [in] (str) "other_transcripts" column info
        unique_det [out] (str)        Transcript ids returned as joined list
    """
    # STEP 1. Check detected by 100% match of primers
    detected = [] # return detected list
    
    for t_id in transcript_l: 
        tseq = cdna_d[t_id]
        
        if search_primers(tseq, forward, reverse): 
            detected.append(t_id)
    
    # STEP 2. Sum with BLAST output
    # clean other_transcripts list
    other_trans_clean = [x.replace("(protein_coding)", "") for x in other_transcripts_str.split(";")]
    
    # build complete list
    c_list = detected + other_trans_clean
    
    # unique detected
    unique_det = ";".join(list(set([x for x in c_list if x != ""])))
 
    return unique_det

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
    
def annotate_detected(final_df, cdna_d, gene_obj): 
    """
    This function annotates the detected and not detected transcripts in the 
    dataframe. Perfect matches to consider detected. 
    Args: 
        final_df [in] (pd.df) Design dataframe
        cdna_d [in] (d)       Dict with trans ids as keys and seqs as values
        gene_obj [in] (Gene obj) Pyensembl gene object
    """
    tlist = ensembl.get_transcript_from_gene(gene_obj, only_id = True)
    
    final_df["detected"] = final_df.apply(lambda row: return_detected_list(cdna_d, 
                                                                           tlist, 
                                                                            row["forward"], 
                                                                            row["reverse"], 
                                                                            row["other_transcripts"]), 
                                          axis = 1)
    
    final_df["not_detected"] = final_df.apply(lambda row: return_not_detected_list(tlist, 
                                                                                   row["detected"]), 
                                              axis = 1)    
    return final_df