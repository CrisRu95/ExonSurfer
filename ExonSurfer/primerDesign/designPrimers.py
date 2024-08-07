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

def extract_error_message(e): 
    """
    This function extracts the reason for the primer failed design from the error
    message. 
    Args: 
        e [in] (error)  Exception as error
        msg [out] (str) Reason for design failure
    """
    msg = getattr(e, 'message', repr(e)).replace("OSError", "")
    msg = msg.replace("(", "")
    msg = msg.replace("'", "")
    msg = msg.replace("(", "")
    return msg

###############################################################################
    
def call_primer3(target_seq, junction_i, design_dict, min_3_overlap, min_5_overlap, 
                 d_option = 1, enum = 2): 
    """
    This function calls primers3 for 2 design options: (1) primers are placed 
    ON the exon junction and (2) primers are placed FLANKING the exon junction. 
    Args: 
        target_seq [in] (str)      Full cDNA sequence of the transcript
        junction_i [in] (int)      Index of the exonic junction on the target_seq
        design_dict [in] (str)     Dict of arguments for primer3
        enum [in] (int)            Exon number, 1 if only 1 exon per transcript 
        d_option [in] (int)        1 if all primers ON junction, "ALL" otherwise
        case1_primers [out] (dict) Dictionary of primers designed with option 1
        case2_primers [out] (dict) Dictionary of primers designed with option 2
    """
    
    if enum == 1: # Only one exon per transcript (no juntion)
        target_dict = {
                'SEQUENCE_ID': 'InternalID',
                'SEQUENCE_TEMPLATE': target_seq,
                'SEQUENCE_TARGET': [junction_i, 1],
                }
        try: 
            case2_primers = primer3.bindings.design_primers(target_dict, design_dict)
        except OSError as e: 
            case2_primers = {"PRIMER_PAIR_NUM_RETURNED":0}
            case2_primers["PRIMER_PAIR_EXPLAIN"] = extract_error_message(e)
            
        return case2_primers
        
    else:
        # Case 1: place forward primer or right primer ON junction
        target_dict = {
                'SEQUENCE_ID': 'InternalID',
                'SEQUENCE_TEMPLATE': target_seq,
                'SEQUENCE_OVERLAP_JUNCTION_LIST': junction_i,
                'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION': min_3_overlap,
                'PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION': min_5_overlap
                }
        
        try: 
            case1_primers = primer3.bindings.design_primers(target_dict, 
                                                            design_dict)
        except OSError as e: 
            case1_primers = {"PRIMER_PAIR_NUM_RETURNED":0}
            case1_primers["PRIMER_PAIR_EXPLAIN"] = extract_error_message(e)
        
        
        if d_option == 1:          
            case2_primers = {}
            case2_primers["PRIMER_PAIR_NUM_RETURNED"] = 0
        else: # if design option != 1, design also case 2 primers
            target_dict = {
                    'SEQUENCE_ID': 'InternalID',
                    'SEQUENCE_TEMPLATE': target_seq,
                    'SEQUENCE_TARGET': [junction_i, 1],
                    }
            
            try: 
                case2_primers = primer3.bindings.design_primers(target_dict, 
                                                                design_dict)
            except OSError as e: 
                case2_primers = {"PRIMER_PAIR_NUM_RETURNED":0}
                case2_primers["PRIMER_PAIR_EXPLAIN"] = extract_error_message(e)
           
            
        return case1_primers, case2_primers

###############################################################################
           
def report_design(c1, c2, exon_len, junction_id, junction_info, df):    
    """
    This function writes a table with the designed primers
    Args: 
        c1 [in] (dict)             Dict of primers designed with option 1
        c2 [in] (dict)             Dict of primers designed with option 2
        exon_len [in] (l)          List of tuples, (exon_id, exon length)
        exon_junction [in] (tuple) Tuple with two strings: (1) "Ensemble exon ID
                                   _ Ensembl exon ID" (2) Design specification
        df [in|out] (ps.df)        Pandas df with the designed primers info
    """
    for pdict in (c1, c2): 
        
        for n in range(pdict["PRIMER_PAIR_NUM_RETURNED"]): 
            
            patt = "_{}".format(n)            
            row = {"option": 1 if pdict == c1 else 2, # design option, 
                   "junction_description": junction_info, 
                   "forward": pdict["PRIMER_LEFT{}_SEQUENCE".format(patt)].upper(), 
                   "reverse": pdict["PRIMER_RIGHT{}_SEQUENCE".format(patt)].upper(), 
                   "amplicon_size": pdict["PRIMER_PAIR{}_PRODUCT_SIZE".format(patt)], 
                   "forward_tm": pdict["PRIMER_LEFT{}_TM".format(patt)], 
                   "reverse_tm": pdict["PRIMER_RIGHT{}_TM".format(patt)], 
                   "forward_gc": pdict["PRIMER_LEFT{}_GC_PERCENT".format(patt)],
                   "reverse_gc": pdict["PRIMER_RIGHT{}_GC_PERCENT".format(patt)], 
                   "amplicon_tm": pdict["PRIMER_PAIR{}_PRODUCT_TM".format(patt)], 
                   "pair_penalty": pdict["PRIMER_PAIR{}_PENALTY".format(patt)]}
            
            # get position of primers
            forpos, forlen = pdict["PRIMER_LEFT{}".format(patt)]
            revpos, revlen = pdict["PRIMER_RIGHT{}".format(patt)]
            if pdict == c1: 
                
                forex1 = [e[0] for e in exon_len if forpos in e[1]][0] # start of primers
                forex2 = [e[0] for e in exon_len if forpos+forlen in e[1]][0] # end of primers   

                revex1 = [e[0] for e in exon_len if revpos-revlen in e[1]][0] # start of primer
                revex2 = [e[0] for e in exon_len if revpos in e[1]][0] # end of primer

                row["for_pos"] = [forex1, forex2]
                row["rev_pos"] = [revex1, revex2]     
                row["for_pos"] = list(set(row["for_pos"]))
                row["rev_pos"] = list(set(row["rev_pos"]))
                  
                first_exon_to_search = forex1 
                last_exon_to_search = revex2
            else: # pdict is c2
                forex = [e[0] for e in exon_len if forpos in e[1]][0] # e1 is range
                revex = [e[0] for e in exon_len if revpos in e[1]][0] # e1 is range
                
                row["for_pos"] = [forex]
                row["rev_pos"] = [revex]
                
                first_exon_to_search = forex   
                last_exon_to_search = revex
            
            # add junction information
            i, first_found, junction = 0, False, ""
            while exon_len[i][0] != last_exon_to_search: 
                if exon_len[i][0] == first_exon_to_search: 
                    first_found = True
                if first_found: 
                    junction += exon_len[i][0] + "_"
                i += 1
            junction += exon_len[i][0] # add last one  
            row["junction"] = junction
            
            # add row to dataframe
            df = pd.concat([df, pd.DataFrame([row])], ignore_index = True)
            
    # convert dataframe type
    df["amplicon_size"] = df["amplicon_size"].astype("int64")
    df["forward_tm"] = df["forward_tm"].astype("float64").round(2)
    df["reverse_tm"] = df["reverse_tm"].astype("float64").round(2)
    df["forward_gc"] = df["forward_gc"].astype("float64").round(2)
    df["reverse_gc"] = df["reverse_gc"].astype("float64").round(2)
    df["amplicon_tm"] = df["amplicon_tm"].astype("float64").round(2)      
    
    # remove possible duplicates
    df = df.drop_duplicates(subset=['forward', 'reverse'])    
    
    return df

###############################################################################

def report_one_exon_design(c2, exon_len, exon_junction, df):    
    """
    This function writes a table with the designed primers, only to use for one
    exon designs. 
    Args: 
        c2 [in] (dict)             Dict of primers
        exon_len [in] (l)          List of tuples, (exon_id, exon length)
        exon_junction [in] (tuple) Tuple with two strings: (1) "Ensemble exon ID
                                   _ Ensembl exon ID" (2) Design specification
        df [in|out] (ps.df)        Pandas df with the designed primers info
    """
    

    for n in range(c2["PRIMER_PAIR_NUM_RETURNED"]): 
        
        patt = "_{}".format(n)
        
        row = {"option": 2, 
               "junction": exon_junction[0], 
               "junction_description": exon_junction[1], 
               "forward": c2["PRIMER_LEFT{}_SEQUENCE".format(patt)].upper(), 
               "reverse": c2["PRIMER_RIGHT{}_SEQUENCE".format(patt)].upper(), 
               "amplicon_size": c2["PRIMER_PAIR{}_PRODUCT_SIZE".format(patt)], 
               "forward_tm": c2["PRIMER_LEFT{}_TM".format(patt)], 
               "reverse_tm": c2["PRIMER_RIGHT{}_TM".format(patt)], 
               "forward_gc": c2["PRIMER_LEFT{}_GC_PERCENT".format(patt)],
               "reverse_gc": c2["PRIMER_RIGHT{}_GC_PERCENT".format(patt)], 
               "amplicon_tm": c2["PRIMER_PAIR{}_PRODUCT_TM".format(patt)], 
               "pair_penalty": c2["PRIMER_PAIR{}_PENALTY".format(patt)]}
        
        # get position of primers
        forpos, forlen = c2["PRIMER_LEFT{}".format(patt)]
        revpos, revlen = c2["PRIMER_RIGHT{}".format(patt)]
        forex = [e[0] for e in exon_len if forpos in e[1]][0] # e1 is range
        revex = [e[0] for e in exon_len if revpos in e[1]][0] # e1 is range
        row["for_pos"] = [(forex, forlen)]
        row["rev_pos"] = [(revex, revlen)]
        
        
        df = pd.concat([df, pd.DataFrame([row])], ignore_index = True)
    
    # convert dataframe type
    df["amplicon_size"] = df["amplicon_size"].astype("int64")
    df["forward_tm"] = df["forward_tm"].astype("float64").round(2)
    df["reverse_tm"] = df["reverse_tm"].astype("float64").round(2)
    df["forward_gc"] = df["forward_gc"].astype("float64").round(2)
    df["reverse_gc"] = df["reverse_gc"].astype("float64").round(2)
    df["amplicon_tm"] = df["amplicon_tm"].astype("float64").round(2)
    
    # remove possible duplicates
    df = df.drop_duplicates(subset=['forward', 'reverse'])
    
    return df