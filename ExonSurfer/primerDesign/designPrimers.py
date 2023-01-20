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

    
def call_primer3(target_seq, junction_i, num_primers): 
    """
    This function calls primers3 for 2 design options: (1) primers are placed 
    ON the exon junction and (2) primers are placed FLANKING the exon junction. 
    Args: 
        target_seq [in] (str)      Full cDNA sequence of the transcript
        junction_i [in] (int)      Index of the exonic junction on the target_seq
        num_primers [in] (int)     Number of primers to design for each option
        case1_primers [out] (dict) Dictionary of primers designed with option 1
        case2_primers [out] (dict) Dictionary of primers designed with option 2
    """
    # annotate number of primers to design
    DESIGN_DICT["PRIMER_NUM_RETURN"] = num_primers
    
    # Case 1: place forward primer or right primer ON junction
    target_dict = {
            'SEQUENCE_ID': 'InternalID',
            'SEQUENCE_TEMPLATE': target_seq,
            'SEQUENCE_OVERLAP_JUNCTION_LIST': junction_i,
            'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION': 4,
            }

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

def annotate_other_transcripts(transcript_list, data): 
    """
    This function annotates the biotype of a given list of transcripts. 
    Args: 
        transcript_list [in] (str) Semicolon separated transcript ID string
        data [in] (Genome obj)     Genome ensembl object
        new_ts [in] (str)          Semicolon separated transcript ID string, with
                                   biotype written in brackets
    """
    
    ts = transcript_list.split(";")
    ts_nov = [t.split(".")[0] for t in ts if t != ""]
    new_ts = ";".join(["{}({})".format(x, data.transcript_by_id(x).biotype) for x in ts_nov])
    
    return new_ts
    
###############################################################################
    
def penalize_final_output(df, transcripts, data, gene_object): 
    """
    This function penalizes the last data design DF (with blast information 
    appended) and returns a list of the best primer pairs for the task
    Args: 
        df [in] (pd.df)          Design df with blast information appended
        transcripts [in] (str)   Target transcript ID (no version) or ALL
        data [in] (Genome obj)   Ensembl Genome object
        gene_obj [in] (Gene obj) Ensembl gene object
        final_df [out] (pd.df)   Filtered df 
    """
    # Annotate other_transcripts and other_genes columns: 
    df["other_transcripts_an"] = df.apply(lambda row: annotate_other_transcripts(row["other_transcripts"], data), axis=1)
    df["other_genes_an"] = df.apply(lambda row: annotate_other_transcripts(row["other_genes"], data), axis=1)
    
    df = df.drop("other_transcripts", axis = 1)
    df = df.drop("other_genes", axis = 1)
    
    df = df.rename({"other_transcripts_an":"other_transcripts", 
                    "other_genes_an": "other_genes"}, axis = "columns")
    
    # annotate number of protein_coding
    df["pcod_transcripts"] = df.apply(lambda row: row["other_transcripts"].count("protein_coding"), axis=1)
    df["pcod_genes"] = df.apply(lambda row: row["other_genes"].count("protein_coding"), axis=1)
    min_pcod_trans = min(df["pcod_transcripts"])
    min_pcod_genes = min(df["pcod_genes"])
    
    if transcripts != "ALL": 
        # prioritize option 1, no other transcripts and no other genes
        if df.loc[(df['other_transcripts'] == "") & (df['other_genes'] == "") & (df["option"] == "1")].shape[0] > 0: 
            final_df = df.loc[(df['other_transcripts'] == "") & (df['other_genes'] == "") & (df["option"] == "1")]
            print("first if")
        # any option, no other transcripts nor genes
        elif df.loc[(df['other_transcripts'] == "") & (df['other_genes'] == "")].shape[0] > 0: 
            final_df = df.loc[(df['other_transcripts'] == "") & (df['other_genes'] == "")]
            print("second if")
        # try with no other genes and only the min protein_coding transcripts
        elif df.loc[(df['other_genes'] == "") & (df['pcod_transcripts'] == min_pcod_trans)].shape[0] > 0: 
            final_df = df.loc[(df['other_genes'] == "") & (df['pcod_transcripts'] == min_pcod_trans)]  
            print("third if")
        # try with no other transcripts at least and the minimum protein coding genes
        elif df.loc[(df['other_transcripts'] == "") & (df['pcod_genes'] == min_pcod_genes)].shape[0] > 0: 
            final_df = df.loc[(df['other_transcripts'] == "") & (df['pcod_genes'] == min_pcod_genes)] 
            print("fourth if")
        # try with the minimum number of protein_coding transcripts and genes
        elif df.loc[(df['pcod_transcripts'] == min_pcod_trans) & (df['pcod_genes'] == min_pcod_genes)].shape[0] > 0:
            final_df = df.loc[(df['pcod_transcripts'] == min_pcod_trans) & (df['pcod_genes'] == min_pcod_genes)]
            print("fifth if")
        else: 
            final_df = df # whatever
            print("whatever")
    
    else: # do not need to prioritize option 1
        # get first transcript id; if 
        vip_trans = gene_object.transcripts[0].transcript_id
        # try without other genes and with first transcript included (not present in junction_description)
        if df.loc[(df['other_genes'] == "") & (df.apply(lambda row: vip_trans not in row["junction_description"], axis=1))].shape[0] > 0: 
            final_df = df.loc[(df['other_genes'] == "") & (df.apply(lambda row: vip_trans not in row["junction_description"], axis=1))]  
        # try with first transcript included and least amount of prot coding genes
        elif df.loc[(df['pcod_genes'] == min_pcod_genes) & (df.apply(lambda row: vip_trans not in row["junction_description"], axis=1))].shape[0] > 0: 
            final_df = df.loc[(df['other_genes'] == "") & (df.apply(lambda row: vip_trans not in row["junction_description"], axis=1))]  
       # try with first transcript incuded
        elif df.loc[(df.apply(lambda row: vip_trans not in row["junction_description"], axis=1))].shape[0] > 0: 
            final_df = df.loc[(df.apply(lambda row: vip_trans not in row["junction_description"], axis=1))]          
        # try with no other genes
        elif df.loc[(df['other_genes'] == "")].shape[0] > 0: 
            final_df = df.loc[(df['other_genes'] == "")]          
        # try with the minimum number of protein_coding transcripts from other genes 
        elif df.loc[(df['pcod_genes'] == min_pcod_genes)].shape[0] > 0:
            final_df = df.loc[(df['pcod_genes'] == min_pcod_genes)]
        else: 
            final_df = df # whatever
            
    return final_df
    
    