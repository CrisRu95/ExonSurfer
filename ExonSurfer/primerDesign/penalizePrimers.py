#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 12:34:41 2023
"""

# imported modules
import re
from itertools import product

###############################################################################
#                penalizePrimers module FUNCTION DEFINITION SITE              #
###############################################################################    

def annotate_other_transcripts(transcript_list, data): 
    """
    This function annotates the biotype of a given list of transcripts. 
    Args: 
        transcript_list [in] (str) Semicolon or hyphen sep transcript ID string
        data [in] (Genome obj)     Genome ensembl object
        new_ts [in] (str)          Semicolon separated transcript ID string, with
                                   biotype written in brackets
    """
    ts = re.split("-|;", transcript_list) # transform to list
    new_ts = ""
    ts_nov = [t.split(".")[0] for t in ts if t != ""] #  no version info
    
    for x in ts_nov: 
        try: 
            new_ts += ";{}({})".format(x, data.transcript_by_id(x).biotype)
        except: 
            new_ts += ";{}({})".format(x, "not_annotated")
        
    return new_ts
    
###############################################################################

def annotate_amp_lens(df): 
    """
    This function transforms the different amplicon lentghs into True or False
    depending if the primer pair amplifies only transcripts of the same size or 
    not. 
    Args: 
        df [in|out] (pd.df)    Design dataframe. Col "diff_amps_lens" required. 
    """
    def check_same(lst):
        if len(set(lst)) == 1: 
            ret = False
        else: 
            ret = True
        return ret

    # Apply the function to each row and create a new column
    df['amp_noteq_len'] = df['diff_lens_amps'].apply(check_same)
    
    #df = df.drop(columns=['diff_amps_lens'])
    #df.rename(columns={'same_sizes': 'diff_amps_lens'}, 
    #inplace=True) 
    
    return df

###############################################################################

def make_df_comparison(df, condition, columns):
    """
    This function subsets the df according to the conditions indicated. 
    Args: 
        df [in] (pd.df)        Design dataframe
        condition [in] (l)     List of 6 elements with the conditions to be met
                               None values indicate no condition for that column
        columns [in] (l)       List of 6 values with the column names for the conds
        final_df [out] (pd.df) Subsetted dataframe (can be shape (0, X)) 
    """
    # keep only conditions that are NOT set to None
    newcols = [columns[i] for i in range(len(columns)) if condition[i] != None]
    newconds = [condition[i] for i in range(len(condition)) if condition[i] != None]
    
    # check conditions length and buid final_df
    if len(newconds) == 9: 
        final_df = df.loc[(df[newcols[0]] <= newconds[0]) &\
                          (df[newcols[1]] <= newconds[1]) &\
                          (df[newcols[2]] <= newconds[2]) &\
                          (df[newcols[3]] <= newconds[3]) &\
                          (df[newcols[4]] <= newconds[4]) &\
                          (df[newcols[5]] <= newconds[5]) &\
                          (df[newcols[6]] <= newconds[6]) &\
                          (df[newcols[7]] <= newconds[7]) &\
                          (df[newcols[8]] <= newconds[8])]
    elif len(newconds) == 8: 
        final_df = df.loc[(df[newcols[0]] <= newconds[0]) &\
                          (df[newcols[1]] <= newconds[1]) &\
                          (df[newcols[2]] <= newconds[2]) &\
                          (df[newcols[3]] <= newconds[3]) &\
                          (df[newcols[4]] <= newconds[4]) &\
                          (df[newcols[5]] <= newconds[5]) &\
                          (df[newcols[6]] <= newconds[6]) &\
                          (df[newcols[7]] <= newconds[7])]
    elif len(newconds) == 7: 
        final_df = df.loc[(df[newcols[0]] <= newconds[0]) &\
                          (df[newcols[1]] <= newconds[1]) &\
                          (df[newcols[2]] <= newconds[2]) &\
                          (df[newcols[3]] <= newconds[3]) &\
                          (df[newcols[4]] <= newconds[4]) &\
                          (df[newcols[5]] <= newconds[5]) &\
                          (df[newcols[6]] <= newconds[6])]
            
    elif len(newconds) == 6: 
        final_df = df.loc[(df[newcols[0]] <= newconds[0]) &\
                          (df[newcols[1]] <= newconds[1]) &\
                          (df[newcols[2]] <= newconds[2]) &\
                          (df[newcols[3]] <= newconds[3]) &\
                          (df[newcols[4]] <= newconds[4]) &\
                          (df[newcols[5]] <= newconds[5])]
            
    elif len(newconds) == 5: 
        final_df = df.loc[(df[newcols[0]] <= newconds[0]) &\
                          (df[newcols[1]] <= newconds[1]) &\
                          (df[newcols[2]] <= newconds[2]) &\
                          (df[newcols[3]] <= newconds[3]) &\
                          (df[newcols[4]] <= newconds[4])]
            
    elif len(newconds) == 4: 
        final_df = df.loc[(df[newcols[0]] <= newconds[0]) &\
                          (df[newcols[1]] <= newconds[1]) &\
                          (df[newcols[2]] <= newconds[2]) &\
                          (df[newcols[3]] <= newconds[3])]
            
    elif len(newconds) == 3: 
        final_df = df.loc[(df[newcols[0]] <= newconds[0]) &\
                          (df[newcols[1]] <= newconds[1]) &\
                          (df[newcols[2]] <= newconds[2])]
            
    elif len(newconds) == 2: 
        final_df = df.loc[(df[newcols[0]] <= newconds[0]) &\
                          (df[newcols[1]] <= newconds[1])]
            
    elif len(newconds) == 1: 
        final_df = df.loc[(df[newcols[0]] <= newconds[0])]
            
    else: 
        final_df = df
                
    return final_df

###############################################################################

def penalize_final_output(df, transcripts, data): 
    """
    This function penalizes the last data design DF (with blast information 
    appended) and returns a list of the best primer pairs for the task
    Args: 
        df [in] (pd.df)          Design df with blast information appended
        transcripts [in] (l|str) List of target transcripts ID (no version) or ALL
        data [in] (Genome obj)   Ensembl Genome object
        final_df [out] (pd.df)   Filtered df 
    """
    # Annotate other_transcripts and other_genes columns: 
    df["other_transcripts_an"] = df.apply(lambda row: annotate_other_transcripts(row["other_transcripts"], 
                                                                                 data), 
                                          axis=1)
    df["other_genes_an"] = df.apply(lambda row: annotate_other_transcripts(row["other_genes"], 
                                                                           data), 
                                    axis=1)
    # remove unnanotated columns
    df = df.drop("other_transcripts", axis = 1)
    df = df.drop("other_genes", axis = 1)
    
    # rename columns 
    df = df.rename({"other_transcripts_an":"other_transcripts", 
                    "other_genes_an": "other_genes"}, axis = "columns")
    
    # annotate number of protein_coding
    df["pcod_trans"] = df.apply(lambda row: row["other_transcripts"].count("protein_coding"),
                                axis=1)
    df["pcod_genes"] = df.apply(lambda row: row["other_genes"].count("protein_coding"), 
                                axis=1)
    # Annotate amplicon lengths
    df = annotate_amp_lens(df)
    
    # Set minimum variables
    min_pcod_trans = min(df["pcod_trans"])
    min_pcod_genes = min(df["pcod_genes"])
    min_indiv = min(df["indiv_als"])
    
    if transcripts != "ALL": 
        
        # combinations of "if statements"
        columns = ["amp_noteq_len", 
                   "other_genes", 
                   "pcod_genes",
                   "other_transcripts", 
                   "pcod_trans", 
                   "indiv_als", 
                   "option"]
        
        conditions = [[False, None], # amp_noteq_len col
                      ["", None], # other_genes col
                      list(set([0, min_pcod_genes, None])), # pcod_genes col
                      ["", None], # other_trans col
                      list(set([0, min_pcod_trans, None])), # pcod_trans col
                      list(set([0, min_indiv, None])), # indiv_als col
                      [1, None]] # option col
        
        all_comb = list(product(*conditions))
        
        # remove ilogical combinations
        all_comb = [comb for comb in all_comb if not (comb[0] == "" and comb[1] == None)]
        all_comb = [comb for comb in all_comb if not (comb[0] == "" and comb[1] > 0)]
        all_comb = [comb for comb in all_comb if not (comb[2] == "" and comb[3] == None)]
        all_comb = [comb for comb in all_comb if not (comb[2] == "" and comb[3] > 0)]
        
        # initialize loop 
        i, shape = 0, 0

        while i < len(all_comb) and shape == 0:
            final_df = make_df_comparison(df, all_comb[i], columns)
            shape = final_df.shape[0]
            i = i + 1 # keep iterating...
    
    else: # transcripts == ALL
        # STEP 2. Get best results according to combinations of "if statements"
        columns = ["amp_noteq_len", 
                   "other_genes", 
                   "pcod_genes",
                   "indiv_als", 
                   "option"]
        
        conditions = [[False, None], # amp_noteq_len col 
                      ["", None], # other_genes col
                      list(set([0, min_pcod_genes, None])), # pcod_genes col
                      list(set([0, min_indiv, None])), # indiv_als col
                      [1, None]] # option col
        
        all_comb = list(product(*conditions))
        
        # remove ilogical combinations
        all_comb = [comb for comb in all_comb if not (comb[0] == "" and comb[1] == None)]
        all_comb = [comb for comb in all_comb if not (comb[0] == "" and comb[1] > 0)]
        
        # initialize loop 
        i, shape = 0, 0

        while i < len(all_comb) and shape == 0:
            final_df = make_df_comparison(df, all_comb[i], columns)
            shape = final_df.shape[0]
            i = i + 1 # keep iterating...    
            
    # Report final conditions for the filter(i-1)
    if transcripts != "ALL": 
        print("FINAL CONDITIONS MET ARE:\nProductive als with other genes\t{}".format(all_comb[i-1][1]))
        print("Amplicons of different lengths: {}".format(all_comb[i-1][0]))
        print("Num of productive als with other pcod genes\t{}".format(all_comb[i-1][2]))
        print("Prod als with other transcripts\t{}".format(all_comb[i-1][3]))
        print("Num of productive als with other pcod trans\t{}".format(all_comb[i-1][4]))
        print("Num of individual als\t{}".format(all_comb[i-1][5]))
        print("Desin option\t{}".format(all_comb[i-1][6]))   
    else: 
        print("FINAL CONDITIONS MET ARE:\nProductive als with other genes\t{}".format(all_comb[i-1][1]))
        print("Amplicons of different lengths: {}".format(all_comb[i-1][0]))
        print("Num of productive als with other pcod genes\t{}".format(all_comb[i-1][2]))
        print("Num of individual als\t{}".format(all_comb[i-1][3]))
        print("Desin option\t{}".format(all_comb[i-1][4]))         
        
    return final_df

###############################################################################

def minmax_norm(v, minv=0, maxv=20): 
    """ This function makes the min max normalization"""
    if maxv != minv: 
        newv = abs((v - maxv) / (maxv - minv))
    else: 
        newv = 1
    return round(newv * 100, 2)

###############################################################################

def make_penalty_score(df, complete_max = 20): 
    """
    This function transforms the penalty value. 
    """
    # Collapse values bigger than 20 at 20
    df['pair_penalty'].values[df['pair_penalty'].values > complete_max] = complete_max
    
    df["pair_score"] = df["pair_penalty"].apply(lambda x: minmax_norm(x))
    
    # remove penalty score
    df = df.drop("pair_penalty", axis = 1)
    
    return df

###############################################################################

def genomic_filter(df): 
    """
    This function filters the primers with genomic amplification (if possible)
    Args: 
        df [in] (pd.df) Final_df with genomic_amp column
    """
    df["genomic_filt"] = df.apply(lambda r: \
                                 len([x for x in r["genomic_amp"].split(";") if x != ""]), 
                                 axis= 1)
        
    minval = min(df["genomic_filt"])
    
    df = df.loc[df["genomic_filt"] == minval]
    
    return df
    
        