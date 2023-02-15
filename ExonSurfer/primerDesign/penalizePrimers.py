#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 12:34:41 2023
"""

# imported modules
from itertools import product


###############################################################################
#                penalizePrimers module FUNCTION DEFINITION SITE              #
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
    if len(newconds) == 8: 
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
    
    # Get number of alignments in transcripts NOT annotated with ensembl
    df["other_genes_rpred_count"] = df.apply(lambda row: row["other_genes_rpred"].count(";"), axis=1)
    df["other_transcripts_rpred_count"] = df.apply(lambda row: row["other_transcripts_rpred"].count(";"), axis=1)
    
    min_pcod_trans = min(df["pcod_trans"])
    min_pcod_genes = min(df["pcod_genes"])
    min_indiv = min(df["indiv_als"])
    min_genes_rpred = min(df["other_genes_rpred_count"])
    min_trans_rpred = min(df["other_transcripts_rpred_count"])    
    
    if transcripts != "ALL": 
        
        # combinations of "if statements"
        columns = ["other_genes", 
                   "pcod_genes",
                   "other_transcripts", 
                   "pcod_trans", 
                   "indiv_als", 
                   "other_genes_rpred_count", 
                   "other_transcripts_rpred_count",
                   "option"]
        
        conditions = [["", None], # other_genes col
                      list(set([0, min_pcod_genes, None])), # pcod_genes col
                      ["", None], # other_trans col
                      list(set([0, min_pcod_trans, None])), # pcod_trans col
                      list(set([0, min_indiv, None])), # indiv_als col
                      list(set([0, min_genes_rpred, None])), # refseq annotated genes
                      list(set([0, min_trans_rpred, None])), # refseq annotated transcripts
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
    
        # STEP 1. Verify all the trans are included or, at least, the most imp
        vip_trans = gene_object.transcripts[0].transcript_id
        
        # without junction description (ALL transcripts included)
        if df.loc[(df["junction_description"] == "")].shape[0] > 0: 
            df = df.loc[(df["junction_description"] == "")]
        
        # first transcript covered by the primers
        elif df.loc[(df.apply(lambda row: vip_trans not in row["junction_description"], 
                              axis=1))].shape[0] > 0: 
            df = df.loc[(df.apply(lambda row: vip_trans not in row["junction_description"], 
                                  axis=1))]
        
        # STEP 2. Get best results according to combinations of "if statements"
        columns = ["other_genes", 
                   "pcod_genes",
                   "indiv_als", 
                   "other_genes_rpred_count", 
                   "other_transcripts_rpred_count",
                   "option"]
        
        conditions = [["", None], # other_genes col
                      list(set([0, min_pcod_genes, None])), # pcod_genes col
                      list(set([0, min_indiv, None])), # indiv_als col
                      list(set([0, min_genes_rpred, None])), # refseq annotated genes
                      list(set([0, min_trans_rpred, None])), # refseq annotated transcripts
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
            print("FINAL CONDITIONS MET ARE:\nProductive als with other genes\t{}".format(all_comb[i-1][0]))
            print("Num of productive als with other pcod genes\t{}".format(all_comb[i-1][1]))
            print("Prod als with other transcripts\t{}".format(all_comb[i-1][2]))
            print("Num of productive als with other pcod trans\t{}".format(all_comb[i-1][3]))
            print("Num of individual als\t{}".format(all_comb[i-1][4]))
            print("Num of other genes refseq preds\t{}".format(all_comb[i-1][5]))
            print("Num of other trans refseq preds\t{}".format(all_comb[i-1][6]))
            print("Desin option\t{}".format(all_comb[i-1][7]))   
        else: 
            print("FINAL CONDITIONS MET ARE:\nProductive als with other genes\t{}".format(all_comb[i-1][0]))
            print("Num of productive als with other pcod genes\t{}".format(all_comb[i-1][1]))
            print("Num of individual als\t{}".format(all_comb[i-1][2]))
            print("Num of other genes refseq preds\t{}".format(all_comb[i-1][3]))
            print("Num of other trans refseq preds\t{}".format(all_comb[i-1][4]))
            print("Desin option\t{}".format(all_comb[i-1][5]))            
    return final_df

###############################################################################

def minmax_norm(v, minv, maxv): 
    """ This function makes the min max normalization"""
    newv = (v - minv) / (maxv - minv) 
    return newv * 100

###############################################################################

def make_penalty_score(df, complete_max = 20): 
    """
    This function transforms the penalty value. 
    """
    # Collapse values bigger than 20 at 20
    df['pair_penalty'].values[df['pair_penalty'].values > complete_max] = complete_max
    
    minv = min(df["pair_penalty"])
    maxv = max(df["pair_penalty"])
    df["pair_score"] = df["pair_penalty"].apply(lambda x: minmax_norm(x, 
                                                                      minv, maxv))
    
    # remove penalty score
    df = df.drop("pair_penalty", axis = 1)
    
    return df