#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 15:03:53 2023

@author: q0791lp
"""


###############################################################################
#                blast for file input FUNCTION DEFINITION SITE                #
###############################################################################

def check_row_spec(blast_df, row, design_df, max_sep): 
    """ USED WITH cDNA BLAST DATABASE
    This function is applied to every row in the blast DF. It extracts the 
    alignment and checks if there are other alignments (from the same primer 
    pair) that could cause off-target amplification and annotates these as: 
        - detected: every transcript from the same gene amplified by the primers
        - other_transcripts: transcripts different than the target ones
        - other_genes: every transcript from different genes
    Args: 
        blast_df      [in] (pd.df)  Blast alignments df
        row           [in] (df row) Row (from blast df) that will be analyzed 
        design_df     [in] (pd.df)  DF with primer information
        max_sep       [in] (int)    Max sep between als to take as off-target
    """
    # get pair name
    pair_name = row["query id"][:-2]
    
    # if not already annotated
    if row["subject id"] not in design_df.loc[pair_name, "detected"] and \
        row["subject id"] not in design_df.loc[pair_name, "other_genes"]: 
            
        # iterate blast df WITHOUT row being inspected
        ppair_num = (row["query id"][:-2] + "_5", 
                     row["query id"][:-2] + "_3")
        
        poss_offts = blast_df.index[(blast_df["query id"].isin(ppair_num)) &\
                                    (blast_df["subject id"] == row["subject id"])]
    
        for i in [ind for ind in poss_offts if ind != row.name]: 
            # if opposite strand
            if row["strand"] != blast_df.iloc[i]["strand"]: 
                # if close enough
                if abs(row["s. start"] - blast_df.iloc[i]["s. start"]) < max_sep: 
                    # if amplifies same gene
                    design_df.loc[pair_name, "other_genes"] += row["subject id"] +";"
                        
###############################################################################

def check_specificity(blast_df, design_df, max_sep): 
    
    # new columns to fill in 
    design_df["other_genes"] = ""
    design_df["other_transcripts"] = ""
    
    # not for filters, just to inform
    design_df["detected"] = ""

    # check off-targets
    blast_df.apply(lambda row: check_row_spec(blast_df, row, design_df, max_sep), 
                   axis = 1)   
    
    # individual alignments with other genes
    design_df["indiv_als"] = design_df.apply(lambda row: blast_df[(blast_df["query id"] == row.name+"_5") | \
                                                                  (blast_df["query id"] == row.name+"_3")].shape[0], 
                                             axis = 1)    
        
    return design_df