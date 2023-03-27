#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 15:04:21 2023

@author: q0791lp
"""

###############################################################################
#                      genomic blast FUNCTION DEFINITION SITE                 #
###############################################################################

def check_row_genomic(blast_df, row, design_df, max_sep): 
    """ USED WITH GENOMIC BLAST DATABASE
    This function is applied to every row in the blast DF. It extracts the 
    alignment and checks if there are other alignments (from the same primer 
    pair) that could cause off-target amplification and annotates these as 
    genomic_amp (chromosome, start and end position). 
    Args: 
        blast_df      [in] (pd.df)  Blast alignments df
        row           [in] (df row) Row (from blast df) that will be analyzed 
        design_df     [in] (pd.df)  DF with primer information
        max_sep       [in] (int)    Max sep between als to take as off-target
    """    
    # get pair name
    pair_name = row["query id"][:-2]
    
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
                
                # get start and end positions
                spos = min(row["s. start"], row["s. end"], 
                           blast_df.iloc[i]["s. start"], blast_df.iloc[i]["s. end"])
                epos = max(row["s. start"], row["s. end"], 
                           blast_df.iloc[i]["s. start"], blast_df.iloc[i]["s. end"]) 
                
                # build chromosome location
                chrom_loc = "({}:{}-{})".format(row["subject id"], spos, epos)
                
                # if not already annotated
                if chrom_loc not in design_df.loc[pair_name, "genomic_amp"]: 
                    design_df.loc[pair_name, "genomic_amp"] += chrom_loc +";"

###############################################################################

def check_genomic_specificity(blast_df, design_df, max_sep): 
    """ USED WITH GENOMIC BLAST DATABASE
    This function annotates genomic_amp (chromosome, start and end position) in
    the final dataframe. 
    Args: 
        blast_df      [in] (pd.df)  Blast alignments df
        row           [in] (df row) Row (from blast df) that will be analyzed 
        design_df     [in] (pd.df)  DF with primer information
        max_sep       [in] (int)    Max sep between als to take as off-target
    """     
    # new columns to fill in 
    design_df["genomic_amp"] = ""

    # check off-targets
    blast_df.apply(lambda row: check_row_genomic(blast_df, row, design_df, 
                                                 max_sep), 
                   axis = 1)    
        
    return design_df