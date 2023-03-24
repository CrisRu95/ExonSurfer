#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import os
import pandas as pd
from subprocess import call, DEVNULL

# own modules
from ExonSurfer.resources import resources

###############################################################################
#                    blast module FUNCTION DEFINITION SITE                    #
###############################################################################
            
def run_blast_list(fastaf, out, db_path, species, tomerge = True, 
                   blastn_path = "blastn", 
                   num_threads = 4):
    """
    This function takes a list of primers and runs blastn on them.
    Args:
        out [in] (path)        Path to the output directory
        num_threads [in] (int) Number of threads to use for blastn
        df [out] (dataframe)   Dataframe with the blast results
    """
    command_line = [blastn_path,
                    '-task', 'blastn-short',
                    '-query', fastaf,
                    '-out', out,
                    '-db', db_path,
                    '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send evalue sstrand',
                    '-num_threads', str(num_threads), 
                    '-strand', 'both']

    # RunBlastDBCommand(command_line)
    call(command_line, stderr = DEVNULL, stdout = DEVNULL)
    
    # Store DF results 
    df_header = ("query id", "subject id","identity", "alignment length", 
                 "q. start", "q. end", "s. start", "s. end", "evalue", "strand")
    
    df = pd.read_csv(out, sep = "\t", names = df_header)
    
    # remove transcript version information
    df["subject id"] = df["subject id"].map(lambda x: x.split(".")[0])
    
    if tomerge == True: 
        table = pd.read_csv(os.path.join(resources.get_blastdb_path(species),
                                         resources.IDS_TABEL), 
                            sep = "\t", names=["id", "gene"], header=None)
        # remove transcript version information
        table["id"] = table["id"].map(lambda x: x.split(".")[0])
        
        df = pd.merge(df, table, left_on="subject id", right_on = "id", how = "left")
        
        # overwrite annotated blast result
        df.to_csv(out, sep = "\t", index = False)
        os.remove(out)
    
    return df

###############################################################################

def filter_3end_al(primer_id, q_start, q_end, strand, design_df): 
    """
    This function assesses whether a blast alignment involves the 3' end of the
    primers. 
    Args: 
        primer_id [in] (str)   Primer identifier. Format: Pair[\d+]_(3|5)
        q_start [in] (int)     Query start
        q_end [in] (int)       Query end
        design_df [in] (pd.df) Design dataframe with the primer sequences
        to_keep [out] (bool)   True if the 3' end of the primer is aligned
    """
    # get complete primer len
    if strand  == "plus": 
        if "_5" in primer_id: 
            plen = len(design_df.loc[primer_id[:-2],"forward"])
        else: 
            plen = len(design_df.loc[primer_id[:-2],"reverse"])
        to_keep = True if q_end == plen else False
    else: 
        to_keep = True if q_start == 1 else False

    return to_keep 

###############################################################################
    
def pre_filter_blast(blast_df, t_transcript, t_gene, design_df, 
                     e_cutoff, i_cutoff, filter2 = True): 
    """
    This function filters the alignments returned by blast in order to keep the
    ones that are unintended (i.e., outside the target transcript), with a good
    3' end alignment, and with good enough e and identity values. 
    Args: 
        blast_df [in] (pd.df)     Alignment dataframe
        t_transcript [in] (l|str) List of target transcripts, w/ version info or ALL
        design_df [in] (pd.df)    Design dataframe
        e_cutoff [in] (float)     Maximum e value to consider an alignment
        i_cutoff [in] (float/int) Minimum identity to consider an alignment
        filter2 [in] (bool)       True if to disregard alignments with mismatches 
                                  at the 3' end'
    """ 
    # keep only alignments where primers 3' end is aligned
    if filter2 == True: 
        blast_df["filter2"] = blast_df.apply(lambda row: filter_3end_al(row["query id"], 
                                                                        row["q. start"], 
                                                                        row["q. end"], 
                                                                        row["strand"],
                                                                        design_df), axis = 1)
        blast_df = blast_df[blast_df["filter2"]]
    
    # apply filters
    blast_df = blast_df[blast_df["evalue"] <= e_cutoff]
    blast_df = blast_df[blast_df["identity"] >= i_cutoff]
    
    # reset index just in case
    blast_df = blast_df.reset_index()
    
    return blast_df

###############################################################################

def filter_big_blast(blast_df, design_df, maxrows = 15000): 
    
    # import median function
    from numpy import median
    
    to_continue = True # initialize
    while blast_df.shape[0] > maxrows and to_continue == True: 
        raw_counts = blast_df["query id"].value_counts()
        clean_counts = []
        for ppair in list(design_df.index): 
            c = raw_counts["{}_3".format(ppair)] + raw_counts["{}_5".format(ppair)]
            clean_counts.append((ppair, c))
    
        # now filter and keep only primer pairs that less appear in the blast 
        cutval = median([x[1] for x in clean_counts])
        to_remove_pairs = [x[0] for x in clean_counts if x[1] >= cutval]
        to_remove_primers = [x+"_5"for x in to_remove_pairs] + [x+"_3"for x in to_remove_pairs]
        
        if len(to_remove_pairs) < design_df.shape[0]: 
            blast_df = blast_df[~blast_df['query id'].isin(to_remove_primers)]
            design_df = design_df.drop(to_remove_pairs)
        else: 
            to_continue = False
            
    # reset index just in case
    blast_df = blast_df.reset_index()
            
    return blast_df, design_df

###############################################################################

def check_row_spec(blast_df, row, design_df, max_sep, t_gene, t_transcripts): 
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
        t_gene        [in] (str)    Target gene name (gene symbol)
        t_transcripts [in] (l|str)  List fo transcripts or "ALL"
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
                    if t_gene == row["gene"]: 
                        design_df.loc[pair_name, "detected"] += row["subject id"] +";"
                        if t_transcripts != "ALL" and row["subject id"] not in t_transcripts: 
                            design_df.loc[pair_name, "other_transcripts"] += row["subject id"] +";"
                    else: 
                        design_df.loc[pair_name, "other_genes"] += row["subject id"] +";"
                        
###############################################################################

def check_specificity(blast_df, design_df, t_gene, t_transcript, max_sep): 
    
    # new columns to fill in 
    design_df["other_genes"] = ""
    design_df["other_transcripts"] = ""
    
    # not for filters, just to inform
    design_df["detected"] = ""

    # check off-targets
    blast_df.apply(lambda row: check_row_spec(blast_df, row, design_df, max_sep, 
                                              t_gene, t_transcript), 
                   axis = 1)   
    
    # individual alignments with other genes
    b_ogenes = blast_df[blast_df["gene"] != t_gene]
    design_df["indiv_als"] = design_df.apply(lambda row: b_ogenes[(b_ogenes["query id"] == row.name+"_5") | \
                                                                  (b_ogenes["query id"] == row.name+"_3")].shape[0], 
                                             axis = 1)    
        
    return design_df

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
    