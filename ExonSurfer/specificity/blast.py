#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import os
import pandas as pd
from subprocess import call , DEVNULL
from copy import deepcopy

# own modules
from ExonSurfer.resources import resources

###############################################################################
#                    blast module FUNCTION DEFINITION SITE                    #
###############################################################################
            
def run_blast_list(fastaf, out, db_path, species, i_cutoff, e_value, 
                   tomerge = True, 
                   genomic = False, 
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
                    '-strand', 'both', 
                    '-perc_identity', str(i_cutoff), 
                    '-evalue', str(e_value)]

    # RunBlastDBCommand(command_line)
    call(command_line, stderr = DEVNULL, stdout = DEVNULL)

    # Store DF results 
    df_header = ("query id", "subject id","identity", "alignment length", 
                 "q. start", "q. end", "s. start", "s. end", "evalue", "strand")
    
    df = pd.read_csv(out, sep = "\t", names = df_header)
    
    if genomic == False: 
        # remove transcript version information
        df["subject id"] = df["subject id"].map(lambda x: x.split(".")[0])
    
    if tomerge == True: 
        table = pd.read_csv(os.path.join(resources.get_blastdb_path(species),
                                         resources.IDS_TABEL), 
                            sep = "\t", names=["id", "gene"], header=None)
        if species=="litomosoides_sigmodontis":
            table["id"] = table["id"].map(lambda x: x.split("-")[0])
        else:
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
    
def pre_filter_blast(blast_df, design_df, e_cutoff, i_cutoff, filter2 = True): 
    """
    This function filters the alignments returned by blast in order to keep the
    ones that are unintended (i.e., outside the target transcript), with a good
    3' end alignment, and with good enough e and identity values. 
    Args: 
        blast_df [in] (pd.df)     Alignment dataframe
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

def filter_big_blast(t_gene, blast_df, design_df, maxrows = 15000): 
    """
    This function controls the maximum blast DF file size. It removes 
    all the alignments of the pairs that have more alignemts to genes different
    than the target, and also removes these pairs from the design df.
    Args: 
        t_gene [in] (str)             Target gene name (as in table.txt)
        blast_df [in|out] (pd.df)     Alignment dataframe
        design_df [in|out] (pd.df)    Design dataframe
    """
    # import median function
    from numpy import median
    
    # ignore all alignments that correspond to target gene
    newblast = deepcopy(blast_df)
    newblast = newblast[newblast["gene"] != t_gene]
    to_continue = True # initialize
    while newblast.shape[0] > maxrows and to_continue == True: 
        raw_counts = newblast["query id"].value_counts()
        clean_counts = []
        for ppair in list(design_df.index): 
            try: 
                c = raw_counts["{}_3".format(ppair)] + raw_counts["{}_5".format(ppair)]
                clean_counts.append((ppair, c))
            except KeyError: 
                pass
            
        if len(clean_counts) == 0: 
            to_continue = False
        else: 
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

def filter_big_gblast(blast_df, design_df, maxrows = 15000): 
    """
    This function controls the maximum blast DF file size. It removes 
    all the alignments of the pairs that have more alignemts, and also removes 
    these pairs from the design df.
    Args: 
        blast_df [in|out] (pd.df)     Alignment dataframe
        design_df [in|out] (pd.df)    Design dataframe
    """
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
