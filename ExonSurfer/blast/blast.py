#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import pandas as pd
from subprocess import call, DEVNULL

# own modules
from ExonSurfer.resources import resources

###############################################################################
#                    blast module FUNCTION DEFINITION SITE                    #
###############################################################################
            
def run_blast_list(fastaf, 
                   out = ".", 
                   blastn_path = "blastn", 
                   db_path = None, 
                   num_threads = 4):
    """
    This function takes a list of primers and runs blastn on them.
    Args:
        lPrimer [in] (list)    List of primers
        out [in] (path)        Path to the output directory
        num_threads [in] (int) Number of threads to use for blastn
        df [out] (dataframe)   Dataframe with the blast results
    """
    #identifier = str(uuid.uuid4())
    command_line = [blastn_path,
                    '-task', 'blastn-short',
                    '-query', fastaf,
                    '-out', out,
                    '-db', db_path,
                    '-outfmt', '6',
                   '-num_threads', str(num_threads)]

    # RunBlastDBCommand(command_line)
    call(command_line, stderr = DEVNULL, stdout = DEVNULL)
    
    # Store DF results 
    df_header = ("query id", "subject id","identity", "alignment length", 
                 "mismatches", "gap opens", "q. start", "q. end", "s. start",
                "s. end", "evalue", "bit score")
    df = pd.read_csv(out, sep = "\t", names = df_header)
    table = pd.read_csv(resources.get_cdna_file(), sep="\t", names=("transcript_id", "gene_symbol"), header=None)
    df = pd.merge(df,table, left_on="subject id", right_on="transcript_id")
    df.to_csv(out, sep = "\t", index = False)
    return df


###############################################################################
                        
def check_blast_result(df, identity = 70, evalue = 0.8):
    """
    This function takes a dataframe with blast results and filters them based on identity and evalue.
    Args:
        df [in] (pd.df)      Dataframe with the blast results
        identity [in] (int)  Minimum identity to keep the blast result
        evalue [in] (float)  Maximum evalue to keep the blast result
        df [out] (pd.df)     Dataframe with the filtered blast results
    """
    df_filter = df.loc[(df.identity >= identity) & (df.evalue <= evalue),"Gene"]
    #print(df.loc[(df.identity >= identity) & (df.evalue <= treshold),:])
    #print(df_filter.unique().tolist())
    
    return len(df_filter.unique().tolist()) == 1

