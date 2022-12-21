#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import uuid
import os
import pandas as pd
from subprocess import call, DEVNULL
from ExonSurfer.resources import resources
###############################################################################
#                    blast module FUNCTION DEFINITION SITE                    #
###############################################################################
            
def run_blast_list(fastaf, out=".", blastn_path="blastn", db_path=None, num_threads = 4):
    """
    This function takes a list of primers and runs blastn on them.
    Args:
        lPrimer [in] (list): List of primers
        out [in] (path): Path to the output directory
        num_threads [in] (int): Number of threads to use for blastn
        df [out] (dataframe): Dataframe with the blast results
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
# Make db using makeblastdb
def make_blast_db(download_path):
    """
    This function makes a blast database from the ensembl cdna fasta file.
        input: download_path Path to the cDna fasta file
        output: None
    """
    blastdb_path = "makeblastdb"
    command = [blastdb_path,
               "-dbtype", "nucl",
               "-in", download_path,
               "-out", download_path]
    call(command, stderr = DEVNULL, stdout = DEVNULL)

###############################################################################
def create_index_table(download_path, realease = "108"):
    """
    This function creates a table with the transcript id and the gene symbol.
    """
    from pyensembl import EnsemblRelease
    data = EnsemblRelease(realease)
    #Save CDNA_GENE.txt in the same directory as the cDNA fasta file
    cdna_path = os.path.dirname(download_path)
    cdna_file = os.path.join(cdna_path, "CDNA_GENE.txt")
    with open(cdna_file, "w") as g:
        with open(download_path, "r") as f:
            for line in f.readlines():
                if line.startswith(">"):
                    lElement = line.replace(">","").split(" ")
                    t_id = lElement[0]
                    try:
                        t = data.transcript_by_id(t_id.split(".")[0])
                    except:
                        try:
                            t = data.transcript_by_id(t_id)
                        except:
                            try:
                                t = None
                                symbol = line.split(" ")[6].split(":")[1]
                            except:
                                t = None
                                symbol = ""


                    if t is not None:
                        symbol = t.gene_name
                        g.write(f"{t_id}\t{symbol}\n")
                    else:
                        g.write(f"{t_id}\t{symbol}\n")

###############################################################################
def check_blast_result(df, identity = 70, evalue = 0.8):
    """
    This function takes a dataframe with blast results and filters them based on identity and evalue.
    Args:
        df [in] (dataframe): Dataframe with the blast results
        identity [in] (int): Minimum identity to keep the blast result
        evalue [in] (float): Maximum evalue to keep the blast result
        df [out] (dataframe): Dataframe with the filtered blast results
    """
    df_filter = df.loc[(df.identity >= identity) & (df.evalue <= evalue),"Gene"]
    #print(df.loc[(df.identity >= identity) & (df.evalue <= treshold),:])
    #print(df_filter.unique().tolist())
    
    return len(df_filter.unique().tolist()) == 1

