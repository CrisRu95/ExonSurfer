#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import pandas as pd
from subprocess import call, DEVNULL

# own modules
from ExonSurfer.resources import resources

# Constants
E_VALUE_CUTOFF = 0.8
IDENTITIY_CUTOFF = 70
MAX_SEP = 700

###############################################################################
#                    blast module FUNCTION DEFINITION SITE                    #
###############################################################################
            
def run_blast_list(fastaf, out, db_path,
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
                 "q. start", "q. end", "s. start",
                "s. end", "evalue", "strand")
    df = pd.read_csv(out, sep = "\t", names = df_header)
    
    table = pd.read_csv(resources.get_cdna_file(), sep="\t", 
                        names=("transcript_id", "gene_symbol"), header=None)
    df = pd.merge(df,table, left_on="subject id", right_on="transcript_id")
    df.to_csv(out, sep = "\t", index = False)
    
    return df


###############################################################################
    
def filter_intended_al(primer_id, identity, transcript_id, gene_id, 
                       t_transcript, t_gene, done_ids) : 
    """
    This function assesses whether a primer alignment is the intended one (i.e, 
    in the target transcript with a high identity) or is a possible off-target.
    Args: 
        primer_id [in] (str)     Primer identifier. Format: Pair[\d+]_(3|5)
        identity [in] (int)      Alignment identity, as returned by blast
        transcript_id [in] (str) Transcript identifier (subject id, from the 
                                 blast df), with version
        gene_id [in] (str)       Gene identifier (form the blast df)
        t_transcript [in] (str)  Target transcript, without version (can be ALL)
        t_gene [in] (str)        Target gene
        done_ids [in] (list)     Already assessed primer_ids, only one alignment
                                 per primer can be removed
        to_keep [out] (bool)     False if the alignment is expected 
    """
    if t_transcript == "ALL": 
        to_keep = True if gene_id != t_gene else False
    else: 
        to_keep = True # default situation
        if primer_id not in done_ids: 
            if transcript_id.split(".")[0] == t_transcript and identity == 100: 
                to_keep = False
                # if there is another instance, it will not be removed
                done_ids.append(primer_id)             

    return to_keep

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
                     e_cutoff = E_VALUE_CUTOFF, i_cutoff = IDENTITIY_CUTOFF): 
    """
    This function filters the alignments returned by blast in order to keep the
    ones that are unintended (i.e., outside the target transcript), with a good
    3' end alignment, and with good enough e and identity values. 
    Args: 
        blast_df [in] (pd.df) Alignment dataframe
        t_transcript [in] (str) Target transcript, without version info
        design_df [in] (pd.df) Design dataframe
        e_cutoff [in] (float) Maximum e value to consider an alignment
        i_cutoff [in] (float/int) Minimum identity to consider an alignment
    """ 
    # keep only alignments outside the intended transcript
    done_ids = []
    blast_df["filter1"] = blast_df.apply(lambda row: filter_intended_al(row["query id"], 
                                                                        row["identity"],
                                                                        row["subject id"],
                                                                        row["gene_symbol"],
                                                                        t_transcript, 
                                                                        t_gene, 
                                                                        done_ids), axis = 1)
    # keep only alignments where primers 3' end is aligned
    blast_df["filter2"] = blast_df.apply(lambda row: filter_3end_al(row["query id"], 
                                                                    row["q. start"], 
                                                                    row["q. end"], 
                                                                    row["strand"],
                                                                    design_df), axis = 1)
    # apply filters
    blast_df = blast_df[blast_df["filter1"]]
    blast_df = blast_df[blast_df["filter2"]]
    blast_df = blast_df[blast_df["evalue"] <= e_cutoff]
    blast_df = blast_df[blast_df["identity"] >= i_cutoff]

    return blast_df


###############################################################################
    
def check_specificity(blast_df, design_df, t_gene, max_sep = MAX_SEP):
    """
    This function check if the specificity of the blast result and annotates: 
    (a) productive blast alignments and (b) unproductive blast alignments as 
    "other_genes" (productive), "other_transcripts" (productive) and "indiv_als"
    (unproductive), in the design_df. 
    Args: 
        blast_df [in] (pd.df)  Filtered alignments dataframe
        design_df [in] (pd.df) Design dataframe
        t_gene [in] (str)      Target gene name
        max_sep [in] (int)     Maximum separation between 2 alignments in order
                               to form an off-target
    """
    # new columns to fill in 
    design_df["other_transcripts"] = ""
    design_df["other_genes"] = ""
    
    # for every forward primer
    for for_id in [x for x in list(set(blast_df["query id"])) if "_5" in x]: 
        # for every alignment for a given forward
        for subj in list(blast_df[blast_df["query id"] == for_id]["subject id"]): 
            forpos = int(blast_df.loc[(blast_df['query id'] == for_id) & \
                                      (blast_df['subject id'] == subj)]["s. start"])
            
            # check if there are any reverse alignments on the same subject id
            rev_id = for_id[:-2] + "_3"
            any_rev = blast_df.loc[(blast_df['query id'] == rev_id) & \
                                   (blast_df['subject id'] == subj)]["s. start"]
            if len(any_rev) > 0: 
                revpos = int(any_rev)
                
                # both alignments are on the same, untargeted, subject
                if abs(revpos-forpos) <= max_sep: 
                    
                    gene = blast_df.loc[(blast_df['query id'] == rev_id) & \
                                        (blast_df['subject id'] == subj)]["gene_symbol"]
                    
                    if t_gene == gene.item(): # same gene, different transcript
                        design_df.loc[for_id[:-2], "other_transcripts"] += subj +";" 
                    else:  # different gene
                        design_df.loc[for_id[:-2], "other_genes"] += subj +";"
    
    
    # individual alignments with other genes
    b_ogenes = blast_df[blast_df["gene_symbol"] != t_gene]
    design_df["indiv_als"] = design_df.apply(lambda row: b_ogenes[(b_ogenes["query id"] == row.name+"_5") | \
                                                                  (b_ogenes["query id"] == row.name+"_3")].shape[0], 
                                             axis = 1)    
    
    return design_df