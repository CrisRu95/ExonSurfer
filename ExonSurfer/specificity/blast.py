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
            
def run_blast_list(fastaf, out, db_path, species,
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
    
    
    table = pd.read_csv(os.path.join(resources.get_blastdb_path(species),
                                     resources.IDS_TABEL), sep = "\t")
    
    df = pd.merge(df, table, left_on="subject id", right_on = "id")
    
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

    return blast_df


###############################################################################
    
def check_specificity(blast_df, design_df, t_gene, t_transcript, max_sep):
    """
    This function check if the specificity of the blast result and annotates: 
    (a) productive blast alignments and (b) unproductive blast alignments as 
    "other_genes" (productive), "other_transcripts" (productive) and "indiv_als"
    (unproductive), in the design_df. 
    Args: 
        blast_df [in] (pd.df)  Filtered alignments dataframe
        design_df [in] (pd.df) Design dataframe
        t_gene [in] (str)      Target gene name
        t_transcript [in] (l|str) List of transcript ids or "ALL"
        max_sep [in] (int)     Maximum separation between 2 alignments in order
                               to form an off-target
    """
    # new columns to fill in 
    design_df["other_transcripts"] = ""
    design_df["other_genes"] = ""
    design_df["other_transcripts_rpred"] = ""
    design_df["other_genes_rpred"] = ""
    
    # not for filters, just to inform

    
    # LOOP 1. only for ensembl nomenclature (not predicte)
    # for every forward primer
    for for_id in [x for x in list(set(blast_df["query id"])) if "_5" in x]: 
        # for every alignment for a given forward
        
        # Subject is annotated in ensembl nomenclature
        for subj in [x for x in list(blast_df[blast_df["query id"] == for_id]["subject id"]) if x != "-"]: 
            
            for forpos in blast_df.loc[(blast_df['query id'] == for_id) & \
                                      (blast_df['subject id'] == subj)]["s. start"]: 
                forpos = int(forpos)

                # check if there are any reverse alignments on the same subject id
                rev_id = for_id[:-2] + "_3"
            
                for revpos in blast_df.loc[(blast_df['query id'] == rev_id) & \
                                       (blast_df['subject id'] == subj)]["s. start"]: 
                    revpos = int(revpos)
                    
                    # both alignments are on the same, untargeted, subject
                    if abs(revpos-forpos) <= max_sep: 
                        
                        gene = blast_df[(blast_df['query id'] == rev_id) & \
                                        (blast_df['subject id'] == subj)].gene.item()
                        
                        ensembl_id = blast_df[(blast_df['query id'] == rev_id) & \
                                              (blast_df['subject id'] == subj)].ensembl_id.item()
                        try: 
                            if ensembl_id != "-": 
                                if t_gene == gene: # same gene, different transcript
                                    # annotate for filters
                                    if ensembl_id not in t_transcript: 
                                        design_df.loc[for_id[:-2], "other_transcripts"] += ensembl_id +";" 
                                else:  # different gene
                                    design_df.loc[for_id[:-2], "other_genes"] += ensembl_id +";"
                            else: 
                                if t_gene == gene: # same gene, different transcript
                                    design_df.loc[for_id[:-2], "other_transcripts_rpred"] += subj +";" 
                                else:  # different gene
                                    design_df.loc[for_id[:-2], "other_genes_rpred"] += subj +";"
                        except: 
                              print("ensembl id: {}".format(ensembl_id))
                              print("rev_id id: {}".format(rev_id))
                              print("subj id: {}".format(subj))
                              
    # individual alignments with other genes
    b_ogenes = blast_df[blast_df["gene"] != t_gene]
    design_df["indiv_als"] = design_df.apply(lambda row: b_ogenes[(b_ogenes["query id"] == row.name+"_5") | \
                                                                  (b_ogenes["query id"] == row.name+"_3")].shape[0], 
                                             axis = 1)    
    
    return design_df

###############################################################################

def show_ot_for_pair(transcripts, other_genes, other_genes_rpred, 
                     other_transcripts, other_transcripts_rpred): 
    """
    This function says if a primer pair has off-targets to return or not. 
    Args: 
        transcripts [in] (str|l)      List of transcripts or ALL
        other_genes [in] (str)        List of other genes or empty
        other_genes_rpred [in] (str)  List other genes refseq annotated
        other_genes [in] (str)        List of other transcripts or empty
        other_genes_rpred [in] (str)  List other transcripts refseq annotated
        to_return [out] (int)         1 if there are off_targets, 0 if not
    """
    to_return = 0 #  default not show 
    
    if transcripts == "ALL": 
        if other_genes != "" or other_genes_rpred != "": 
            to_return = 1
    else: 
        if other_genes != "" or other_genes_rpred != "" or other_transcripts != "" or other_transcripts_rpred != "":
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
                                                                          row["other_genes_rpred"],
                                                                          row["other_transcripts"], 
                                                                          row["other_transcripts_rpred"]), 
                                             axis = 1)
    
    return final_df
    