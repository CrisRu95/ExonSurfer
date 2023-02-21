#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import os
import pandas as pd

# own modules
from ExonSurfer.ensembl import ensembl
from ExonSurfer.blast import blast
from ExonSurfer.resources import resources
from ExonSurfer.primerDesign import designPrimers, chooseTarget, designConfig
from ExonSurfer.primerDesign import penalizePrimers

# Constants
NPRIMERS = 100

def CreatePrimers(gene, transcripts = "ALL", species = "homo_sapiens_masked",
                  release = 108, design_dict = designConfig.design_dict, 
                  path_out = ".", save_files = True, e_value = 0.8,
                  i_cutoff = 70, max_sep = 700):
    """
    This function is the main function of the pipeline. It takes a gene name and
    a transcript name and returns a list of primers.
    Args:
        gene [in] (str):        Gene name
        transcripts [in] (str): Transcript name or ALL
        species [in] (str):     Species for the data (human, mus_musculus or 
                                                      rattus_novegicus)
        release [in] (str):     Ensembl release number. Change with caution
        design_dict [in] (d):   Dict with primer3 parameters
        path_out [in] (str):    Path to output direct
        BLAST_OUT [out] (df):   Dataframe with blast results
        DESIGN_OUT [out] (df):  Ddataframe with primer design results
    """ 
        
    # Construct transcripts dictionary
    print("Extracting ensemble info")
    data = ensembl.create_ensembl_data(release, 
                                       species.replace("_masked", ""))
    gene_obj = ensembl.get_gene_by_symbol(gene, data)
        
    d = ensembl.get_transcripts_dict(gene_obj, exclude_noncoding = True)
    
    # Get best exonic junction
    print("Getting exon junction")
    junctions_d = chooseTarget.format_junctions(d)
    junction = chooseTarget.choose_target(d, junctions_d, transcripts)
    print("Exon junctions: {}".format(junction))
    
    # Get sequence and junction index
    print("Starting primer design")
    
    # Create dataframe for design
    cols = ("option", "junction", "junction_description", "forward", "reverse", 
            "amplicon_size", "forward_tm", "reverse_tm", "forward_gc", "reverse_gc", 
            "amplicon_tm", "pair_penalty")
    df = pd.DataFrame(columns = cols)
    
    if junction == None: # only one exon
        design_dict["PRIMER_NUM_RETURN"] = NPRIMERS
        target, index, elen = ensembl.construct_one_exon_cdna(resources.MASKED_SEQS(species), 
                                                              data, transcripts)        
        # Design primers
        c2 = designPrimers.call_primer3(target, index, design_dict, enum = 1)
        item = [data.transcript_by_id(transcripts).exons[0].exon_id, 
                "one exon"]
        df = designPrimers.report_one_exon_design(c2, item, df)
        
    else: # Normal design (more than 1 exon)
        # number of primers to design for each junction and option
        num_primers = int(NPRIMERS / (len(junction)*2))
        design_dict["PRIMER_NUM_RETURN"] = num_primers

        for item in junction: 
            target, index, elen = ensembl.construct_target_cdna(resources.MASKED_SEQS(species), 
                                                                gene_obj,
                                                                data, 
                                                                transcripts, 
                                                                item)
            # Design primers
            c1, c2 = designPrimers.call_primer3(target, index, design_dict)
            df = designPrimers.report_design(c1, c2,elen, item, df)
    

    df["pair_num"] = ["Pair{}".format(x) for x in range(0, df.shape[0])]
    df = df.set_index('pair_num')
    
    if df.shape[0] > 1: # we designed primers
        
        # Define three output files, with the same name as the gene and transcript
        DESIGN_OUT = os.path.join(path_out, 
                                  "{}_{}_design.txt".format(gene, transcripts))
        BLAST_OUT = os.path.join(path_out, 
                                 "{}_{}_blast.txt".format(gene, transcripts))
        
        # remove in case remaining from previous design
        for file in (DESIGN_OUT, BLAST_OUT): 
            if os.path.exists(file):
                os.remove(file)
    
        # Write blast input
        FASTA_F = resources.create_temp_fasta()
        resources.fillin_temp_fasta(df, FASTA_F)
        
        # Call blast
        blast_df = blast.run_blast_list(FASTA_F, BLAST_OUT, 
                                        resources.BLAST_DB(species), 
                                        species)
        
        # Delete fasta file
        os.remove(FASTA_F)
        
        # Filter blast results
        blast_df = blast.pre_filter_blast(blast_df, transcripts, gene, df, 
                                          e_value, i_cutoff, False)
        
        # Check blast results positions
        df = blast.check_specificity(blast_df, df, gene, max_sep)
        
        # Filter final DF 
        final_df = penalizePrimers.penalize_final_output(df, transcripts, data, 
                                                         gene_obj)
        final_df = penalizePrimers.make_penalty_score(final_df)
        
        if save_files == True:
            final_df.to_csv(DESIGN_OUT, sep = "\t")
        else: 
            os.remove(BLAST_OUT)
        
        logresult = False # not empty
    
    else: 
        # extract c2 reason
        msg = c2["PRIMER_PAIR_EXPLAIN"].split(",")
        msg_f = [x for x in msg if "ok" not in x and "considered" not in x]
        
        if len(msg_f) > 1: 
            msg_f = [x for x in msg_f if "product size" not in x][0]
        else: 
            msg_f = msg_f[0]
        
        final_df = None
        blast_df = None
        logresult = "Main reason for design failure: {}".format(msg_f)
        
    return blast_df, final_df, logresult
    
    
    
if __name__ == "__main__":
    CreatePrimers()
