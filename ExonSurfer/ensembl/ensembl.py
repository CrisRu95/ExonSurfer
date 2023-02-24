#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
from pyensembl import EnsemblRelease

###############################################################################
#                  ensemble module FUNCTION DEFINITION SITE                   #
###############################################################################

def create_ensembl_data(release, species): 
    return EnsemblRelease(release = release, species = species)

###############################################################################
    
def get_gene_by_symbol(gene_symbol, data):
    # release 77 uses human reference genome GRCh38
    """
    This function takes a gene symbol and returns the gene object.
    Args:
        gene_symbol [in] (str)   Gene symbol
        data [in] (Genome obj)   Genome pyensmebl object 
        gene [out] (gene object) Gene object
    """
    gene = data.genes_by_name(gene_symbol)
    
    try: 
        return gene[0]
    except: 
        print("Invalid name")
    
###############################################################################

def get_transcript_from_gene(gene):
    """
    This function takes a gene object and returns a list of transcript objects.
    Args:
        gene [in] (gene object)  Gene object
        transcripts [out] (list) List of transcript objects
    """
    return gene.transcripts

###############################################################################
    
def get_exons_from_transcript(transcript):
    """
    This function takes a transcript object and returns a list of exon objects.
    Args:
        transcripts [in] (transcript object): Transcript object
        exons [out] (list):                   List of exon objects
        """
    return [ex.exon_id for ex in transcript.exons]

###############################################################################
    
def get_coding_transcript(transcripts):
    """
    This function takes a list of transcript objects and keep only the
    proteing coding transcript.
    Args:
        transcripts [in] (list):        List of transcript objects
        coding_transcript [out] (list): List of protein coding transcript objects
    """
    return [x for x in transcripts if x.biotype == "protein_coding"]
    
###############################################################################
    
def get_transcripts_dict(gene, exclude_noncoding = True):
    """
    This function takes a gene object and returns a dictionary of transcript
    objects, with transcript ID as keys, and exon objects as values.
    Args:
        gene [in] (gene object)   Gene object
        exclude_noncoding [in] (bool) False if all transcripts, True to exclude non
                          coding
        dTranscripts [out] (dict) Dictionary of transcript objects, with
                     transcript ID as keys, and exon objects as values
    """
    d = {}
    
    # get list of transcripts to iterate
    all_transcripts = get_transcript_from_gene(gene)
    if exclude_noncoding == False: 
        tcripts = all_transcripts
    else: 
        tcripts = get_coding_transcript(all_transcripts)
    

    for tcript in tcripts:
        d[tcript.id] = "-".join(get_exons_from_transcript(tcript))
     
    return d


  
    