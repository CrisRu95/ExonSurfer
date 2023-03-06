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

def get_transcript_from_gene(gene_obj, only_id = False):
    """
    This function takes a gene object and returns a list of transcript objects.
    Args:
        gene_obj [in] (gene object)  Gene object
        only_id [in] (bool)      True if to return only tcript ids, not objects
        to_return [out] (l)      List of transcript object of transcript ids
    """
    if only_id == False: 
        to_return = gene_obj.transcripts
    else: 
        to_return = [t.id for t in gene_obj.transcripts]
        
    return to_return

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
    
def get_cdna_seq(data, transcript_id, masked_chr): 
    """
    This function returns the CDNA of a specific transcript
    Args: 
        data [in] (Genome obj)    EnsemblRelease object
        transcript_it [in] (str)  Ensembl transcript id without versoin info
        masked_chr  [in] (str)    Full path to the masked chromosome files
    """
    # transcript object
    t_obj = data.transcript_by_id(transcript_id)
    
    # read chromosome 
    chrom_open = open(masked_chr.format(t_obj.contig), "r")
    tt = chrom_open.read() # full chromosome sequence
    chrom_open.close()
    
    # initialize all
    cdna = ""
    tosum = 0
    
    # check strand
    if t_obj.strand == "+": 
        list_of_exons = t_obj.exons
    else: 
        list_of_exons = reversed(t_obj.exons)
    
    # iterate ALL exons in the transcript and build cdna
    for exon in list_of_exons: # exon is pyensembl object
        cdna += tt[exon.start-1:exon.end]
        tosum += exon.end - exon.start + 1
        
    return cdna         

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
    
def get_transcripts_dict(gene, exclude_noncoding):
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


  
    