#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
from pyensembl import EnsemblRelease

MASKED_SEQS = "/home/ecrisru/Desktop/refgenome/Chrs_masked_001/Chr_{}_masked001.fa"

###############################################################################
#                  ensemble module FUNCTION DEFINITION SITE                   #
###############################################################################

    
def get_gene_by_symbol(gene_symbol, release = 108):
    # release 77 uses human reference genome GRCh38
    """
    This function takes a gene symbol and returns the gene object.
    Args:
        gene_symbol [in] (str):   Gene symbol
        release [in] (int):       Ensembl release
        gene [out] (gene object): Gene object
    """
    data = EnsemblRelease(release)
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
        gene [in] (gene object):  Gene object
        transcripts [out] (list): List of transcript objects
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
    
def get_coding_sequence(transcripts):
    """
    This function takes a list of transcript objects and returns a list of
    coding sequences.
    Args:
        transcripts [in] (list):       List of transcript objects
        coding_sequences [out] (list): List of coding sequences
        """
    return [x.coding_sequence for x in transcripts]

###############################################################################
    
def get_exon_locations(exon_id, release = 108): 
    """
    This function return exon start and end of a given exonic identifer. 
    """
    data = EnsemblRelease(release)
    exonlocus = data.locus_of_exon_id(exon_id)
    
    return exonlocus
    

###############################################################################
    
def get_transcripts_dict(gene, exclude_coding = False):
    """
    This function takes a gene object and returns a dictionary of transcript
    objects, with transcript ID as keys, and exon objects as values.
    Args:
        gene [in] (gene object):   Gene object
        dTranscripts [out] (dict): Dictionary of transcript objects, with
        transcript ID as keys, and exon objects as values
    """
    d = {}
    
    # get list of transcripts to iterate
    all_transcripts = get_transcript_from_gene(gene)
    if exclude_coding == False: 
        tcripts = all_transcripts
    else: 
        tcripts = get_coding_transcript(all_transcripts)
    

    for tcript in tcripts:
        d[tcript.id] = "-".join(get_exons_from_transcript(tcript))
     
    return d

###############################################################################

def constructu_target_cdna(masked_chr, gene_obj, transcript, exon_junction): 
    """
    This function takes a transcript and an exon junction inside this transcript
    and returns the complete transcript cDNA + the index of the junction on the
    sequence. 
    Args: 
        masked_chr [in] (str)    Full path to the masked chromosome files
        gene_obj [in] (Gene obj) Gene object returned by ensembl
        transcript [in] (str)    Ensembl transcript ID
        exon_junction [in] (str) Ensembl exon IDs (e.g. ENS001-ENS002)
        cdna [out] (str)         Complete cDNA of the transcript 
        junction_i [out] (int)   exon_junction location on the cdna
    """
    # t is a transcript object, and transcript is a string
    t = [x for x in gene_obj.transcripts if x.transcript_id == transcript]
    
    if len(t) > 0: 
        t = t[0]

        cdna, junction_i, found_junction = "", 0, False # initialize
        
        chrom_open = open(masked_chr.format(gene_obj.contig), "r")
        tt = chrom_open.read() # full chromosome sequence
        chrom_open.close()
        
        if gene_obj.on_positive_strand: 
            for exon in t.exons: 
                cdna += tt[exon.start:exon.end+1]
                
                if found_junction == False: 
                    junction_i += exon.end - exon.start + 1
                
                if exon.exon_id in exon_junction[0]: 
                    found_junction = True # stop summing on junction index
        
        else: # reverse strand genes
            for exon in reversed(t.exons): 
                print("exon: {} - {}".format(exon.start, exon.end))
                cdna += tt[exon.start:exon.end+1]
                
                if found_junction == False: 
                    junction_i += exon.end - exon.start + 1
                
                if exon.exon_id in exon_junction[0]: 
                    found_junction = True # stop summing on junction index            
                    
    
    else: 
        print("Transcript not found in gene")
        cdna, junction_i = None, None
        
    return cdna, junction_i
    
    
    
    
    
    
    
    
    
    
    
    