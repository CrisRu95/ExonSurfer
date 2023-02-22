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

###############################################################################

def construct_target_cdna(masked_chr, gene_obj, data, transcript, exon_junction): 
    """
    This function takes a transcript and an exon junction inside this transcript
    and returns the complete transcript cDNA + the index of the junction on the
    sequence. 
    Args: 
        masked_chr [in] (str)    Full path to the masked chromosome files
        gene_obj [in] (Gene obj) Gene object returned by ensembl
        data [in] (Genome obj)   Genome object returned by ensembl
        transcript [in] (str)    Ensembl transcript ID
        exon_junction [in] (str) Ensembl exon IDs (e.g. ENS001-ENS002)
        cdna [out] (str)         Complete cDNA of the transcript 
        junction_i [out] (int)   exon_junction location on the cdna
        exon_len [out] (l)       List of exon ids and pos range from the junction
    """
    # t is a transcript object, and transcript is a string
    #l = [str(x.transcript_id) for x in gene_obj.transcripts ]
    #print(set(l))
    
    # read chromosome 
    chrom_open = open(masked_chr.format(gene_obj.contig), "r")
    tt = chrom_open.read() # full chromosome sequence
    chrom_open.close()
            
    if transcript == "ALL": 
        # exon_junction[1] contains information on the transcripts NOT detected by that junction
        e1 = exon_junction[0].split("-")[0]
        e2 = exon_junction[0].split("-")[1]
        
        exon_obj1 = data.exon_by_id(e1)
        exon_obj2 = data.exon_by_id(e2)
        
        if gene_obj.on_positive_strand: 
            cdna = tt[exon_obj1.start-1:exon_obj1.end] + tt[exon_obj2.start-1:exon_obj2.end]
            junction_i = exon_obj1.end - exon_obj1.start
            e1_len = exon_obj1.end - exon_obj1.start + 1
            e2_len = exon_obj2.end - exon_obj2.start + 1
            exon_len = [(exon_obj1.exon_id, range(0, e1_len)), 
                        (exon_obj2.exon_id, range(e1_len, e1_len + e2_len))]
        else: 
            cdna = tt[exon_obj2.start-1:exon_obj2.end] + tt[exon_obj1.start-1:exon_obj1.end]
            junction_i = exon_obj2.end - exon_obj2.start     
            e1_len = exon_obj1.end - exon_obj1.start + 1
            e2_len = exon_obj2.end - exon_obj2.start + 1
            
            exon_len = [(exon_obj2.exon_id, range(0, e2_len)), 
                        (exon_obj1.exon_id, range(e2_len, e2_len + e1_len))]
    else:
        # t is a target transcript object
        t_filt = [x for x in gene_obj.transcripts if x.transcript_id == transcript]
        
        if len(t_filt) > 0: 
            t = t_filt[0]
    
            cdna, junction_i, found_junction, exon_len = "", 0, False, [] # initialize
            
            if gene_obj.on_positive_strand: 
                tosum = 0
                for exon in t.exons: 
                    cdna += tt[exon.start-1:exon.end]
                    
                    elen = (exon.exon_id, 
                            range(tosum, exon.end - exon.start + 1))
                    exon_len.append(elen)
                    tosum += exon.end - exon.start + 1
                    
                    if found_junction == False: 
                        junction_i += exon.end - exon.start + 1
                    
                    if exon.exon_id in exon_junction[0]: 
                        found_junction = True # stop summing on junction index
            
            else: # reverse strand genes
                for exon in reversed(t.exons): 
                    tosum = 0
                    cdna += tt[exon.start-1:exon.end]
                    
                    elen = (exon.exon_id, 
                            range(tosum, exon.end - exon.start + 1))
                    exon_len.append(elen)
                    tosum += exon.end - exon.start + 1
                    
                    if found_junction == False: 
                        junction_i += exon.end - exon.start + 1
                    
                    if exon.exon_id in exon_junction[0]: 
                        found_junction = True # stop summing on junction index            

        else: 
            print("Transcript not found in gene")
            cdna, junction_i, exon_len = None, None, None
    
    return cdna, junction_i, exon_len

###############################################################################

def construct_one_exon_cdna(masked_chr, gene_obj, data, transcript): 
    """

    """
    # get only transcript object
    if transcript == "ALL": 
        t_obj = gene_obj.transcripts[0]
    else: 
        t_obj = data.transcript_by_id(transcript)
    
    # read chromosome 
    chrom_open = open(masked_chr.format(t_obj.contig), "r")
    tt = chrom_open.read() # full chromosome sequence
    chrom_open.close()
    

    cdna = tt[t_obj.exons[0].start-1:t_obj.exons[0].end]
    
    junction_i = int(len(cdna) / 2)
    exon_len = [(t_obj.exons[0].exon_id, 
                 range(0, t_obj.exons[0].end - t_obj.exons[0].start + 1))]
    
    return cdna, junction_i, exon_len
        
    
    
    