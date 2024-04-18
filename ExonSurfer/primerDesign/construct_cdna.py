#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 15:16:19 2023

@author: q0791lp
"""
###############################################################################
#                   construct DNA FUNCTION DEFINITION SITE                    #
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
        transcript [in] (l|str)  List of ensembl transcript IDs or ALL
        exon_junction [in] (str) Ensembl exon IDs (e.g. ENS001-ENS002)
        to_return [out] (l)      List of tuples, each tuple formed by: 
            - cdna (str)       Complete cDNA of the transcript 
            - junction_i (int) exon_junction location on the cdna
            - exon_len (l)     List of exon ids and pos range from the junction
    """
    to_return = [] 
    
    # read chromosome 
    chrom_open = open(masked_chr.format(gene_obj.contig), "r")
    tt = chrom_open.read() # full chromosome sequence
    chrom_open.close()
            
    if transcript == "ALL" or len(transcript) > 1: 
        
        # initialize all 
        cdna, exon_len, tosum = "", [], 0 
        
        # check strand
        if gene_obj.on_positive_strand: 
            list_of_exons = exon_junction[0].split("-")
        else: 
            list_of_exons = exon_junction[0].split("-")[::-1]
         
        # iterate all exons in the junction and build cdna
        for exon in list_of_exons: # exon is string id
            exon_obj = data.exon_by_id(exon)
            cdna += tt[exon_obj.start-1:exon_obj.end]
            
            length = exon_obj.end - exon_obj.start + 1
            elen = (exon_obj.exon_id, range(tosum, tosum+length))
            exon_len.append(elen)
            tosum += length
        
        # construct distance-info objects
        e_obj1 = data.exon_by_id(list_of_exons[0])
        e_objF = data.exon_by_id(list_of_exons[-1])
        ji1 = e_obj1.end - e_obj1.start + 1
        ji2 = len(cdna) - (e_objF.end - e_objF.start + 1)
        
        # build list to return
        if ji1 != ji2: 
            to_return = [(cdna, ji1, exon_len, exon_junction[0], exon_junction[1]), 
                         (cdna, ji2, exon_len, exon_junction[0], exon_junction[1])]
        else: 
            to_return = [(cdna, ji1, exon_len, exon_junction[0], exon_junction[1])]
        
    else:
        # t is a target transcript object
        t_filt = [x for x in gene_obj.transcripts if x.transcript_id == transcript[0]]
        
        if len(t_filt) > 0: 
            t = t_filt[0]
            
            # initialize all
            cdna, junction_i, found_junction, exon_len, tosum = "", 0, False, [] , 0
            
            # check strand
            if gene_obj.on_positive_strand: 
                list_of_exons = t.exons
            else: 
                list_of_exons = reversed(t.exons)
            
            # iterate ALL exons in the transcript and build cdna
            for exon in list_of_exons: # exon is pyensembl object
                cdna += tt[exon.start-1:exon.end]
                
                elen = (exon.exon_id, 
                        range(tosum, tosum+exon.end - exon.start + 1))
                exon_len.append(elen)
                tosum += exon.end - exon.start + 1
                
                if found_junction == False: 
                    junction_i += exon.end - exon.start + 1
                
                if exon.exon_id in exon_junction[0]: 
                    found_junction = True # stop summing on junction index
                           
            to_return = [(cdna, junction_i, exon_len, exon_junction[0], exon_junction[1])]
            
        else: 
            print("Transcript not found in gene")
            to_return = [(None, None, None, exon_junction[0], exon_junction[1])]
    
    return to_return

###############################################################################

def construct_one_exon_cdna(masked_chr, gene_obj, data, transcript, i, window): 
    """
    This function returns the CDNA for the transcripts that only have one exon. 
    Args: 
        masked_chr [in] (str)    Full path to the masked chromosome files
        gene_obj [in] (Gene obj) Gene object returned by ensembl
        data [in] (Genome obj)   Genome object returned by ensembl
        transcript [in] (str)    Ensembl transcript ID
        i [in] (int)             Iteration round (to decide index)
        window [in] (int)        An index is defined every window bases
        cdna [out] (str)         Complete cDNA of the transcript 
        junction_i [out] (int)   Exon_junction location on the cdna
        exon_len [out] (l)       List of exon ids and pos range from the junction
    """
    
    # get only transcript object
    if transcript == "ALL": 
        t_obj = gene_obj.transcripts[0]
    else: 
        t_obj = data.transcript_by_id(transcript[0])
    
    # read chromosome 
    chrom_open = open(masked_chr.format(t_obj.contig), "r")
    tt = chrom_open.read() # full chromosome sequence
    chrom_open.close()
    
    cdna = tt[t_obj.exons[0].start-1:t_obj.exons[0].end]
    
    junction_i = int(window * i) # define "target" position
    exon_len = [(t_obj.exons[0].exon_id, 
                 range(0, t_obj.exons[0].end - t_obj.exons[0].start + 1))]
    
    return cdna, junction_i, exon_len
        
    