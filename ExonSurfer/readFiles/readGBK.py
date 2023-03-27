#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 13:27:28 2023

@author: q0791lp
"""
# imported modules
from Bio import SeqIO

###############################################################################
#                  GenBank information FUNCTION DEFINITION SITE               #
############################################################################### 

def read_genbank_file(gb_file): 
    """
    This function reads GenBank files. 
    Args: 
        gb_file [in] (str)          Full path to the GenBank file
        gb_record [out] (SeqRecord) GenBank record
    """
    gb_record = SeqIO.read(gb_file, "genbank")
    return gb_record

###############################################################################

def extract_junctions(gb_record): 
    """
    This function takes a SeqRecord object representing a GenBank format and 
    extracts the exon junctions locations (NOT exon length)
    Args: 
        gb_record [in] (SeqRecord) GenBank file content
        exon_junctions [out] (l)   List of integers, with one int per junction
    """
    exon_junctions = []
    for feature in gb_record.features:
        if feature.type == "exon":
            end = feature.location.end
            exon_junctions.append(int(end))
            
    return exon_junctions[:-1]

###############################################################################

def extract_cdna(gb_record): 
    """
    This function takes a SeqRecord object representing a GenBank format and 
    extract the cDNA sequence from the exons in the GenBank record. 
    Args: 
        gb_record [in] (SeqRecord) GenBank file content
        exon_sequences [out] (str) Exon sequences.
    """
    exon_sequences = ""
    for feature in gb_record.features:
        if feature.type == "exon":
            exon_location = feature.location
            exon_sequence = gb_record.seq[exon_location.start:exon_location.end]
            if exon_location.strand == -1:
                exon_sequence = exon_sequence.reverse_complement()
            exon_sequences += str(exon_sequence)
    
    return exon_sequences
