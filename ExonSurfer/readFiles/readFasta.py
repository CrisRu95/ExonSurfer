#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 15:49:37 2023

@author: q0791lp
"""

###############################################################################
#                   Fasta information FUNCTION DEFINITION SITE                #
############################################################################### 

def read_fasta_file(f_file): 
    """
    This function reads fasta files. 
    Args: 
        f_file [in] (str)  Full path to the GenBank file
        header [out] (str) Header information (separated by blank spaces)
        seq [out] (str)    cDNA sequence
    """
    seq = ""
    with open(f_file, "r") as f_open: 
        for line in f_open.readlines(): 
            if ">" in line: 
                header = line.rstrip().replace(">", "")
            else: 
                seq += line.rstrip()

    return header, seq

###############################################################################

def extract_junctions(header): 
    """
    This function takes a header from a fasta file and extracts the exon 
    junctions locations (NOT exon length)
    Args: 
        header [in] (str)        Header information (separated by blank spaces)
        exon_junctions [out] (l) List of integers, with one int per junction
    """
    exon_junctions = header.split(" ")[1:]
    exon_junctions = [int(x) for x in exon_junctions]
            
    return exon_junctions

