#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 16:38:33 2023

@author: ecrisru
"""

from ExonSurfer.resources import resources

def extract_canonical(gene_obj): 
    """
    This function takes a gene object and searches the canonical transcript in 
    the ensembl database. 
    Args: 
        gene_obj [in] (gene obj)  Pyensembl gene object
        canonical_t [out] (str|l) Transcript identifier (str) or empty list
    """
    cfile = open(resources.CANONICAL(), "r")
    lines = cfile.read().split("\n")
    canonical_t = [  # list comprehension here
        l.split("\t")[1] 
        for l 
        in lines 
        if gene_obj.gene_id in l and "Ensembl Canonical" in l
        ]
    
    if len(canonical_t) > 0: 
        canonical_t = canonical_t[0].split(".")[0]
        
    return canonical_t