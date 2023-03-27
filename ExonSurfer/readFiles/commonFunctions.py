#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 15:59:03 2023

@author: q0791lp
"""
###############################################################################
#               Files information common FUNCTION DEFINITION SITE             #
############################################################################### 

def get_junc_dict(exon_junctions): 
    """
    This function transforms the junctions list into a dictionary with exon ids
    Args: 
        exon_junctions [in] (l) List of integers, with one int per junction
        junctions_d [out] (d)   Dict of junctions id and locations
    """
    junctions_d = {}
    for i in range(0, len(exon_junctions)-1): 
        key = "EXON{}-EXON{}".format(i+1, i+2)        
        junctions_d[key] = exon_junctions[i]
        
    return junctions_d    

###############################################################################

def get_elen(exon_junctions, exon_sequence): 
    """
    This function return a exon_length object, a list with tuples with 1 tuple
    per exon and 2 elements per tuple: exon identifier and a range object. 
    Args: 
        exon_junctions [in] (l)  List of integers
        exon_sequence [in] (str) cDNA sequence
        elen [out] (l)           E.g.: [("EXON1", range(0, 10)), 
                                        ("EXON2", range(10, 20))]
    """
    elen = []
    
    # Just one exon
    if len(exon_junctions) == 0: 
        t = ("EXON1", range(0, len(exon_sequence)))
        elen.append(t)
    else: 
        first_exon = ["EXON1", range(0, exon_junctions[0])]
        last_exon = ["EXON{}".format(len(exon_junctions)+1), 
                     range(exon_junctions[-1], len(exon_sequence))]
        
        elen.append(first_exon)
        for i in range(0, len(exon_junctions)-1): 
            name = "EXON{}".format(i+2) # EXON 1 is already done
            
            rng = range(exon_junctions[i], exon_junctions[i+1])
            
            t = (name, rng)
            elen.append(t)
        
        elen.append(last_exon)
    
    return elen