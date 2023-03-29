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

###############################################################################

def get_junction_i(jdict, junction, nm_dna, ji): 
    """
    This function returns only the sequence for a specific junction. 
    Args: 
        jdict [in] (d)         Dict of junctions as keys and indexes as values
        junction [in] (str)    Junction identifier (Ex "EXON1-EXON2")
        nm_dna [in] (str)      Complete cDNA
        ji [in] (l)            List of junction indixes
        cut_nm_dna [out] (str) Junction sequence
        cut_ji [out] (l)       Indexes corresponding to the cut_nm_dna
    """
    fe = junction.split("-")[0] # first exon
    le = junction.split("-")[-1] # last exon
    
    c1 = [jdict[k] for k in jdict.keys() if k.split("-")[1] == fe][0]
    c2 = [jdict[k] for k in jdict.keys() if k.split("-")[0] == le][0]
    
    cut_nm_dna = nm_dna[c1:c2]
    cut_ji = [x for x in ji if x in range(c1+1, c2)]
    cut_ji = [x-c1 for x in cut_ji]
    
    return cut_nm_dna, cut_ji
    
    