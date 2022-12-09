#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 13:10:49 2022

@author: ecrisru
"""

# imported modules
import re
import primer3

# my constants
from .designConfig import design_dict as DESIGN_DICT

# Constants
FLANK = 400
CASE2_PRIMERS = 100

###############################################################################
#                 primerDesign module FUNCTION DEFINITION SITE                #
###############################################################################

    
def call_primer3(target_seq, junction_i): 
    """
    This function calls primers3 for 2 design options: (1) primers are placed 
    ON the exon junction and (2) primers are placed FLANKING the exon junction. 
    Args: 
        target_seq [in] (str)      Full cDNA sequence of the transcript
        junction_i [in] (int)      Index of the exonic junction on the target_seq
        case1_primers [out] (dict) Dictionary of primers designed with option 1
        case2_primers [out] (dict) Dictionary of primers designed with option 2
    """
    
    # Case 1: place forward primer or right primer ON junction
    target_dict = {
            'SEQUENCE_ID': 'InternalID',
            'SEQUENCE_TEMPLATE': target_seq,
            'SEQUENCE_OVERLAP_JUNCTION_LIST': junction_i,
            'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION': 4,
            }

    # PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION=4]
    case1_primers = primer3.bindings.designPrimers(target_dict, DESIGN_DICT)
    
    # Case 2: place each primer on one exon
    target_dict = {
            'SEQUENCE_ID': 'InternalID',
            'SEQUENCE_TEMPLATE': target_seq,
            'SEQUENCE_TARGET': [junction_i, 1],
            }
    case2_primers = primer3.bindings.designPrimers(target_dict, DESIGN_DICT)
    
    return case1_primers, case2_primers

###############################################################################

def write_blast_fasta(case1_primers, case2_primers, fastaf): 
    """
    This function writes a fasta file with the designed primers. 
    Args: 
        case1_primers [in] (dict) Dictionary of primers designed with option 1
        case2_primers [in] (dict) Dictionary of primers designed with option 2
        fastaf [in] (str)         Full path to the fasta file to write
    """
    with open(fastaf, "w") as f_open: 
        for k in case1_primers: 
            if "PRIMER" in k and "SEQUENCE" in k: 
                string = ">{}\n{}\n".format(k, case1_primers[k])
                f_open.write(string)
                
        for k in case2_primers: 
            if "PRIMER" in k and "SEQUENCE" in k: 
                
                # this is done to differenciate between case1 and case2
                pnum_toreplace = re.search("_\d+_", k).group()
                new_pnum = int(pnum_toreplace.replace("_", "")) + CASE2_PRIMERS
                new_pnum = "_{}_".format(new_pnum)
                
                newk = k.replace(pnum_toreplace, new_pnum)
                
                string = ">{}\n{}\n".format(newk, case1_primers[k])
                f_open.write(string)
                
###############################################################################
                   
def report_design(c1, c2, exon_junction, ofile):    
    """
    This function writes a table with the designed primers
    Args: 
        c1 [in] (dict)             Dict of primers designed with option 1
        c2 [in] (dict)             Dict of primers designed with option 2
        exon_junction [in] (tuple) Tuple with two strings: (1) "Ensemble exon ID
                                   - Ensembl exon ID" (2) Design specification
        ofile [in] (str)           Path to the output file
    """
    out_open = open(ofile, "a+")
    
    for pdict in (c1, c2): 
        for n in range(DESIGN_DICT["PRIMER_NUM_RETURN"]): 
            patt = "_{}_".format(n)
            forseq = [pdict[k] for k in pdict if "PRIMER_LEFT{}SEQUENCE".format(patt) in k][0]
            revseq = [pdict[k] for k in pdict if "PRIMER_RIGHT{}SEQUENCE".format(patt) in k][0]
            
            ampsize = [pdict[k] for k in pdict if "{}PRODUCT_SIZE".format(patt) in k][0]
            if pdict == c1: 
                num = str(n)
                opt = "1"
            else: 
                num = str(n + CASE2_PRIMERS)
                opt = "2"
            
            string = "\t".join((opt, num, exon_junction[0], forseq, revseq, 
                               str(ampsize), exon_junction[1])) + "\n"
            
            out_open.write(string)
            
    out_open.close()
    
    