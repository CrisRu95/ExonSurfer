#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 17:31:44 2023

@author: Elena Cristina Rusu 
"""
# imported modules
import re

# own modules
from ExonSurfer.ensembl import ensembl
from ExonSurfer.resources import resources

###############################################################################
#                        CONSTANTS FOR THE HTML FILE                          #
###############################################################################

UPPER_LINES =['<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">', 
              '<html>', 
              '<head>', 
              
              '<style>', 
              '@page { size: 8.27in 11.69in; margin: 0.79in }',
              'div {text-align:justify; word-break: break-all; max-width:40em; line-height: 20%;}',   
              'span{ width:100%; height:1em;}', 
              
              '.ex1 {color: #7b3294; display: inline;}', 
              '.ex1H {color: #7b3294; background-color: #c2a5cf; display: inline;}', 
              '.ex2 {color: #008837; display: inline;}', 
              '.ex2H {color: #008837; background-color: #a6dba0; display: inline;}', 
              
              '.divS { width: 30em; max-width: 100%; margin-bottom: 0in; line-height: 100%}', 
              '.two p:nth-child(1) { float:left; }', 
              '.two p:nth-child(2) { float:right; }', 
              
              '.header {font-family:"Liberation Sans", sans-serif; size: 14pt; }',
              '</style>',
              
              '<meta http-equiv="content-type" content="text/html; charset=utf-8"/>', 
              '<title>{}</title>',
              '<meta name="generator" content="LibreOffice 6.4.7.2 (Linux)"/>', 
              '<meta name="created" content="2022-03-31T10:34:31.650656446"/>', 
              '<meta name="changed" content="2022-05-25T14:09:40.676889449"/>', 
              '</head>', 
              
              '<body lang="en-US" link="#000080" vlink="#800000" dir="ltr">', 
              '<h3 class="header"> {} {}</h3>',
              '<div class="divS">', 
              '<font face="Liberation Sans, sans-serif" size="2" style="font-size: 10pt">']


FINAL_LINES = ['</font>', '</div>', '</body>', '</html>']


###############################################################################
#                          Function definition site                           #
###############################################################################

def reverse_complement(seq):
    """
    This function returns a sequence read from right to left. 
    Args: 
        seq [in] (str) Sequence to be inverted. 
    """
    seq = seq.upper()
    compdict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'} 
    bases = [compdict[base] for base in seq]
    bases = "".join(bases)
    return bases[::-1]     
                
###############################################################################

def get_junction_seqs(junction, masked_chr, data): 
    """
    This func obtains the exon sequences and prints them in different colors. 
    Args: 
        junction [in] (str)    Exon IDs separated by "-"
        masked_chr [in] (str)  Path to the chr files, should have "{}" 
        data [in] (Genome obj) Genome object returned by ensembl
        target_s [out] (str)   Exon seqs marked in diff colors (html style)
    """
    exon1 = junction.split("-")[0]
    exon2 = junction.split("-")[1]
    
    e_obj1 = data.exon_by_id(exon1)
    e_obj2 = data.exon_by_id(exon2)
    
    # read chromosome 
    chrom_open = open(masked_chr.format(e_obj1.contig), "r")
    tt = chrom_open.read() # full chromosome sequence
    chrom_open.close()
    
    if e_obj1.on_positive_strand: 
        nm_dna = tt[e_obj1.start-1:e_obj1.end] + tt[e_obj2.start-1:e_obj2.end]
        ji = e_obj1.end - e_obj1.start
    else: 
        nm_dna = tt[e_obj2.start-1:e_obj2.end] + tt[e_obj1.start-1:e_obj1.end]
        ji = e_obj2.end - e_obj2.start  

    return nm_dna, ji
    
    
###############################################################################

def get_primers_i(junction_dna, forward, reverse): 
    """
    This func returns the start and end of forward and reverse primers on the 
    junction_dna
    Args: 
        junction_dna [in] (str) Junction seq
        forward [in] (str)      Forward primer seq
        reverse [in] (str)      Reverse primer seq
        f1, f2 [out] (int)      Forward start and end indices
        r1, r2 [out] (int)      Reverse start and end indices
    """
    
    f1 = re.search(forward, junction_dna, re.I).start()
    f2 =  re.search(forward, junction_dna, re.I).end()

    r1 = re.search(reverse_complement(reverse), junction_dna, re.I).start()
    r2 =  re.search(reverse_complement(reverse), junction_dna, re.I).end()
    
    return f1, f2, r1, r2

###############################################################################

def create_par_string(nm_dna, option, ji, f1, f2, r1, r2): 
    
    if option == 1: 
        if f2 > ji: # the forward is ON the junction
            
            string = '<p class="ex1">' + nm_dna[:f1] + '</p>'
            string += '<p class="ex1H">' + nm_dna[f1:ji] + '</p>' # forward ini
            string += '<p class="ex2H">' + nm_dna[ji:f2] + "➛"+ '</p>' # forward final
            string += '<p class="ex2">' + nm_dna[f2:r1] + '</p>'
            string += '<p class="ex2H">' + "￩"+nm_dna[r1:r2] + '</p>' # full reverse
            string += '<p class="ex2">' + nm_dna[r2:] + '</p>'
        
        else: # the reverse is ON the junction
        
            string = '<p class="ex1">' + nm_dna[:f1] + '</p>'
            string += '<p class="ex1H">' + nm_dna[f1:f2] + "➛" + '</p>' # full forward
            string += '<p class="ex1">' + nm_dna[f2:r1] + '</p>' 
            string += '<p class="ex1H">' + "￩"+nm_dna[r1:ji] + '</p>' # reverse on ex1
            string += '<p class="ex2H">' + nm_dna[ji:r2] + '</p>' # reverse on ex2
            string += '<p class="ex2">' + nm_dna[r2:] + '</p>'
    else: 
        string = '<p class="ex1">' + nm_dna[:f1] + '</p>'
        string += '<p class="ex1H">' + nm_dna[f1:f2] + "➛"+ '</p>' # full forward
        string += '<p class="ex1">' + nm_dna[f2:ji] + '</p>' 
        string += '<p class="ex2">' + nm_dna[ji:r1] + '</p>'
        string += '<p class="ex2H">' + "￩"+nm_dna[r1:r2] + '</p>' # full reverse
        string += '<p class="ex2">' + nm_dna[r2:] + '</p>'      
        
    return string

###############################################################################

def highlight_primers(pair_id, final_df, species, release, outfile): 
    """
    """
    junction = final_df.loc[pair_id]["junction"]
    masked_chr = resources.MASKED_SEQS(species)
    data = ensembl.create_ensembl_data(release, 
                                       species.replace("_masked", ""))
    
    nm_dna, ji = get_junction_seqs(junction, masked_chr, data)
    f1, f2, r1, r2 = get_primers_i(nm_dna, final_df.loc[pair_id]["forward"], 
                                   final_df.loc[pair_id]["reverse"])
    

    string = create_par_string(nm_dna, final_df.loc[pair_id]["option"],
                               ji, f1, f2, r1, r2)    
    
    with open(outfile, "w") as file_open: 
        for l in UPPER_LINES: 
            if "title" in l: 
                file_open.write(l.format(pair_id))
            elif "h3" in l: 
                file_open.write(l.format(pair_id, junction))
            else: 
                file_open.write(l)

        file_open.write(string)
        for l in FINAL_LINES: 
            file_open.write(l)


