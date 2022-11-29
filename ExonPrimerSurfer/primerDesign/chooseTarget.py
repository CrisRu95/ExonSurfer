#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 16:35:07 2022

@author: Elena Cristina Rusu 
Exon Primer Surfer
"""


###############################################################################
#           read input and choose target FUNCTION DEFINITION SITE             #
###############################################################################

def read_input(file): 
    """
    This function takes a file where transcript ids are separated from the exon
    ids with tabs and returns a dictionary. 
    Args: 
        file [in] (str): Complete path to the transcript list file
        d [out] (dict): Transcript IDs as keys, exon list (ORDERED) as values 
    """
    d = {}
    with open(file, "r") as f_open: 
        for line in f_open.readlines(): 
            if line != "": 
                line = line.rstrip().split("\t")
                d[line[0]] = line[1]
    
    return d

###############################################################################
    
def format_junctions(d): 
    """
    This function creates the dictionary junctions. 
    Args: 
        d [in] (dict): Transcript IDs as keys, exon list (ORDERED) as values 
        junctions [out] (dict): Exon junctions as keys, transcript IDs as values
    """
    junctions = {}
    
    for transcript in d: 
        exons = d[transcript].split("-")
        for i in range(0, len(exons)-1): 
            key = "{}-{}".format(exons[i], exons[i+1])
            if key not in junctions: 
                junctions[key] = [transcript]
            else: 
                 junctions[key].append(transcript)
    
    return junctions

###############################################################################
    
def choose_target(d, junctions, to_detect): 
    """
    This function chooses the best exon junctions to target with primer design, 
    depending whether we want to target 1 transcript in particular or ALL. If 
    there are no specific / universal junctions, the function returns the next 
    most acceptable solution. 
    Args: 
        d [in] (dict): Keys are transcripts, values are exons separated by "-"
        junctions [in] (dict): Keys are exon junctions, values are transcript ids
        to_detect [in] (str): transcript ID or "ALL"
        toreturn [out] (list): list of tuples, each tuple with a junction and 
                               a explanation message
    """
    toreturn, p_sols = [], [] # initialize
    
    # Case 1: Detect only one transcript ID
    if to_detect in d.keys(): 
        perf_j = [j for j in junctions if junctions[j] == [to_detect]]
        
        if len(perf_j) >= 1:  # solution found, no need for loop 
            p_sols = perf_j
        
        else: 
            i, found = 2, False    # loop initialization
            while i <= len(d) and found == False: 
                less_perf_j = [j for j in junctions if to_detect in junctions[j] and len(junctions[j]) == i]
                if len(less_perf_j) >= 1: 
                    found = True
                    p_sols = less_perf_j 
                i += 1
        
        for sol in p_sols: 
            if junctions[sol] == [to_detect]: 
                t = (sol, "unique to target")
            else: 
                string = ", ".join([x for x in junctions[sol] if x != to_detect])
                t = (sol, "also detects: {}".format(string))
            toreturn.append(t)
        
    # Case 2: try to detect all transcript
    elif to_detect == "ALL": 
        perf_j = [x for x in junctions if len(junctions[x]) == len(d)]

        if len(perf_j) >= 1: # solution found, no need for loop 
            p_sols = perf_j
        
        else: 
            i, found = len(d), False    # loop initialization
            while i >= 1 and found == False: 
                less_perf_j = [j for j in junctions if len(junctions[j]) == i]
                if len(less_perf_j) >= 1: 
                    found = True
                    p_sols = less_perf_j 
                i -= 1
        
        for sol in p_sols: 
            if len(junctions[sol]) == len(d): 
                t = (sol, "targets all transcripts")
            else: 
                string = ", ".join([x for x in d.keys() if x not in junctions[sol]])
                t = (sol, "does NOT detect: {}".format(string))
            toreturn.append(t)
    
    # Error: exit function            
    else:
        print("The tanscript identifier is NOT valid")
        return
    
    
    return toreturn

