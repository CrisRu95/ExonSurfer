#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 16:35:07 2022

@author: Elena Cristina Rusu 
Exon Primer Surfer


mock inputfile: 
T1	ENSE1-ENSE2-ENSE4
T2	ENSE1-ENSE2-ENSE3-ENSE4
T3	ENSE1-ENSE2
T4	ENSE1-ENSE2-ENSE3
"""


###############################################################################
#           read input and choose target FUNCTION DEFINITION SITE             #
###############################################################################

def get_junction_len(key, data):
    """
    This function returns the length of a single or multi exon junction, where
    exon ids are separated by a hyphen. 
    Args: 
        key [in] (str)         Two or more exon ids separated by a hyphen
        data [in] (Genome obj) Genome object returned by ensembl
        j_len [out] (int)      Combined length of the exons
    """
    j_len = 0
    for exon in key.split("-"): 
        e_len = data.exon_by_id(exon).end - data.exon_by_id(exon).start
        j_len += e_len
    
    return j_len

###############################################################################

def format_junctions(d, to_detect, opt_amp_len, data): 
    """
    This function creates the dictionary junctions. 
    Args: 
        d [in] (dict)          Trans IDs as keys, exon list (ORDERED) as vals 
        to_detect [in] (str)   Transcript ID or "ALL"
        opt_amp_len [in] (int) Optimum amplicon length 
        data [in] (Genome obj) Genome object returned by ensembl
        junctions [out] (dict) Exon junctions as keys, transcript IDs as values
    """
    junctions = {}
    
    # exon junctions of 2 exons, in ensembl.construct_target_cdna we use all cdna
    for transcript in d: 
        exons = d[transcript].split("-")
        for i in range(0, len(exons)-1): 
            key = "{}-{}".format(exons[i], exons[i+1])
            
            # we will NOT use all the cdna but only the junction
            # so we append multiple junctions in case opt_amp_len is longer      
            if to_detect == "ALL": 
                j_len = get_junction_len(key, data)
                
                j = i + 2 # initialize
                while j_len < opt_amp_len and j < len(exons):
                    key += "-{}".format(exons[j]) # append exon to key
                    j_len = get_junction_len(key, data)
                    j += 1  
                    
            if key not in junctions: 
                junctions[key] = [transcript]
            else: 
                 junctions[key].append(transcript)
                     
    return junctions

###############################################################################
    
def choose_target(d, junctions, to_detect, canonical_t): 
    """
    This function chooses the best exon junctions to target with primer design, 
    depending whether we want to target 1 transcript in particular or ALL. If 
    there are no specific / universal junctions, the function returns the next 
    most acceptable solution. 
    Args: 
        d [in] (dict)         Keys are trans, values are exons separated by "-"
        junctions [in] (dict) Keys are exon junctions, values are transcript ids
        to_detect [in] (str)  Transcript ID or "ALL"
        canonical_t [in] (l)  List with the canonical transcript (can be empty)
        toreturn [out] (l)    List of tuples, each tuple with a junction and 
                              a explanation message
    """
    toreturn, p_sols = [], [] # initialize
    
    # Case 1: Detect only one transcript ID
    if to_detect in d.keys():
        print("Detecting only {}".format(to_detect)) 
        perf_j = [j for j in junctions if junctions[j] == [to_detect]]
        
        if len(perf_j) >= 1:  # solution found, no need for loop 
            p_sols = perf_j
        
        else: 
            i, found = 2, False    # loop initialization
            while i <= len(d) and found == False: 
                less_perf_j = [j 
                               for j 
                               in junctions 
                               if to_detect in junctions[j] and len(junctions[j]) == i]
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
        print("Detecting all transcripts")
        perf_j = [x for x in junctions if len(junctions[x]) == len(d)]

        if len(perf_j) >= 1: # solution found, no need for loop 
            p_sols = perf_j
        
        else: # no solution found, need for filters
        
            # keep solutions where the canonical is present
            if len(canonical_t) > 0: 
                if len([j 
                        for j 
                        in junctions 
                        if canonical_t[0].split(".")[0] in junctions[j]]) > 0: 
                    todel = [j 
                             for j 
                             in junctions 
                             if canonical_t[0].split(".")[0] not in junctions[j]]
                    for key in todel: 
                        del junctions[key]
                            
            # keep solutions where the most other transcripts are present
            i, found = len(d), False    # loop initialization
            while i >= 1 and found == False: 
                less_perf_j = [j for j in junctions if len(junctions[j]) == i]
                if len(less_perf_j) >= 1: 
                    found = True
                    p_sols = less_perf_j 
                i -= 1
                
        # annotate junctions
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
    
    return toreturn


