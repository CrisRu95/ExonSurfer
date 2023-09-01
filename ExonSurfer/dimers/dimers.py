#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 15:32:42 2021

@author: Elena Cristina Rusu Hutu
"""

# Constants
MATCH_DICT = {"A": "T", "C": "G", "G": "C", "T": "A"}

###############################################################################
#                     FUNCTIONS DEFINITION SITE                               #
###############################################################################

def find_max_complementary_bases(seq1, seq2):    
    """
    Returns maximum number of matches between two sequences. 
    """
    max_matches = 0
    
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            matches = 0
            while (i + matches < len(seq1)) and (j - matches >= 0) and \
            (seq1[i + matches] == MATCH_DICT[seq2[j - matches]]):
                matches += 1
            max_matches = max(max_matches, matches)
    
    return max_matches

###############################################################################

def get_max_comp(forward, reverse): 
    """
    Call find_max_complementary_bases 3 times for each combination. 
    Returns maximum number of matches. 
    """
    n1 = find_max_complementary_bases(forward, reverse)
    n2 = find_max_complementary_bases(forward, forward)
    n3 = find_max_complementary_bases(reverse, reverse)

    return max((n1, n2, n3))





