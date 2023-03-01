#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 17:31:44 2023

@author: Elena Cristina Rusu 
"""
# imported modules
import os
import re
import regex

# own modules
from ExonSurfer.ensembl import ensembl
from ExonSurfer.resources import resources

###############################################################################
#                        CONSTANTS FOR THE HTML FILE                          #
###############################################################################

HEADER_LINE ='<h3 class="header"> {} {}</h3>'

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
        nm_dna [out] (str)     CDNA seqs
        ji_l [out] (list)      List of ints, marking exon junctions in the nm_dna
    """
    # get pyensembl object
    e_obj = data.exon_by_id(junction.split("-")[0])
    
    # read chromosome 
    chrom_open = open(masked_chr.format(e_obj.contig), "r")
    tt = chrom_open.read() # full chromosome sequence
    chrom_open.close()
    
    # get strand
    if e_obj.on_positive_strand: 
        list_of_exons = junction.split("-")
    else: 
        list_of_exons = junction.split("-")[::-1]
    
    # build cdna and junction index
    ji_l = [] #  list of indices
    tosum = 0 # value to sum to the indices
    nm_dna = "" # cdna
    
    for exon in list_of_exons: # exon is string id
        exon_obj = data.exon_by_id(exon)
        nm_dna += tt[exon_obj.start-1:exon_obj.end]
        exon_length = exon_obj.end - exon_obj.start + 1
        ji_l.append(tosum + exon_length) # list of indices
        tosum += exon_length
        
    return nm_dna, ji_l
    
###############################################################################

def get_primers_i(dna, forward, reverse, e = 0): 
    """
    This func returns the start and end of forward and reverse primers on a DNA
    Args: 
        dna [in] (str)          DNA seq where to search
        forward [in] (str)      Forward primer seq
        reverse [in] (str)      Reverse primer seq
        e [in] (int)            Mismatches allowed (mismatches)
        f1, f2 [out] (int)      Forward start and end indices
        r1, r2 [out] (int)      Reverse start and end indices
    """
    if e == 0: 
        f1 = re.search(forward, dna, re.I).start()
        f2 =  re.search(forward, dna, re.I).end()
    
        r1 = re.search(reverse_complement(reverse), dna, re.I).start()
        r2 =  re.search(reverse_complement(reverse), dna, re.I).end()
    
    else: 
        err_patt = "){s<=" + str(e) + "}" # only substitutions allowed
        try: 
            f1 = regex.search("(?:"+forward + err_patt, dna).span()[0]
            f2 = regex.search("(?:"+forward + err_patt, dna).span()[1]
            
        except: 
            print("not found: dna")
        
        try: 
            r1 = regex.search("(?:"+reverse_complement(reverse) + err_patt, dna).span()[0]
            r2 = regex.search("(?:"+reverse_complement(reverse) + err_patt, dna).span()[1]        
        except: 
            print("not found: dna")    
    return f1, f2, r1, r2

###############################################################################

def format_indices(ji, f1, f2, r1, r2): 
    """
    This func returns the start and end of forward and reverse primers on a DNA
    Args: 
        ji [in] (l)        List of ints, marking exon junctions in the nm_dna
        f1, f2 [in] (int)  Forward start and end indices
        r1, r2 [in] (int)  Reverse start and end indices
        indices [out] (l)  List of ordered tuples (index, description)
    """    
    indices = [(i, "j") for i in ji]
    indices += [(f1, "f1")]
    indices += [(f2, "f2")]
    indices += [(r1, "r1")]
    indices += [(r2, "r2")]
    
    # order
    print("indices preorder is: {}".format(indices))
    indices = sorted(indices, key = lambda x: x[0])
    
    return indices
    
###############################################################################

def create_par_string(nm_dna, indices): 
    """
    This func creates the html string with the marked primers
    Args: 
        nm_dna [in] (str)  Junction seq
        indices [out] (l)  List of ordered tuples (index, description)
        string [out] (str) Annotated html string
    """
    # initialize 
    string = ""
    ex = 1
    h = False
    
    # start string until first index
    string = '<p class="ex{}">'.format(ex) + nm_dna[0:indices[0][0]] + '</p>'
    
    for i in range(0, len(indices)-1): 
        
        # exon - exon junction
        if indices[i][1] == "j": # here change the pointer of the exon 
            ex = 1 if ex == 2 else 2 # change exon color
            if h == True: 
                flag = '<p class="ex{}H">'.format(ex) # highlighted
            else: 
                flag = '<p class="ex{}">'.format(ex)
            string += flag + nm_dna[indices[i][0]:indices[i+1][0]] + '</p>'
        
        # start of primer (highlight)
        elif indices[i][1] == "f1" or indices[i][1] == "r1": 
            h = True
            flag = '<p class="ex{}H">'.format(ex)
            string += flag + nm_dna[indices[i][0]:indices[i+1][0]] + '</p>'
        
        # end of highlight
        elif indices[i][1] == "f2" or indices[i][1] == "r2": 
            h = False
            flag = '<p class="ex{}">'.format(ex)
            string += flag + nm_dna[indices[i][0]:indices[i+1][0]] + '</p>'

    # end sequence
    if indices[-1][1] == "j": 
        ex = 1 if ex == 2 else 2 # change exon color
        flag = '<p class="ex{}">'.format(ex)  
    elif indices[-1][1] == "f2" or indices[-1][1] == "r2": 
        flag = '<p class="ex{}">'.format(ex)   
        
    string += flag + nm_dna[indices[-1][0]:] + '</p>'
    
    return string

###############################################################################

def create_offt_string(seq, f1, f2, r1, r2, mm_i): 
    """
    This func creates the html string with the marked primers
    Args: 
        seq    [in] (str)  Offtarget DNA seq
        f1, f2 [in] (int)  Forward start and end indices
        r1, r2 [in] (int)  Reverse start and end indices   
        mm_i   [in] (l)    List of indices where mismatches are found
        string [out] (str) Sequence with marked primers and mismatches
    """

    string = '<p class="ex">' + seq[:f1] + '</p>'
    # all possible mismatches in the forward primer
    if any([x for x in mm_i if x in range(f1, f2)]):
        for_mmi = [x for x in mm_i if x in range(f1, f2)]
        start = f1
        for i in range(0, len(for_mmi)): 
            string += '<p class="exH">' + seq[start:for_mmi[i]] + '</p>'
            string += '<p class="mismatch">' + seq[for_mmi[i]] + '</p>'
            start = for_mmi[i] + 1
        string += '<p class="exH">' + seq[start:f2] + '</p>'
    else: 
        string += '<p class="exH">' + seq[f1:f2] + '</p>' # full forward
    
    string += '<p class="ex">' + seq[f2:r1] + '</p>' # in between primers
    
    if any([x for x in mm_i if x in range(r1, r2)]):
        rev_mmi = [x for x in mm_i if x in range(r1, r2)]
        start = r1
        for mi in range(0, len(rev_mmi)): 
            string += '<p class="exH">' + seq[start:rev_mmi[i]] + '</p>'
            string += '<p class="mismatch">' + seq[rev_mmi[i]] + '</p>'
            start = rev_mmi[i] + 1
        string += '<p class="exH">' + seq[start:r2] + '</p>'
    else: 
        string += '<p class="exH">' + seq[r1:r2] + '</p>' # full reverse   
        
    string += '<p class="ex">' + seq[r2:] + '</p>'
        
    # cut start and end of sequence if too long: 
    if len(seq[:f1]) > 200: 
        string = '<p class="ex">' + "..." + string[100:]
        if len(seq[r2:]) > 350: 
            string = string[:r2+200] + "..." + '</p>'
    else: 
        if len(seq[r2:]) > 350: 
            string = string[:r2+200] + "..." + '</p>'
        
    return string

###############################################################################

def obtain_offtarget_list(pair_id, final_df, species): 
    """
    This function returns a list of refseq identifiers, corresponding to the
    possible off-target amplification. 
    Args: 
        pair_id [in] (str)   Primer pair identifier (ex "Pair1")
        final_df [in] (str)  Final design DF returned by exon surfer
        species [in] (str)   Organism
        refseq_ids [out] (l) List of refseq identifiers
    """
    refseq_ids = [] # to return
    
    # Obtain ensembl ids
    ensembl_ids = final_df.loc[pair_id]["other_genes"].split(";")
    ensembl_ids += final_df.loc[pair_id]["other_transcripts"].split(";")    
    
    # Remove empty ones
    ensembl_ids = [x for x in ensembl_ids if x != ""]
    
    # Remove annotation
    ensembl_ids = [x.replace("(protein_coding)", "") for x in ensembl_ids]
    
    # Transform to refseq ids
    if len(ensembl_ids) > 0: 
        table_file = os.path.join(resources.get_blastdb_path(species), 
                                  resources.IDS_TABEL)
        with open(table_file, "r") as op: 
            lines = op.readlines()
        refseq_ids += [l.split("\t")[0] for l in lines if any([x for x in ensembl_ids if x in l])]

    # Obtain refseq_ids
    refseq_ids += final_df.loc[pair_id]["other_genes_rpred"].split(";")
    refseq_ids += final_df.loc[pair_id]["other_transcripts_rpred"].split(";")    

    # Remove empty ones
    refseq_ids = [x for x in refseq_ids if x != ""]
    
    return refseq_ids
  
###############################################################################

def get_offtarget_sequence(refseq_id, species): 
    """
    This function returns the complete sequence of an off-target location. 
    Args: 
        refseq_id [in] (str) Refseq identifier
        species [in] (str)   Organism
        seq [out] (str)      Complete sequence of the refseq identifier

    """
    blast_file = os.path.join(resources.get_blastdb_path(species), 
                              resources.BLAST_DB_NAMES[species])
    
    with open(blast_file, "r") as op: 
        cseq = op.read() # complete sequence
        cseq = cseq.split(">") # list
        
    # add "." to ensure complete match
    string = [x for x in cseq if refseq_id+"." in x][0]
    string = string.split("\n") # list
    string = string[1:] # remove header
    string = "".join(string).upper()
        
    return string

###############################################################################

def get_mismatch_indices(seq, f1, r1, forward, reverse): 
    """
    This function returns the mismatch indices between the off-target 
    alignment and the primers. 
    Args: 
        seq [in] (str)     Off-target DNA seq
        f1, r1 [in] (int)  Start indices for forward and reverse primer
        forward [in] (str) Forward sequence
        reverse [in] (str) Reverse sequence
        mm_i [out] (l)     List of indices where mismatches are found
    """
    mm_i = [] # to return
    
    for i in range(0, len(forward)): 
        # there is a mismatch
        if forward[i] != seq[i + f1]: 
            mm_i.append(i + f1)
            
    for i in range(0, len(reverse)): 
        # there is a mismatch
        if reverse_complement(reverse)[i] != seq[i + r1]: 
            mm_i.append(i + r1)        
        
    return mm_i

###############################################################################
#                MAIN FUNCTION FOR ON TARGET HIGHLIGHTING                     #
###############################################################################

def highlight_ontarget(pair_id, final_df, species, release): 
    """
    MAIN FUNCTION: highlights the ON target alignment of the primers. 
    Args: 
        pair_id [in] (str)  Primer pair identifier (ex "Pair1")
        final_df [in] (str) Final design DF returned by exon surfer
        species [in] (str)  Organism
        release [in] (int)  Ensembl release
        string [out] (str)  String to write to the html file
    """
    # obtain exon junction information
    junction = final_df.loc[pair_id]["junction"]
    masked_chr = resources.MASKED_SEQS(species)
    data = ensembl.create_ensembl_data(release, 
                                       species.replace("_masked", ""))
    # obtain sequence and indices
    nm_dna, ji = get_junction_seqs(junction, masked_chr, data)
    f1, f2, r1, r2 = get_primers_i(nm_dna, final_df.loc[pair_id]["forward"], 
                                   final_df.loc[pair_id]["reverse"])
    # format indices
    indices = format_indices(ji, f1, f2, r1, r2)
    
    # marked string
    string = create_par_string(nm_dna, indices)    
    
    return HEADER_LINE.format(pair_id, junction) + "<br>" + string

###############################################################################
#                MAIN FUNCTION FOR off TARGET HIGHLIGHTING                     #
###############################################################################
        
def highlight_offtarget(pair_id, final_df, species): 
    """
    MAIN FUNCTION: highlights the OFF target alignment of the primers. 
    Args: 
        pair_id [in] (str)  Primer pair identifier (ex "Pair1")
        final_df [in] (str) Final design DF returned by exon surfer
        species [in] (str)  Organism
        strings [out] (l)   List of strings to write to the html file
    """
    strings = [] # list to return
    
    # Obtain all off-target ids
    refseq_ids = obtain_offtarget_list(pair_id, final_df, species)
    
    for item in refseq_ids: 
        seq = get_offtarget_sequence(item, species) # refseq sequence
        # get primer match indices
        f1, f2, r1, r2 = get_primers_i(seq, 
                                       final_df.loc[pair_id]["forward"], 
                                       final_df.loc[pair_id]["reverse"], 
                                       e = 3)
        
        mm_i = get_mismatch_indices(seq, f1, r1, 
                                    final_df.loc[pair_id]["forward"], 
                                    final_df.loc[pair_id]["reverse"])
        
    
        string = create_offt_string(seq, f1, f2, r1, r2, mm_i)
        
        fullstr = HEADER_LINE.format(pair_id, item) + "<br>" + string 
        
        strings.append(fullstr)
    
    return strings
        
