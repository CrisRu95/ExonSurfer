#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 11:09:05 2023

@author: q0791lp
"""
# imported modules
import plotly.graph_objects as go
import plotly.offline as opy

# own modules
from ExonSurfer.ensembl import ensembl

# configuration
config = {
    'toImageButtonOptions': {
        'format': 'svg', # one of png, svg, jpeg, webp
    }}

# colors
COL_UNDETECTED = "#999999" #grey
COL_DETECTED = "#4daf4a" # green
COL_EXON = "#E23F3F" # red
COLS = ['#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', 
        '#a65628','#f781bf', '#999999', '#1b9e77', '#d95f02', '#7570b3', 
        '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666', '#b3e2cd', 
        '#fdb462'] * 200

###############################################################################
#                  ensemble module FUNCTION DEFINITION SITE                   #
###############################################################################

def get_transcripts_exons_dict(gene, exclude_noncoding = True):
    """
    This function takes a gene object and returns a dictionary of transcript
    objects, with transcript ID as keys, and exon objects as values.
    Args:
        gene [in] (gene object)       Gene object
        exclude_noncoding [in] (bool) False if all transcripts, True to exclude 
                                      non-coding
        dT [out] (dict)               Dict of transcript objects, with trans ID 
                                      as keys, and exon objects as values
    """
    dT = {}
    dE = {}
    
    # get list of transcripts to iterate
    all_transcripts = ensembl.get_transcript_from_gene(gene)
    
    if exclude_noncoding == True:
        tcripts = ensembl.get_coding_transcript(all_transcripts)
    else: 
        tcripts = all_transcripts

    for tcript in tcripts:
        dT[tcript.id] = ensembl.get_exons_from_transcript(tcript)

        for exon in tcript.exons:
            dE[exon.id] = (exon.start, exon.end)

    return dT, dE

###############################################################################

def get_exon_transcript_information(species, symbol, transcript, release = 108):
    """
    Function that obtain the information of the transcripts and exons positions 
    of a gene
    Args:
        species [in] (str)  Species name
        symbol  [in] (str)  Gene symbol
        release  [in] (int) Ensembl release
        dt [out] (dict)     Dict with the exons of each transcript
        de  [out] (dict)    Dict with exons positions
    """
    data = ensembl.create_ensembl_data(release, species.replace("_masked", ""))
    gene_obj = ensembl.get_gene_by_symbol(symbol, data)
        
    dT, dE = get_transcripts_exons_dict(gene_obj)

    # If transcript are provided, only return the information of that transcript
    if transcript != "ALL":
        dT = {transcript: dT[transcript]}

    return dT, dE

###############################################################################

def transform_primers_pos(for_pos, rev_pos, de): 
    """
    This function transforms the positions returned in the design DF to a dict
    with genomic positions. 
    Args: 
        for_pos [in] (l)   List of 1|2 tuples. Each tuple contains: (ENSE, plen)
        rev_pos [in] (l)   List of 1|2 tuples. Each tuple contains: (ENSE, plen)
        de [in] (dict)     Dict with exons positions
        primd [out] (dict) Dict with 2 keys and genomic positions
    """
    primd = {} # to return
    
    # only one exon design
    if len(for_pos) == 1 and len(rev_pos) == 1 and for_pos[0][0] == rev_pos[0][0]: 
        exon = for_pos[0][0]
        midpoint1 = int(de[exon][0] + (de[exon][1] - de[exon][0])*0.25)
        midpoint2 = int(de[exon][0] + (de[exon][1] - de[exon][0])*0.75)
        primd["forward1"] = (midpoint1 - 10, midpoint1 + 10)
        primd["forward2"] = (midpoint1 - 10, midpoint1 + 10) # same pos
        primd["reverse1"] = (midpoint2 - 10, midpoint2 + 10)
        primd["reverse2"] = (midpoint2 - 10, midpoint2 + 10) # same pos
        
    # typical multi exon design
    else: 
        if len(for_pos) == 1: 
            exon = for_pos[0][0]
            midpoint = int(de[exon][0] + (de[exon][1] - de[exon][0])/2)
            primd["forward1"] = (midpoint - 10, midpoint + 10)
            primd["forward2"] = (midpoint - 10, midpoint + 10) # same pos
        else: 
            exon1 = for_pos[0][0]
            exon2 = for_pos[1][0]
            primd["forward1"] = (de[exon1][1] - 10, de[exon1][1])
            primd["forward2"] = (de[exon2][0], de[exon2][0] + 10)
    
        # reverse position
        if len(rev_pos) == 1: 
            exon = rev_pos[0][0]
            midpoint = int(de[exon][0] + (de[exon][1] - de[exon][0])/2)
            primd["reverse1"] = (midpoint - 10, midpoint + 10)
            primd["reverse2"] = (midpoint - 10, midpoint + 10) # same pos
        else: 
            exon1 = rev_pos[0][0]
            exon2 = rev_pos[1][0]
            primd["reverse1"] = (de[exon1][1] - 10, de[exon1][1])
            primd["reverse2"] = (de[exon2][0], de[exon2][0] + 10)
    
    return primd

###############################################################################

def get_transcript_color(transd, detected): 
    
    transcol = {}
    
    for key in transd: 
        if key in detected: 
            transcol[key] = COL_DETECTED
        else: 
            transcol[key] = COL_UNDETECTED
    return transcol

###############################################################################

def plot_transcripts_alone(species, gene, transcripts, release): 
    """
    This function plots the list of transcripts a gene has. 
    """
    transd, exd = get_exon_transcript_information(species, gene, transcripts, 
                                                  release)
    # create colors dict
    colors = dict(zip(transd.keys(), COLS[:len(transd)]))
    
    # define spacing between exon boxes
    box_spacing = 0.5

    # create the figure
    fig = go.Figure()    
    
    # loop over transcripts
    for i, transcript in enumerate(transd):
        # create the transcript line
        fig.add_shape(type = 'line',
                    x0 = min(exd[e][0] for e in transd[transcript]),
                    y0 = i + 0.25,
                    x1 = max(exd[e][1] for e in transd[transcript]),
                    y1 = i + 0.25,
                    line = dict(color = 'black', width = 1))

        # loop over exons in transcript
        for j, exon in enumerate(transd[transcript]):
            # determine x-coordinates of exon box
            x0 = exd[exon][0]
            x1 = exd[exon][1]
            width = x1 - x0
            
            fig.add_shape(type = 'rect',
                        x0 = x0,
                        y0 = i,
                        x1 = x1,
                        y1 = i + 0.5,
                        fillcolor = colors[transcript],
                        #line=dict(color='black'),
                        opacity = 1)
            # add hover with exon_id
            # Adding a trace with a fill, setting opacity to 0
            fig.add_trace(
                go.Scatter(
                    x = [x0 + width/2],
                    y = [i + 0.25],
                    mode = 'markers',
                    marker = dict(size = 0.1,
                                  color = colors[transcript],
                                  opacity = 0),
                    hovertext = exon,
                    hoverinfo = 'text')
            )
    # set x-axis range
    x_range = [min(exd[e][0] for t in transd.values() for e in t) - 1,
               max(exd[e][1] for t in transd.values() for e in t) + len(transd) * (box_spacing + 1)]
    fig.update_xaxes(range = x_range)


    # set layout properties
    fig.update_layout(
        #title='Exon Plot',
        xaxis_title = 'Position',
        #yaxis_title='Transcript',
        showlegend = False,
        #height=500,
        #width=800,
    )

    fig.update_layout(template = "plotly_white")
    fig.update_layout(yaxis = dict(tickmode = 'array',
                                  tickvals = [x + 0.25 for x in list(range(len(transd))) ],
                                  ticktext = ["<b> %s </b>"%x for x in list(transd.keys())],
                                  range = [-0.6, len(transd)-0.4]))

    div = opy.plot(fig, 
                   auto_open = False, 
                   config = config,
                   output_type = 'div')

    return div    
    
###############################################################################  
  
def plot_transcripts_marked(species, gene, transcripts, release, pair_id, final_df): 
    """
    This function plots the list of transcripts a gene has and the primers locations. 
    """
    transd, exd = get_exon_transcript_information(species, gene, transcripts, 
                                                  release)
    for_pos = final_df.loc[pair_id]["for_pos"]
    rev_pos = final_df.loc[pair_id]["rev_pos"]
    
    mex = []
    for ex in exd: 
        if ex in [exon[0] for exon in for_pos] + [exon[0] for exon in rev_pos]: 
            mex.append(ex)

    
    colors = get_transcript_color(transd, final_df.loc[pair_id]["detected"])
    
    # define spacing between exon boxes
    box_spacing = 0.5
    
    # create the figure
    fig = go.Figure()

    # loop over transcripts
    for i, transcript in enumerate(transd):
        # create the transcript line
        fig.add_shape(type = 'line',
                    x0 = min(exd[e][0] for e in transd[transcript]),
                    y0 = i + 0.25,
                    x1 = max(exd[e][1] for e in transd[transcript]),
                    y1 = i + 0.25,
                    line = dict(color = 'black', width = 2))

        # loop over exons in transcript
        for j, exon in enumerate(transd[transcript]):
            # determine x-coordinates of exon box
            x0 = exd[exon][0]
            x1 = exd[exon][1]
            width = x1 - x0
            
            # # exon is targeted by that primer pair
            if exon in mex and transcript in final_df.loc[pair_id]["detected"]: 
                fig.add_shape(type = 'rect',
                            x0 = x0,
                            y0 = i,
                            x1 = x1,
                            y1 = i + 0.5,
                            fillcolor = COL_EXON,
                            line = dict(color = COL_EXON),
                            opacity = 1)  

            else: # exon is not targeted
                fig.add_shape(type = 'rect',
                            x0 = x0,
                            y0 = i,
                            x1 = x1,
                            y1 = i + 0.5,
                            fillcolor = colors[transcript],
                            #line=dict(color='black'),
                            opacity = 1)
            # add hover with exon_id
            # Adding a trace with a fill, setting opacity to 0
            fig.add_trace(
                go.Scatter(
                    x = [x0 + width/2],
                    y = [i + 0.25],
                    mode = 'markers',
                    marker = dict(size = 0.1,
                                  color = colors[transcript],
                                  opacity = 0),
                    hovertext = exon,
                    hoverinfo = 'text')
            )
    # set x-axis range
    x_range = [min(exd[e][0] for t in transd.values() for e in t) - 1,
               max(exd[e][1] for t in transd.values() for e in t) + len(transd) * (box_spacing + 1)]
    fig.update_xaxes(range = x_range)


    # set layout properties
    fig.update_layout(
        #title='Exon Plot',
        xaxis_title = 'Position',
        #yaxis_title='Transcript',
        showlegend = False,
        #height=500,
        #width=800,
    )

    fig.update_layout(template = "plotly_white")
    fig.update_layout(yaxis = dict(tickmode = 'array',
                                  tickvals = [x + 0.25 for x in list(range(len(transd))) ],
                                  ticktext = ["<b> %s </b>"%x for x in list(transd.keys())],
                                  range = [-0.6, len(transd)-0.4]))
    div = opy.plot(fig, 
                   auto_open = False, 
                   config = config,
                   output_type = 'div')

    return div

###############################################################################
        
def plot_transcripts_with_primers(species, gene, transcripts, release, pair_id, 
                                  final_df): 
    """
    This function plots the list of transcripts a gene has and the primers locations. 
    """
    transd, exd = get_exon_transcript_information(species, gene, transcripts, 
                                                  release)
    for_pos = final_df.loc[pair_id]["for_pos"]
    rev_pos = final_df.loc[pair_id]["rev_pos"]
    primd = transform_primers_pos(for_pos, rev_pos, exd)
    
    mex = []
    for ex in exd: 
        if ex in [exon[0] for exon in for_pos] + [exon[0] for exon in rev_pos]: 
            mex.append(ex)

    
    colors = get_transcript_color(transd, final_df.loc[pair_id]["detected"])
    
    # define spacing between exon boxes
    box_spacing = 0.5
    aline = 2
    pline = 4
    
    # define spacing according to number of transcripts: 
    if len(transd) < 8: 
        updist = 0.55
        updist_t = 0.60  
        fsize = 10
    else: 
        updist = 0.62
        updist_t = 0.75  
        fsize = 8       
    
    # create the figure
    fig = go.Figure()

    # loop over transcripts
    for i, transcript in enumerate(transd):
        # create the transcript line
        fig.add_shape(type = 'line',
                    x0 = min(exd[e][0] for e in transd[transcript]),
                    y0 = i + 0.25,
                    x1 = max(exd[e][1] for e in transd[transcript]),
                    y1 = i + 0.25,
                    line = dict(color = 'black', width = 2))

        # loop over exons in transcript
        for j, exon in enumerate(transd[transcript]):
            # determine x-coordinates of exon box
            x0 = exd[exon][0]
            x1 = exd[exon][1]
            width = x1 - x0
            
            # # exon is targeted by that primer pair
            if exon in mex and transcript in final_df.loc[pair_id]["detected"]: 
                fig.add_shape(type = 'rect',
                            x0 = x0,
                            y0 = i,
                            x1 = x1,
                            y1 = i + 0.5,
                            fillcolor = COL_EXON,
                            line = dict(color = COL_EXON),
                            opacity = 1)  
                

                # add amplicon lines
                fig.add_shape(type = 'line',
                            x0 = primd["forward1"][1],
                            y0 = i + updist,
                            x1 = primd["forward2"][0],
                            y1 = i + updist,
                            line = dict(color = 'grey', width = aline, dash='dot'))   
                # add amplicon lines
                fig.add_shape(type = 'line',
                            x0 = primd["forward2"][1],
                            y0 = i + updist,
                            x1 = primd["reverse1"][0],
                            y1 = i + updist,
                            line = dict(color = 'grey', width = aline))   
                # add amplicon lines
                fig.add_shape(type = 'line',
                            x0 = primd["reverse1"][1],
                            y0 = i + updist,
                            x1 = primd["reverse2"][0],
                            y1 = i + updist,
                            line = dict(color = 'grey', width = aline, dash='dot'))   
                
                
                for primer in primd:
                    fig.add_annotation(text = primer[:1].upper(),
                                    x = (primd[primer][0] + primd[primer][1])/2,
                                    y = i + updist_t,
                                    showarrow = False, 
                                    font = dict(size = fsize, color = "black"))
                    fig.add_shape(type='line',
                                x0 = primd[primer][0],
                                y0 = i + updist,
                                x1 = primd[primer][1],
                                y1 = i + updist,
                                line = dict(color = 'black', width = pline))  

            else: # exon is not targeted
                fig.add_shape(type = 'rect',
                            x0 = x0,
                            y0 = i,
                            x1 = x1,
                            y1 = i + 0.5,
                            fillcolor = colors[transcript],
                            #line=dict(color='black'),
                            opacity = 1)
            # add hover with exon_id
            # Adding a trace with a fill, setting opacity to 0
            fig.add_trace(
                go.Scatter(
                    x = [x0 + width/2],
                    y = [i + 0.25],
                    mode = 'markers',
                    marker = dict(size = 0.1,
                                  color = colors[transcript],
                                  opacity = 0),
                    hovertext = exon,
                    hoverinfo = 'text')
            )
    # set x-axis range
    x_range = [min(exd[e][0] for t in transd.values() for e in t) - 1,
               max(exd[e][1] for t in transd.values() for e in t) + len(transd) * (box_spacing + 1)]
    fig.update_xaxes(range = x_range)


    # set layout properties
    fig.update_layout(
        #title='Exon Plot',
        xaxis_title = 'Position',
        #yaxis_title='Transcript',
        showlegend = False,
        #height=500,
        #width=800,
    )

    fig.update_layout(template = "plotly_white")
    fig.update_layout(yaxis = dict(tickmode = 'array',
                                  tickvals = [x + 0.25 for x in list(range(len(transd))) ],
                                  ticktext = ["<b> %s </b>"%x for x in list(transd.keys())],
                                  range = [-0.6, len(transd)-0.4]))

    div = opy.plot(fig, 
                   auto_open = False, 
                   config = config,
                   output_type = 'div')

    return div
    
    