#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import os
import pandas as pd


# Temporary constants for db locations
MASKED_SEQS = "/home/ecrisru/Desktop/refgenome/Chrs_masked_001/Chr_{}_masked001.fa"
BLAST_DB = "/home/ecrisru/Desktop/refgenome/ExonPrimerSurfer/blast/Homo_sapiens.GRCh38.cdna.all.fa"
FASTA_F = "/home/ecrisru/Desktop/refgenome/ExonPrimerSurfer/blast/primers.fa"
BLAST_OUT = "/home/ecrisru/Desktop/refgenome/ExonPrimerSurfer/blast/primers.out"
BLAST_PATH = "blastn" 
DESIGN_OUT = "/home/ecrisru/Desktop/ExonPrimerSurfer/pruebas/OFILE.txt"


def _get_path_data():
    """Return the path to the data directory."""
    return os.path.dirname(__file__)

def get_blastn_path():
    """
    This function returns the path to the blastn executable.
    """
    blast_path = os.path.join(str(_get_path_data()),"ncbi-blast-2.13.0+/bin/blastn")
    return blast_path

def get_db_path():
    """
    This function returns the path to the blastn database.
    """
    db_path = os.path.join(str(_get_path_data()),"db/Homo_sapiens.GRCh38.cdna.all.fa")
    return db_path

def load_mapping():
    import pkg_resources
    """Return a dataframe about Gene/Gene synthetic lehal

    Contains the following fields:
        col1        gene symbol
        col2        gene symbol
    """
    stream = pkg_resources.resource_stream(__name__, 'mapping.csv')
    return pd.read_csv(stream, sep = ",",header=None, index_col=0, names=["","tx_id","name"])
    
    
    
def get_maskedseq_path(): 
    db_path = os.path.join(str(_get_path_data()),"db/Homo_sapiens.GRCh38.cdna.all.fa")
    return db_path    


def test():
    return _get_path_data()