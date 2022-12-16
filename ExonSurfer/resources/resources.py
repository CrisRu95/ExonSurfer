#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import os
import pandas as pd
from ExonSurfer.blast import blast

def _get_path_data():
    """Return the path to the data directory."""
    return str(os.path.dirname(__file__))


def get_blastn_path():
    """
    This function returns the path to the blastn executable.
    """
    blast_path = os.path.join(str(_get_path_data()),"ncbi-blast-2.13.0+/bin/blastn")
    return blast_path


def download_ensembl_cdna():
    """
    This function downloads the ensembl cdna fasta file.
    """
    import datacache
    
    url = "https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    print("Downloading ensembl cdna fasta file...")
    download_path = os.path.join(str(get_db_path()),"Homo_sapiens.GRCh38.cdna.all.fa")
    datacache.download._download_and_decompress_if_necessary(
                    full_path=download_path,
                    download_url=url,
                    timeout=3600)

    blast.make_blast_db(download_path)


def get_db_path():
    """
    This function returns the path to the blastn database.
    """
    db_path = os.path.join(str(_get_path_data()),"db/")
    #check if db folder exists, if not create it
    if not os.path.exists(db_path):
        os.makedirs(db_path)
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
    db_path = os.path.join(str(_get_path_data()),"Homo_sapiens")
    #check if db folder exists, if not create it
    if not os.path.exists(db_path):
        os.makedirs(db_path)
    return db_path


def test():
    return _get_path_data()


# Function to download the Human Masked Sequence
def download_maskedseq():
    import datacache
    
    for chr in list(range(1,24)) + ["X","MT"]:
        url = f"https://sandbox.zenodo.org/record/1136239/files/Chr_{chr}_masked001.fa?download=1"
        download_path = os.path.join(str(get_maskedseq_path()),f"Chr_{chr}_masked001.fa")
        print("Downloading masked sequence for chromosome",chr)
        datacache.download._download_and_decompress_if_necessary(
                        full_path=download_path,
                        download_url=url,
                        timeout=3600)

## Function to create a temporary file to save the primers fasta sequence with uuid
def create_temp_fasta():
    import uuid
    import os

    # Create a temporary file
    temp_file = os.path.join(str(_get_path_data()),"temp",str(uuid.uuid4())+".fa")
    #check if temp folder exists, if not create it
    if not os.path.exists(os.path.join(str(_get_path_data()),"temp")):
        os.makedirs(os.path.join(str(_get_path_data()),"temp"))
    return temp_file

# Temporary constants for db locations

def MASKED_SEQS():
    mask_path = os.path.join(str(get_maskedseq_path()),f"Chr_1_masked001.fa")
    if not os.path.exists(mask_path):
        download_maskedseq()
    return get_maskedseq_path()+"/Chr_{}_masked001.fa"

def BLAST_DB():
    blast_db_path = os.path.join(get_db_path()+"Homo_sapiens.GRCh38.cdna.all.fa")

    ## Check if the Homo_sapiens.GRCh38.cdna.all.fa file exists, if not download it
    if not os.path.exists(blast_db_path):
        download_ensembl_cdna()

    return blast_db_path
    
