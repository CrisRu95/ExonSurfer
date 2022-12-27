#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import os
import pandas as pd
from subprocess import call, DEVNULL


###############################################################################
#                   resources module FUNCTION DEFINITION SITE                 #
###############################################################################

def _get_path_data():
    """Return the path to the data directory."""
    return str(os.path.dirname(__file__))

###############################################################################
    
def get_blastn_path():
    """
    This function returns the path to the blastn executable.
    """
    blast_path = os.path.join(str(_get_path_data()),
                              "ncbi-blast-2.13.0+/bin/blastn")
    return blast_path

###############################################################################
    
def make_blast_db(download_path):
    """
    This function makes a blast database from the ensembl cdna fasta file.
        - input: download_path Path to the cDna fasta file
        - output: None
    """
    blastdb_path = "makeblastdb"
    command = [blastdb_path,
               "-dbtype", "nucl",
               "-in", download_path,
               "-out", download_path]
    
    call(command, stderr = DEVNULL, stdout = DEVNULL)

###############################################################################
    
def create_index_table(download_path, realease = "108"):
    """
    This function creates a table with the transcript id and the gene symbol.
    """
    # import necessary module (only inside this function)
    from pyensembl import EnsemblRelease
    data = EnsemblRelease(realease)
    
    # Save CDNA_GENE.txt in the same directory as the cDNA fasta file
    cdna_file = os.path.join(os.path.dirname(download_path), "CDNA_GENE25.txt")
    
    with open(cdna_file, "w") as table_out: # output file
        fasta_open = open(download_path, "r")    
        line = fasta_open.readline()
        while line != "":     
        
            if line.startswith(">"):
                lElement = line.replace(">","").split(" ") # line list
                t_id = lElement[0] # transcript id (with version info)
                try:
                    t_obj = data.transcript_by_id(t_id.split(".")[0])
                    symbol = t_obj.gene_name
                except ValueError: # no gene information in EnsemblRelease
                    try:
                        symbol = [x for x in lElement if "gene_symbol" in x][0].split(":")[1]
                    except:
                        symbol = [x for x in lElement if "gene:" in x][0].split(":")[1]

                table_out.write(f"{t_id}\t{symbol}\n")
            line = fasta_open.readline()
            
    fasta_open.close()

###############################################################################
                        
def download_ensembl_cdna():
    """
    This function downloads the ensembl cdna fasta file.
    """
    import datacache
    
    url = "https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    print("Downloading ensembl cdna fasta file...")
    download_path = os.path.join(str(get_db_path()),"Homo_sapiens.GRCh38.cdna.all.fa")
    datacache.download._download_and_decompress_if_necessary(full_path = download_path,
                                                             download_url = url,
                                                             timeout = 3600)
    
    print("Download complete.")
    print("Creating index table...")
    create_index_table(download_path)
    print("Index table created.")
    print("Creating blast database...")
    make_blast_db(download_path)
    print("Blast database created.")

###############################################################################
    
def get_cdna_file():
    """
    This function returns the path to the ensembl cdna table file.
    """
    return str(os.path.join(str(get_db_path()),"CDNA_GENE.txt"))

###############################################################################
    
def get_db_path():
    """
    This function returns the path to the blastn database.
    """
    db_path = os.path.join(str(_get_path_data()),"db/")
    # check if db folder exists
    if not os.path.exists(db_path):
        os.makedirs(db_path) # if not, create it
    return db_path

###############################################################################
    
def load_mapping():
    """
    Return a dataframe about Gene/Gene synthetic lehal

    Contains the following fields:
        col1        gene symbol
        col2        gene symbol
    """
    import pkg_resources
    
    stream = pkg_resources.resource_stream(__name__, 'mapping.csv')
    df = pd.read_csv(stream, sep = ",", header = None, index_col = 0, 
                     names = ["","tx_id","name"])
    
    return df
    
###############################################################################    
    
def get_maskedseq_path(): 
    db_path = os.path.join(str(_get_path_data()),"Homo_sapiens")
    #check if db folder exists, if not create it
    if not os.path.exists(db_path):
        os.makedirs(db_path)
    return db_path

###############################################################################
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

###############################################################################
        
## Function to create a temporary file to save the primers fasta sequence with uuid
def create_temp_fasta():
    import uuid

    # Create a temporary file
    temp_file = os.path.join(str(_get_path_data()),"temp",str(uuid.uuid4())+".fa")
    #check if temp folder exists, if not create it
    if not os.path.exists(os.path.join(str(_get_path_data()),"temp")):
        os.makedirs(os.path.join(str(_get_path_data()),"temp"))
    return temp_file

# Temporary constants for db locations
###############################################################################
    
def MASKED_SEQS():
    mask_path = os.path.join(str(get_maskedseq_path()),f"Chr_1_masked001.fa")
    if not os.path.exists(mask_path):
        download_maskedseq()
    return get_maskedseq_path()+"/Chr_{}_masked001.fa"

###############################################################################
    
def BLAST_DB():
    blast_db_path = os.path.join(get_db_path()+"Homo_sapiens.GRCh38.cdna.all.fa")

    ## Check if the Homo_sapiens.GRCh38.cdna.all.fa file exists, if not download it
    if not os.path.exists(blast_db_path):
        download_ensembl_cdna()

    return blast_db_path
    
