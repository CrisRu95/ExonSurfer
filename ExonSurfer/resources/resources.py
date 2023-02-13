#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import os
import uuid
import pandas as pd
from subprocess import call, DEVNULL

# Constants
CDNA_LINKS = {
    "homo_sapiens": "https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz", 
    "mus_musculus": "https://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz", 
    "rattus_norvegicus": "https://ftp.ensembl.org/pub/release-108/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz"
    }


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
    command = ["makeblastdb",
               "-dbtype", "nucl",
               "-in", download_path,
               "-out", download_path]
    
    call(command, stderr = DEVNULL, stdout = DEVNULL)

###############################################################################
    
def create_index_table(download_path, cdna_file, realease, species):
    """
    This function creates a table with the transcript id and the gene symbol.
    """
    # import necessary module (only inside this function)
    from pyensembl import EnsemblRelease
    data = EnsemblRelease(realease, species)
    
    with open(cdna_file, "w") as table_out: # output file
        fasta_open = open(download_path, "r")    
        line = fasta_open.readline()
        while line != "":     
        
            if line.startswith(">"):
                lElement = line[1:].split(" ") # line list, remove ">"
                t_id = lElement[0] # transcript id (with version info)
                try:
                    t_obj = data.transcript_by_id(t_id.split(".")[0])
                    symbol = t_obj.gene_name
                except ValueError: # no gene information in EnsemblRelease
                    try:
                        symbol = [x for x in lElement if "gene_symbol" in x][0].split(":")[1]
                    except: # at least keep gene ensemble id
                        symbol = [x for x in lElement if "gene:" in x][0].split(":")[1]

                table_out.write(f"{t_id}\t{symbol}\n")
            line = fasta_open.readline()
            
    fasta_open.close()
    
    return cdna_file

###############################################################################
                        
def download_ensembl_cdna(release, species):
    """
    This function downloads the ensembl cdna fasta file.
    """
    from datacache import download
    print("Downloading ensembl cdna fasta file...")
    download_path = os.path.join(str(get_db_path(species)), CDNA_LINKS[species].split("/")[-1].replace(".gz", ""))
    download._download_and_decompress_if_necessary(full_path = download_path,
                                                   download_url = CDNA_LINKS[species],
                                                   timeout = 3600)
    
    print("Download complete.")
    print("Creating regular index table...")
    
    create_index_table(download_path, get_cdna_file(species), release, species)
    print("Index table created.")
    print("Creating blast database...")
    make_blast_db(download_path)

###############################################################################
    
def get_cdna_file(species):
    """
    This function returns the path to the ensembl cdna table file.
    """
    return str(os.path.join(str(get_db_path(species)),"CDNA_GENE.txt"))

###############################################################################
    
def get_filtered_cdna_file(species):
    """
    This function returns the path to the ensembl cdna table file.
    """
    return str(os.path.join(str(get_db_path(species)),"CDNA_GENE_FILT.txt"))

###############################################################################
    
def get_db_path(species):
    """
    This function returns the path to the blastn database.
    """
    db_path = os.path.join(str(_get_path_data()),"db/", species)
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
    
def get_maskedseq_path(species): 
    db_path = os.path.join(str(_get_path_data()),species)
    #check if db folder exists, if not create it
    if not os.path.exists(db_path):
        os.makedirs(db_path)
    return db_path

###############################################################################
# Function to download the Human Masked Sequence
def download_maskedseq(species):
    import datacache
    
    for chr in list(range(1,24)) + ["X","MT"]:
        url = f"https://sandbox.zenodo.org/record/1136239/files/Chr_{chr}_masked001.fa?download=1"
        download_path = os.path.join(str(get_maskedseq_path(species)),f"Chr_{chr}_masked001.fa")
        print("Downloading masked sequence for chromosome",chr)
        datacache.download._download_and_decompress_if_necessary(
                        full_path=download_path,
                        download_url=url,
                        timeout = 3600)

###############################################################################
        
## Function to create a temporary file to save the primers fasta sequence with uuid
def create_temp_fasta():
    """This function creates a temporary file name inside a temp folder."""
    # Create a temporary file
    temp_file = os.path.join(str(_get_path_data()),"temp",str(uuid.uuid4())+".fa")
    # check if temp folder exists, if not create it
    if not os.path.exists(os.path.join(str(_get_path_data()),"temp")):
        os.makedirs(os.path.join(str(_get_path_data()),"temp"))
        
    return temp_file

###############################################################################

def fillin_temp_fasta(df, temp_file): 
    """
    This function writes a temporary fasta file with primer sequences for BLAST. 
    Args: 
        df [in|out] (pd.df)  Dataframe with designed primers
        temp_file [in] (str) Temporary fasta file
    """
    with open(temp_file, "w") as t_open: 
        for p in list(df.index): 
               
            string1 = ">{}_5\n{}\n".format(p, df.loc[p]["forward"]) 
            string2 = ">{}_3\n{}\n".format(p, df.loc[p]["reverse"]) 
            
            t_open.write(string1 + string2)
            
###############################################################################
    
def MASKED_SEQS(species):
    mask_path = os.path.join(str(get_maskedseq_path(species)),"Chr_1_masked001.fa")
    if not os.path.exists(mask_path):
        download_maskedseq(species)
    return get_maskedseq_path(species)+"/Chr_{}_masked001.fa"

###############################################################################
    
def BLAST_DB(release, species):
    """
    This function returns the blast DB path AND constructs it if necessary. 
    Args: 
        blast_db_path [out] (str) Path to the BLAST db
    """
    blast_db_path = os.path.join(get_db_path(species) + "/"+CDNA_LINKS[species].split("/")[-1].replace(".gz", ""))

    # Check if Homo_sapiens.GRCh38.cdna.all.fa file exists, if not download it
    if not os.path.exists(blast_db_path):
        print("this path does not exist: {}".format(blast_db_path))
        download_ensembl_cdna(release, species)

    return blast_db_path
    
###############################################################################
    
def FILT_BLAST_DB(release, species): 
    """
    This function returns the path to the blast DB composed ONLY by protein
    coding transcripts, AND constructs it if necessary
    """
    # download regular DB if it does not exist
    blast_db_path = BLAST_DB(release, species)
    
    # Open filtered DB file
    filt_db_path = blast_db_path.replace(".cdna.", ".filtered_cdna.")
    
    # If filtered DB does not exist
    if not os.path.exists(filt_db_path): 
        
        to_write = open(filt_db_path, "a+")
        
        # import module
        from pyensembl import EnsemblRelease
        data = EnsemblRelease(release, species)
        
        all_cdna = open(blast_db_path, "r").read()
        all_cdna = all_cdna.split(">")
        all_cdna = [x for x in all_cdna if len(x) > 0]
        
        for item in all_cdna: 
            transcript_id = item.split(" ")[0].split(".")[0] # remove version info
            try: # if transcript is protein coding
                if data.transcript_by_id(transcript_id).biotype == "protein_coding": 
                    to_write.write(">" + item)
            except: 
                pass # transcript is not in the DB because it is not prot cod
        print("Created filtered PCOD fasta")            
        create_index_table(filt_db_path, get_filtered_cdna_file(species), release)
        print("Created filtered index table")   
        make_blast_db(filt_db_path)     
        print("created filtered blast db")   
        
    return filt_db_path
