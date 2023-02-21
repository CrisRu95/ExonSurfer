#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# imported modules
import os
import uuid

# Constants
BLAST_DB_LINKS = {
    "homo_sapiens": "https://zenodo.org/record/7638757/files/human.zip?download=1", 
    "homo_sapiens_masked": "https://zenodo.org/record/7638757/files/human.zip?download=1", 
    "mus_musculus": "https://zenodo.org/record/7638757/files/mouse.zip?download=1", 
    "rattus_norvegicus": "https://zenodo.org/record/7638757/files/rat.zip?download=1"
    }

BLAST_DB_NAMES = {
    "homo_sapiens": "human.rna.fna", 
    "homo_sapiens_masked": "human.rna.fna", 
    "mus_musculus": "mouse.rna.fna", 
    "rattus_norvegicus": "rat.rna.fna"
    }

CHROM_LINKS = {
    "homo_sapiens": "https://zenodo.org/record/7638757/files/human_seqs.zip?download=1", 
    "homo_sapiens_masked": "https://zenodo.org/record/7638757/files/human_masked_seqs.zip?download=1", 
    "mus_musculus": "https://zenodo.org/record/7638757/files/mouse_seqs.zip?download=1", 
    "rattus_norvegicus": "https://zenodo.org/record/7638757/files/rat_seqs.zip?download=1"}

FIRST_CHR = {
    "homo_sapiens": "human{}.txt", 
    "homo_sapiens_masked": "human{}_masked.txt", 
    "mus_musculus": "mouse{}.txt", 
    "rattus_norvegicus": "rat{}.txt"    
    }

IDS_TABEL = "table.txt"

###############################################################################
#                   resources module FUNCTION DEFINITION SITE                 #
###############################################################################

def _get_path_data():
    """Return the path to the data directory."""
    custom_path = os.environ.get("EXONSURFER_CACHE_DIR", None)

    if custom_path is not None:
        path_data = custom_path
    else:
        path_data = str(os.path.dirname(__file__))
    return path_data

###############################################################################
    
def get_blastn_path():
    """
    This function returns the path to the blastn executable.
    """
    blast_path = os.path.join(str(_get_path_data()),
                              "ncbi-blast-2.13.0+/bin/blastn")
    return blast_path

###############################################################################
    
def get_blastdb_path(species):
    """
    This function returns the path to the BLAST database.
    """
    if species == "homo_sapiens_masked":
        species = "homo_sapiens"
        
    db_path = os.path.join(str(_get_path_data()),"db/", species)
    # check if db folder exists
    if not os.path.exists(db_path):
        os.makedirs(db_path) # if not, create it
        
    return db_path
    
###############################################################################    
    
def get_maskedseq_path(species): 
    """
    This function returns the path to the masked sequences.
    """
    db_path = os.path.join(str(_get_path_data()), species)
    
    #check if db folder exists, if not create it
    if not os.path.exists(db_path):
        os.makedirs(db_path)
    return db_path

###############################################################################
# Function to download the Human Masked Sequence
def download_maskedseq(species):
    """
    This function downloads and extracts the zenodo records for the sequences.
    Args: 
        species [in] (str) homo_sapiens, mus_musculus or rattus_norvegicus
    """
    # imported modules
    import urllib.request
    import zipfile
    
    url = CHROM_LINKS[species]
    download_path = os.path.join(str(get_maskedseq_path(species)), 
                                 species + ".gz")
    
    print("Starting download for {} sequences.\nThis will take a while...".format(species))
    urllib.request.urlretrieve(url, download_path)
    print("Download completed.")

    # unzip 
    with zipfile.ZipFile(download_path, 'r') as zip_ref:
        zip_ref.extractall(str(get_maskedseq_path(species)))
    
    # remove extraction file
    os.remove(download_path)


###############################################################################
                        
def download_blast_db(species):
    """
    This function downloads the BLAST database from the zenodo records.
    """
    # imported modules
    import urllib.request
    import zipfile
    
    url = BLAST_DB_LINKS[species]
    download_path = os.path.join(str(get_blastdb_path(species)), species + ".gz")
    
    print("Starting download for BLAST DB. This will take a while...")
    urllib.request.urlretrieve(url, download_path)
    print("Download completed.")    
    
    # unzip 
    with zipfile.ZipFile(download_path, 'r') as zip_ref:
        zip_ref.extractall(str(get_blastdb_path(species)))
    
    # remove extraction file
    os.remove(download_path)

          
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
    """
    Thi function acts as a constant; downloads data if necessary and returns 
    its location. 
    """
    # filename to search for
    mask_path = os.path.join(str(get_maskedseq_path(species)),
                             FIRST_CHR[species].format(1))
    # check if file exists
    if not os.path.exists(mask_path):
        download_maskedseq(species)
        
    return os.path.join(str(get_maskedseq_path(species)), FIRST_CHR[species])

###############################################################################
    
def BLAST_DB(species):
    """
    This function returns the blast DB path AND constructs it if necessary. 
    Args: 
        species [in] (str)        homo_sapiens, mus_musculus or rattus_norvegicus
        blast_db_path [out] (str) Path to the BLAST db
    """
    # filename to search for
    blast_db_path = os.path.join(get_blastdb_path(species), 
                                 BLAST_DB_NAMES[species])

    # check if file exists
    if not os.path.exists(blast_db_path):
        download_blast_db(species)

    return blast_db_path
###############################################################################

def download_all_db():
    """
    This function downloads all the databases.
    Fasta sequences and BLAST databases.
    """
    for species in BLAST_DB_LINKS.keys():
        print("[+] Downloading {} BLAST DB...".format(species), flush=True)
        BLAST_DB(species)
        print("[+] Downloading {} masked sequences...".format(species), flush=True)
        MASKED_SEQS(species)
        print("[+] Done.", flush=True)
