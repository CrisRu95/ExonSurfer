import re
from pyensembl import EnsemblRelease
import multiprocessing

def request_fasta(id_ensembl):
    import requests, sys
 
    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{id_ensembl}?"

    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()
    return r.text

def get_gene_by_symbol(gene_symbol = None, release = 108):
    # release 77 uses human reference genome GRCh38
    """
    This function takes a gene symbol and returns the gene object.
    Args:
        gene_symbol [in] (str): Gene symbol
        release [in] (int): Ensembl release
        gene [out] (gene object): Gene object
    """
    data = EnsemblRelease(release)
    gene = data.genes_by_name(gene_symbol)
    
    return gene

def get_transcript_from_gene(gene = None):
    """
    This function takes a gene object and returns a list of transcript objects.
    Args:
        gene [in] (gene object): Gene object
        transcripts [out] (list): List of transcript objects
    """
    return gene.transcripts

def get_exons_from_transcript(transcripts = None):
    """
    This function takes a transcript object and returns a list of exon objects.
    Args:
        transcripts [in] (transcript object): Transcript object
        exons [out] (list): List of exon objects
        """
    return transcripts.exons

def get_coding_transcript(transcripts):
    """
    This function takes a list of transcript objects and keep only the
    proteing coding transcript.
    Args:
        transcripts [in] (list): List of transcript objects
        coding_transcript [out] (list): List of protein coding transcript objects
    """
    return [x for x in transcripts if x.biotype == "protein_coding"]

def get_coding_sequence(transcripts):
    """
    This function takes a list of transcript objects and returns a list of
    coding sequences.
    Args:
        transcripts [in] (list): List of transcript objects
        coding_sequences [out] (list): List of coding sequences
        """
    return [x.coding_sequence for x in transcripts]

def get_transcripts_dict(gene):
    """
    This function takes a gene object and returns a dictionary of transcript
    objects, with transcript ID as keys, and exon objects as values.
    Args:
        gene [in] (gene object): Gene object
        dTranscripts [out] (dict): Dictionary of transcript objects, with
        transcript ID as keys, and exon objects as values
    """
    dTranscripts = {}

    for transcripts in get_transcript_from_gene(gene):
        dTranscripts[transcripts.id] = get_exons_from_transcript(transcripts)

    return dTranscripts

if __name__ == '__main__':
    print(get_transcripts_dict(get_gene_by_symbol("BRCA1")[0]))