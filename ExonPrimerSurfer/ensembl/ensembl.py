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

def get_gene_by_symbol(gene_symbol = None):
    # release 77 uses human reference genome GRCh38
    data = EnsemblRelease()
    gene = data.genes_by_name(gene_symbol)
    
    return gene

def get_transcript_from_gene(gene = None):
    return gene.transcripts

def get_exons_from_transcript(transcripts = None):
    return transcripts.exons

def get_coding_transcript(transcripts):
    return [x for x in transcripts if x.biotype == "protein_coding"]

def get_coding_sequence(transcripts):
    return [x.coding_sequence for x in transcripts]

def get_transcripts_dict(gene):
    dTranscripts = {}

    for transcripts in get_transcript_from_gene(gene):
        dTranscripts[transcripts.id] = get_exons_from_transcript(transcripts)

    return dTranscripts

if __name__ == '__main__':
    print(get_transcripts_dict(get_gene_by_symbol("BRCA1")[0]))