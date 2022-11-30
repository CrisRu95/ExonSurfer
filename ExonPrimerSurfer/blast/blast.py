from ..resources import get_blastn_path, get_db_path
import uuid
import os
import pandas as pd

def writing_fasta_from_list(lPrimer, out, identifier):
    """
    This function takes a list of primers and writes a fasta file with the primers. 
    Args: 
        lPrimer [in] (list): List of primers
        out [in] (path): Path to the output directory
        identifier [in] (str): Identifier for the fasta file
    """
    path_out = f"{out}/{identifier}.fa"
    with open(path_out, "w") as f:
        for primer in lPrimer:
            f.write(f">{primer}")
            f.write("\n")
            f.write(primer)
            f.write("\n")       


def run_blast_list(lPrimer, out, num_threads = 4):
    """
    This function takes a list of primers and runs blastn on them.
    Args:
        lPrimer [in] (list): List of primers
        out [in] (path): Path to the output directory
        num_threads [in] (int): Number of threads to use for blastn
        df [out] (dataframe): Dataframe with the blast results
    """
    import subprocess
    
    blastn_path = get_blastn_path()
    db_path = get_db_path()

    identifier = str(uuid.uuid4())
    writing_fasta_from_list(lPrimer, out, identifier = identifier)
    command_line = [blastn_path,
                    '-task', 'blastn-short',
                    '-query', f'{out}/{identifier}.fa',
                    '-out',f'{out}/{identifier}.out',
                    '-db',db_path,
                    '-outfmt','6',
                   '-num_threads', str(num_threads)]

    #RunBlastDBCommand(command_line)
    subprocess.call(command_line, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    df = pd.read_csv(f'{out}/{identifier}.out', sep = "\t", names=("query id", "subject id","identity", "alignment length", "mismatches", 
                                                              "gap opens", "q. start", "q. end", "s. start"," s. end", "evalue", "bit score"))
    os.remove(f'{out}/{identifier}.fa')
    os.remove(f'{out}/{identifier}.out')
    return df
    

def check_blast_result(df, identity = 90, evalue = 0.5):
    """
    This function takes a dataframe with blast results and filters them based on identity and evalue.
    Args:
        df [in] (dataframe): Dataframe with the blast results
        identity [in] (int): Minimum identity to keep the blast result
        evalue [in] (float): Maximum evalue to keep the blast result
        df [out] (dataframe): Dataframe with the filtered blast results
    """
    df_filter = df.loc[(df.identity >= identity) & (df.evalue <= evalue),"Gene"]
    #print(df.loc[(df.identity >= identity) & (df.evalue <= treshold),:])
    #print(df_filter.unique().tolist())
    
    return len(df_filter.unique().tolist()) == 1

