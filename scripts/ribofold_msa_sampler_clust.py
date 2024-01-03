"""
The sequence (de)serialization part is modified from Sailer's Phylopandas package.
Users can call `python3 rhofold_msa_sampler_clust.py -i MSA_PATH -o OUT_DIR -n NUM_CLUST`
liang.hong@link.cuhk.edu.hk
"""

import argparse
import os
import sys
import torch
from sklearn.cluster import KMeans
import numpy as np
import pandas as pd
import math
import rhofold.model.rna_fm as rna_esm

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord

# Load ESM model
model, alphabet = rna_esm.pretrained.esm1b_rna_t12()
batch_converter = alphabet.get_batch_converter()
model.eval()

def _read(
    filename,
    use_uids=False,
    **kwargs):
    """Use BioPython's sequence parsing module to convert any file format to
    a Pandas DataFrame.

    The resulting DataFrame has the following columns:
        - name
        - id
        - description
        - sequence
    """

    # Prepare DataFrame fields.
    data = {
        'id': [],
        'sequence': [],
        'description': [],
        'label': []
    }
    if use_uids:
        data['uid'] = []

    # Parse Fasta file.
    for i, s in enumerate(SeqIO.parse(filename, format='fasta', **kwargs)):
        data['id'].append(s.id)
        data['sequence'].append(str(s.seq))
        data['description'].append(s.description)
        data['label'].append(s.name)

        if use_uids:
            data['uid'].append(get_random_id(10))

    # Port to DataFrame.
    return pd.DataFrame(data)

def pandas_df_to_biopython_seqrecord(
    df,
    id_col='id',
    sequence_col='sequence',
    extra_data=None,
    mtype=None,
    ):
    """Convert pandas dataframe to biopython seqrecord for easy writing.

    Parameters
    ----------
    df : Dataframe
        Pandas dataframe to convert

    id_col : str
        column in dataframe to use as sequence label

    sequence_col : str
        column in dataframe to use as sequence data

    extra_data : list
        extra columns to use in sequence description line

    mtype : str
        molecule type that should be set in annotation

    Returns
    -------
    seq_records :
        List of biopython seqrecords.
    """
    seq_records = []

    for i, row in df.iterrows():
        # Tries getting sequence data. If a TypeError at the seqrecord
        # creation is thrown, it is assumed that this row does not contain
        # sequence data and therefore the row is ignored.
        try:
            # Get sequence
            seq = Seq(row[sequence_col])

            # Get id
            id = row[id_col]

            # Build a description
            description = ""
            if extra_data is not None:
                description = " ".join([row[key] for key in extra_data])

            # Build a record
            annotations = {"molecule_type": mtype} if mtype else None
            record = SeqRecord(
                seq=seq,
                id=id,
                description=description,
                annotations=annotations,
            )
            seq_records.append(record)
        except TypeError:
            pass

    return seq_records

def _write(
    data,
    filename=None,
    id_col='id',
    sequence_col='sequence',
    extra_data=['description', 'label'],
    mtype=None,
    **kwargs):
    """General write function. Write phylopanda data to biopython format.

    Parameters
    ----------
    filename : str
        File to write string to. If no filename is given, a string
        will be returned.

    sequence_col : str (default='sequence')
        Sequence column name in DataFrame.

    id_col : str (default='id')
        ID column name in DataFrame

    id_only : bool (default=False)
        If True, use only the ID column to label sequences in fasta.

    mtype : str
        molecule type that should be set in annotation
    """

    # Build a list of records from a pandas DataFrame
    seq_records = pandas_df_to_biopython_seqrecord(
        data,
        id_col=id_col,
        sequence_col=sequence_col,
        extra_data=extra_data,
        mtype=mtype,
    )

    SeqIO.write(seq_records, filename, format='fasta-2line', **kwargs)

def batch_sequences(sequences, max_tokens):
    """Yield batches of sequences, each with no more than max_tokens tokens."""
    batch = []
    n_tokens = 0
    idx = 0
    for sequence in sequences:
        n_tokens_in_sequence = len(sequence)
        if n_tokens + n_tokens_in_sequence > max_tokens and batch:
            yield batch
            batch = []
            n_tokens = 0
        batch.append((str(idx), sequence))
        idx += 1
        n_tokens += n_tokens_in_sequence
    if batch:
        yield batch

def compute_embeddings(sequences):
    # Initialize list to store embeddings
    embeddings = []
    
    # Compute embeddings for each sequence in the alignment
    for batch in batch_sequences(sequences, max_tokens=200000):
        batch_labels, batch_strs, data = batch_converter(batch)
        with torch.no_grad():
            results = model(data, repr_layers=[12])
            token_embeddings = results['representations'][12]
            embeddings.extend([embed[0].numpy() for embed in token_embeddings])
    
    return np.array(embeddings)

def cluster_sequences(embeddings, n_clusters):
    # Perform clustering
    kmeans = KMeans(n_clusters=n_clusters)
    labels = kmeans.fit_predict(embeddings)
    
    return labels

def generate_clusters(input_file, output_dir, n_clusters=None):
    # Read the input file into a DataFrame
    df = _read(input_file)

    # Extract sequences from the DataFrame
    sequences = df['sequence'].tolist()

    # If the number of sequences is less than 256, write the file directly to the output directory
    if len(sequences) < 256:
        output_file = os.path.join(output_dir, os.path.basename(input_file))
        _write(df, output_file)
        return

    # Compute embeddings for each sequence
    embeddings = compute_embeddings(sequences)

    # If n_clusters is not provided, set it to num_seq / 64
    if n_clusters is None:
        n_clusters = max(1, math.ceil(len(sequences) / 64.0))

    # Perform clustering on the embeddings
    labels = cluster_sequences(embeddings, n_clusters)

    # Create output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)

    # Generate new MSAs for each cluster
    for i in range(n_clusters):
        cluster_df = df[labels == i]
        output_file = os.path.join(output_dir, f"cluster_{i+1:03}.msa")
        _write(cluster_df, output_file)

def main():
    parser = argparse.ArgumentParser(description='Generate clustered MSAs from a single MSA file.')
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input MSA file.')
    parser.add_argument('-o', '--output_dir', required=True, help='Path to the output directory.')
    parser.add_argument('-n', '--n_clusters', type=int, default=None, help='Number of clusters to generate. If not provided, it will be set to num_seq / 64.')
    
    args = parser.parse_args()

    generate_clusters(args.input_file, args.output_dir, args.n_clusters)

if __name__ == "__main__":
    main()

