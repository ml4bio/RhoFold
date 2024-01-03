"""
The sequence (de)serialization part is modified from Sailer's Phylopandas package.
Users can call `python3 rhofold_msa_sampler_rnd.py -i MSA_PATH -o OUT_DIR -n NUM_SAMPLE`
liang.hong@link.cuhk.edu.hk
"""

import argparse
import random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os
import shutil

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

def sample_sequences(input_file, output_directory, sample_times):
    # Constants
    SAMPLE_SIZE = 128
    MIN_SEQUENCES = 256

    # Load the MSA file
    alignment_df = _read(input_file) 

    # If the alignment has fewer than the minimum number of sequences, just copy the input file to the output directory
    if len(alignment_df) < MIN_SEQUENCES:
        shutil.copy(input_file, os.path.join(output_directory, "original.fasta"))
        return

    # Create output directory if it does not exist
    os.makedirs(output_directory, exist_ok=True)

    # Randomly sample sequences and save them
    for i in range(sample_times):
        sampled_df = alignment_df.sample(n=SAMPLE_SIZE)

        # Save the sampled alignment
        output_file = os.path.join(output_directory, f"sample_{i+1:03d}.afa")
        _write(sampled_df, output_file)

def main():
    parser = argparse.ArgumentParser(description='Sample sequences from an MSA file.')
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input MSA file.')
    parser.add_argument('-o', '--output_dir', required=True, help='Path to the output directory.')
    parser.add_argument('-t', '--sample_times', type=int, default=200, help='Number of times to sample.')
    
    args = parser.parse_args()

    sample_sequences(args.input_file, args.output_dir, args.sample_times)

if __name__ == "__main__":
    main()

