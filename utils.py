import csv

import pandas as pd
from Bio import SeqIO

import csv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def csv_to_fasta(df: pd.DataFrame, fasta_path: str):
    """
    Converts a CSV file with 'sequence_id' and 'sequence' columns to a FASTA file.

    Parameters:
    - csv_path: Path to the input CSV file.
    - fasta_path: Path to the output FASTA file.
    """
    records = []

    for _, row in df.iterrows():
        seq_id = str(row["Sequence ID"])
        sequence = str(row["Sequence"])
        record = SeqRecord(Seq(sequence), id=seq_id, description="")
        records.append(record)

    with open(fasta_path, "w") as fasta_file:
        SeqIO.write(records, fasta_file, "fasta")

    print(f"FASTA file written to {fasta_path}")


def fasta_to_csv(fasta_path: str, csv_path: str):
    """
    Converts a FASTA file to a CSV file with sequence ID and sequence.

    Parameters:
    - fasta_path: Path to the input FASTA file.
    - csv_path: Path to output CSV file.
    """
    with open(fasta_path, "r") as fasta_file, open(csv_path, "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["sequence_id", "sequence"])  # Header row

        for record in SeqIO.parse(fasta_file, "fasta"):
            writer.writerow([record.id, str(record.seq)])

    print(f"CSV file written to {csv_path}")