import subprocess

import pandas as pd
from Bio import AlignIO
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from utils import fasta_to_csv, csv_to_fasta

import pymsaviz
from smolagents import tool

import subprocess


@tool
def align_and_visualize_fasta(data: pd.DataFrame) -> Path:
    """
    Aligns sequences from a FASTA file using Clustal Omega and visualizes the alignment.

    Args:
        data: A pandas DataFrame containing the dataset to analyze. The DataFrame
            should contain at least two numerical columns for correlation analysis.

    Returns:
        path: The path of the alignment image file.
    """

    fasta_path = Path(Path(__file__).parent, "seqs.fasta")
    align_path = Path(Path(__file__).parent, "seqs.aln.fasta")
    png_path = Path(Path(__file__).parent, "seqs.aln.png")

    csv_to_fasta(data, fasta_path)

    # Step 1: Align using Clustal Omega
    clustalo_cmd = [
        "clustalo",
        "-i", fasta_path,
        "-o", align_path,
        "--force",          # Overwrite existing output
        "--outfmt=fasta"    # FASTA format
    ]

    try:
        subprocess.run(clustalo_cmd, check=True)
        print(f"Alignment written to {align_path}")
    except subprocess.CalledProcessError as e:
        print("Error during Clustal Omega execution:", e)
        return

    from pymsaviz import MsaViz
    mv = MsaViz(align_path)
    mv.savefig(png_path)

    return png_path

    # import base64
    # with open(png_path, "rb") as image_file:
    #     return base64.b64encode(image_file.read())