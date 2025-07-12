import re
from typing import List, Dict, Any
from Bio import AlignIO, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import draw
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import io
import base64
import tempfile
import os
from io import StringIO

def parse_aligned_sequences(alignment_data: str) -> List[Dict[str, str]]:
    """Parse aligned sequences from string format."""
    sequences = []
    current_name = ""
    current_sequence = ""
    
    for line in alignment_data.strip().split('\n'):
        line = line.strip()
        if line.startswith('>'):
            # Save previous sequence if exists
            if current_name and current_sequence:
                sequences.append({
                    "name": current_name,
                    "sequence": current_sequence
                })
            # Start new sequence
            current_name = line[1:].strip()
            current_sequence = ""
        else:
            # Add to current sequence
            current_sequence += line
    
    # Add the last sequence
    if current_name and current_sequence:
        sequences.append({
            "name": current_name,
            "sequence": current_sequence
        })
    
    # Error handling: if no sequences parsed, return empty list
    if not sequences:
        print("[ERROR] No sequences parsed from input. Check input format.")
    return sequences


def create_phylogenetic_tree(aligned_sequences: List[Dict[str, str]]) -> Dict[str, Any]:
    """Create phylogenetic tree from aligned sequences and return Newick string."""
    # Error handling: check for empty or insufficient sequences
    if not aligned_sequences or len(aligned_sequences) < 2:
        return {
            "text": "Error: At least 2 aligned sequences in valid FASTA format are required for phylogenetic tree construction. Example: >seq1\\nATCG...\\n>seq2\\nATCG...",
            "error": "Insufficient or invalid sequences"
        }
    try:
        # Create temporary FASTA file with aligned sequences
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
            for seq_data in aligned_sequences:
                temp_fasta.write(f">{seq_data['name']}\n{seq_data['sequence']}\n")
            temp_fasta_path = temp_fasta.name
        # Read the alignment
        alignment = AlignIO.read(temp_fasta_path, "fasta")
        # Calculate distance matrix
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        # Construct tree using UPGMA method
        constructor = DistanceTreeConstructor(calculator, 'upgma')
        tree = constructor.build_tree(alignment)
        # Export tree to Newick format
        from Bio import Phylo
        newick_io = io.StringIO()
        Phylo.write(tree, newick_io, "newick")
        newick_str = newick_io.getvalue().strip()
        # Calculate tree statistics
        num_sequences = len(aligned_sequences)
        alignment_length = len(aligned_sequences[0]["sequence"])
        tree_height = tree.root.branch_length if tree.root.branch_length else 0
        num_clades = len(list(tree.get_terminals()))
        result_text = f"""Phylogenetic tree constructed successfully from {num_sequences} aligned sequences.\n\nTree Statistics:\n- Number of sequences: {num_sequences}\n- Alignment length: {alignment_length}\n- Tree height: {tree_height:.4f}\n- Number of clades: {num_clades}\n- Distance method: Identity\n- Tree construction method: UPGMA\n\nThe tree shows the evolutionary relationships between the sequences based on sequence identity."""
        os.unlink(temp_fasta_path)
        return {
            "text": result_text,
            "tree_newick": newick_str,
            "statistics": {
                "num_sequences": num_sequences,
                "alignment_length": alignment_length,
                "tree_height": tree_height,
                "num_clades": num_clades,
                "distance_method": "identity",
                "tree_method": "upgma"
            }
        }
    except Exception as e:
        return {
            "text": f"Error creating phylogenetic tree: {str(e)}",
            "error": str(e)
        }

def run_phylogenetic_tree_raw(aligned_sequences: str):
    """Create phylogenetic tree from aligned sequences."""
    
    print(f"🔧 Phylogenetic tree tool called with aligned_sequences: '{aligned_sequences}'")
    
    if not aligned_sequences or not aligned_sequences.strip():
        return {
            "text": "Error: No aligned sequences provided",
            "error": "No sequences"
        }
    
    # Parse sequences
    try:
        parsed_sequences = parse_aligned_sequences(aligned_sequences)
        print(f"🔧 Parsed {len(parsed_sequences)} sequences: {parsed_sequences}")
    except Exception as e:
        return {
            "text": f"Error parsing sequences: {str(e)}",
            "error": str(e)
        }
    
    # Check if sequences need alignment (different lengths)
    sequence_lengths = [len(seq["sequence"]) for seq in parsed_sequences]
    if len(set(sequence_lengths)) > 1:
        print(f"🔧 Sequences have different lengths: {sequence_lengths}, aligning first...")
        
        # Import and use the alignment tool to align sequences
        import sys
        import os
        sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
        from alignment import run_alignment_tool
        
        # Convert sequences back to FASTA format for alignment
        fasta_input = "\n".join([f">{seq['name']}\n{seq['sequence']}" for seq in parsed_sequences])
        
        # Align the sequences
        alignment_result = run_alignment_tool(fasta_input)
        
        if "alignment" in alignment_result and alignment_result["alignment"]:
            # Use the aligned sequences
            aligned_parsed_sequences = alignment_result["alignment"]
            print(f"🔧 Aligned sequences: {aligned_parsed_sequences}")
        else:
            return {
                "text": f"Error aligning sequences: {alignment_result.get('text', 'Unknown error')}",
                "error": "Alignment failed"
            }
    else:
        print(f"🔧 All sequences have same length: {sequence_lengths[0]}")
        aligned_parsed_sequences = parsed_sequences
    
    # Create phylogenetic tree
    return create_phylogenetic_tree(aligned_parsed_sequences)

from langchain.agents import tool

@tool
def run_phylogenetic_tree(aligned_sequences: str):
    """Create phylogenetic tree visualization from aligned sequences."""
    return run_phylogenetic_tree_raw(aligned_sequences) 