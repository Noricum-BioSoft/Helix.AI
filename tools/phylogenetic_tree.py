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
    
    return sequences

def create_phylogenetic_tree(aligned_sequences: List[Dict[str, str]]) -> Dict[str, Any]:
    """Create phylogenetic tree from aligned sequences."""
    
    if len(aligned_sequences) < 2:
        return {
            "text": "Error: At least 2 sequences are required for phylogenetic tree construction",
            "error": "Insufficient sequences"
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
        
        # Create visualization
        plt.figure(figsize=(12, 8))
        draw(tree, axes=plt.gca(), do_show=False)
        plt.title("Phylogenetic Tree from Aligned Sequences")
        plt.tight_layout()
        
        # Save plot to base64 string
        img_buffer = io.BytesIO()
        plt.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
        img_buffer.seek(0)
        img_base64 = base64.b64encode(img_buffer.getvalue()).decode()
        plt.close()
        
        # Calculate tree statistics
        num_sequences = len(aligned_sequences)
        alignment_length = len(aligned_sequences[0]["sequence"])
        
        # Get tree topology information
        tree_height = tree.root.branch_length if tree.root.branch_length else 0
        num_clades = len(list(tree.get_terminals()))
        
        # Create result text
        result_text = f"""Phylogenetic tree constructed successfully from {num_sequences} aligned sequences.

Tree Statistics:
- Number of sequences: {num_sequences}
- Alignment length: {alignment_length}
- Tree height: {tree_height:.4f}
- Number of clades: {num_clades}
- Distance method: Identity
- Tree construction method: UPGMA

The tree shows the evolutionary relationships between the sequences based on sequence identity."""
        
        # Clean up temporary file
        os.unlink(temp_fasta_path)
        
        return {
            "text": result_text,
            "plot": {
                "data": [{
                    "type": "image",
                    "source": f"data:image/png;base64,{img_base64}",
                    "format": "png"
                }],
                "layout": {
                    "title": "Phylogenetic Tree",
                    "width": 800,
                    "height": 600
                }
            },
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
    
    if not aligned_sequences or not aligned_sequences.strip():
        return {
            "text": "Error: No aligned sequences provided",
            "error": "No sequences"
        }
    
    # Parse aligned sequences
    try:
        parsed_sequences = parse_aligned_sequences(aligned_sequences)
    except Exception as e:
        return {
            "text": f"Error parsing aligned sequences: {str(e)}",
            "error": str(e)
        }
    
    # Create phylogenetic tree
    return create_phylogenetic_tree(parsed_sequences)

from langchain.agents import tool

@tool
def run_phylogenetic_tree(aligned_sequences: str):
    """Create phylogenetic tree visualization from aligned sequences."""
    return run_phylogenetic_tree_raw(aligned_sequences) 