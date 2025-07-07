# /backend/tools/alignment.py

import re
from typing import List, Dict, Any
from Bio import AlignIO, SeqIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline
import subprocess
import tempfile
import os
from io import StringIO

def parse_fasta_string(fasta_string: str) -> List[Dict[str, str]]:
    """Parse FASTA format string into list of sequences."""
    print(f"🔍 parse_fasta_string input: '{fasta_string}'")
    sequences = []
    current_name = ""
    current_sequence = ""
    
    for line in fasta_string.strip().split('\n'):
        line = line.strip()
        print(f"🔍 Processing line: '{line}'")
        if line.startswith('>'):
            # Save previous sequence if exists
            if current_name and current_sequence:
                print(f"🔍 Adding sequence: name='{current_name}', sequence='{current_sequence}'")
                sequences.append({
                    "name": current_name,
                    "sequence": current_sequence
                })
            # Start new sequence
            current_name = line[1:].strip()
            current_sequence = ""
            print(f"🔍 Starting new sequence with name: '{current_name}'")
        else:
            # Add to current sequence
            current_sequence += line
            print(f"🔍 Added to current sequence: '{line}' -> current_sequence='{current_sequence}'")
    
    # Add the last sequence
    if current_name and current_sequence:
        print(f"🔍 Adding final sequence: name='{current_name}', sequence='{current_sequence}'")
        sequences.append({
            "name": current_name,
            "sequence": current_sequence
        })
    
    print(f"🔍 Final parsed sequences: {sequences}")
    return sequences

def validate_sequences(sequences: List[Dict[str, str]]) -> str:
    """Validate that sequences are valid DNA/RNA."""
    valid_bases = set("ATCGU")
    
    print(f"🔍 Validating {len(sequences)} sequences")
    for seq_data in sequences:
        sequence = seq_data["sequence"].upper()
        print(f"🔍 Validating sequence '{sequence}' from '{seq_data['name']}'")
        invalid_bases = [base for base in sequence if base not in valid_bases]
        if invalid_bases:
            print(f"🔍 Found invalid bases: {invalid_bases}")
            return f"Invalid sequence '{sequence}' in {seq_data['name']}. Only A, T, C, G, U allowed."
    
    print(f"🔍 All sequences validated successfully")
    return ""

def run_alignment(sequences: str):
    """Perform sequence alignment on a given set of sequences."""
    
    if not sequences or not sequences.strip():
        return {
            "text": "Error: No sequences provided",
            "alignment": [],
            "statistics": {}
        }
    
    # Parse FASTA input
    try:
        parsed_sequences = parse_fasta_string(sequences)
    except Exception as e:
        return {
            "text": f"Error parsing FASTA format: {str(e)}",
            "alignment": [],
            "statistics": {}
        }
    
    if len(parsed_sequences) < 2:
        return {
            "text": "Error: At least 2 sequences are required for alignment",
            "alignment": [],
            "statistics": {}
        }
    
    # Validate sequences
    validation_error = validate_sequences(parsed_sequences)
    if validation_error:
        return {
            "text": f"Error: {validation_error}",
            "alignment": [],
            "statistics": {}
        }
    
    # Create temporary FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
        for seq_data in parsed_sequences:
            temp_fasta.write(f">{seq_data['name']}\n{seq_data['sequence']}\n")
        temp_fasta_path = temp_fasta.name
    
    # Create temporary output file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.aln', delete=False) as temp_output:
        temp_output_path = temp_output.name
    
    try:
        # Try to use ClustalW for alignment
        try:
            clustalw_cline = ClustalwCommandline("clustalw", infile=temp_fasta_path, outfile=temp_output_path)
            stdout, stderr = clustalw_cline()
            
            # Read the alignment
            alignment = AlignIO.read(temp_output_path, "clustal")
            alignment_method = "ClustalW"
            
        except (subprocess.CalledProcessError, FileNotFoundError):
            # Fallback to Muscle if ClustalW is not available
            try:
                muscle_cline = MuscleCommandline(input=temp_fasta_path, out=temp_output_path)
                stdout, stderr = muscle_cline()
                
                # Read the alignment
                alignment = AlignIO.read(temp_output_path, "fasta")
                alignment_method = "Muscle"
                
            except (subprocess.CalledProcessError, FileNotFoundError):
                # Manual alignment as fallback
                alignment = perform_manual_alignment(parsed_sequences)
                alignment_method = "Manual"
        
        # Convert alignment to list format
        aligned_sequences = []
        for record in alignment:
            aligned_sequences.append({
                "name": record.id,
                "sequence": str(record.seq)
            })
        
        # Calculate alignment statistics
        alignment_length = len(aligned_sequences[0]["sequence"])
        total_sequences = len(aligned_sequences)
        
        # Calculate identity matrix
        identity_matrix = []
        for i in range(total_sequences):
            row = []
            for j in range(total_sequences):
                if i == j:
                    row.append(100.0)
                else:
                    seq1 = aligned_sequences[i]["sequence"]
                    seq2 = aligned_sequences[j]["sequence"]
                    matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-')
                    total = sum(1 for a, b in zip(seq1, seq2) if a != '-' or b != '-')
                    identity = (matches / total * 100) if total > 0 else 0
                    row.append(round(identity, 2))
            identity_matrix.append(row)
        
        # Calculate average identity
        total_identity = 0
        count = 0
        for i in range(total_sequences):
            for j in range(i + 1, total_sequences):
                total_identity += identity_matrix[i][j]
                count += 1
        average_identity = total_identity / count if count > 0 else 0
        
        # Create result text
        result_text = f"""Sequences aligned successfully using {alignment_method}.

Alignment Statistics:
- Number of sequences: {total_sequences}
- Alignment length: {alignment_length}
- Average pairwise identity: {average_identity:.2f}%

Aligned sequences:
"""
        for seq in aligned_sequences:
            result_text += f"{seq['name']}: {seq['sequence']}\n"
        
        return {
            "text": result_text,
            "alignment": aligned_sequences,
            "statistics": {
                "method": alignment_method,
                "num_sequences": total_sequences,
                "alignment_length": alignment_length,
                "average_identity": round(average_identity, 2),
                "identity_matrix": identity_matrix
            }
        }
        
    except Exception as e:
        return {
            "text": f"Error during alignment: {str(e)}",
            "alignment": [],
            "statistics": {}
        }
    
    finally:
        # Clean up temporary files
        try:
            os.unlink(temp_fasta_path)
            os.unlink(temp_output_path)
        except:
            pass

def perform_manual_alignment(sequences: List[Dict[str, str]]) -> Any:
    """Perform a simple manual alignment as fallback."""
    from Bio.Align import MultipleSeqAlignment
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    # This is a very basic alignment - in practice, you'd want a more sophisticated algorithm
    max_length = max(len(seq["sequence"]) for seq in sequences)
    
    # Create SeqRecord objects for the alignment
    seq_records = []
    for seq_data in sequences:
        sequence = seq_data["sequence"]
        # Pad with gaps to match the longest sequence
        padded_sequence = sequence + "-" * (max_length - len(sequence))
        
        # Create SeqRecord with proper id and seq attributes
        seq_record = SeqRecord(
            Seq(padded_sequence),
            id=seq_data["name"],
            description=""
        )
        seq_records.append(seq_record)
    
    # Create MultipleSeqAlignment object
    alignment = MultipleSeqAlignment(seq_records)
    return alignment

from langchain.agents import tool

@tool
def run_alignment_tool(sequences: str):
    """Perform sequence alignment on a given set of sequences."""
    return run_alignment(sequences)

