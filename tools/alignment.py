# /backend/tools/alignment.py

import re
from typing import List, Dict, Any
from Bio import AlignIO, SeqIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import ClustalwCommandline
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
    
    # Check if this is inline FASTA format (multiple >seq on same line)
    # Pattern: >seq1 ATGC >seq2 ATGC >seq3 ATGC
    if fasta_string.count('>') > 1 and '\n' not in fasta_string.strip():
        # Inline FASTA format - split by > and parse each
        print(f"🔍 Detected inline FASTA format")
        parts = fasta_string.strip().split('>')
        for part in parts:
            part = part.strip()
            if not part:
                continue
            # Split on whitespace - first part is name, rest is sequence
            parts_split = part.split(None, 1)
            if len(parts_split) == 2:
                name, sequence = parts_split
                sequences.append({
                    "name": name,
                    "sequence": sequence
                })
                print(f"🔍 Parsed inline sequence: name='{name}', sequence='{sequence}'")
            elif len(parts_split) == 1:
                # Just a name, sequence might be on next iteration or empty
                current_name = parts_split[0]
                current_sequence = ""
        if sequences:
            print(f"🔍 Final parsed sequences from inline format: {sequences}")
            return sequences
    
    # Standard FASTA format - split by lines
    lines = fasta_string.strip().split('\n')
    
    for line in lines:
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
            
            # Parse the header line - it might contain sequence data
            header_parts = line[1:].split(None, 1)  # Split on whitespace, max 1 split
            if len(header_parts) == 2:
                # Header contains both name and sequence
                current_name = header_parts[0]
                current_sequence = header_parts[1]
                print(f"🔍 Header contains sequence: name='{current_name}', sequence='{current_sequence}'")
            else:
                # Header contains only name
                current_name = header_parts[0]
                current_sequence = ""
                print(f"🔍 Starting new sequence with name: '{current_name}'")
        else:
            # Add to current sequence
            if current_name:  # Only add if we have a name (we're in a sequence)
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
                # Use Muscle directly with subprocess (matching command line: muscle -align input.fasta -output output.fasta)
                # BioPython's MuscleCommandline is deprecated and doesn't support -align/-output flags
                result = subprocess.run(
                    ["muscle", "-align", temp_fasta_path, "-output", temp_output_path],
                    capture_output=True,
                    text=True,
                    check=True
                )
                
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
    """Perform alignment using BioPython's PairwiseAligner with star alignment."""
    from Bio.Align import MultipleSeqAlignment, PairwiseAligner
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    if len(sequences) < 2:
        # If only one sequence, just return it
        seq_data = sequences[0]
        seq_record = SeqRecord(Seq(seq_data["sequence"]), id=seq_data["name"], description="")
        return MultipleSeqAlignment([seq_record])
    
    # Use BioPython's PairwiseAligner for better alignment
    aligner = PairwiseAligner()
    aligner.mode = 'global'  # Global alignment
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.5
    
    # Convert sequences to Seq objects
    seq_objects = [(Seq(seq_data["sequence"]), seq_data["name"]) for seq_data in sequences]
    
    # Star alignment: use the longest sequence as reference
    ref_idx = max(range(len(seq_objects)), key=lambda i: len(seq_objects[i][0]))
    ref_seq, ref_name = seq_objects[ref_idx]
    
    # Align all sequences to the reference and collect alignments
    alignments_to_ref = []
    for i, (seq, name) in enumerate(seq_objects):
        if i == ref_idx:
            # Reference sequence - no gaps
            alignments_to_ref.append((str(ref_seq), name, []))
            continue
        
        alignments = aligner.align(ref_seq, seq)
        
        if alignments:
            best = alignments[0]
            # best[0] is aligned reference (with gaps), best[1] is aligned query
            aligned_ref = str(best[0])
            aligned_query = str(best[1])
            
            # Find gap positions in the aligned reference
            gap_positions = [j for j, char in enumerate(aligned_ref) if char == '-']
            alignments_to_ref.append((aligned_query, name, gap_positions))
        else:
            # Fallback: pad to reference length
            padded = str(seq) + "-" * (len(ref_seq) - len(seq))
            alignments_to_ref.append((padded, name, []))
    
    # Merge all gap positions - we need to insert gaps into all sequences at all gap positions
    all_gap_positions = set()
    for _, _, gaps in alignments_to_ref:
        all_gap_positions.update(gaps)
    
    # Build final aligned sequences
    final_seqs = []
    for aligned_seq, name, gaps in alignments_to_ref:
        # Insert gaps at all positions where any sequence has gaps
        seq_list = list(aligned_seq)
        for gap_pos in sorted(all_gap_positions, reverse=True):
            if gap_pos < len(seq_list):
                if seq_list[gap_pos] != '-':  # Don't double-insert
                    seq_list.insert(gap_pos, '-')
            elif gap_pos >= len(seq_list):
                # Gap is beyond current sequence - add at end
                seq_list.append('-')
        final_seqs.append(SeqRecord(Seq(''.join(seq_list)), id=name, description=""))
    
    return MultipleSeqAlignment(final_seqs)

from langchain.agents import tool

@tool
def run_alignment_tool(sequences: str):
    """Perform sequence alignment on a given set of sequences."""
    return run_alignment(sequences)

