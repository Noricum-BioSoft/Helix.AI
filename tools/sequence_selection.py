import random
import re
from typing import List, Dict, Any
from collections import Counter

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

def calculate_sequence_statistics(sequences: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    """Calculate statistics for each sequence."""
    enriched_sequences = []
    
    for seq_data in sequences:
        sequence = seq_data["sequence"]
        
        # Calculate basic statistics
        length = len(sequence)
        gc_content = (sequence.count('G') + sequence.count('C')) / length * 100 if length > 0 else 0
        gap_count = sequence.count('-')
        gap_percentage = (gap_count / length * 100) if length > 0 else 0
        
        # Calculate conservation score (how many positions match the consensus)
        if len(sequences) > 1:
            # Find consensus at each position
            consensus = ""
            for pos in range(length):
                bases = [seq["sequence"][pos] for seq in sequences if pos < len(seq["sequence"])]
                base_counts = Counter(bases)
                most_common_base = base_counts.most_common(1)[0][0]
                consensus += most_common_base
            
            # Calculate how well this sequence matches consensus
            matches = sum(1 for a, b in zip(sequence, consensus) if a == b)
            conservation_score = (matches / length * 100) if length > 0 else 0
        else:
            conservation_score = 100.0
        
        enriched_sequences.append({
            **seq_data,
            "statistics": {
                "length": length,
                "gc_content": round(gc_content, 2),
                "gap_count": gap_count,
                "gap_percentage": round(gap_percentage, 2),
                "conservation_score": round(conservation_score, 2)
            }
        })
    
    return enriched_sequences

def select_sequences_by_criteria(sequences: List[Dict[str, str]], 
                               selection_type: str = "random", 
                               num_sequences: int = 1,
                               criteria: str = None) -> Dict[str, Any]:
    """Select sequences based on specified criteria."""
    
    if not sequences:
        return {
            "text": "Error: No sequences provided for selection",
            "selected_sequences": [],
            "error": "No sequences"
        }
    
    # Calculate statistics for all sequences
    enriched_sequences = calculate_sequence_statistics(sequences)
    
    selected_sequences = []
    selection_explanation = ""
    
    if selection_type == "random":
        # Random selection
        if num_sequences > len(sequences):
            num_sequences = len(sequences)
        
        selected_indices = random.sample(range(len(sequences)), num_sequences)
        selected_sequences = [enriched_sequences[i] for i in selected_indices]
        selection_explanation = f"Randomly selected {len(selected_sequences)} sequence(s) from {len(sequences)} total sequences."
        
    elif selection_type == "best_conservation":
        # Select sequences with highest conservation scores
        sorted_sequences = sorted(enriched_sequences, 
                                key=lambda x: x["statistics"]["conservation_score"], 
                                reverse=True)
        selected_sequences = sorted_sequences[:num_sequences]
        selection_explanation = f"Selected {len(selected_sequences)} sequence(s) with highest conservation scores."
        
    elif selection_type == "lowest_gaps":
        # Select sequences with fewest gaps
        sorted_sequences = sorted(enriched_sequences, 
                                key=lambda x: x["statistics"]["gap_count"])
        selected_sequences = sorted_sequences[:num_sequences]
        selection_explanation = f"Selected {len(selected_sequences)} sequence(s) with fewest gaps."
        
    elif selection_type == "highest_gc":
        # Select sequences with highest GC content
        sorted_sequences = sorted(enriched_sequences, 
                                key=lambda x: x["statistics"]["gc_content"], 
                                reverse=True)
        selected_sequences = sorted_sequences[:num_sequences]
        selection_explanation = f"Selected {len(selected_sequences)} sequence(s) with highest GC content."
        
    elif selection_type == "longest":
        # Select longest sequences
        sorted_sequences = sorted(enriched_sequences, 
                                key=lambda x: x["statistics"]["length"], 
                                reverse=True)
        selected_sequences = sorted_sequences[:num_sequences]
        selection_explanation = f"Selected {len(selected_sequences)} longest sequence(s)."
        
    elif selection_type == "shortest":
        # Select shortest sequences
        sorted_sequences = sorted(enriched_sequences, 
                                key=lambda x: x["statistics"]["length"])
        selected_sequences = sorted_sequences[:num_sequences]
        selection_explanation = f"Selected {len(selected_sequences)} shortest sequence(s)."
        
    else:
        # Default to random selection
        if num_sequences > len(sequences):
            num_sequences = len(sequences)
        selected_indices = random.sample(range(len(sequences)), num_sequences)
        selected_sequences = [enriched_sequences[i] for i in selected_indices]
        selection_explanation = f"Randomly selected {len(selected_sequences)} sequence(s) from {len(sequences)} total sequences."
    
    # Create result text
    result_text = f"""Sequence selection completed successfully.

{selection_explanation}

Selected sequences:"""
    
    for i, seq in enumerate(selected_sequences, 1):
        stats = seq["statistics"]
        result_text += f"""

{i}. {seq['name']}
   Sequence: {seq['sequence']}
   Length: {stats['length']} bp
   GC Content: {stats['gc_content']}%
   Gaps: {stats['gap_count']} ({stats['gap_percentage']}%)
   Conservation Score: {stats['conservation_score']}%"""
    
    return {
        "text": result_text,
        "selected_sequences": selected_sequences,
        "selection_criteria": {
            "type": selection_type,
            "num_requested": num_sequences,
            "num_selected": len(selected_sequences),
            "total_available": len(sequences),
            "explanation": selection_explanation
        }
    }

def run_sequence_selection_raw(aligned_sequences: str, 
                              selection_type: str = "random", 
                              num_sequences: int = 1):
    """Select sequences from aligned sequences based on criteria."""
    
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
    
    if len(parsed_sequences) < 1:
        return {
            "text": "Error: At least 1 sequence is required for selection",
            "error": "Insufficient sequences"
        }
    
    # Select sequences
    return select_sequences_by_criteria(parsed_sequences, selection_type, num_sequences)

from langchain.agents import tool

@tool
def run_sequence_selection(aligned_sequences: str, selection_type: str = "random", num_sequences: int = 1):
    """Select sequences from aligned sequences based on various criteria (random, best_conservation, lowest_gaps, highest_gc, longest, shortest)."""
    return run_sequence_selection_raw(aligned_sequences, selection_type, num_sequences) 