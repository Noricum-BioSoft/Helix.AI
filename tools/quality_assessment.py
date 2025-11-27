"""
Quality assessment utilities for sequencing data analysis.

Provides quality metrics and statistics for merged sequences,
including length distribution, GC content, base composition, and quality scores.
"""

from __future__ import annotations

from typing import Dict, List
from collections import Counter
import statistics


def _parse_fasta(content: str) -> List[tuple[str, str]]:
    """Parse FASTA content into (header, sequence) tuples."""
    sequences: List[tuple[str, str]] = []
    lines = [line.strip() for line in content.strip().splitlines() if line.strip()]
    
    current_header = None
    current_sequence = []
    
    for line in lines:
        if line.startswith(">"):
            # Save previous sequence if exists
            if current_header and current_sequence:
                sequences.append((current_header, "".join(current_sequence)))
            # Start new sequence
            current_header = line[1:].strip()
            current_sequence = []
        else:
            current_sequence.append(line)
    
    # Don't forget the last sequence
    if current_header and current_sequence:
        sequences.append((current_header, "".join(current_sequence)))
    
    return sequences


def _calculate_gc_content(sequence: str) -> float:
    """Calculate GC content percentage."""
    if not sequence:
        return 0.0
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100


def _calculate_base_composition(sequence: str) -> Dict[str, int]:
    """Calculate base composition."""
    seq_upper = sequence.upper()
    return {
        'A': seq_upper.count('A'),
        'T': seq_upper.count('T'),
        'G': seq_upper.count('G'),
        'C': seq_upper.count('C'),
        'N': seq_upper.count('N'),
        'other': len(sequence) - (seq_upper.count('A') + seq_upper.count('T') + 
                                  seq_upper.count('G') + seq_upper.count('C') + 
                                  seq_upper.count('N'))
    }


def run_quality_assessment_raw(sequences: str) -> Dict[str, object]:
    """
    Perform quality assessment on merged sequences.
    
    Args:
        sequences: FASTA format sequences (merged reads).
    
    Returns:
        Dictionary containing quality metrics, statistics, and visualization data.
    """
    parsed_sequences = _parse_fasta(sequences)
    
    if not parsed_sequences:
        return {
            "text": "No sequences found for quality assessment.",
            "status": "error",
            "message": "No sequences provided"
        }
    
    # Calculate metrics for each sequence
    lengths: List[int] = []
    gc_contents: List[float] = []
    base_compositions: List[Dict[str, int]] = []
    
    for header, sequence in parsed_sequences:
        seq_len = len(sequence)
        lengths.append(seq_len)
        gc_contents.append(_calculate_gc_content(sequence))
        base_compositions.append(_calculate_base_composition(sequence))
    
    # Aggregate statistics
    total_bases = sum(lengths)
    total_sequences = len(parsed_sequences)
    
    # Base composition totals
    total_composition = {
        'A': sum(comp['A'] for comp in base_compositions),
        'T': sum(comp['T'] for comp in base_compositions),
        'G': sum(comp['G'] for comp in base_compositions),
        'C': sum(comp['C'] for comp in base_compositions),
        'N': sum(comp['N'] for comp in base_compositions),
        'other': sum(comp['other'] for comp in base_compositions)
    }
    
    # Calculate statistics
    metrics = {
        "total_sequences": total_sequences,
        "total_bases": total_bases,
        "average_length": statistics.mean(lengths) if lengths else 0,
        "median_length": statistics.median(lengths) if lengths else 0,
        "min_length": min(lengths) if lengths else 0,
        "max_length": max(lengths) if lengths else 0,
        "length_std_dev": statistics.stdev(lengths) if len(lengths) > 1 else 0,
        "average_gc_content": statistics.mean(gc_contents) if gc_contents else 0,
        "gc_content_std_dev": statistics.stdev(gc_contents) if len(gc_contents) > 1 else 0,
        "base_composition": total_composition,
        "base_percentages": {
            base: (count / total_bases * 100) if total_bases > 0 else 0
            for base, count in total_composition.items()
        }
    }
    
    # Create length distribution for visualization
    length_distribution = Counter(lengths)
    sorted_lengths = sorted(length_distribution.keys())
    length_counts = [length_distribution[length] for length in sorted_lengths]
    
    # Create GC content distribution
    gc_bins = [0, 20, 30, 40, 50, 60, 70, 100]
    gc_distribution = [0] * (len(gc_bins) - 1)
    for gc in gc_contents:
        for i in range(len(gc_bins) - 1):
            if gc_bins[i] <= gc < gc_bins[i + 1]:
                gc_distribution[i] += 1
                break
    
    # Create visualization data for Plotly
    plot_data = {
        "length_distribution": {
            "x": sorted_lengths,
            "y": length_counts,
            "type": "bar",
            "name": "Sequence Length Distribution"
        },
        "gc_content_distribution": {
            "x": [f"{gc_bins[i]}-{gc_bins[i+1]}%" for i in range(len(gc_bins) - 1)],
            "y": gc_distribution,
            "type": "bar",
            "name": "GC Content Distribution"
        },
        "base_composition": {
            "x": list(total_composition.keys()),
            "y": [total_composition[base] for base in total_composition.keys()],
            "type": "bar",
            "name": "Base Composition"
        }
    }
    
    return {
        "text": f"Quality assessment completed successfully. Analyzed {total_sequences} sequences.",
        "metrics": metrics,
        "plot_data": plot_data,
        "summary": {
            "total_sequences": total_sequences,
            "total_bases": total_bases,
            "average_length": round(metrics["average_length"], 2),
            "average_gc_content": round(metrics["average_gc_content"], 2),
            "length_range": f"{metrics['min_length']}-{metrics['max_length']} bp"
        }
    }


__all__ = ["run_quality_assessment_raw"]


