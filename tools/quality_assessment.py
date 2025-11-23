"""
Quality assessment utilities for sequencing data.

Provides quality metrics and reports for merged sequences from RNA-seq preprocessing.
"""

from __future__ import annotations

from typing import Dict, List


def _parse_fasta(sequences: str) -> List[tuple[str, str]]:
    """Parse FASTA format sequences into (header, sequence) tuples."""
    entries = []
    current_header = None
    current_seq = []
    
    for line in sequences.strip().split('\n'):
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current_header and current_seq:
                entries.append((current_header, ''.join(current_seq)))
            current_header = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    
    if current_header and current_seq:
        entries.append((current_header, ''.join(current_seq)))
    
    return entries


def _calculate_quality_metrics(sequences: List[tuple[str, str]]) -> Dict[str, any]:
    """Calculate quality metrics for a list of sequences."""
    if not sequences:
        return {
            "total_sequences": 0,
            "average_length": 0,
            "min_length": 0,
            "max_length": 0,
            "total_bases": 0,
            "gc_content": 0.0,
            "n_content": 0.0
        }
    
    lengths = [len(seq) for _, seq in sequences]
    total_bases = sum(lengths)
    
    # Calculate GC and N content
    all_bases = ''.join(seq for _, seq in sequences)
    gc_count = all_bases.upper().count('G') + all_bases.upper().count('C')
    n_count = all_bases.upper().count('N')
    
    return {
        "total_sequences": len(sequences),
        "average_length": sum(lengths) / len(lengths) if lengths else 0,
        "min_length": min(lengths) if lengths else 0,
        "max_length": max(lengths) if lengths else 0,
        "total_bases": total_bases,
        "gc_content": (gc_count / total_bases * 100) if total_bases > 0 else 0.0,
        "n_content": (n_count / total_bases * 100) if total_bases > 0 else 0.0
    }


def run_quality_assessment_raw(sequences: str) -> Dict[str, object]:
    """
    Generate a quality assessment report for merged sequences.
    
    Args:
        sequences: FASTA-formatted sequences (merged reads)
    
    Returns:
        Dictionary containing quality metrics and report text.
    """
    if not sequences or len(sequences.strip()) == 0:
        return {
            "status": "error",
            "message": "No sequences provided for quality assessment"
        }
    
    # Parse sequences
    parsed_sequences = _parse_fasta(sequences)
    
    if not parsed_sequences:
        return {
            "status": "error",
            "message": "No valid sequences found in input"
        }
    
    # Calculate metrics
    metrics = _calculate_quality_metrics(parsed_sequences)
    
    # Generate report text
    report_lines = [
        "Quality Assessment Report",
        "=" * 50,
        f"Total Sequences: {metrics['total_sequences']}",
        f"Total Bases: {metrics['total_bases']:,}",
        f"Average Length: {metrics['average_length']:.1f} bp",
        f"Length Range: {metrics['min_length']} - {metrics['max_length']} bp",
        f"GC Content: {metrics['gc_content']:.2f}%",
        f"N Content: {metrics['n_content']:.2f}%",
        "",
        "Quality Assessment:",
    ]
    
    # Add quality assessment
    if metrics['total_sequences'] > 0:
        if metrics['average_length'] > 100:
            report_lines.append("✓ Good average read length")
        else:
            report_lines.append("⚠ Low average read length")
        
        if metrics['gc_content'] > 30 and metrics['gc_content'] < 70:
            report_lines.append("✓ GC content within normal range")
        else:
            report_lines.append("⚠ GC content outside normal range (30-70%)")
        
        if metrics['n_content'] < 5:
            report_lines.append("✓ Low N content (good quality)")
        else:
            report_lines.append("⚠ High N content (may indicate low quality)")
    else:
        report_lines.append("⚠ No sequences to assess")
    
    report_text = "\n".join(report_lines)
    
    return {
        "text": report_text,
        "metrics": metrics,
        "status": "success"
    }


__all__ = ["run_quality_assessment_raw"]

