"""
Lightweight read merging utilities for paired-end sequencing data.

The implementation provides a simple consensus merge for demonstration
purposes. It identifies an overlap between forward reads and the reverse
complement of reverse reads, then stitches them together. This avoids
external dependencies while supporting MCP workflows.
"""

from __future__ import annotations

from typing import Dict, List, Tuple


COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def _parse_fastq_pairs(forward: str, reverse: str) -> List[Tuple[str, str]]:
    def chunk_fastq(content: str) -> List[str]:
        lines = [line.strip() for line in content.strip().splitlines() if line.strip()]
        return [lines[i + 1] for i in range(0, len(lines) - 3, 4)]

    forward_reads = chunk_fastq(forward)
    reverse_reads = chunk_fastq(reverse)
    return list(zip(forward_reads, reverse_reads))


def _reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def _merge_pair(
    forward: str,
    reverse: str,
    min_overlap: int,
) -> Tuple[str, int]:
    """Merge a forward/reverse pair using the longest overlap ≥ min_overlap."""
    rc_reverse = _reverse_complement(reverse)
    max_overlap = 0
    merged_sequence = forward + rc_reverse

    for overlap in range(min(len(forward), len(rc_reverse)), min_overlap - 1, -1):
        if forward.endswith(rc_reverse[:overlap]):
            merged_sequence = forward + rc_reverse[overlap:]
            max_overlap = overlap
            break

    return merged_sequence, max_overlap


def run_read_merging_raw(
    forward_reads: str,
    reverse_reads: str,
    min_overlap: int = 12,
) -> Dict[str, object]:
    """
    Merge paired-end reads into consensus sequences.

    Args:
        forward_reads: FASTQ string containing forward reads.
        reverse_reads: FASTQ string containing reverse reads.
        min_overlap: Minimum overlap length to merge reads.

    Returns:
        Dictionary containing merged FASTA content and summary metrics.
    """
    pairs = _parse_fastq_pairs(forward_reads, reverse_reads)

    merged_sequences: List[str] = []
    overlaps: List[int] = []

    for idx, (fwd, rev) in enumerate(pairs, start=1):
        merged, overlap = _merge_pair(fwd, rev, min_overlap)
        merged_sequences.append(f">merged_{idx}\n{merged}")
        overlaps.append(overlap)

    merged_fasta = "\n".join(merged_sequences)
    summary = {
        "total_pairs": len(pairs),
        "merged_pairs": sum(1 for o in overlaps if o >= min_overlap),
        "min_overlap": min_overlap,
        "average_overlap": sum(overlaps) / len(overlaps) if overlaps else 0,
    }

    return {
        "text": "Read merging completed successfully.",
        "merged_sequences": merged_fasta,
        "summary": summary,
    }


__all__ = ["run_read_merging_raw"]


