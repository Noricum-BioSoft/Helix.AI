"""
Utility functions for basic read trimming operations.

These helpers provide a lightweight, pure-Python fallback so the MCP
tooling can demonstrate read trimming behaviour without requiring external
executables such as cutadapt or trimmomatic. They operate on FASTQ-formatted
strings and support optional adapter removal and quality-based tail trimming.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple


PHRED_OFFSET = 33


@dataclass
class FastqRecord:
    header: str
    sequence: str
    separator: str
    quality: str

    def to_fastq(self) -> str:
        return "\n".join([self.header, self.sequence, self.separator, self.quality])


def _parse_fastq(content: str) -> List[FastqRecord]:
    """Parse FASTQ content into records. Ignores incomplete trailing records."""
    lines = [line.strip() for line in content.strip().splitlines() if line.strip()]
    records: List[FastqRecord] = []
    for i in range(0, len(lines) - 3, 4):
        header, sequence, separator, quality = lines[i : i + 4]
        if not header.startswith("@") or not separator.startswith("+"):
            continue
        records.append(FastqRecord(header, sequence, separator, quality))
    return records


def _quality_trim(sequence: str, quality: str, threshold: int) -> Tuple[str, str]:
    """Trim low-quality bases from the end of a read."""
    # Ensure quality and sequence have the same length
    min_length = min(len(sequence), len(quality))
    if min_length == 0:
        return "", ""
    
    # Trim both to the same length if they differ
    sequence = sequence[:min_length]
    quality = quality[:min_length]
    
    trim_index = min_length
    for idx in range(min_length - 1, -1, -1):
        if idx >= len(quality):
            break
        q_score = ord(quality[idx]) - PHRED_OFFSET
        if q_score >= threshold:
            break
        trim_index = idx
    return sequence[:trim_index], quality[:trim_index]


def _adapter_trim(sequence: str, adapter: str) -> str:
    """Remove adapter sequence from the end of the sequence."""
    if not adapter:
        return sequence
    
    # Check if adapter is present (exact match)
    if adapter in sequence:
        return sequence.replace(adapter, "")
    
    # Also check for adapter at the end (common case)
    # This handles cases where adapter might be partially present
    if sequence.endswith(adapter):
        return sequence[:-len(adapter)]
    
    # Check if sequence ends with a prefix of the adapter (partial adapter contamination)
    # Remove the longest matching suffix
    for i in range(len(adapter), 0, -1):
        if sequence.endswith(adapter[:i]):
            return sequence[:-i]
    
    return sequence


def run_read_trimming_raw(
    reads: str,
    adapter: str | None = None,
    quality_threshold: int = 20,
) -> Dict[str, object]:
    """
    Perform simple read trimming on FASTQ-formatted reads.

    Args:
        reads: FASTQ data as a single string.
        adapter: Optional adapter sequence to remove.
        quality_threshold: Phred score threshold for right-tail trimming.

    Returns:
        Dictionary with trimmed FASTQ content and summary statistics.
    """
    adapter = adapter or ""
    records = _parse_fastq(reads)

    trimmed_records: List[FastqRecord] = []
    total_bases = 0
    trimmed_bases = 0

    for record in records:
        total_bases += len(record.sequence)
        seq = _adapter_trim(record.sequence, adapter)
        # Ensure quality string matches sequence length
        # If adapter was removed, we need to adjust quality accordingly
        if len(seq) < len(record.sequence):
            # Adapter was removed, take quality corresponding to remaining sequence
            # This assumes adapter was at the end (most common case)
            qual = record.quality[:len(seq)]
        else:
            qual = record.quality[:len(seq)]
        
        # Ensure quality and sequence are same length before quality trimming
        min_len = min(len(seq), len(qual))
        seq = seq[:min_len]
        qual = qual[:min_len]
        
        seq, qual = _quality_trim(seq, qual, quality_threshold)
        trimmed_bases += len(record.sequence) - len(seq)
        trimmed_records.append(
            FastqRecord(record.header, seq, record.separator, qual)
        )

    trimmed_fastq = "\n".join(rec.to_fastq() for rec in trimmed_records)
    summary = {
        "total_reads": len(records),
        "total_bases": total_bases,
        "trimmed_bases": trimmed_bases,
        "adapter_removed": bool(adapter),
        "quality_threshold": quality_threshold,
    }

    return {
        "text": "Read trimming completed successfully.",
        "trimmed_reads": trimmed_fastq,
        "summary": summary,
    }


__all__ = ["run_read_trimming_raw"]


