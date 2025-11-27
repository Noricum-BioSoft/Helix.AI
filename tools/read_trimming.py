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
    trim_index = len(sequence)
    for idx in range(len(sequence) - 1, -1, -1):
        q_score = ord(quality[idx]) - PHRED_OFFSET
        if q_score >= threshold:
            break
        trim_index = idx
    return sequence[:trim_index], quality[:trim_index]


def _adapter_trim(sequence: str, adapter: str) -> str:
    """
    Remove adapter sequence from read.
    Tries multiple strategies:
    1. Exact match anywhere in sequence
    2. Adapter at the start (prefix match) - common in reverse reads
    3. Adapter at the end (suffix match) - common in forward reads
    4. Partial match at start (first N bases matching adapter start)
    5. Partial match at end (last N bases matching adapter start)
    """
    if not adapter:
        return sequence
    
    adapter_len = len(adapter)
    
    # Strategy 1: Exact match anywhere
    if adapter in sequence:
        return sequence.replace(adapter, "", 1)  # Replace only first occurrence
    
    # Strategy 2: Adapter at the start (prefix) - common in reverse reads
    if sequence.startswith(adapter):
        return sequence[adapter_len:]
    
    # Strategy 3: Adapter at the end (suffix) - common in forward reads
    if sequence.endswith(adapter):
        return sequence[:-adapter_len]
    
    # Strategy 4: Partial match at start - check if first part of sequence matches start of adapter
    # This handles cases where the adapter is partially sequenced at the beginning
    for overlap_len in range(min(adapter_len, len(sequence)), max(3, adapter_len // 2), -1):
        sequence_prefix = sequence[:overlap_len]
        adapter_prefix = adapter[:overlap_len]
        # Allow 1-2 mismatches for quality issues
        if _sequences_match(sequence_prefix, adapter_prefix, max_mismatches=2):
            return sequence[overlap_len:]
    
    # Strategy 5: Partial match at end - check if last part of sequence matches start of adapter
    # This handles cases where the adapter is partially sequenced at the end
    for overlap_len in range(min(adapter_len, len(sequence)), max(3, adapter_len // 2), -1):
        sequence_suffix = sequence[-overlap_len:]
        adapter_prefix = adapter[:overlap_len]
        # Allow 1-2 mismatches for quality issues
        if _sequences_match(sequence_suffix, adapter_prefix, max_mismatches=2):
            return sequence[:-overlap_len]
    
    return sequence


def _sequences_match(seq1: str, seq2: str, max_mismatches: int = 0) -> bool:
    """Check if two sequences match with allowed mismatches."""
    if len(seq1) != len(seq2):
        return False
    mismatches = sum(1 for a, b in zip(seq1, seq2) if a != b)
    return mismatches <= max_mismatches


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
    adapter_bases_removed = 0

    for record in records:
        original_len = len(record.sequence)
        total_bases += original_len
        
        # First remove adapter
        seq_after_adapter = _adapter_trim(record.sequence, adapter)
        adapter_removed = original_len - len(seq_after_adapter)
        adapter_bases_removed += adapter_removed
        
        # Then trim quality
        qual = record.quality[: len(seq_after_adapter)]
        seq, qual = _quality_trim(seq_after_adapter, qual, quality_threshold)
        
        # Total trimmed bases = adapter removal + quality trimming
        trimmed_bases += original_len - len(seq)
        
        trimmed_records.append(
            FastqRecord(record.header, seq, record.separator, qual)
        )
    
    print(f"🔧 [DEBUG] Adapter removal stats: {adapter_bases_removed} bases removed from {len(records)} reads")

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


