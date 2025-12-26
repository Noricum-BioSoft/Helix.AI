"""
Lightweight read merging utilities for paired-end sequencing data.

The implementation provides a simple consensus merge for demonstration
purposes. It identifies an overlap between forward reads and the reverse
complement of reverse reads, then stitches them together. This avoids
external dependencies while supporting MCP workflows.
"""

from __future__ import annotations

from typing import Dict, List, Tuple, Optional
import boto3
import tempfile
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

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


def _parse_s3_path(s3_path: str) -> Tuple[str, str]:
    """Parse S3 path into bucket and key."""
    if not s3_path.startswith("s3://"):
        raise ValueError(f"Invalid S3 path: {s3_path}")
    path_without_prefix = s3_path[5:]
    parts = path_without_prefix.split("/", 1)
    bucket = parts[0]
    key = parts[1] if len(parts) > 1 else ""
    return bucket, key


def _parse_fastq_with_quality(content: str) -> List[Tuple[str, str, str, str]]:
    """Parse FASTQ content into (header, sequence, plus, quality) tuples."""
    lines = [line.rstrip() for line in content.splitlines()]
    records = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("@"):
            header = lines[i]
            if i + 3 < len(lines):
                sequence = lines[i + 1]
                plus = lines[i + 2]
                quality = lines[i + 3]
                records.append((header, sequence, plus, quality))
                i += 4
            else:
                break
        else:
            i += 1
    return records


def _reverse_complement_quality(quality: str) -> str:
    """Reverse complement of quality string (just reverse, no complement needed)."""
    return quality[::-1]


def _merge_pair_with_quality(
    forward_seq: str,
    forward_qual: str,
    reverse_seq: str,
    reverse_qual: str,
    min_overlap: int,
) -> Tuple[str, str, int]:
    """Merge a forward/reverse pair with quality scores."""
    rc_reverse_seq = _reverse_complement(reverse_seq)
    rc_reverse_qual = _reverse_complement_quality(reverse_qual)
    
    max_overlap = 0
    merged_sequence = forward_seq + rc_reverse_seq
    merged_quality = forward_qual + rc_reverse_qual
    
    for overlap in range(min(len(forward_seq), len(rc_reverse_seq)), min_overlap - 1, -1):
        if forward_seq.endswith(rc_reverse_seq[:overlap]):
            merged_sequence = forward_seq + rc_reverse_seq[overlap:]
            # For quality in overlap region, take maximum
            overlap_qual_fwd = forward_qual[-overlap:]
            overlap_qual_rev = rc_reverse_qual[:overlap]
            merged_overlap_qual = "".join(
                max(q1, q2) for q1, q2 in zip(overlap_qual_fwd, overlap_qual_rev)
            )
            merged_quality = forward_qual[:-overlap] + merged_overlap_qual + rc_reverse_qual[overlap:]
            max_overlap = overlap
            break
    
    return merged_sequence, merged_quality, max_overlap


def merge_reads_from_s3(
    r1_path: str,
    r2_path: str,
    output_path: str,
    min_overlap: int = 12,
) -> Dict[str, object]:
    """
    Merge paired-end reads from S3 files and output merged FASTQ to S3.
    
    Args:
        r1_path: S3 path to R1 reads (e.g., "s3://bucket/key/R1.fq")
        r2_path: S3 path to R2 reads (e.g., "s3://bucket/key/R2.fq")
        output_path: S3 path for output merged reads (e.g., "s3://bucket/key/merged.fq")
        min_overlap: Minimum overlap length to merge reads.
    
    Returns:
        Dictionary containing status and summary metrics.
    """
    try:
        s3_client = boto3.client('s3')
        
        # Parse S3 paths
        r1_bucket, r1_key = _parse_s3_path(r1_path)
        r2_bucket, r2_key = _parse_s3_path(r2_path)
        output_bucket, output_key = _parse_s3_path(output_path)
        
        # Download files to temporary locations
        logger.info(f"Downloading R1 from s3://{r1_bucket}/{r1_key}...")
        tmp_r1_file = tempfile.NamedTemporaryFile(suffix='.fq', delete=False)
        tmp_r1_path = tmp_r1_file.name
        tmp_r1_file.close()
        s3_client.download_file(r1_bucket, r1_key, tmp_r1_path)
        
        logger.info(f"Downloading R2 from s3://{r2_bucket}/{r2_key}...")
        tmp_r2_file = tempfile.NamedTemporaryFile(suffix='.fq', delete=False)
        tmp_r2_path = tmp_r2_file.name
        tmp_r2_file.close()
        s3_client.download_file(r2_bucket, r2_key, tmp_r2_path)
        
        # Read FASTQ files
        logger.info("Reading FASTQ files...")
        with open(tmp_r1_path, 'r') as f:
            r1_content = f.read()
        with open(tmp_r2_path, 'r') as f:
            r2_content = f.read()
        
        # Parse FASTQ records
        r1_records = _parse_fastq_with_quality(r1_content)
        r2_records = _parse_fastq_with_quality(r2_content)
        
        if len(r1_records) != len(r2_records):
            logger.warning(f"R1 has {len(r1_records)} reads, R2 has {len(r2_records)} reads")
        
        # Merge reads
        logger.info(f"Merging {min(len(r1_records), len(r2_records))} read pairs...")
        merged_records = []
        overlaps = []
        
        for idx, ((r1_header, r1_seq, r1_plus, r1_qual), (r2_header, r2_seq, r2_plus, r2_qual)) in enumerate(
            zip(r1_records, r2_records)
        ):
            merged_seq, merged_qual, overlap = _merge_pair_with_quality(
                r1_seq, r1_qual, r2_seq, r2_qual, min_overlap
            )
            # Use R1 header but update to indicate merged read
            merged_header = r1_header.replace(" ", "_merged_", 1) if " " in r1_header else f"{r1_header}_merged"
            merged_records.append((merged_header, merged_seq, "+", merged_qual))
            overlaps.append(overlap)
        
        # Write merged FASTQ to temporary file
        logger.info("Writing merged FASTQ...")
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fq', delete=False) as tmp_output:
            tmp_output_path = tmp_output.name
            for header, seq, plus, qual in merged_records:
                tmp_output.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
        
        # Upload to S3
        logger.info(f"Uploading merged reads to s3://{output_bucket}/{output_key}...")
        s3_client.upload_file(tmp_output_path, output_bucket, output_key)
        
        # Clean up temporary files
        try:
            Path(tmp_r1_path).unlink()
            Path(tmp_r2_path).unlink()
            Path(tmp_output_path).unlink()
        except Exception as e:
            logger.warning(f"Failed to clean up temporary files: {e}")
        
        summary = {
            "total_pairs": len(merged_records),
            "merged_pairs": sum(1 for o in overlaps if o >= min_overlap),
            "min_overlap": min_overlap,
            "average_overlap": sum(overlaps) / len(overlaps) if overlaps else 0,
            "output_path": output_path,
        }
        
        logger.info(f"Successfully merged {len(merged_records)} read pairs")
        
        return {
            "status": "success",
            "text": f"Read merging completed successfully. Merged {len(merged_records)} read pairs.",
            "summary": summary,
            "output_path": output_path,
        }
        
    except Exception as e:
        logger.error(f"Error merging reads from S3: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e),
            "text": f"Failed to merge reads: {str(e)}",
        }


__all__ = ["run_read_merging_raw", "merge_reads_from_s3"]


