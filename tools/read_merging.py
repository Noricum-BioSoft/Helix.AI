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
    if not s3_path or not s3_path.startswith("s3://"):
        raise ValueError(f"Invalid S3 path: {s3_path!r}")
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


def _get_s3_file_size(s3_path: str) -> int:
    """Get file size in bytes from S3."""
    try:
        s3_client = boto3.client('s3')
        bucket, key = _parse_s3_path(s3_path)
        response = s3_client.head_object(Bucket=bucket, Key=key)
        return response['ContentLength']
    except Exception as e:
        logger.warning(f"Could not get file size for {s3_path}: {e}")
        return 0


def merge_reads_pyspark(
    r1_path: str,
    r2_path: str,
    output_path: str,
    min_overlap: int = 12,
) -> Dict[str, object]:
    """
    Merge paired-end reads using PySpark for distributed processing.
    
    This implementation reads directly from S3 (no local downloads), processes data
    in partitions across cluster nodes, and writes output directly to S3. This avoids
    disk space issues for large files.
    
    Args:
        r1_path: S3 path to R1 reads (e.g., "s3://bucket/key/R1.fq")
        r2_path: S3 path to R2 reads (e.g., "s3://bucket/key/R2.fq")
        output_path: S3 path for output merged reads (e.g., "s3://bucket/key/merged.fq")
        min_overlap: Minimum overlap length to merge reads.
    
    Returns:
        Dictionary containing status and summary metrics.
    """
    try:
        from pyspark.sql import SparkSession
        from pyspark import SparkContext
    except ImportError:
        return {
            "status": "error",
            "error": "PySpark is not available. Cannot use distributed processing.",
            "text": "PySpark is required for distributed read merging but is not installed.",
        }
    
    try:
        logger.info("Initializing PySpark session for distributed read merging...")
        
        # Convert s3:// to s3a:// for Spark's S3 filesystem
        r1_s3a_path = r1_path.replace("s3://", "s3a://")
        r2_s3a_path = r2_path.replace("s3://", "s3a://")
        output_s3a_path = output_path.replace("s3://", "s3a://")
        
        # Initialize Spark with S3A configuration
        spark = SparkSession.builder \
            .appName("Read Merging - Distributed") \
            .config("spark.sql.adaptive.enabled", "true") \
            .config("spark.sql.adaptive.coalescePartitions.enabled", "true") \
            .config("spark.hadoop.fs.s3a.impl", "org.apache.hadoop.fs.s3a.S3AFileSystem") \
            .config("spark.hadoop.fs.s3a.aws.credentials.provider", 
                    "com.amazonaws.auth.InstanceProfileCredentialsProvider") \
            .config("spark.hadoop.fs.s3a.path.style.access", "true") \
            .getOrCreate()
        
        sc = spark.sparkContext
        sc.setLogLevel("WARN")
        
        logger.info(f"Reading R1 from {r1_s3a_path}")
        logger.info(f"Reading R2 from {r2_s3a_path}")
        
        # Read FASTQ files line-by-line from S3
        # Spark automatically partitions large files
        r1_lines = sc.textFile(r1_s3a_path)
        r2_lines = sc.textFile(r2_s3a_path)
        
        # Add line numbers (index) to each line
        r1_indexed = r1_lines.zipWithIndex()
        r2_indexed = r2_lines.zipWithIndex()
        
        # Group lines by record index (4 lines per FASTQ record)
        # Record index = line_number // 4
        # Line type = line_number % 4 (0=header, 1=seq, 2=+, 3=qual)
        def group_by_record(line_with_index):
            """Convert (line, index) to (record_index, (line_type, line))."""
            line, line_idx = line_with_index
            record_idx = line_idx // 4
            line_type = line_idx % 4
            return (record_idx, (line_type, line.strip()))
        
        r1_grouped = r1_indexed.map(group_by_record)
        r2_grouped = r2_indexed.map(group_by_record)
        
        # Aggregate lines by record index to form complete FASTQ records
        def parse_record(group):
            """Parse grouped lines into a FASTQ record tuple."""
            record_idx, lines = group
            # Sort by line_type to ensure correct order (header, seq, plus, qual)
            lines_dict = dict(lines)
            if len(lines_dict) == 4:
                header = lines_dict[0].lstrip('@')  # Remove @ from header
                sequence = lines_dict[1]
                quality = lines_dict[3]
                return (record_idx, (header, sequence, quality))
            return None
        
        logger.info("Parsing FASTQ records...")
        r1_records = r1_grouped.groupByKey().map(parse_record).filter(lambda x: x is not None)
        r2_records = r2_grouped.groupByKey().map(parse_record).filter(lambda x: x is not None)
        
        # Join R1 and R2 records by index
        logger.info("Pairing R1 and R2 records...")
        paired_records = r1_records.join(r2_records)
        
        # Merge each pair
        def merge_record_pair(pair):
            """Merge a single R1/R2 pair."""
            (r1_header, r1_seq, r1_plus, r1_qual), (r2_header, r2_seq, r2_plus, r2_qual) = pair
            merged_seq, merged_qual, overlap = _merge_pair_with_quality(
                r1_seq, r1_qual, r2_seq, r2_qual, min_overlap
            )
            # Update header to indicate merged read
            merged_header = r1_header.replace(" ", "_merged_", 1) if " " in r1_header else f"{r1_header}_merged"
            return (merged_header, merged_seq, merged_qual, overlap)
        
        logger.info("Merging read pairs...")
        merged_records = paired_records.map(merge_record_pair)
        
        # Format as FASTQ
        def format_fastq(record):
            """Format record as FASTQ string."""
            header, seq, qual, overlap = record
            return f"@{header}\n{seq}\n+\n{qual}"
        
        formatted_output = merged_records.map(format_fastq)
        
        # Collect overlap statistics (need to collect before stopping Spark)
        overlaps_rdd = merged_records.map(lambda x: x[3])  # Extract overlap values
        overlaps_list = overlaps_rdd.collect()
        
        # Write output to S3 as partitioned files
        # Note: saveAsTextFile creates a directory with part files
        # We'll combine them into a single file after
        logger.info(f"Writing merged reads to {output_s3a_path}...")
        # Create a temporary directory path for Spark output
        output_bucket, output_key_base = _parse_s3_path(output_path)
        output_key_dir = output_key_base.rsplit('/', 1)[0] if '/' in output_key_base else ''
        output_dir = f"s3a://{output_bucket}/{output_key_dir}/_spark_temp_merged" if output_key_dir else f"s3a://{output_bucket}/_spark_temp_merged"
        
        formatted_output.saveAsTextFile(output_dir)
        
        # Combine part files into single output file
        # Use boto3 to read all part files and combine
        logger.info("Combining partition files into single output...")
        s3_client = boto3.client('s3')
        output_bucket, output_key = _parse_s3_path(output_path)
        
        # Extract prefix from output_dir (remove s3a:// or s3:// prefix)
        prefix = output_dir
        if prefix.startswith('s3a://'):
            prefix = prefix[6:]
        elif prefix.startswith('s3://'):
            prefix = prefix[5:]
        
        # Remove bucket name if present
        if prefix.startswith(f"{output_bucket}/"):
            prefix = prefix[len(f"{output_bucket}/"):]
        
        if not prefix.endswith('/'):
            prefix += '/'
        
        # List all part files
        paginator = s3_client.get_paginator('list_objects_v2')
        pages = paginator.paginate(Bucket=output_bucket, Prefix=prefix)
        
        part_files = []
        for page in pages:
            if 'Contents' in page:
                for obj in page['Contents']:
                    key = obj['Key']
                    if key.endswith('/_SUCCESS'):
                        continue
                    if 'part-' in key or key.endswith('.fq') or key.endswith('.fastq'):
                        part_files.append(key)
        
        # Sort part files by name (they should be ordered)
        part_files.sort()
        
        # Combine into single file
        if part_files:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fq', delete=False) as tmp_output:
                tmp_output_path = tmp_output.name
                for part_key in part_files:
                    try:
                        obj = s3_client.get_object(Bucket=output_bucket, Key=part_key)
                        content = obj['Body'].read().decode('utf-8')
                        tmp_output.write(content)
                    except Exception as e:
                        logger.warning(f"Error reading part file {part_key}: {e}")
            
            # Upload combined file
            logger.info(f"Uploading combined file to s3://{output_bucket}/{output_key}...")
            s3_client.upload_file(tmp_output_path, output_bucket, output_key)
            
            # Clean up part files and directory
            logger.info("Cleaning up partition files...")
            for part_key in part_files:
                try:
                    s3_client.delete_object(Bucket=output_bucket, Key=part_key)
                except Exception as e:
                    logger.warning(f"Error deleting part file {part_key}: {e}")
            
            # Clean up directory marker and _SUCCESS file if they exist
            try:
                s3_client.delete_object(Bucket=output_bucket, Key=f"{prefix}_SUCCESS")
                s3_client.delete_object(Bucket=output_bucket, Key=prefix.rstrip('/'))
            except:
                pass
            
            # Clean up temp file
            Path(tmp_output_path).unlink()
        
        total_pairs = len(overlaps_list)
        merged_pairs = sum(1 for o in overlaps_list if o >= min_overlap)
        average_overlap = sum(overlaps_list) / len(overlaps_list) if overlaps_list else 0
        
        summary = {
            "total_pairs": total_pairs,
            "merged_pairs": merged_pairs,
            "min_overlap": min_overlap,
            "average_overlap": average_overlap,
            "output_path": output_path,
            "processing_mode": "pyspark_distributed",
        }
        
        logger.info(f"Successfully merged {total_pairs} read pairs using PySpark")
        
        spark.stop()
        
        return {
            "status": "success",
            "text": f"Read merging completed successfully using distributed PySpark processing. Merged {total_pairs} read pairs.",
            "summary": summary,
            "output_path": output_path,
        }
        
    except Exception as e:
        logger.error(f"Error in PySpark read merging: {e}", exc_info=True)
        try:
            spark.stop()
        except:
            pass
        return {
            "status": "error",
            "error": str(e),
            "text": f"Failed to merge reads using PySpark: {str(e)}",
        }


def merge_reads_from_s3(
    r1_path: str,
    r2_path: str,
    output_path: str,
    min_overlap: int = 12,
    use_pyspark_threshold_gb: float = 5.0,
) -> Dict[str, object]:
    """
    Merge paired-end reads from S3 files and output merged FASTQ to S3.
    
    For large files (>5GB by default), automatically uses PySpark distributed processing
    to avoid disk space issues. For smaller files, uses single-node processing.
    
    Args:
        r1_path: S3 path to R1 reads (e.g., "s3://bucket/key/R1.fq")
        r2_path: S3 path to R2 reads (e.g., "s3://bucket/key/R2.fq")
        output_path: S3 path for output merged reads (e.g., "s3://bucket/key/merged.fq")
        min_overlap: Minimum overlap length to merge reads.
        use_pyspark_threshold_gb: Use PySpark if total input size exceeds this threshold (default: 5.0 GB)
    
    Returns:
        Dictionary containing status and summary metrics.
    """
    # Check file sizes to decide on processing method
    try:
        r1_size = _get_s3_file_size(r1_path)
        r2_size = _get_s3_file_size(r2_path)
        total_size_gb = (r1_size + r2_size) / (1024 ** 3)
        
        logger.info(f"R1 size: {r1_size / (1024**2):.2f} MB, R2 size: {r2_size / (1024**2):.2f} MB, Total: {total_size_gb:.2f} GB")
        
        # Use PySpark for large files
        if total_size_gb > use_pyspark_threshold_gb:
            logger.info(f"Total file size ({total_size_gb:.2f} GB) exceeds threshold ({use_pyspark_threshold_gb} GB). Using PySpark distributed processing.")
            try:
                from pyspark.sql import SparkSession
                # PySpark is available, use it
                return merge_reads_pyspark(r1_path, r2_path, output_path, min_overlap)
            except ImportError:
                logger.warning("PySpark not available, falling back to single-node processing despite large file size")
                # Fall through to single-node processing
        
        # Use single-node processing for smaller files or if PySpark unavailable
        logger.info("Using single-node processing (file size below threshold or PySpark unavailable)")
        
    except Exception as e:
        logger.warning(f"Could not determine file sizes, using single-node processing: {e}")
        # Fall through to single-node processing
    
    # Single-node processing (original implementation)
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


__all__ = ["run_read_merging_raw", "merge_reads_from_s3", "merge_reads_pyspark"]


