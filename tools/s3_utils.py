"""
Utility functions for accessing files from S3 in tools.

This module provides helpers for tools to access files that are either:
1. Uploaded directly to session bucket (small files)
2. Referenced from dataset bucket (large files)
"""

import boto3
import tempfile
from pathlib import Path
from typing import Dict, Any, List, Optional
import logging

logger = logging.getLogger(__name__)

def get_s3_client():
    """Get or create S3 client."""
    try:
        return boto3.client('s3')
    except Exception as e:
        logger.error(f"Failed to create S3 client: {e}")
        return None

def download_from_s3(s3_bucket: str, s3_key: str, local_path: Optional[Path] = None) -> Path:
    """Download a file from S3 to a local path.
    
    Args:
        s3_bucket: S3 bucket name
        s3_key: S3 object key
        local_path: Optional local path. If None, creates a temp file.
    
    Returns:
        Path to the downloaded file
    """
    s3_client = get_s3_client()
    if s3_client is None:
        raise RuntimeError("S3 client not available")
    
    if local_path is None:
        # Create temp file with same extension
        suffix = Path(s3_key).suffix
        local_path = Path(tempfile.mktemp(suffix=suffix))
    
    local_path = Path(local_path)
    local_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        s3_client.download_file(s3_bucket, s3_key, str(local_path))
        logger.info(f"Downloaded {s3_key} from {s3_bucket} to {local_path}")
        return local_path
    except Exception as e:
        logger.error(f"Failed to download {s3_key} from {s3_bucket}: {e}")
        raise

def get_session_files_from_context(session_context: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Get all files available to a session from session context.
    
    This function checks both:
    1. uploaded_files (small files in session bucket)
    2. dataset_references (large files in dataset bucket)
    
    Args:
        session_context: Session context dictionary
    
    Returns:
        List of file dictionaries with s3_bucket, s3_key, filename, size, type
    """
    files = []
    metadata = session_context.get("metadata", {})
    
    # Add small uploaded files
    for file_info in metadata.get("uploaded_files", []):
        files.append({
            "type": "uploaded",
            "s3_bucket": file_info.get("s3_bucket", ""),
            "s3_key": file_info.get("s3_key", ""),
            "filename": file_info.get("filename", file_info.get("name", "")),
            "size": file_info.get("size", 0)
        })
    
    # Add large dataset references
    for ref in metadata.get("dataset_references", []):
        files.append({
            "type": "dataset_reference",
            "s3_bucket": ref.get("s3_bucket", ""),
            "s3_key": ref.get("s3_key", ""),
            "filename": ref.get("filename", ""),
            "size": ref.get("size", 0),
            "dataset_id": ref.get("dataset_id", "")
        })
    
    return files

def find_files_by_pattern(files: List[Dict[str, Any]], pattern: str) -> List[Dict[str, Any]]:
    """Find files matching a pattern in filename.
    
    Args:
        files: List of file dictionaries
        pattern: Pattern to search for (e.g., "R1", ".fastq", "mate")
    
    Returns:
        List of matching file dictionaries
    """
    pattern_lower = pattern.lower()
    matches = []
    for file_info in files:
        filename_lower = file_info.get("filename", "").lower()
        if pattern_lower in filename_lower:
            matches.append(file_info)
    return matches

def download_session_file(session_context: Dict[str, Any], filename_pattern: str, 
                         local_path: Optional[Path] = None) -> Optional[Path]:
    """Download a file from session by filename pattern.
    
    Args:
        session_context: Session context dictionary
        filename_pattern: Pattern to match filename (e.g., "R1", "mate_R1.fq")
        local_path: Optional local path for downloaded file
    
    Returns:
        Path to downloaded file, or None if not found
    """
    files = get_session_files_from_context(session_context)
    matches = find_files_by_pattern(files, filename_pattern)
    
    if not matches:
        logger.warning(f"No files found matching pattern: {filename_pattern}")
        return None
    
    # Use first match
    file_info = matches[0]
    
    try:
        return download_from_s3(
            file_info["s3_bucket"],
            file_info["s3_key"],
            local_path
        )
    except Exception as e:
        logger.error(f"Failed to download file {file_info['filename']}: {e}")
        return None

def get_paired_end_reads(session_context: Dict[str, Any], 
                        r1_pattern: str = "R1", 
                        r2_pattern: str = "R2") -> Optional[Dict[str, Path]]:
    """Get paired-end read files (R1 and R2) from session.
    
    Args:
        session_context: Session context dictionary
        r1_pattern: Pattern to identify R1 file (default: "R1")
        r2_pattern: Pattern to identify R2 file (default: "R2")
    
    Returns:
        Dictionary with "r1" and "r2" keys pointing to local Path objects,
        or None if files not found
    """
    files = get_session_files_from_context(session_context)
    
    r1_files = find_files_by_pattern(files, r1_pattern)
    r2_files = find_files_by_pattern(files, r2_pattern)
    
    if not r1_files or not r2_files:
        logger.warning(f"Could not find paired-end reads: R1={len(r1_files)}, R2={len(r2_files)}")
        return None
    
    try:
        r1_path = download_from_s3(r1_files[0]["s3_bucket"], r1_files[0]["s3_key"])
        r2_path = download_from_s3(r2_files[0]["s3_bucket"], r2_files[0]["s3_key"])
        
        return {
            "r1": r1_path,
            "r2": r2_path,
            "r1_info": r1_files[0],
            "r2_info": r2_files[0]
        }
    except Exception as e:
        logger.error(f"Failed to download paired-end reads: {e}")
        return None