"""
File Metadata Inspector - Read-only tool for file size and location discovery.

This tool provides factual file metadata to the Infrastructure Decision Agent,
reducing reliance on LLM hallucination or user-provided estimates.

Capabilities:
- S3 object sizes via head_object (respects bucket regions)
- Local file sizes via stat
- Mock catalog mode for testing (no actual S3/filesystem access)
- Batch processing of multiple files
- Graceful handling of permissions errors

Key principle: Always return confidence level with size estimates.
"""

import os
import json
import logging
from pathlib import Path
from typing import Optional, List, Dict, Any, Literal
from dataclasses import dataclass

from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)


class FileMetadata(BaseModel):
    """Metadata for a single file."""
    
    uri: str = Field(..., description="File URI (s3://bucket/key or local path)")
    
    size_bytes: Optional[int] = Field(
        default=None,
        ge=0,
        description="Size in bytes (None if unknown)"
    )
    
    size_confidence: float = Field(
        default=0.0,
        ge=0.0,
        le=1.0,
        description="Confidence in size (1.0=known, 0.0=unknown)"
    )
    
    source: Literal["head_object", "stat", "mock_catalog", "unknown"] = Field(
        default="unknown",
        description="How the size was determined"
    )
    
    location_type: Literal["S3", "Local", "Unknown"] = Field(
        default="Unknown",
        description="Storage location type"
    )
    
    accessible: bool = Field(
        default=True,
        description="Whether file is accessible (permissions, existence)"
    )
    
    error: Optional[str] = Field(
        default=None,
        description="Error message if inaccessible"
    )
    
    region: Optional[str] = Field(
        default=None,
        description="AWS region (for S3 files only)"
    )


class FileMetadataInspector:
    """
    Read-only tool for inspecting file metadata (sizes, locations, accessibility).
    
    Usage:
        inspector = FileMetadataInspector(mock_mode=False)
        metadata = inspector.inspect_files([
            "s3://my-bucket/file1.fq",
            "/local/path/file2.fq"
        ])
    
    Mock mode:
        inspector = FileMetadataInspector(
            mock_mode=True,
            mock_catalog={"s3://my-bucket/file1.fq": 1000000}
        )
        # Uses mock catalog instead of real S3/filesystem calls
    """
    
    def __init__(
        self,
        mock_mode: bool = False,
        mock_catalog: Optional[Dict[str, int]] = None,
        default_region: str = "us-east-1"
    ):
        """
        Initialize FileMetadataInspector.
        
        Args:
            mock_mode: If True, use mock_catalog instead of real filesystem/S3
            mock_catalog: Dict mapping URIs to sizes (for testing)
            default_region: Default AWS region for S3 operations
        """
        self.mock_mode = mock_mode
        self.mock_catalog = mock_catalog or {}
        self.default_region = os.getenv("AWS_REGION", os.getenv("AWS_DEFAULT_REGION", default_region))
        
        # Cache for S3 clients per region (avoid recreating)
        self._s3_clients: Dict[str, Any] = {}
    
    def inspect_files(self, uris: List[str]) -> List[FileMetadata]:
        """
        Inspect multiple files and return metadata.
        
        Args:
            uris: List of file URIs (s3://bucket/key or local paths)
        
        Returns:
            List of FileMetadata objects with sizes and confidence
        """
        results = []
        for uri in uris:
            metadata = self.inspect_file(uri)
            results.append(metadata)
        return results
    
    def inspect_file(self, uri: str) -> FileMetadata:
        """
        Inspect a single file and return metadata.
        
        Args:
            uri: File URI (s3://bucket/key or local path)
        
        Returns:
            FileMetadata with size and confidence
        """
        # Mock mode: use catalog
        if self.mock_mode:
            return self._inspect_mock(uri)
        
        # Determine location type
        if uri.startswith("s3://"):
            return self._inspect_s3(uri)
        elif uri.startswith("/") or uri.startswith("./") or uri.startswith("../"):
            return self._inspect_local(uri)
        else:
            return FileMetadata(
                uri=uri,
                location_type="Unknown",
                accessible=False,
                error="Unknown URI format (not S3 or local path)"
            )
    
    def _inspect_mock(self, uri: str) -> FileMetadata:
        """Inspect file using mock catalog."""
        if uri in self.mock_catalog:
            size = self.mock_catalog[uri]
            location_type = "S3" if uri.startswith("s3://") else "Local"
            return FileMetadata(
                uri=uri,
                size_bytes=size,
                size_confidence=1.0,
                source="mock_catalog",
                location_type=location_type,
                accessible=True
            )
        else:
            return FileMetadata(
                uri=uri,
                location_type="Unknown",
                accessible=False,
                error="File not in mock catalog"
            )
    
    def _inspect_s3(self, uri: str) -> FileMetadata:
        """
        Inspect S3 object using head_object.
        
        Handles:
        - Bucket region detection
        - Permission errors (403)
        - Missing objects (404)
        - Invalid URIs
        """
        try:
            import boto3
            from botocore.exceptions import ClientError
        except ImportError:
            return FileMetadata(
                uri=uri,
                location_type="S3",
                accessible=False,
                error="boto3 not available"
            )
        
        # Parse S3 URI
        if not uri.startswith("s3://"):
            return FileMetadata(
                uri=uri,
                location_type="S3",
                accessible=False,
                error="Invalid S3 URI"
            )
        
        parts = uri[5:].split("/", 1)
        if len(parts) != 2:
            return FileMetadata(
                uri=uri,
                location_type="S3",
                accessible=False,
                error="Invalid S3 URI format (missing key)"
            )
        
        bucket, key = parts
        
        try:
            # Try default region first
            region = self.default_region
            s3_client = self._get_s3_client(region)
            
            try:
                response = s3_client.head_object(Bucket=bucket, Key=key)
                size = response.get("ContentLength")
                
                if size is not None:
                    return FileMetadata(
                        uri=uri,
                        size_bytes=size,
                        size_confidence=1.0,
                        source="head_object",
                        location_type="S3",
                        accessible=True,
                        region=region
                    )
            
            except ClientError as e:
                error_code = e.response.get("Error", {}).get("Code", "")
                
                # Try to detect bucket region and retry
                if error_code in ["PermanentRedirect", "301"]:
                    try:
                        bucket_region = self._get_bucket_region(bucket, s3_client)
                        if bucket_region and bucket_region != region:
                            logger.debug(f"Bucket {bucket} is in {bucket_region}, retrying")
                            s3_client = self._get_s3_client(bucket_region)
                            response = s3_client.head_object(Bucket=bucket, Key=key)
                            size = response.get("ContentLength")
                            
                            if size is not None:
                                return FileMetadata(
                                    uri=uri,
                                    size_bytes=size,
                                    size_confidence=1.0,
                                    source="head_object",
                                    location_type="S3",
                                    accessible=True,
                                    region=bucket_region
                                )
                    except Exception:
                        pass  # Fall through to error handling
                
                # Handle specific error codes
                if error_code in ["403", "AccessDenied", "Forbidden"]:
                    return FileMetadata(
                        uri=uri,
                        location_type="S3",
                        accessible=False,
                        error=f"Access denied (check IAM permissions)",
                        region=region
                    )
                elif error_code in ["404", "NoSuchKey"]:
                    return FileMetadata(
                        uri=uri,
                        location_type="S3",
                        accessible=False,
                        error="Object not found",
                        region=region
                    )
                else:
                    return FileMetadata(
                        uri=uri,
                        location_type="S3",
                        accessible=False,
                        error=f"S3 error: {error_code}",
                        region=region
                    )
        
        except Exception as e:
            error_type = type(e).__name__
            if "NoCredentials" in error_type or "Credentials" in str(e):
                error_msg = "AWS credentials not configured"
            else:
                error_msg = f"{error_type}: {str(e)[:100]}"
            
            return FileMetadata(
                uri=uri,
                location_type="S3",
                accessible=False,
                error=error_msg
            )
    
    def _inspect_local(self, uri: str) -> FileMetadata:
        """
        Inspect local file using stat.
        
        Handles:
        - File not found
        - Permission errors
        - Invalid paths
        """
        try:
            path = Path(uri).expanduser().resolve()
            
            if not path.exists():
                return FileMetadata(
                    uri=uri,
                    location_type="Local",
                    accessible=False,
                    error="File not found"
                )
            
            size = path.stat().st_size
            return FileMetadata(
                uri=uri,
                size_bytes=size,
                size_confidence=1.0,
                source="stat",
                location_type="Local",
                accessible=True
            )
        
        except PermissionError:
            return FileMetadata(
                uri=uri,
                location_type="Local",
                accessible=False,
                error="Permission denied"
            )
        except Exception as e:
            return FileMetadata(
                uri=uri,
                location_type="Local",
                accessible=False,
                error=f"{type(e).__name__}: {str(e)[:100]}"
            )
    
    def _get_s3_client(self, region: str):
        """Get or create S3 client for region (cached)."""
        if region not in self._s3_clients:
            import boto3
            self._s3_clients[region] = boto3.client("s3", region_name=region)
        return self._s3_clients[region]
    
    def _get_bucket_region(self, bucket: str, s3_client) -> Optional[str]:
        """Get bucket region."""
        try:
            response = s3_client.get_bucket_location(Bucket=bucket)
            region = response.get("LocationConstraint")
            # get_bucket_location returns None for us-east-1
            return region or "us-east-1"
        except Exception:
            return None
    
    @classmethod
    def load_mock_catalog_from_file(cls, filepath: str) -> "FileMetadataInspector":
        """
        Load mock catalog from JSON file and create inspector in mock mode.
        
        JSON format:
        {
            "s3://bucket/file1.fq": 1000000,
            "/local/file2.fq": 500000,
            ...
        }
        
        Args:
            filepath: Path to JSON catalog file
        
        Returns:
            FileMetadataInspector in mock mode
        """
        with open(filepath, "r") as f:
            catalog = json.load(f)
        
        return cls(mock_mode=True, mock_catalog=catalog)
