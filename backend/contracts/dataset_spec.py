"""
Dataset Specification Contract - Describes data inputs with uncertainty tracking.

This contract is used to pass dataset information to the Infrastructure Decision Agent,
with explicit tracking of:
- Known file sizes
- Unknown file sizes (due to permissions, missing files, etc.)
- Location types (S3, Local, URL)
- Confidence in size estimates

Key principle: Uncertainty is tracked explicitly, not hidden.
"""

from typing import List, Optional, Literal
from pydantic import BaseModel, Field, field_validator, ConfigDict


class FileSpec(BaseModel):
    """Specification of a single file with size and location uncertainty."""
    
    uri: str = Field(
        ...,
        description="URI of the file (s3://bucket/key, /local/path, etc.)"
    )
    
    # Size information with uncertainty
    size_bytes: Optional[int] = Field(
        default=None,
        ge=0,
        description="Size in bytes (None if unknown)"
    )
    size_confidence: float = Field(
        default=1.0,
        ge=0.0,
        le=1.0,
        description="Confidence in size estimate (1.0 = known exactly, 0.0 = unknown)"
    )
    size_source: Literal["head_object", "stat", "user_provided", "estimated", "unknown"] = Field(
        default="unknown",
        description="How the size was determined"
    )
    
    # Location information
    location_type: Literal["S3", "Local", "URL", "Unknown"] = Field(
        default="Unknown",
        description="Type of storage location"
    )
    
    # Access information
    accessible: bool = Field(
        default=True,
        description="Whether the file is accessible (permissions, existence)"
    )
    access_error: Optional[str] = Field(
        default=None,
        description="Error message if file is not accessible"
    )
    
    # Metadata
    file_format: Optional[str] = Field(
        default=None,
        description="File format (e.g., 'fastq', 'bam', 'vcf')"
    )
    compressed: bool = Field(
        default=False,
        description="Whether the file is compressed (affects processing)"
    )
    
    @field_validator('location_type', mode='before')
    @classmethod
    def infer_location_type(cls, v, info):
        """Infer location_type from URI if not provided."""
        if v == "Unknown" and 'uri' in info.data:
            uri = info.data['uri']
            if uri.startswith('s3://'):
                return "S3"
            elif uri.startswith('http://') or uri.startswith('https://'):
                return "URL"
            elif uri.startswith('/') or uri.startswith('./') or uri.startswith('../'):
                return "Local"
        return v
    
    @field_validator('size_confidence', mode='before')
    @classmethod
    def compute_size_confidence(cls, v, info):
        """Compute size_confidence from size_source if not explicitly provided."""
        if v == 1.0 and 'size_source' in info.data:
            source = info.data['size_source']
            if source == "unknown":
                return 0.0
            elif source == "estimated":
                return 0.5
            elif source == "user_provided":
                return 0.7
            elif source in ["head_object", "stat"]:
                return 1.0
        return v
    
    @field_validator('accessible', mode='before')
    @classmethod
    def check_accessible(cls, v, info):
        """Mark as not accessible if size is unknown and error present."""
        if 'size_bytes' in info.data and 'access_error' in info.data:
            if info.data['size_bytes'] is None and info.data['access_error'] is not None:
                return False
        return v


class DatasetSpec(BaseModel):
    """
    Dataset Specification - describes a set of input files with uncertainty.
    
    Used by Infrastructure Decision Agent to make recommendations based on:
    - Total data volume (sum of known file sizes)
    - Number of unknown files (lowers confidence)
    - Location distribution (S3 vs Local vs mixed)
    - Access issues (permissions, missing files)
    
    Key principle: Uncertainty is first-class - we track what we don't know.
    """
    
    # File specifications
    files: List[FileSpec] = Field(
        default_factory=list,
        description="List of file specifications"
    )
    
    # Computed statistics (from files)
    total_size_bytes: int = Field(
        default=0,
        ge=0,
        description="Total size of all known files in bytes"
    )
    total_size_mb: float = Field(
        default=0.0,
        ge=0.0,
        description="Total size in MB"
    )
    file_count: int = Field(
        default=0,
        ge=0,
        description="Total number of files"
    )
    known_size_count: int = Field(
        default=0,
        ge=0,
        description="Number of files with known sizes"
    )
    unknown_size_count: int = Field(
        default=0,
        ge=0,
        description="Number of files with unknown sizes"
    )
    inaccessible_count: int = Field(
        default=0,
        ge=0,
        description="Number of inaccessible files"
    )
    
    # Location distribution
    locations: List[str] = Field(
        default_factory=list,
        description="Unique location types in this dataset"
    )
    
    # Confidence
    overall_confidence: float = Field(
        default=1.0,
        ge=0.0,
        le=1.0,
        description="Overall confidence in dataset specs (based on unknown counts)"
    )
    
    @field_validator('file_count', mode='after')
    @classmethod
    def compute_file_count(cls, v, info):
        """Compute file_count from files list."""
        if 'files' in info.data:
            return len(info.data['files'])
        return v
    
    @field_validator('total_size_bytes', mode='after')
    @classmethod
    def compute_total_size_bytes(cls, v, info):
        """Compute total_size_bytes from files list."""
        if 'files' in info.data:
            return sum(f.size_bytes for f in info.data['files'] if f.size_bytes is not None)
        return v
    
    @field_validator('total_size_mb', mode='after')
    @classmethod
    def compute_total_size_mb(cls, v, info):
        """Compute total_size_mb from files list."""
        if 'files' in info.data:
            total_bytes = sum(f.size_bytes for f in info.data['files'] if f.size_bytes is not None)
            return total_bytes / (1024 * 1024) if total_bytes > 0 else 0.0
        return v
    
    @field_validator('known_size_count', mode='after')
    @classmethod
    def compute_known_size_count(cls, v, info):
        """Compute known_size_count from files list."""
        if 'files' in info.data:
            return sum(1 for f in info.data['files'] if f.size_bytes is not None)
        return v
    
    @field_validator('unknown_size_count', mode='after')
    @classmethod
    def compute_unknown_size_count(cls, v, info):
        """Compute unknown_size_count from files list."""
        if 'files' in info.data:
            return sum(1 for f in info.data['files'] if f.size_bytes is None)
        return v
    
    @field_validator('inaccessible_count', mode='after')
    @classmethod
    def compute_inaccessible_count(cls, v, info):
        """Compute inaccessible_count from files list."""
        if 'files' in info.data:
            return sum(1 for f in info.data['files'] if not f.accessible)
        return v
    
    @field_validator('locations', mode='after')
    @classmethod
    def compute_locations(cls, v, info):
        """Compute locations from files list."""
        if 'files' in info.data:
            return list(set(f.location_type for f in info.data['files']))
        return v
    
    @field_validator('overall_confidence', mode='after')
    @classmethod
    def compute_overall_confidence(cls, v, info):
        """Compute overall_confidence from files list."""
        if 'files' not in info.data:
            return v
        
        files = info.data['files']
        if not files:
            return 0.0
        
        # Confidence is average of file confidences
        avg_confidence = sum(f.size_confidence for f in files) / len(files)
        # Penalize for inaccessible files
        inaccessible = sum(1 for f in files if not f.accessible)
        if inaccessible > 0:
            penalty = inaccessible / len(files) * 0.5  # Up to 50% penalty
            avg_confidence = max(0.0, avg_confidence - penalty)
        return round(avg_confidence, 2)
    
    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "files": [
                    {
                        "uri": "s3://my-bucket/sample_R1.fq",
                        "size_bytes": 250000000,
                        "size_confidence": 1.0,
                        "size_source": "head_object",
                        "location_type": "S3",
                        "accessible": True,
                        "file_format": "fastq",
                        "compressed": False
                    },
                    {
                        "uri": "s3://my-bucket/sample_R2.fq",
                        "size_bytes": 240000000,
                        "size_confidence": 1.0,
                        "size_source": "head_object",
                        "location_type": "S3",
                        "accessible": True,
                        "file_format": "fastq",
                        "compressed": False
                    },
                    {
                        "uri": "s3://my-bucket/reference.fa",
                        "size_bytes": None,
                        "size_confidence": 0.0,
                        "size_source": "unknown",
                        "location_type": "S3",
                        "accessible": False,
                        "access_error": "Access denied (403)",
                        "file_format": "fasta",
                        "compressed": False
                    }
                ],
                "total_size_bytes": 490000000,
                "total_size_mb": 467.3,
                "file_count": 3,
                "known_size_count": 2,
                "unknown_size_count": 1,
                "inaccessible_count": 1,
                "locations": ["S3"],
                "overall_confidence": 0.5
            }
        }
    )
