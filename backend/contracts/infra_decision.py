"""
Infrastructure Decision Contract - Pydantic models for infrastructure recommendations.

This contract defines the output of the Infrastructure Decision Agent, which
recommends optimal execution environments (Local/EC2/EMR/Batch/Lambda) based on
data characteristics and computational requirements.

Key improvements over previous dataclass implementation:
- Strict Pydantic validation
- Confidence scoring (0-1)
- Cost ranges instead of exact values
- Uncertainty tracking
- Automatic JSON serialization
"""

from typing import List, Optional, Tuple, Literal
from pydantic import BaseModel, Field, field_validator, model_validator, ConfigDict
from pydantic_core import ValidationError


class FileAnalysis(BaseModel):
    """Analysis of input files for infrastructure decision."""
    
    total_size_bytes: int = Field(
        ...,
        ge=0,
        description="Total size of all known input files in bytes"
    )
    total_size_mb: float = Field(
        ...,
        ge=0.0,
        description="Total size in megabytes (computed from bytes)"
    )
    file_count: int = Field(
        ...,
        ge=0,
        description="Number of input files"
    )
    unknown_sizes: int = Field(
        default=0,
        ge=0,
        description="Number of files with unknown sizes (permissions/not found)"
    )
    locations: List[str] = Field(
        default_factory=list,
        description="File locations: S3, Local, Unknown, etc."
    )
    largest_file_mb: float = Field(
        default=0.0,
        ge=0.0,
        description="Size of largest file in megabytes"
    )

    @property
    def all_local(self) -> bool:
        """True if all analyzed inputs are local paths."""
        return bool(self.locations) and all(loc == "Local" for loc in self.locations)

    @property
    def all_s3(self) -> bool:
        """True if all analyzed inputs are S3 URIs."""
        return bool(self.locations) and all(loc == "S3" for loc in self.locations)

    @property
    def all_in_s3(self) -> bool:
        """Back-compat alias."""
        return self.all_s3

    @property
    def mixed_locations(self) -> bool:
        """True if inputs span multiple location types (e.g., Local + S3)."""
        return len(set(self.locations or [])) > 1

    @property
    def has_unknown_sizes(self) -> bool:
        return self.unknown_sizes > 0
    
    @field_validator('total_size_mb', mode='before')
    @classmethod
    def compute_total_size_mb(cls, v, info):
        """Compute total_size_mb from total_size_bytes if not provided."""
        # In Pydantic V2, use info.data to access other field values
        if v is None or v == 0.0:
            if 'total_size_bytes' in info.data:
                return info.data['total_size_bytes'] / (1024 * 1024)
        return v


class ComputationalRequirements(BaseModel):
    """Estimated computational requirements for the operation."""
    
    estimated_cpu_cores: int = Field(
        default=1,
        ge=1,
        description="Estimated CPU cores needed"
    )
    estimated_memory_gb: float = Field(
        default=2.0,
        ge=0.1,
        description="Estimated memory in GB"
    )
    estimated_runtime_minutes: float = Field(
        default=5.0,
        ge=0.1,
        description="Estimated runtime in minutes"
    )
    parallelizable: bool = Field(
        default=False,
        description="Whether the operation can be parallelized/distributed"
    )
    gpu_required: bool = Field(
        default=False,
        description="Whether GPU acceleration is required"
    )


class CostAnalysis(BaseModel):
    """Cost analysis with ranges and confidence, not exact values.
    
    Important: Costs should be ranges with explicit assumptions, not fake precision.
    """
    
    estimated_cost_range_usd: Tuple[float, float] = Field(
        ...,
        description="Cost range in USD (min, max) based on assumptions"
    )
    cost_assumptions: str = Field(
        ...,
        description="Explicit assumptions used for cost estimation (region, instance type, etc.)"
    )
    cost_confidence: float = Field(
        ...,
        ge=0.0,
        le=1.0,
        description="Confidence in cost estimate (0-1). Low if many unknowns."
    )
    data_transfer_cost_usd: Optional[float] = Field(
        default=None,
        ge=0.0,
        description="Estimated data transfer cost (if applicable)"
    )
    breakdown: Optional[dict] = Field(
        default=None,
        description="Optional cost breakdown by component (compute, storage, transfer)"
    )

    @property
    def cost_class(self) -> Literal["Free", "Low", "Medium", "High"]:
        """Coarse cost class derived from the estimated max cost."""
        _, max_cost = self.estimated_cost_range_usd
        if max_cost <= 0:
            return "Free"
        if max_cost <= 2:
            return "Low"
        if max_cost <= 10:
            return "Medium"
        return "High"
    
    @field_validator('estimated_cost_range_usd')
    @classmethod
    def validate_cost_range(cls, v):
        """Ensure cost range is valid (min <= max)."""
        min_cost, max_cost = v
        if min_cost < 0 or max_cost < 0:
            raise ValueError("Costs cannot be negative")
        if min_cost > max_cost:
            raise ValueError(f"Cost range invalid: min ({min_cost}) > max ({max_cost})")
        return v


class InfraAlternative(BaseModel):
    """An alternative infrastructure option with tradeoffs."""
    
    infrastructure: Literal["Local", "EC2", "EMR", "Batch", "Lambda"] = Field(
        ...,
        description="Alternative infrastructure option"
    )
    reasoning: str = Field(
        ...,
        min_length=10,
        description="Why this alternative could work"
    )
    tradeoffs: str = Field(
        ...,
        min_length=10,
        description="What you'd gain/lose compared to primary recommendation"
    )
    confidence: float = Field(
        ...,
        ge=0.0,
        le=1.0,
        description="Confidence in this alternative (0-1)"
    )


class InfraDecision(BaseModel):
    """
    Infrastructure Decision - output of Infrastructure Decision Agent.
    
    This contract replaces the previous dataclass implementation with strict Pydantic
    validation, confidence scoring, and cost range estimation.
    
    Key principles:
    - confidence_score: Always present, 0-1, reflects uncertainty
    - warnings: Populated when confidence is low or unknowns exist
    - Cost ranges (not exact values) with explicit assumptions
    - decision_summary: High-level justification (1-2 sentences)
    """
    
    # Core decision
    infrastructure: Literal["Local", "EC2", "EMR", "Batch", "Lambda"] = Field(
        ...,
        description="Recommended execution infrastructure"
    )
    
    # Confidence and summary
    confidence_score: float = Field(
        ...,
        ge=0.0,
        le=1.0,
        description="Confidence in this recommendation (0-1). <0.5 suggests human review."
    )
    decision_summary: str = Field(
        ...,
        min_length=20,
        description="High-level justification for this decision (1-2 sentences)"
    )
    
    # Detailed reasoning
    reasoning: str = Field(
        ...,
        min_length=50,
        description="Detailed explanation of why this infrastructure was chosen"
    )
    
    # Supporting analysis
    file_analysis: FileAnalysis = Field(
        ...,
        description="Analysis of input files (sizes, locations, counts)"
    )
    computational_requirements: ComputationalRequirements = Field(
        ...,
        description="Estimated computational needs (CPU, memory, runtime)"
    )
    cost_analysis: CostAnalysis = Field(
        ...,
        description="Cost analysis with ranges and assumptions"
    )
    
    # Alternatives and warnings
    alternatives: List[InfraAlternative] = Field(
        default_factory=list,
        description="Alternative infrastructure options with tradeoffs"
    )
    warnings: List[str] = Field(
        default_factory=list,
        description="Warnings about decision (unknown sizes, missing data, etc.)"
    )
    
    # Metadata
    inputs_analyzed: int = Field(
        default=0,
        ge=0,
        description="Number of input files analyzed"
    )
    
    @model_validator(mode='after')
    def check_warnings_for_low_confidence(self):
        """Ensure warnings are populated when confidence is low."""
        if self.confidence_score < 0.5 and not self.warnings:
            # Low confidence without warnings - add a default warning
            self.warnings = ["Low confidence in recommendation - consider human review"]
        return self
    
    @field_validator('decision_summary', 'reasoning')
    @classmethod
    def validate_text_not_empty(cls, v):
        """Ensure text fields are not just whitespace."""
        if not v or not v.strip():
            raise ValueError("Field cannot be empty or whitespace")
        return v.strip()
    
    model_config = ConfigDict(
        # Enable JSON schema generation with examples
        json_schema_extra={
            "example": {
                "infrastructure": "EMR",
                "confidence_score": 0.85,
                "decision_summary": "Large S3 files (500MB) exceed threshold; EMR recommended to avoid costly data transfers and process where data lives.",
                "reasoning": "Input files totaling 500MB are stored on S3. Downloading to local/EC2 would incur significant transfer costs (~$0.045) and time. EMR processes data in-place on S3, providing cost-effective distributed processing. Confidence is high (0.85) because file sizes are known and exceed the 100MB threshold significantly.",
                "file_analysis": {
                    "total_size_bytes": 524288000,
                    "total_size_mb": 500.0,
                    "file_count": 2,
                    "unknown_sizes": 0,
                    "locations": ["S3"],
                    "largest_file_mb": 250.0
                },
                "computational_requirements": {
                    "estimated_cpu_cores": 4,
                    "estimated_memory_gb": 8.0,
                    "estimated_runtime_minutes": 15.0,
                    "parallelizable": True,
                    "gpu_required": False
                },
                "cost_analysis": {
                    "estimated_cost_range_usd": (1.5, 3.5),
                    "cost_assumptions": "us-east-1, m5.xlarge nodes, 15min runtime, standard S3 storage",
                    "cost_confidence": 0.7,
                    "data_transfer_cost_usd": 0.0,
                    "breakdown": {
                        "compute": (1.5, 3.0),
                        "storage": (0.0, 0.5),
                        "transfer": 0.0
                    }
                },
                "alternatives": [
                    {
                        "infrastructure": "EC2",
                        "reasoning": "Could download files and process on EC2 with pre-installed tools",
                        "tradeoffs": "Faster setup but higher transfer cost ($0.045) and longer total time",
                        "confidence": 0.4
                    }
                ],
                "warnings": [],
                "inputs_analyzed": 2
            }
        }
    )


# -----------------------------------------------------------------------------
# Backward-compatible names (older code/tests)
# -----------------------------------------------------------------------------

class AlternativeRecommendation(InfraAlternative):
    """Back-compat wrapper with default confidence."""

    confidence: float = Field(
        default=0.4,
        ge=0.0,
        le=1.0,
        description="Confidence in this alternative (0-1)"
    )


InfrastructureDecision = InfraDecision
