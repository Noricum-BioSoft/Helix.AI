"""
Workflow Plan Contract - Input to Infrastructure Decision Agent.

This contract defines the structure of workflow plans that are passed to the
Infrastructure Decision Agent for analysis and infrastructure recommendation.

A WorkflowPlan describes:
- Data inputs (files, locations, sizes)
- Operations to perform (tool names, parameters)
- Constraints (time limits, cost limits, reproducibility needs)
- Expected scale (data volume, computation intensity)
"""

from typing import List, Optional, Dict, Any, Literal
from pydantic import BaseModel, Field, field_validator, ConfigDict


class DataInput(BaseModel):
    """A data input for the workflow (file, sequence, etc.)."""
    
    uri: str = Field(
        ...,
        description="URI of the input (s3://bucket/key, /local/path, etc.)"
    )
    size_bytes: Optional[int] = Field(
        default=None,
        ge=0,
        description="Size in bytes (None if unknown)"
    )
    location_type: Literal["S3", "Local", "URL", "Unknown"] = Field(
        default="Unknown",
        description="Type of storage location"
    )
    description: Optional[str] = Field(
        default=None,
        description="Human-readable description of this input"
    )
    metadata: Dict[str, Any] = Field(
        default_factory=dict,
        description="Additional metadata (format, checksum, etc.)"
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


class OperationSpec(BaseModel):
    """Specification of an operation to perform in the workflow."""
    
    operation_name: str = Field(
        ...,
        min_length=1,
        description="Name of the operation (e.g., 'read_merging', 'alignment')"
    )
    tool_name: Optional[str] = Field(
        default=None,
        description="Specific tool to use (e.g., 'bbmerge', 'bwa'), if known"
    )
    parameters: Dict[str, Any] = Field(
        default_factory=dict,
        description="Operation parameters"
    )
    expected_output_size_mb: Optional[float] = Field(
        default=None,
        ge=0.0,
        description="Expected output size in MB (if estimable)"
    )
    parallelizable: bool = Field(
        default=False,
        description="Whether this operation can be parallelized"
    )


class ConstraintSpec(BaseModel):
    """Constraints on workflow execution."""
    
    max_runtime_minutes: Optional[float] = Field(
        default=None,
        ge=0.1,
        description="Maximum acceptable runtime in minutes"
    )
    max_cost_usd: Optional[float] = Field(
        default=None,
        ge=0.0,
        description="Maximum acceptable cost in USD"
    )
    reproducibility_required: bool = Field(
        default=False,
        description="Whether reproducible results are critical (affects infrastructure choice)"
    )
    containerization_preferred: bool = Field(
        default=False,
        description="Whether containerized execution is preferred (Batch/Lambda over EC2)"
    )
    region_restriction: Optional[str] = Field(
        default=None,
        description="AWS region restriction (e.g., 'us-east-1')"
    )


class WorkflowPlan(BaseModel):
    """
    Workflow Plan - input to Infrastructure Decision Agent.
    
    Describes a bioinformatics workflow that needs infrastructure recommendation.
    The Infrastructure Decision Agent analyzes this plan and recommends the optimal
    execution environment (Local/EC2/EMR/Batch/Lambda).
    
    This contract is intentionally flexible to accommodate:
    - Simple single-operation workflows (e.g., "merge these 2 files")
    - Complex multi-step workflows (e.g., "align, call variants, annotate")
    - Workflows with unknown inputs (e.g., "user will upload files later")
    """
    
    # Core workflow description
    workflow_id: Optional[str] = Field(
        default=None,
        description="Unique workflow identifier (for tracking)"
    )
    description: str = Field(
        ...,
        min_length=10,
        description="Human-readable description of the workflow"
    )
    
    # Data inputs
    data_inputs: List[DataInput] = Field(
        default_factory=list,
        description="Input files/data for the workflow"
    )
    
    # Operations
    operations: List[OperationSpec] = Field(
        default_factory=list,
        description="Operations to perform (in order, if multi-step)"
    )
    
    # Constraints
    constraints: ConstraintSpec = Field(
        default_factory=ConstraintSpec,
        description="Execution constraints (time, cost, reproducibility)"
    )
    
    # Expected scale
    expected_data_volume_mb: Optional[float] = Field(
        default=None,
        ge=0.0,
        description="Expected total data volume in MB (sum of inputs + intermediates + outputs)"
    )
    expected_compute_intensity: Literal["Low", "Medium", "High", "Unknown"] = Field(
        default="Unknown",
        description="Expected computational intensity"
    )
    
    # Metadata
    session_id: Optional[str] = Field(
        default=None,
        description="User session ID (for tracking)"
    )
    user_id: Optional[str] = Field(
        default=None,
        description="User ID (for cost tracking)"
    )
    
    @field_validator('expected_data_volume_mb', mode='before')
    @classmethod
    def compute_expected_data_volume(cls, v, info):
        """Compute expected data volume from inputs if not provided."""
        if v is None and 'data_inputs' in info.data:
            total_bytes = 0
            for inp in info.data['data_inputs']:
                if inp.size_bytes is not None:
                    total_bytes += inp.size_bytes
            if total_bytes > 0:
                return total_bytes / (1024 * 1024)
        return v
    
    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "workflow_id": "workflow_123",
                "description": "Merge paired-end reads from RNA-seq experiment",
                "data_inputs": [
                    {
                        "uri": "s3://my-bucket/sample_R1.fq",
                        "size_bytes": 250000000,
                        "location_type": "S3",
                        "description": "Forward reads (R1)"
                    },
                    {
                        "uri": "s3://my-bucket/sample_R2.fq",
                        "size_bytes": 240000000,
                        "location_type": "S3",
                        "description": "Reverse reads (R2)"
                    }
                ],
                "operations": [
                    {
                        "operation_name": "read_merging",
                        "tool_name": "bbmerge",
                        "parameters": {
                            "min_overlap": 12,
                            "max_mismatch_rate": 0.05
                        },
                        "expected_output_size_mb": 200.0,
                        "parallelizable": False
                    }
                ],
                "constraints": {
                    "max_runtime_minutes": 30.0,
                    "max_cost_usd": 5.0,
                    "reproducibility_required": True,
                    "containerization_preferred": False
                },
                "expected_data_volume_mb": 466.15,
                "expected_compute_intensity": "Medium",
                "session_id": "session_456"
            }
        }
    )
