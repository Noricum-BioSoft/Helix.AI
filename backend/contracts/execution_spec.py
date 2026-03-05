"""
Execution Specification Contract - Output from ImplementationAgent.

The ImplementationAgent consumes WorkflowPlan + InfraDecision and produces
an ExecutionToolSpec that describes HOW to execute (but does not execute).

This contract defines:
- Container/environment setup
- Execution commands
- Retry policies
- Resource requirements
- Expected outputs

Key principle: Separation of planning vs execution.
ImplementationAgent plans, external runner executes.
"""

from typing import List, Optional, Dict, Any, Literal
from pydantic import BaseModel, Field, field_validator, model_validator, ConfigDict


class ContainerSpec(BaseModel):
    """Container/environment specification for execution."""
    
    image: str = Field(
        ...,
        description="Container image (e.g., docker.io/biocontainers/fastqc:0.11.9)"
    )
    
    image_type: Literal["docker", "singularity", "conda", "native"] = Field(
        ...,
        description="Container/environment type"
    )
    
    pull_policy: Literal["Always", "IfNotPresent", "Never"] = Field(
        default="IfNotPresent",
        description="Image pull policy"
    )
    
    env_vars: Dict[str, str] = Field(
        default_factory=dict,
        description="Environment variables to set"
    )
    
    mount_paths: List[str] = Field(
        default_factory=list,
        description="Paths to mount into container"
    )
    
    working_dir: Optional[str] = Field(
        default=None,
        description="Working directory inside container"
    )


class CommandSpec(BaseModel):
    """Single command/step specification."""
    
    name: str = Field(
        ...,
        min_length=1,
        description="Human-readable command name"
    )
    
    command: str = Field(
        ...,
        min_length=1,
        description="Shell command to execute"
    )
    
    inputs: List[str] = Field(
        default_factory=list,
        description="Input file URIs required by this command"
    )
    
    outputs: List[str] = Field(
        default_factory=list,
        description="Output file URIs produced by this command"
    )
    
    success_criteria: Optional[str] = Field(
        default=None,
        description="How to determine success (e.g., 'exit_code==0', 'output_file_exists')"
    )
    
    timeout_minutes: Optional[float] = Field(
        default=None,
        ge=0.0,
        description="Command timeout in minutes (None = no timeout)"
    )


class RetryPolicy(BaseModel):
    """Retry policy for execution failures."""
    
    max_retries: int = Field(
        default=0,
        ge=0,
        le=5,
        description="Maximum retry attempts (0 = no retries)"
    )
    
    retry_on: List[Literal["exit_code", "timeout", "oom", "network", "all"]] = Field(
        default_factory=lambda: ["all"],
        description="Conditions that trigger retry"
    )
    
    backoff_multiplier: float = Field(
        default=2.0,
        ge=1.0,
        description="Backoff multiplier for retries (2.0 = double wait time)"
    )
    
    initial_delay_seconds: float = Field(
        default=10.0,
        ge=0.0,
        description="Initial delay before first retry"
    )


class ResourceRequirements(BaseModel):
    """Resource requirements for execution."""
    
    min_cpu_cores: Optional[float] = Field(
        default=None,
        ge=0.0,
        description="Minimum CPU cores"
    )
    
    min_memory_gb: Optional[float] = Field(
        default=None,
        ge=0.0,
        description="Minimum memory in GB"
    )
    
    min_disk_gb: Optional[float] = Field(
        default=None,
        ge=0.0,
        description="Minimum disk space in GB"
    )
    
    gpu_required: bool = Field(
        default=False,
        description="Whether GPU is required"
    )
    
    gpu_count: Optional[int] = Field(
        default=None,
        ge=0,
        description="Number of GPUs required (if gpu_required=True)"
    )


class OutputSpec(BaseModel):
    """Specification for expected outputs."""
    
    uri: str = Field(
        ...,
        description="Output file URI (s3://bucket/key or local path)"
    )
    
    format: Optional[str] = Field(
        default=None,
        description="Expected file format (e.g., 'fastq', 'bam', 'vcf')"
    )
    
    required: bool = Field(
        default=True,
        description="Whether this output is required for success"
    )
    
    size_estimate_mb: Optional[float] = Field(
        default=None,
        ge=0.0,
        description="Expected output size in MB (for validation)"
    )


class ExecutionToolSpec(BaseModel):
    """
    Complete specification for tool execution.
    
    This is the output from ImplementationAgent.
    It describes HOW to execute a workflow step, but does NOT execute it.
    
    External runner consumes this spec to perform actual execution.
    """
    
    tool_name: str = Field(
        ...,
        min_length=1,
        description="Tool/workflow name"
    )
    
    infrastructure: Literal["Local", "EC2", "EMR", "Batch", "Lambda"] = Field(
        ...,
        description="Target execution environment (from InfraDecision)"
    )
    
    container_spec: Optional[ContainerSpec] = Field(
        default=None,
        description="Container specification (None for native execution)"
    )
    
    commands: List[CommandSpec] = Field(
        ...,
        min_length=1,
        description="Commands to execute (in order)"
    )
    
    retry_policy: RetryPolicy = Field(
        default_factory=RetryPolicy,
        description="Retry policy for failures"
    )
    
    resource_requirements: ResourceRequirements = Field(
        default_factory=ResourceRequirements,
        description="Resource requirements"
    )
    
    expected_outputs: List[OutputSpec] = Field(
        default_factory=list,
        description="Expected output files"
    )
    
    # Metadata
    confidence_score: float = Field(
        ...,
        ge=0.0,
        le=1.0,
        description="Confidence in this execution plan (0-1)"
    )
    
    reasoning: str = Field(
        ...,
        min_length=50,
        description="Why this execution plan was chosen"
    )
    
    warnings: List[str] = Field(
        default_factory=list,
        description="Warnings or caveats about this execution plan"
    )
    
    estimated_runtime_minutes: Optional[float] = Field(
        default=None,
        ge=0.0,
        description="Estimated runtime in minutes"
    )
    
    # Traceability
    request_id: Optional[str] = Field(
        default=None,
        description="Request ID for tracing"
    )
    
    workflow_plan_hash: Optional[str] = Field(
        default=None,
        description="Hash of input WorkflowPlan (for reproducibility)"
    )
    
    infra_decision_hash: Optional[str] = Field(
        default=None,
        description="Hash of input InfraDecision (for reproducibility)"
    )
    
    @field_validator('commands')
    @classmethod
    def validate_commands_not_empty(cls, v):
        """Ensure at least one command is specified."""
        if not v:
            raise ValueError("At least one command must be specified")
        return v
    
    @model_validator(mode='after')
    def check_gpu_consistency(self):
        """If GPU required, ensure gpu_count is set."""
        if self.resource_requirements.gpu_required:
            if self.resource_requirements.gpu_count is None or self.resource_requirements.gpu_count == 0:
                self.warnings.append("GPU required but gpu_count not specified - defaulting to 1")
                self.resource_requirements.gpu_count = 1
        return self
    
    @model_validator(mode='after')
    def check_container_for_batch_lambda(self):
        """Batch and Lambda require container specs."""
        if self.infrastructure in ["Batch", "Lambda"]:
            if self.container_spec is None:
                self.warnings.append(
                    f"{self.infrastructure} requires container spec - execution may fail"
                )
        return self
    
    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "tool_name": "fastqc",
                "infrastructure": "EMR",
                "container_spec": {
                    "image": "docker.io/biocontainers/fastqc:0.11.9",
                    "image_type": "docker",
                    "pull_policy": "IfNotPresent",
                    "env_vars": {},
                    "mount_paths": ["/data"],
                    "working_dir": "/data"
                },
                "commands": [
                    {
                        "name": "run_fastqc",
                        "command": "fastqc -o /data/output /data/input_R1.fq /data/input_R2.fq",
                        "inputs": ["s3://bucket/input_R1.fq", "s3://bucket/input_R2.fq"],
                        "outputs": ["s3://bucket/output/input_R1_fastqc.html", "s3://bucket/output/input_R2_fastqc.html"],
                        "success_criteria": "exit_code==0",
                        "timeout_minutes": 30.0
                    }
                ],
                "retry_policy": {
                    "max_retries": 2,
                    "retry_on": ["exit_code", "timeout"],
                    "backoff_multiplier": 2.0,
                    "initial_delay_seconds": 10.0
                },
                "resource_requirements": {
                    "min_cpu_cores": 4.0,
                    "min_memory_gb": 8.0,
                    "min_disk_gb": 50.0,
                    "gpu_required": False,
                    "gpu_count": None
                },
                "expected_outputs": [
                    {
                        "uri": "s3://bucket/output/input_R1_fastqc.html",
                        "format": "html",
                        "required": True,
                        "size_estimate_mb": 5.0
                    }
                ],
                "confidence_score": 0.85,
                "reasoning": "FastQC is containerized and well-tested on EMR with Spark. Files are in S3, EMR can process in-place.",
                "warnings": [],
                "estimated_runtime_minutes": 15.0,
                "request_id": "req_12345",
                "workflow_plan_hash": "abc123...",
                "infra_decision_hash": "def456..."
            }
        }
    )
