# =========================
# Agent Contracts (v1 draft)
# =========================
from enum import Enum
from typing import Any, Dict, List, Optional, Literal
from pydantic import BaseModel, Field


class AgentIntent(str, Enum):
    ASK = "ask"
    EXECUTE = "execute"


class IntentResult(BaseModel):
    """
    OWNER: Intent Detector
    Purpose: Decide whether the user wants an answer ("ask") or to run a workflow ("execute").
    """
    intent: AgentIntent
    confidence: float = Field(..., ge=0.0, le=1.0)
    clarifying_questions: List[str] = Field(default_factory=list)
    rationale: Optional[str] = None


class PlanScale(str, Enum):
    SMALL = "small"
    MEDIUM = "medium"
    LARGE = "large"


class PlanParallelism(str, Enum):
    NONE = "none"
    BY_SAMPLE = "by_sample"
    DISTRIBUTED = "distributed"


class ReproducibilityLevel(str, Enum):
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"


class WorkflowInput(BaseModel):
    """
    OWNER: Bioinformatics Executor (Planner)
    Represents a logical input to the workflow.
    """
    uri: str
    description: Optional[str] = None
    # Optional metadata. Prefer grounding via tools when possible.
    size_bytes: Optional[int] = None
    location_hint: Optional[Literal["Local", "S3", "Unknown"]] = None


class WorkflowOutput(BaseModel):
    """
    OWNER: Bioinformatics Executor (Planner)
    Represents an expected output artifact.
    """
    uri: Optional[str] = None
    description: Optional[str] = None
    artifact_type: Optional[str] = None  # e.g., "bam", "vcf", "qc_report", "plot"


class PlanMetadata(BaseModel):
    """
    OWNER: Bioinformatics Executor (Planner)
    High-signal annotations that downstream agents use (infra, codegen, viz).
    """
    scale: PlanScale = PlanScale.MEDIUM
    parallelism: PlanParallelism = PlanParallelism.NONE
    reproducibility: ReproducibilityLevel = ReproducibilityLevel.MEDIUM
    expected_outputs: List[WorkflowOutput] = Field(default_factory=list)
    risk_flags: List[str] = Field(default_factory=list)
    missing_info: List[str] = Field(default_factory=list)


class WorkflowPlan(BaseModel):
    """
    OWNER: Bioinformatics Executor (Planner)
    Wrapper around the Plan IR (plan_ir.Plan) with additional metadata and I/O declarations.

    Note: Keep this wrapper even if you later move Plan into shared/contracts.
    """
    version: str = Field(default="v1")
    metadata: PlanMetadata = Field(default_factory=PlanMetadata)
    inputs: List[WorkflowInput] = Field(default_factory=list)
    # Steps are your Plan IR steps. Keep as dicts if you don't want an import dependency here.
    # Later: replace with `from backend.plan_ir import PlanStep` or similar.
    steps: List[Dict[str, Any]] = Field(default_factory=list)


class InfrastructureEnvironment(str, Enum):
    LOCAL = "Local"
    EC2 = "EC2"
    AWS_BATCH = "AWS Batch"
    EMR = "EMR"
    LAMBDA = "Lambda"


class InfraConstraint(BaseModel):
    """
    OWNER: Infrastructure Expert
    Non-invasive requirements the infra selection implies (without rewriting the plan).
    """
    key: str  # e.g. "container_required", "data_locality", "split_by_sample"
    value: Any
    rationale: Optional[str] = None


class CostConfidence(str, Enum):
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"


class CostAnalysisContract(BaseModel):
    """
    OWNER: Infrastructure Expert
    Cost reasoning should prefer ranges and explicit assumptions.
    """
    estimated_cost_range_usd: str = Field(..., description="e.g., '$5-$20', '$50-$200', or 'unknown'")
    assumptions: List[str] = Field(default_factory=list)
    confidence: CostConfidence = CostConfidence.LOW


class InfraAlternativeContract(BaseModel):
    environment: InfrastructureEnvironment
    tradeoffs: str


class InfraDecisionContract(BaseModel):
    """
    OWNER: Infrastructure Expert
    Output of infrastructure selection. This is a thin, shared contract.
    (You may already have a richer InfraDecision model elsewhere; consider aliasing to it.)
    """
    recommended_environment: InfrastructureEnvironment
    confidence_score: float = Field(..., ge=0.0, le=1.0)
    decision_summary: str
    constraints: List[InfraConstraint] = Field(default_factory=list)
    cost_analysis: Optional[CostAnalysisContract] = None
    alternatives: List[InfraAlternativeContract] = Field(default_factory=list)
    warnings: List[str] = Field(default_factory=list)
    assumptions: List[str] = Field(default_factory=list)


class ExecutionTarget(str, Enum):
    """
    OWNER: Code Generator (Tooling/Packaging Agent)
    """
    LOCAL = "local"
    EC2 = "ec2"
    EMR = "emr"
    AWS_BATCH = "aws_batch"
    LAMBDA = "lambda"


class ExecutionSpec(BaseModel):
    """
    OWNER: Code Generator (Tooling/Packaging Agent)
    Runnable description of how to execute the WorkflowPlan on the selected infrastructure.
    Does NOT execute anything by itself.
    """
    target: ExecutionTarget
    # container/image is optional but recommended when reproducibility >= medium
    container_image: Optional[str] = None
    entrypoint: Optional[str] = None
    command: List[str] = Field(default_factory=list)
    env: Dict[str, str] = Field(default_factory=dict)
    inputs: List[WorkflowInput] = Field(default_factory=list)
    outputs: List[WorkflowOutput] = Field(default_factory=list)
    resources: Dict[str, Any] = Field(default_factory=dict)  # cpu/mem/runtime hints
    retry_strategy: Dict[str, Any] = Field(default_factory=dict)
    notes: List[str] = Field(default_factory=list)


class ExecutionStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    CANCELLED = "cancelled"


class ExecutionResult(BaseModel):
    """
    OWNER: Execution Broker (service)
    Captures side-effectful execution outcomes for downstream visualization/reporting.
    """
    status: ExecutionStatus
    job_id: Optional[str] = None
    environment: Optional[InfrastructureEnvironment] = None
    logs_uri: Optional[str] = None
    artifacts: List[WorkflowOutput] = Field(default_factory=list)
    metrics: Dict[str, Any] = Field(default_factory=dict)
    error: Optional[str] = None


class VisualizationArtifact(BaseModel):
    """
    OWNER: Data Visualizer
    """
    artifact_type: str  # e.g., "plot", "table", "report"
    uri: Optional[str] = None
    description: Optional[str] = None
    metadata: Dict[str, Any] = Field(default_factory=dict)


class VisualizationArtifacts(BaseModel):
    """
    OWNER: Data Visualizer
    """
    summary: str
    artifacts: List[VisualizationArtifact] = Field(default_factory=list)
    warnings: List[str] = Field(default_factory=list)
