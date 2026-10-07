"""
ExecutionProvider — the generic contract every execution backend implements.

Lifecycle: describe_capabilities → validate_request → estimate →
prepare_execution → execute → get_status* → retrieve_results | cancel →
provenance. Adapters receive a ``ProviderConfig`` (execution profile) and
never read environment variables themselves.

Idempotency: ``ExecutionRequest.idempotency_key`` is persisted by the fabric
before ``execute``. An adapter declares via ``idempotency_support`` whether
it honours the key natively, by client-side de-duplication (look up the key
before submitting again), or not at all; its docstring must say how.
"""

from __future__ import annotations

from datetime import datetime
from typing import Any, Dict, List, Optional, Protocol, runtime_checkable

from pydantic import BaseModel, ConfigDict, Field

from backend.contracts.execution_request import ExecutionRequest, RunStatus
from backend.contracts.provenance import ProvenanceRecord
from backend.contracts.scientific_objective import ArtifactRef
from shared.capability_registry import CapabilityDescriptor, IdempotencySupport


class ProviderNotConfigured(RuntimeError):
    """Raised by stub providers and by adapters whose ProviderConfig is incomplete."""


class ValidationResult(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: bool
    errors: List[str] = Field(default_factory=list)
    warnings: List[str] = Field(default_factory=list)


class Estimate(BaseModel):
    model_config = ConfigDict(extra="forbid")

    cost_range_usd: List[float] = Field(..., min_length=2, max_length=2)
    turnaround_minutes: Optional[float] = Field(default=None, ge=0)
    confidence: float = Field(..., ge=0, le=1)
    assumptions: List[str] = Field(default_factory=list)


class PreparedExecution(BaseModel):
    model_config = ConfigDict(extra="forbid")

    execution_request_id: str
    idempotency_key: str
    provider_id: str
    capability_id: str
    payload: Dict[str, Any] = Field(default_factory=dict, description="Provider-specific submission payload (no secrets).")
    staged_inputs: List[str] = Field(default_factory=list)


class ExecutionHandle(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    provider_id: str
    provider_handle: str
    idempotency_key: str
    execution_request_id: str
    submitted_at: datetime


class ExecutionStatus(BaseModel):
    model_config = ConfigDict(extra="forbid")

    status: RunStatus
    progress: Optional[float] = Field(default=None, ge=0, le=1)
    message: Optional[str] = None
    cost_so_far_usd: Optional[float] = Field(default=None, ge=0)


class ExecutionResult(BaseModel):
    model_config = ConfigDict(extra="forbid")

    status: RunStatus
    outputs: List[ArtifactRef] = Field(default_factory=list)
    metrics: Dict[str, Any] = Field(default_factory=dict)
    cost_actual_usd: Optional[float] = Field(default=None, ge=0)
    error: Optional[str] = None
    failure_report: Optional[str] = Field(default=None, description="Provider debug report, attached as a failure Observation.")


@runtime_checkable
class ExecutionProvider(Protocol):
    provider_id: str

    @property
    def idempotency_support(self) -> IdempotencySupport: ...

    def describe_capabilities(self) -> List[CapabilityDescriptor]: ...

    def validate_request(self, request: ExecutionRequest) -> ValidationResult: ...

    def estimate(self, request: ExecutionRequest) -> Estimate: ...

    def prepare_execution(self, request: ExecutionRequest) -> PreparedExecution: ...

    def execute(self, prepared: PreparedExecution) -> ExecutionHandle: ...

    def get_status(self, handle: ExecutionHandle) -> ExecutionStatus: ...

    def retrieve_results(self, handle: ExecutionHandle) -> ExecutionResult: ...

    def cancel(self, handle: ExecutionHandle) -> None: ...

    def provenance(self, handle: ExecutionHandle) -> ProvenanceRecord: ...


PROVIDER_METHODS = (
    "describe_capabilities",
    "validate_request",
    "estimate",
    "prepare_execution",
    "execute",
    "get_status",
    "retrieve_results",
    "cancel",
    "provenance",
)
