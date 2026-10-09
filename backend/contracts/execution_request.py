"""
ExecutionRequest / ExecutionRun — idempotent submission to a provider.

``idempotency_key`` is derived from the execution intent hash plus a nonce,
generated once and persisted *before* any provider call, so a retry after a
network timeout re-attaches to the existing submission instead of creating
duplicate work (money for compute, physical duplicates for wet-lab).
"""

from __future__ import annotations

import hashlib
import secrets
from datetime import datetime
from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field

from backend.contracts.base import TracedModel
from backend.contracts.scientific_objective import ArtifactRef

RunStatus = Literal["pending", "submitted", "running", "succeeded", "failed", "cancelled", "unknown"]


def make_idempotency_key(execution_intent_hash: str, nonce: Optional[str] = None) -> str:
    nonce = nonce or secrets.token_hex(16)
    return hashlib.sha256(f"{execution_intent_hash}:{nonce}".encode("utf-8")).hexdigest()


class ExecutionInput(BaseModel):
    model_config = ConfigDict(extra="forbid")

    uri: str
    size_bytes: Optional[int] = Field(default=None, ge=0)
    content_hash: Optional[str] = None
    role: Optional[str] = Field(default=None, description="e.g. reads_r1, samplesheet, reference")


class ExecutionRequest(TracedModel):
    execution_request_id: str
    execution_intent_id: str
    execution_intent_hash: str = Field(..., min_length=1)
    approval_id: Optional[str] = Field(default=None, description="Required for anything that is not read-only.")
    idempotency_key: str = Field(..., min_length=16)
    provider_id: str
    capability_id: str
    inputs: List[ExecutionInput] = Field(default_factory=list)
    parameters: Dict[str, Any] = Field(default_factory=dict)
    submitted_at: Optional[datetime] = None
    provider_handle: Optional[str] = Field(default=None, description="Provider-side identifier once submitted.")


class ExecutionRun(TracedModel):
    execution_run_id: str
    execution_request_id: str
    status: RunStatus = "pending"
    provider_id: str
    provider_handle: Optional[str] = None
    started_at: Optional[datetime] = None
    finished_at: Optional[datetime] = None
    cost_actual_usd: Optional[float] = Field(default=None, ge=0)
    outputs: List[ArtifactRef] = Field(default_factory=list)
    error: Optional[str] = None
