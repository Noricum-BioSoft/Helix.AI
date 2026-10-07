"""ScientificObjective — the structured intent a scientific loop starts from."""

from __future__ import annotations

from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field

from backend.contracts.base import TracedModel

ObjectiveStatus = Literal["proposed", "active", "answered", "abandoned"]


class ArtifactRef(BaseModel):
    model_config = ConfigDict(extra="forbid")

    artifact_id: Optional[str] = None
    uri: Optional[str] = None
    kind: Optional[str] = Field(default=None, description="dataset | file | sequence | report | …")
    content_hash: Optional[str] = None


class ExecutionConstraints(BaseModel):
    """Per-objective narrowing of where execution may happen (profile policy still wins)."""

    model_config = ConfigDict(extra="forbid")

    allowed_providers: List[str] = Field(default_factory=list)
    preferred_provider: Optional[str] = None
    data_residency: Optional[str] = None
    max_cost_usd: Optional[float] = Field(default=None, ge=0)


class ObjectiveConstraints(BaseModel):
    model_config = ConfigDict(extra="forbid")

    execution: Optional[ExecutionConstraints] = None
    time_budget_days: Optional[float] = Field(default=None, ge=0)
    other: List[str] = Field(default_factory=list, description="Free-form user constraints, kept as given.")


class ScientificObjective(TracedModel):
    objective_id: str
    version: int = Field(default=1, ge=1)
    objective: str = Field(..., min_length=1, description="One-sentence objective as the user framed it.")
    question: str = Field(..., min_length=1, description="The scientific question the loop should answer.")
    biological_system: Optional[str] = None
    constraints: ObjectiveConstraints = Field(default_factory=ObjectiveConstraints)
    desired_evidence: List[str] = Field(default_factory=list)
    known_inputs: List[ArtifactRef] = Field(default_factory=list)
    expected_outputs: List[str] = Field(default_factory=list)
    success_criteria: List[str] = Field(default_factory=list)
    user_context: Dict[str, Any] = Field(default_factory=dict)
    status: ObjectiveStatus = "proposed"
    parent_decision_id: Optional[str] = Field(
        default=None, description="Set when this objective was proposed by a NextDecision; stays 'proposed' until a user activates it."
    )
    supersedes_version: Optional[int] = None
