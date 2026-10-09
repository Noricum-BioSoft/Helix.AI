"""ExecutionRecommendation — explainable provider/capability choice over the registry."""

from __future__ import annotations

from typing import Dict, List, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, model_validator

from backend.contracts.base import TracedModel
from backend.contracts.rationale import RationaleItem

LegacyInfrastructure = Literal["Local", "EC2", "EMR", "Batch", "Lambda"]

CRITERIA = (
    "suitability",
    "capability_match",
    "inputs_ready",
    "cost",
    "turnaround",
    "throughput",
    "security",
    "privacy",
    "availability",
    "prior_performance",
    "provenance",
)


class CriteriaScores(BaseModel):
    """Each criterion scored 0–1 (higher is better). Missing criteria are not allowed."""

    model_config = ConfigDict(extra="forbid")

    suitability: float = Field(..., ge=0, le=1)
    capability_match: float = Field(..., ge=0, le=1)
    inputs_ready: float = Field(..., ge=0, le=1)
    cost: float = Field(..., ge=0, le=1)
    turnaround: float = Field(..., ge=0, le=1)
    throughput: float = Field(..., ge=0, le=1)
    security: float = Field(..., ge=0, le=1)
    privacy: float = Field(..., ge=0, le=1)
    availability: float = Field(..., ge=0, le=1)
    prior_performance: float = Field(..., ge=0, le=1)
    provenance: float = Field(..., ge=0, le=1)

    def weighted(self, weights: Dict[str, float]) -> float:
        total_w = sum(weights.get(c, 0.0) for c in CRITERIA) or 1.0
        return sum(getattr(self, c) * weights.get(c, 0.0) for c in CRITERIA) / total_w


class RecommendationAlternative(BaseModel):
    model_config = ConfigDict(extra="forbid")

    provider_id: str
    capability_id: str
    score: float = Field(..., ge=0, le=1)
    criteria_scores: CriteriaScores
    tradeoffs: List[str] = Field(default_factory=list)


class ExecutionRecommendation(TracedModel):
    recommendation_id: str
    plan_id: str
    step_id: Optional[str] = None
    provider_id: str
    capability_id: str
    candidate_set: List[str] = Field(..., description="All provider_ids that were eligible; approval may only pick from these.")
    criteria_scores: CriteriaScores
    weights: Dict[str, float] = Field(default_factory=dict)
    score: float = Field(..., ge=0, le=1)
    rationale: List[RationaleItem] = Field(default_factory=list, description="Items should reference criteria via evidence_refs 'criterion:<name>'.")
    alternatives: List[RecommendationAlternative] = Field(default_factory=list)
    warnings: List[str] = Field(default_factory=list)
    confidence: float = Field(..., ge=0, le=1)
    legacy_infrastructure: Optional[LegacyInfrastructure] = Field(
        default=None, description="Back-compat mapping for InfraDecision consumers; None for non-compute providers."
    )

    @model_validator(mode="after")
    def _provider_in_candidates(self) -> "ExecutionRecommendation":
        if self.provider_id not in self.candidate_set:
            raise ValueError("recommended provider_id must be in candidate_set")
        for alt in self.alternatives:
            if alt.provider_id not in self.candidate_set:
                raise ValueError(f"alternative provider {alt.provider_id} not in candidate_set")
        return self
