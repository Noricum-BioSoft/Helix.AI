"""
Secure Science pre-routing assessment.

Answers "can this class of action proceed?" *before* a provider is chosen and
emits a ``PolicyEnvelope`` that provider selection and provider-specific
authorization (P4) must obey. Keyword-based dual-use logic is triage only and
may contribute ALLOW or REQUIRE_REVIEW, never a safe/unsafe verdict.
"""

from __future__ import annotations

from typing import List, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, model_validator

from backend.contracts.base import TracedModel
from backend.contracts.ids import stable_hash
from backend.contracts.rationale import RationaleItem

SecurityOutcome = Literal["ALLOW", "ALLOW_WITH_APPROVAL", "REQUIRE_REVIEW", "DENY"]
OUTCOME_SEVERITY = {"ALLOW": 0, "ALLOW_WITH_APPROVAL": 1, "REQUIRE_REVIEW": 2, "DENY": 3}

DataClass = Literal["public", "internal", "sensitive", "phi", "restricted"]
ProviderCategory = Literal["computational", "experimental", "data", "advisory"]


def combine_outcomes(outcomes: List[SecurityOutcome]) -> SecurityOutcome:
    """DENY > REQUIRE_REVIEW > ALLOW_WITH_APPROVAL > ALLOW."""
    if not outcomes:
        return "ALLOW"
    return max(outcomes, key=lambda o: OUTCOME_SEVERITY[o])


class PolicyEnvelope(BaseModel):
    """Constraints later provider selection must satisfy. Intersection semantics when combined."""

    model_config = ConfigDict(extra="forbid")

    allowed_data_classes: List[DataClass] = Field(default_factory=lambda: ["public", "internal"])
    external_execution_allowed: bool = True
    screening_required: bool = False
    human_approval_required: bool = False
    allowed_regions: List[str] = Field(default_factory=list, description="Empty = no regional restriction.")
    allowed_provider_categories: List[ProviderCategory] = Field(
        default_factory=lambda: ["computational", "experimental", "data", "advisory"]
    )
    max_cost_usd: Optional[float] = Field(default=None, ge=0)

    def intersect(self, other: "PolicyEnvelope") -> "PolicyEnvelope":
        regions: List[str]
        if not self.allowed_regions:
            regions = list(other.allowed_regions)
        elif not other.allowed_regions:
            regions = list(self.allowed_regions)
        else:
            regions = [r for r in self.allowed_regions if r in other.allowed_regions]
        costs = [c for c in (self.max_cost_usd, other.max_cost_usd) if c is not None]
        return PolicyEnvelope(
            allowed_data_classes=[d for d in self.allowed_data_classes if d in other.allowed_data_classes],
            external_execution_allowed=self.external_execution_allowed and other.external_execution_allowed,
            screening_required=self.screening_required or other.screening_required,
            human_approval_required=self.human_approval_required or other.human_approval_required,
            allowed_regions=regions,
            allowed_provider_categories=[
                c for c in self.allowed_provider_categories if c in other.allowed_provider_categories
            ],
            max_cost_usd=min(costs) if costs else None,
        )


class AppliedPolicy(BaseModel):
    model_config = ConfigDict(extra="forbid")

    policy_id: str
    version: str
    result: SecurityOutcome


class ScreeningResult(BaseModel):
    model_config = ConfigDict(extra="forbid")

    provider: str = Field(..., description="Screening adapter id; 'mock' or 'not_configured' are explicit.")
    status: Literal["accepted", "flagged", "rejected", "not_run", "error"]
    summary: Optional[str] = None
    reference_id: Optional[str] = None


class RequiredApproval(BaseModel):
    model_config = ConfigDict(extra="forbid")

    role: str
    reason: str


class AssessmentAudit(BaseModel):
    model_config = ConfigDict(extra="forbid")

    actor: str = Field(..., description="agent:<id> | system | user:<subject>")
    gate_version: str
    profile: Optional[str] = None


def compute_assessment_hash(outcome: SecurityOutcome, envelope: PolicyEnvelope, policies: List[AppliedPolicy]) -> str:
    return stable_hash(
        {
            "outcome": outcome,
            "envelope": envelope.model_dump(mode="json"),
            "policies": sorted((p.model_dump(mode="json") for p in policies), key=lambda p: (p["policy_id"], p["version"])),
        }
    )


class SecurityAssessment(TracedModel):
    assessment_id: str
    plan_id: str
    plan_hash: str
    assessment_hash: str = Field(default="", description="Computed over outcome + envelope + applied policies.")
    outcome: SecurityOutcome
    rationale: List[RationaleItem] = Field(default_factory=list)
    applied_policies: List[AppliedPolicy] = Field(default_factory=list)
    screening_results: List[ScreeningResult] = Field(default_factory=list)
    required_approvals: List[RequiredApproval] = Field(default_factory=list)
    policy_envelope: PolicyEnvelope = Field(default_factory=PolicyEnvelope)
    audit: AssessmentAudit

    @model_validator(mode="after")
    def _fill_hash(self) -> "SecurityAssessment":
        expected = compute_assessment_hash(self.outcome, self.policy_envelope, self.applied_policies)
        if self.assessment_hash and self.assessment_hash != expected:
            raise ValueError("assessment_hash does not match assessment content")
        if not self.assessment_hash:
            object.__setattr__(self, "assessment_hash", expected)
        return self
