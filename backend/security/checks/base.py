"""
SecurityCheck protocol and shared types.

A check never sees a provider. It receives the plan, the objective and an
``AssessmentContext`` (session facts: uploads, principal, flags) and returns a
``CheckResult`` whose outcome and envelope the assessor combines.

Policies are declarative YAML under ``backend/config/security/policies/``;
``load_policy`` reads them once per process.
"""

from __future__ import annotations

import functools
from pathlib import Path
from typing import Any, Dict, List, Optional, Protocol, runtime_checkable

import yaml
from pydantic import BaseModel, ConfigDict, Field

from backend.contracts.human_approval import Principal
from backend.contracts.rationale import RationaleItem
from backend.contracts.scientific_objective import ScientificObjective
from backend.contracts.scientific_plan import ScientificPlan
from backend.contracts.security_assessment import (
    AppliedPolicy,
    PolicyEnvelope,
    RequiredApproval,
    ScreeningResult,
    SecurityOutcome,
)

POLICIES_DIR = Path(__file__).resolve().parents[2] / "config" / "security" / "policies"


class PolicyError(RuntimeError):
    pass


@functools.lru_cache(maxsize=None)
def load_policy(name: str) -> Dict[str, Any]:
    path = POLICIES_DIR / f"{name}.yaml"
    if not path.exists():
        raise PolicyError(f"security policy not found: {path}")
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict) or "policy_id" not in data or "version" not in data:
        raise PolicyError(f"security policy {name} must define policy_id and version")
    return data


class AssessmentContext(BaseModel):
    """Session facts a check may consult. Built by the caller; checks never read global state."""

    model_config = ConfigDict(extra="forbid", arbitrary_types_allowed=True)

    session_id: Optional[str] = None
    principal: Optional[Principal] = None
    uploaded_files: List[Dict[str, Any]] = Field(default_factory=list, description="Session upload metadata incl. intake_policy.")
    sequences: List[str] = Field(default_factory=list, description="Explicit sequences to screen, beyond those found in the plan.")
    requires_manual_review: bool = False
    profile: Optional[str] = None
    extra: Dict[str, Any] = Field(default_factory=dict)


class CheckResult(BaseModel):
    model_config = ConfigDict(extra="forbid")

    check_id: str
    version: str
    outcome: SecurityOutcome
    rationale: List[RationaleItem] = Field(default_factory=list)
    policies: List[AppliedPolicy] = Field(default_factory=list)
    screening_results: List[ScreeningResult] = Field(default_factory=list)
    required_approvals: List[RequiredApproval] = Field(default_factory=list)
    envelope: Optional[PolicyEnvelope] = Field(default=None, description="None = no constraint contributed.")
    risk_categories: List[str] = Field(default_factory=list, description="Informational; set by action classification.")


@runtime_checkable
class SecurityCheck(Protocol):
    check_id: str
    version: str

    def applies_to(self, plan: ScientificPlan) -> bool: ...

    def evaluate(self, plan: ScientificPlan, objective: ScientificObjective, context: AssessmentContext) -> CheckResult: ...


def applied(policy: Dict[str, Any], outcome: SecurityOutcome) -> AppliedPolicy:
    return AppliedPolicy(policy_id=str(policy["policy_id"]), version=str(policy["version"]), result=outcome)


def plan_text(plan: ScientificPlan, objective: ScientificObjective) -> str:
    """All human-authored text of a loop, lower-cased, for triage-style checks."""
    parts: List[str] = [objective.objective, objective.question, objective.biological_system or ""]
    parts += objective.desired_evidence + objective.expected_outputs + objective.success_criteria + objective.constraints.other
    for step in plan.steps:
        parts.append(step.description or "")
        parts += step.assumptions
        for value in step.arguments.values():
            if isinstance(value, str):
                parts.append(value)
    return "\n".join(p for p in parts if p).lower()
