"""
ScientificPlan — versioned, hashed wrapper around the existing plan IR.

Wraps, does not replace, ``backend.plan_ir.Plan``. ``plan_hash`` covers the
executable content (steps, IR, required data/capabilities) and deliberately
excludes rationale so that re-phrasing an explanation never invalidates an
approval.
"""

from __future__ import annotations

from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, model_validator

from backend.contracts.base import TracedModel
from backend.contracts.ids import stable_hash
from backend.contracts.rationale import RationaleItem
from backend.plan_ir import Plan as PlanIR
from backend.plan_ir import PlanStep as IRStep

StepKind = Literal["reasoning", "computational", "experimental"]
PlanStatus = Literal[
    "draft", "assessed", "recommended", "authorized", "approved", "rejected", "executing", "completed", "failed"
]


class ScientificPlanStep(BaseModel):
    """Plan-IR step plus the scientific metadata the gate/recommender need."""

    model_config = ConfigDict(extra="forbid")

    id: str
    action_type: str = "execute_tool"
    tool_name: Optional[str] = None
    arguments: Dict[str, Any] = Field(default_factory=dict)
    description: Optional[str] = None
    kind: StepKind = "computational"
    required_data: List[str] = Field(default_factory=list)
    required_capabilities: List[str] = Field(default_factory=list, description="capability_ids, e.g. local_compute:fastqc")
    assumptions: List[str] = Field(default_factory=list)
    depends_on: List[str] = Field(default_factory=list)

    @classmethod
    def from_ir(cls, step: IRStep, **extra: Any) -> "ScientificPlanStep":
        return cls(
            id=step.id,
            action_type=step.action_type,
            tool_name=step.tool_name,
            arguments=dict(step.arguments),
            description=step.description,
            **extra,
        )


def compute_plan_hash(steps: List[ScientificPlanStep], ir: PlanIR, security_requirements: List[str]) -> str:
    return stable_hash(
        {
            "steps": [s.model_dump(mode="json") for s in steps],
            "ir": ir.model_dump(mode="json"),
            "security_requirements": sorted(security_requirements),
        }
    )


class ScientificPlan(TracedModel):
    plan_id: str
    version: int = Field(default=1, ge=1)
    objective_id: str
    plan_hash: str = Field(default="", description="Computed; covers steps+ir+security_requirements, not rationale.")
    plan_rationale: List[RationaleItem] = Field(default_factory=list)
    steps: List[ScientificPlanStep]
    security_requirements: List[str] = Field(default_factory=list)
    expected_outputs: List[str] = Field(default_factory=list)
    status: PlanStatus = "draft"
    supersedes_plan_id: Optional[str] = None
    ir: PlanIR

    @model_validator(mode="after")
    def _fill_hash(self) -> "ScientificPlan":
        expected = compute_plan_hash(self.steps, self.ir, self.security_requirements)
        if self.plan_hash and self.plan_hash != expected:
            raise ValueError("plan_hash does not match plan content")
        if not self.plan_hash:
            object.__setattr__(self, "plan_hash", expected)
        return self

    @classmethod
    def from_plan_ir(
        cls,
        ir: PlanIR,
        *,
        plan_id: str,
        trace_id: str,
        objective_id: str,
        rationale: Optional[List[RationaleItem]] = None,
        step_kinds: Optional[Dict[str, StepKind]] = None,
        step_capabilities: Optional[Dict[str, List[str]]] = None,
        **kwargs: Any,
    ) -> "ScientificPlan":
        step_kinds = step_kinds or {}
        step_capabilities = step_capabilities or {}
        steps = [
            ScientificPlanStep.from_ir(
                s,
                kind=step_kinds.get(s.id, "computational"),
                required_capabilities=step_capabilities.get(s.id, [f"local_compute:{s.tool_name}"] if s.tool_name else []),
            )
            for s in ir.steps
        ]
        return cls(
            plan_id=plan_id,
            trace_id=trace_id,
            objective_id=objective_id,
            plan_rationale=rationale or [],
            steps=steps,
            ir=ir,
            **kwargs,
        )
