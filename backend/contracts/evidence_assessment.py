"""
Learning layer: Observation ≠ Interpretation ≠ Decision.

Observations are measured facts and carry no claims. Interpretation lists what
the evidence supports/contradicts and the open uncertainties. NextDecision
proposes the next step and may carry a proposed objective, which stays
``proposed`` until a user explicitly activates it — never auto-executed.
"""

from __future__ import annotations

from typing import List, Literal, Optional, Union

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from backend.contracts.base import TracedModel
from backend.contracts.rationale import RationaleItem
from backend.contracts.scientific_objective import ScientificObjective

DecisionKind = Literal[
    "stop", "rerun_analysis", "acquire_data", "repeat_experiment", "new_condition", "switch_provider", "escalate"
]


class Observation(TracedModel):
    observation_id: str
    run_ids: List[str] = Field(default_factory=list)
    metric: str = Field(..., min_length=1)
    value: Union[float, int, str, bool]
    unit: Optional[str] = None
    artifact_refs: List[str] = Field(default_factory=list)
    method: Optional[str] = Field(default=None, description="How the value was obtained (tool/version), not what it means.")


class Interpretation(BaseModel):
    model_config = ConfigDict(extra="forbid")

    supports: List[str] = Field(default_factory=list)
    contradicts: List[str] = Field(default_factory=list)
    uncertainties: List[str] = Field(default_factory=list)
    rationale: List[RationaleItem] = Field(default_factory=list)


class NextDecision(TracedModel):
    decision_id: str
    kind: DecisionKind
    rationale: List[RationaleItem] = Field(..., min_length=1)
    evidence_ids: List[str] = Field(..., min_length=1, description="Observation/artifact ids this decision rests on.")
    proposed_objective: Optional[ScientificObjective] = None

    @model_validator(mode="after")
    def _proposed_stays_proposed(self) -> "NextDecision":
        if self.proposed_objective is not None:
            if self.proposed_objective.status != "proposed":
                raise ValueError("proposed_objective.status must be 'proposed'")
            if self.proposed_objective.parent_decision_id != self.decision_id:
                raise ValueError("proposed_objective.parent_decision_id must reference this decision")
        return self


class EvidenceAssessment(TracedModel):
    objective_id: str
    run_ids: List[str] = Field(default_factory=list)
    observations: List[Observation] = Field(default_factory=list)
    interpretation: Interpretation
    decision: NextDecision

    @field_validator("observations")
    @classmethod
    def _observations_share_runs(cls, v: List[Observation]) -> List[Observation]:
        return v
