"""
ExecutionIntent — the immutable digest a human approves.

An approval bound only to ``plan_hash`` would survive a change of provider
("run locally" → "send to external provider Y") even though the risk changed
materially. The intent therefore covers plan + assessment + recommendation +
authorization + provider + capability + inputs + parameters, and approval
binds to ``execution_intent_hash``. Execution fails closed on any mismatch.

Assessment/recommendation/authorization fields are nullable only while the
corresponding feature flags are off (Phase 1 populates ``provider_id`` from
the legacy InfraDecision).
"""

from __future__ import annotations

from typing import Any, Dict, Optional

from pydantic import ConfigDict, Field, model_validator

from backend.contracts.base import TracedModel
from backend.contracts.ids import stable_hash

INTENT_HASH_FIELDS = (
    "plan_id",
    "plan_hash",
    "assessment_id",
    "assessment_hash",
    "recommendation_id",
    "authorization_id",
    "provider_id",
    "capability_id",
    "input_manifest_hash",
    "execution_parameters_hash",
)


def compute_intent_hash(values: Dict[str, Any]) -> str:
    return stable_hash({k: values.get(k) for k in INTENT_HASH_FIELDS})


class ExecutionIntent(TracedModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    execution_intent_id: str
    plan_id: str
    plan_hash: str
    assessment_id: Optional[str] = None
    assessment_hash: Optional[str] = None
    recommendation_id: Optional[str] = None
    authorization_id: Optional[str] = None
    provider_id: str = Field(..., description="Registry provider id, or legacy:<Local|EC2|EMR|Batch|Lambda> in Phase 1.")
    capability_id: str
    input_manifest_hash: str = Field(..., description="stable_hash over resolved input URIs + sizes + content hashes where known.")
    execution_parameters_hash: str = Field(..., description="stable_hash over step arguments after reference resolution.")
    execution_intent_hash: str = Field(default="", description="Computed over INTENT_HASH_FIELDS.")

    @model_validator(mode="before")
    @classmethod
    def _fill_hash(cls, data: Any) -> Any:
        if not isinstance(data, dict):
            return data
        expected = compute_intent_hash(data)
        given = data.get("execution_intent_hash") or ""
        if given and given != expected:
            raise ValueError("execution_intent_hash does not match intent content")
        return {**data, "execution_intent_hash": expected}

    def matches(self, other_hash: str) -> bool:
        return bool(other_hash) and other_hash == self.execution_intent_hash
