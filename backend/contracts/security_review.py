"""SecurityReview — a reviewer's decision on a REQUIRE_REVIEW assessment (Phase 2)."""

from __future__ import annotations

from datetime import datetime
from typing import Literal, Optional

from pydantic import ConfigDict, Field

from backend.contracts.base import TracedModel, utcnow
from backend.contracts.human_approval import Principal
from backend.contracts.security_assessment import SecurityOutcome

ReviewDecision = Literal["approved", "denied"]


class SecurityReview(TracedModel):
    """Bound to the assessment hash so a review cannot be replayed against a changed assessment."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    review_id: str
    assessment_id: str
    assessment_hash: str = Field(..., min_length=1)
    plan_id: str
    plan_hash: str = Field(..., min_length=1)
    decision: ReviewDecision
    resulting_outcome: SecurityOutcome = Field(..., description="ALLOW_WITH_APPROVAL when approved (human approval still required), DENY when denied.")
    principal: Principal
    timestamp: datetime = Field(default_factory=utcnow)
    note: Optional[str] = None

    @property
    def approved(self) -> bool:
        return self.decision == "approved"
