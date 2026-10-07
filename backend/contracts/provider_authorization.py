"""
Provider-specific authorization (post-routing).

Answers "may this proceed through THIS provider under THESE conditions?" after
a recommendation exists. Re-checks the recommendation against the
PolicyEnvelope from the pre-routing assessment, which stays authoritative.
"""

from __future__ import annotations

from typing import List, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field

from backend.contracts.base import TracedModel
from backend.contracts.rationale import RationaleItem

AuthorizationResult = Literal["AUTHORIZED", "DENIED", "REQUIRES_SCREENING", "REQUIRES_REVIEW"]


class AuthorizationCheck(BaseModel):
    model_config = ConfigDict(extra="forbid")

    check_id: str = Field(
        ...,
        description="data_class | capability_supported | residency | screening_satisfied | human_approval_flag | envelope_compliance | provider_constraint:<name>",
    )
    result: Literal["pass", "fail", "not_applicable"]
    detail: Optional[str] = None


class ProviderAuthorization(TracedModel):
    authorization_id: str
    recommendation_id: str
    assessment_id: str
    provider_id: str
    capability_id: str
    result: AuthorizationResult
    checks: List[AuthorizationCheck] = Field(default_factory=list)
    rationale: List[RationaleItem] = Field(default_factory=list)

    @property
    def authorized(self) -> bool:
        return self.result == "AUTHORIZED"
