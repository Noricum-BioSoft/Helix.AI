"""HumanApproval — persisted, auditable decision bound to an immutable ExecutionIntent."""

from __future__ import annotations

from datetime import datetime
from typing import List, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field

from backend.contracts.base import TracedModel, utcnow

ApprovalDecision = Literal["approved", "rejected", "changes_requested"]
IdentityProvider = Literal["dev_header", "oidc", "natural_language", "system"]


class Principal(BaseModel):
    """Who approved. ``identity_provider`` records how much the identity can be trusted."""

    model_config = ConfigDict(extra="forbid")

    subject_id: str = Field(..., min_length=1, description="Stable subject identifier from the identity provider.")
    display_name: Optional[str] = None
    identity_provider: IdentityProvider
    auth_method: str = Field(..., description="header | oidc_bearer | chat_message | …")
    roles: List[str] = Field(default_factory=list)
    tenant_id: Optional[str] = None
    project_id: Optional[str] = None
    audit_signature: Optional[str] = Field(default=None, description="Reserved for signed audit metadata in production.")


class HumanApproval(TracedModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    approval_id: str
    execution_intent_id: str
    execution_intent_hash: str = Field(..., min_length=1)
    plan_id: str
    plan_hash: str = Field(..., min_length=1)
    decision: ApprovalDecision
    principal: Principal
    timestamp: datetime = Field(default_factory=utcnow)
    note: Optional[str] = None
    scope: Literal["plan", "steps"] = "plan"
    step_ids: List[str] = Field(default_factory=list, description="Populated when scope == 'steps'.")
    selected_recommendation_id: Optional[str] = Field(
        default=None, description="Approver chose an alternative; must be in the recommendation's candidate_set."
    )

    @property
    def approved(self) -> bool:
        return self.decision == "approved"
