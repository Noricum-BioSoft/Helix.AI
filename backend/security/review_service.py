"""
Security review service (Phase 2): resolve a WAITING_FOR_SECURITY_REVIEW assessment.

A reviewer (role in ``policies/identity.yaml: review_roles``) approves or
denies. Approve → an ExecutionIntent is built (bound to the original
assessment hash + the review) and the session returns to
WAITING_FOR_APPROVAL: review is **not** approval, the normal human approval
still follows. Deny → DENIED (terminal for this plan).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Dict, Optional

from backend.contracts.human_approval import Principal
from backend.contracts.ids import new_id
from backend.contracts.security_review import ReviewDecision, SecurityReview
from backend.orchestration.ledger import LocalLedger, get_ledger
from backend.orchestration.plan_staging import build_intent_for
from backend.security.checks.identity import review_roles
from backend.workflow_checkpoint import WorkflowCheckpoint, WorkflowState

logger = logging.getLogger(__name__)


class ReviewError(Exception):
    status_code = 400

    def __init__(self, detail: str, *, code: str):
        super().__init__(detail)
        self.detail = detail
        self.code = code


class NothingToReview(ReviewError):
    status_code = 404


class ReviewConflict(ReviewError):
    status_code = 409


class ReviewForbidden(ReviewError):
    status_code = 403


@dataclass(frozen=True)
class ReviewOutcome:
    review: SecurityReview
    checkpoint: WorkflowCheckpoint
    intent_id: Optional[str]


def review_assessment(
    session_id: str,
    assessment_id: str,
    decision: ReviewDecision,
    principal: Optional[Principal],
    *,
    note: Optional[str] = None,
    checkpoint: Optional[WorkflowCheckpoint] = None,
    ledger: Optional[LocalLedger] = None,
    save: bool = True,
) -> ReviewOutcome:
    from backend.history_manager import history_manager

    if principal is None:
        raise ReviewForbidden("a principal is required to review", code="principal_required")
    allowed = set(review_roles())
    if allowed and not allowed & set(principal.roles):
        raise ReviewForbidden(f"reviewer role required ({', '.join(sorted(allowed))})", code="reviewer_role_required")

    ledger = ledger or get_ledger()
    checkpoint = checkpoint or history_manager.load_checkpoint(session_id)
    if checkpoint.state != WorkflowState.WAITING_FOR_SECURITY_REVIEW or not checkpoint.assessment_id:
        raise NothingToReview("no assessment is awaiting review for this session", code="no_pending_review")
    if checkpoint.assessment_id != assessment_id:
        raise ReviewConflict(f"assessment {assessment_id} is not the one awaiting review", code="assessment_not_pending")
    assessment = ledger.load_assessment(session_id, assessment_id)
    if assessment is None:
        raise NothingToReview("assessment not found in ledger", code="assessment_missing")
    if ledger.reviews_for_assessment(session_id, assessment_id):
        raise ReviewConflict("this assessment has already been reviewed", code="already_reviewed")
    if assessment.outcome != "REQUIRE_REVIEW":
        raise ReviewConflict(f"assessment outcome is {assessment.outcome}; only REQUIRE_REVIEW can be reviewed", code="not_reviewable")

    review = SecurityReview(
        review_id=new_id(),
        trace_id=assessment.trace_id,
        assessment_id=assessment.assessment_id,
        assessment_hash=assessment.assessment_hash,
        plan_id=assessment.plan_id,
        plan_hash=assessment.plan_hash,
        decision=decision,
        resulting_outcome="ALLOW_WITH_APPROVAL" if decision == "approved" else "DENY",
        principal=principal,
        note=note,
    )
    ledger.record_review(session_id, review)

    intent_id: Optional[str] = None
    if decision == "approved":
        plan = next(
            (p for p in ledger.list_plans(session_id) if p.plan_id == assessment.plan_id and p.plan_hash == assessment.plan_hash),
            None,
        )
        if plan is None:
            raise NothingToReview("assessed plan not found in ledger", code="plan_missing")
        intent = build_intent_for(plan, assessment, None, review=review)
        ledger.record_intent(session_id, intent)
        intent_id = intent.execution_intent_id
        new_cp = WorkflowCheckpoint.waiting_for_approval(pending_plan=checkpoint.pending_plan or {}).with_platform_records(
            **{**checkpoint.platform_records(),
               "pending_execution_intent_id": intent.execution_intent_id,
               "pending_execution_intent_hash": intent.execution_intent_hash,
               "security_outcome": review.resulting_outcome,
               "approval_id": None}
        )
        if save:
            _restore_legacy_pending_plan(session_id, checkpoint.pending_plan)
    else:
        new_cp = WorkflowCheckpoint.denied(pending_plan=checkpoint.pending_plan).with_platform_records(
            **{**checkpoint.platform_records(), "security_outcome": "DENY"}
        )
        if save:
            _clear_legacy_pending_plan(session_id)

    if save:
        history_manager.save_checkpoint(session_id, new_cp)
        ledger.record_audit_event(
            session_id,
            {
                "event_type": "security_review",
                "trace_id": review.trace_id,
                "assessment_id": assessment_id,
                "review_id": review.review_id,
                "decision": decision,
                "principal": principal.subject_id,
                "timestamp": datetime.now(timezone.utc).isoformat(),
            },
        )
    logger.info("[security] review trace=%s assessment=%s decision=%s by=%s", review.trace_id, assessment_id, decision, principal.subject_id)
    return ReviewOutcome(review=review, checkpoint=new_cp, intent_id=intent_id)


def _restore_legacy_pending_plan(session_id: str, pending: Optional[Dict[str, Any]]) -> None:
    """The legacy approval path reads ``session["pending_plan"]``; put it back after a review approval."""
    from backend.history_manager import history_manager

    if not isinstance(pending, dict) or not isinstance(pending.get("plan"), dict):
        return
    if hasattr(history_manager, "sessions") and session_id in history_manager.sessions:
        history_manager.sessions[session_id]["pending_plan"] = {
            "command": pending.get("command") or "",
            "plan": pending["plan"],
            "created_at": datetime.now(timezone.utc).isoformat(),
            "restored_by": "security_review",
        }


def _clear_legacy_pending_plan(session_id: str) -> None:
    from backend.history_manager import history_manager

    if hasattr(history_manager, "sessions") and session_id in history_manager.sessions:
        history_manager.sessions[session_id].pop("pending_plan", None)


def review_summary(outcome: ReviewOutcome) -> Dict[str, Any]:
    r = outcome.review
    return {
        "review_id": r.review_id,
        "trace_id": r.trace_id,
        "assessment_id": r.assessment_id,
        "assessment_hash": r.assessment_hash,
        "decision": r.decision,
        "resulting_outcome": r.resulting_outcome,
        "principal": r.principal.model_dump(mode="json"),
        "timestamp": r.timestamp.isoformat(),
        "execution_intent_id": outcome.intent_id,
        "workflow_state": outcome.checkpoint.state.value,
    }
