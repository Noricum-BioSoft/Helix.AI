"""
Security endpoints (Phase 2).

    GET  /session/{session_id}/assessments/pending
    GET  /session/{session_id}/assessments/{assessment_id}
    POST /session/{session_id}/assessments/{assessment_id}/review   {decision: approve|deny, note?}

Reviewing requires a principal with a review role (``policies/identity.yaml``).
Approving a review does not execute and does not approve: it builds the
ExecutionIntent and returns the session to WAITING_FOR_APPROVAL.
"""

from __future__ import annotations

from typing import Literal, Optional

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel, Field

from backend.config.auth_mode import AuthConfigurationError, principal_from_headers
from backend.config.feature_flags import science_gate_enabled
from backend.orchestration.ledger import get_ledger
from backend.security.review_service import ReviewError, review_assessment, review_summary

router = APIRouter(tags=["security"])


class ReviewBody(BaseModel):
    decision: Literal["approve", "deny"]
    note: Optional[str] = Field(default=None, max_length=2000)


def _require_gate() -> None:
    if not science_gate_enabled():
        raise HTTPException(status_code=404, detail="security endpoints require HELIX_SCIENCE_GATE_V1=1")


def _session_or_404(session_id: str):
    from backend.history_manager import history_manager

    if not history_manager.get_session(session_id):
        raise HTTPException(status_code=404, detail=f"Session {session_id} not found")
    return history_manager


def _principal(request: Request):
    try:
        principal = principal_from_headers(request.headers)
    except AuthConfigurationError as exc:
        raise HTTPException(status_code=500, detail=str(exc))
    if principal is None:
        raise HTTPException(status_code=401, detail="missing principal (X-Helix-User header in dev_header mode)")
    return principal


@router.get("/session/{session_id}/assessments/pending")
async def pending_assessment(session_id: str):
    _require_gate()
    hm = _session_or_404(session_id)
    cp = hm.load_checkpoint(session_id)
    if not cp.assessment_id:
        return {"success": True, "session_id": session_id, "assessment": None, "workflow_state": cp.state.value}
    assessment = get_ledger().load_assessment(session_id, cp.assessment_id)
    return {
        "success": True,
        "session_id": session_id,
        "workflow_state": cp.state.value,
        "security_outcome": cp.security_outcome,
        "assessment": assessment.model_dump(mode="json") if assessment else None,
    }


@router.get("/session/{session_id}/assessments/{assessment_id}")
async def get_assessment(session_id: str, assessment_id: str):
    _require_gate()
    _session_or_404(session_id)
    ledger = get_ledger()
    assessment = ledger.load_assessment(session_id, assessment_id)
    if assessment is None:
        raise HTTPException(status_code=404, detail="assessment not found")
    return {
        "success": True,
        "session_id": session_id,
        "assessment": assessment.model_dump(mode="json"),
        "reviews": [r.model_dump(mode="json") for r in ledger.reviews_for_assessment(session_id, assessment_id)],
    }


@router.post("/session/{session_id}/assessments/{assessment_id}/review")
async def review(session_id: str, assessment_id: str, body: ReviewBody, request: Request):
    _require_gate()
    _session_or_404(session_id)
    principal = _principal(request)
    try:
        outcome = review_assessment(
            session_id,
            assessment_id,
            "approved" if body.decision == "approve" else "denied",
            principal,
            note=body.note,
        )
    except ReviewError as exc:
        raise HTTPException(status_code=exc.status_code, detail={"code": exc.code, "message": exc.detail})
    return {"success": True, "session_id": session_id, **review_summary(outcome)}
