"""
Approval endpoints (Phase 1.4).

    POST /session/{session_id}/intents/{intent_id}/approve
    POST /session/{session_id}/intents/{intent_id}/reject
    POST /session/{session_id}/intents/{intent_id}/request-changes
    GET  /session/{session_id}/intents/pending

The principal comes from ``backend.config.auth_mode`` (dev_header mode:
``X-Helix-User`` is required → 401 otherwise). The body may pin the
``execution_intent_hash`` and ``plan_hash`` the client reviewed; a mismatch
with the staged records is a 409 and nothing is recorded.

Approving does not execute. It records the ``HumanApproval`` and moves the
checkpoint to ``READY_TO_EXECUTE``; execution still goes through ``/execute``
(and, from Phase 3A, the Execution Fabric), which re-verifies the approval.
"""

from __future__ import annotations

from typing import Optional

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel, Field

from backend.config.auth_mode import AuthConfigurationError, principal_from_headers
from backend.config.feature_flags import science_gate_enabled
from backend.orchestration.approval_service import (
    ApprovalError,
    approval_summary,
    decide_intent,
)
from backend.orchestration.ledger import get_ledger

router = APIRouter(tags=["approvals"])


class ApprovalBody(BaseModel):
    execution_intent_hash: Optional[str] = Field(default=None, description="Hash the client reviewed; 409 on mismatch.")
    plan_hash: Optional[str] = Field(default=None, description="Plan hash the client reviewed; 409 on mismatch.")
    note: Optional[str] = Field(default=None, max_length=2000)
    selected_recommendation_id: Optional[str] = None


def _require_gate() -> None:
    if not science_gate_enabled():
        raise HTTPException(status_code=404, detail="approval endpoints require HELIX_SCIENCE_GATE_V1=1")


def _principal(request: Request):
    try:
        principal = principal_from_headers(request.headers)
    except AuthConfigurationError as exc:
        raise HTTPException(status_code=500, detail=str(exc))
    if principal is None:
        raise HTTPException(status_code=401, detail="missing principal (X-Helix-User header in dev_header mode)")
    return principal


def _decide(session_id: str, intent_id: str, decision: str, body: Optional[ApprovalBody], request: Request):
    _require_gate()
    from backend.history_manager import history_manager

    if not history_manager.get_session(session_id):
        raise HTTPException(status_code=404, detail=f"Session {session_id} not found")
    principal = _principal(request)
    body = body or ApprovalBody()
    try:
        outcome = decide_intent(
            session_id,
            intent_id,
            decision,  # type: ignore[arg-type]
            principal,
            expected_intent_hash=body.execution_intent_hash,
            expected_plan_hash=body.plan_hash,
            note=body.note,
            selected_recommendation_id=body.selected_recommendation_id,
        )
    except ApprovalError as exc:
        raise HTTPException(status_code=exc.status_code, detail={"code": exc.code, "message": exc.detail})
    return {"success": True, "session_id": session_id, **approval_summary(outcome)}


@router.post("/session/{session_id}/intents/{intent_id}/approve")
async def approve_intent(session_id: str, intent_id: str, request: Request, body: Optional[ApprovalBody] = None):
    return _decide(session_id, intent_id, "approved", body, request)


@router.post("/session/{session_id}/intents/{intent_id}/reject")
async def reject_intent(session_id: str, intent_id: str, request: Request, body: Optional[ApprovalBody] = None):
    return _decide(session_id, intent_id, "rejected", body, request)


@router.post("/session/{session_id}/intents/{intent_id}/request-changes")
async def request_changes(session_id: str, intent_id: str, request: Request, body: Optional[ApprovalBody] = None):
    return _decide(session_id, intent_id, "changes_requested", body, request)


@router.get("/session/{session_id}/intents/pending")
async def pending_intent(session_id: str):
    """What a client needs to render an approval card: the staged intent, plan hash and objective."""
    _require_gate()
    from backend.history_manager import history_manager

    if not history_manager.get_session(session_id):
        raise HTTPException(status_code=404, detail=f"Session {session_id} not found")
    checkpoint = history_manager.load_checkpoint(session_id)
    if not checkpoint.pending_execution_intent_id:
        return {"success": True, "session_id": session_id, "pending": None, "workflow_state": checkpoint.state.value}
    ledger = get_ledger()
    intent = ledger.load_intent(session_id, checkpoint.pending_execution_intent_id)
    plan = next(
        (p for p in ledger.list_plans(session_id) if p.plan_id == checkpoint.pending_plan_id and p.plan_hash == checkpoint.pending_plan_hash),
        None,
    )
    objective = ledger.load_objective(session_id, checkpoint.objective_id) if checkpoint.objective_id else None
    return {
        "success": True,
        "session_id": session_id,
        "workflow_state": checkpoint.state.value,
        "approval_id": checkpoint.approval_id,
        "pending": {
            "intent": intent.model_dump(mode="json") if intent else None,
            "plan": plan.model_dump(mode="json", exclude={"ir"}) if plan else None,
            "objective": objective.model_dump(mode="json") if objective else None,
        },
    }
