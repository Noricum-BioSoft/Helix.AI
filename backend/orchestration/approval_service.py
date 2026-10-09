"""
Approval service (Phase 1.4 / 1.5).

One code path for every way a human approves: the HTTP router
(``backend/api/approvals.py``, principal from headers) and the natural-language
"yes"/"approve" path in ``/execute`` (principal ``identity_provider=
"natural_language"``) both call ``decide_intent``. Both perform the same
checks:

* the intent being approved is the one currently staged on the checkpoint,
* ``execution_intent_hash`` matches the staged intent,
* ``plan_hash`` matches the staged plan **and** the pending plan dict still
  hashes to it (a plan that changed since staging cannot be approved),
* a decision has not already been recorded for the intent.

``verify_approved_for_execution`` is the pre-execution re-check the broker
path runs (invariant: ``approval.execution_intent_hash == intent.execution_intent_hash``).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Dict, Optional

from backend.contracts.human_approval import ApprovalDecision, HumanApproval, Principal
from backend.contracts.ids import new_id
from backend.contracts.scientific_plan import ScientificPlan
from backend.orchestration import invariants
from backend.orchestration.invariants import InvariantViolation
from backend.orchestration.ledger import LocalLedger, get_ledger
from backend.orchestration.plan_staging import plan_ir_from_dict
from backend.contracts.scientific_plan import compute_plan_hash
from backend.workflow_checkpoint import WorkflowCheckpoint, WorkflowState

logger = logging.getLogger(__name__)


class ApprovalError(Exception):
    """Base class; ``status_code`` maps to HTTP in the router."""

    status_code = 400

    def __init__(self, detail: str, *, code: str):
        super().__init__(detail)
        self.detail = detail
        self.code = code


class NothingToApprove(ApprovalError):
    status_code = 404


class StaleApproval(ApprovalError):
    """Hash mismatch: plan or intent changed since staging."""

    status_code = 409


class AlreadyDecided(ApprovalError):
    status_code = 409


class MissingPrincipal(ApprovalError):
    status_code = 401


@dataclass(frozen=True)
class ApprovalOutcome:
    approval: HumanApproval
    checkpoint: WorkflowCheckpoint
    plan: ScientificPlan


def natural_language_principal(session_id: str, command: str) -> Principal:
    """Principal for the chat approval path. Lowest-trust identity, recorded as such."""
    return Principal(
        subject_id=f"session:{session_id}",
        display_name="chat user",
        identity_provider="natural_language",
        auth_method="chat_message",
        roles=[],
        audit_signature=None,
    )


def current_plan_hash_of_pending(checkpoint: WorkflowCheckpoint, staged_plan: ScientificPlan) -> Optional[str]:
    """Hash of the plan dict as it is *now* on the checkpoint, computed the same way as at staging."""
    pending = checkpoint.pending_plan
    plan_dict = pending.get("plan") if isinstance(pending, dict) else None
    if not isinstance(plan_dict, dict):
        return None
    try:
        ir = plan_ir_from_dict(plan_dict)
    except Exception:
        return None
    # Rebuild steps the way ScientificPlan.from_plan_ir did, preserving kinds/capabilities from the staged plan.
    staged_by_id = {s.id: s for s in staged_plan.steps}
    steps = []
    for ir_step in ir.steps:
        staged = staged_by_id.get(ir_step.id)
        steps.append(
            staged.model_copy(
                update={
                    "action_type": ir_step.action_type,
                    "tool_name": ir_step.tool_name,
                    "arguments": dict(ir_step.arguments),
                    "description": ir_step.description,
                }
            )
            if staged
            else None
        )
    if any(s is None for s in steps) or len(steps) != len(staged_plan.steps):
        return "__structure_changed__"
    return compute_plan_hash(steps, ir, staged_plan.security_requirements)


def _load_staged(session_id: str, checkpoint: WorkflowCheckpoint, ledger: LocalLedger):
    if not checkpoint.pending_execution_intent_id or not checkpoint.pending_plan_id:
        raise NothingToApprove("no execution intent is staged for this session", code="no_pending_intent")
    intent = ledger.load_intent(session_id, checkpoint.pending_execution_intent_id)
    if intent is None:
        raise NothingToApprove("staged execution intent not found in ledger", code="intent_missing")
    plan = None
    for candidate in ledger.list_plans(session_id):
        if candidate.plan_id == checkpoint.pending_plan_id and candidate.plan_hash == checkpoint.pending_plan_hash:
            plan = candidate
            break
    if plan is None:
        raise NothingToApprove("staged plan not found in ledger", code="plan_missing")
    return intent, plan


def decide_intent(
    session_id: str,
    execution_intent_id: str,
    decision: ApprovalDecision,
    principal: Optional[Principal],
    *,
    expected_intent_hash: Optional[str] = None,
    expected_plan_hash: Optional[str] = None,
    note: Optional[str] = None,
    selected_recommendation_id: Optional[str] = None,
    checkpoint: Optional[WorkflowCheckpoint] = None,
    ledger: Optional[LocalLedger] = None,
    save_checkpoint: bool = True,
) -> ApprovalOutcome:
    """Record a human decision for the staged intent and transition the checkpoint.

    ``expected_*_hash`` are what the caller believes it is approving (HTTP
    body); when omitted (NL path) the staged values are used, but the pending
    plan dict is still re-hashed so a mutated plan cannot be approved.
    """
    from backend.history_manager import history_manager

    if principal is None:
        raise MissingPrincipal("a principal is required to record an approval decision", code="principal_required")
    ledger = ledger or get_ledger()
    checkpoint = checkpoint or history_manager.load_checkpoint(session_id)

    intent, plan = _load_staged(session_id, checkpoint, ledger)
    if intent.execution_intent_id != execution_intent_id:
        raise StaleApproval(
            f"intent {execution_intent_id} is not the staged intent ({intent.execution_intent_id})", code="intent_not_staged"
        )
    if ledger.approvals_for_intent(session_id, intent.execution_intent_id):
        raise AlreadyDecided("a decision has already been recorded for this execution intent", code="already_decided")

    intent_hash = expected_intent_hash or intent.execution_intent_hash
    plan_hash = expected_plan_hash or plan.plan_hash
    if not intent.matches(intent_hash):
        raise StaleApproval("execution_intent_hash does not match the staged intent", code="intent_hash_mismatch")
    if plan_hash != plan.plan_hash or plan.plan_hash != checkpoint.pending_plan_hash:
        raise StaleApproval("plan_hash does not match the staged plan", code="plan_hash_mismatch")
    live_hash = current_plan_hash_of_pending(checkpoint, plan)
    if live_hash is not None and live_hash != plan.plan_hash:
        raise StaleApproval("the pending plan changed after it was staged; re-plan before approving", code="plan_changed")

    approval = HumanApproval(
        approval_id=new_id(),
        trace_id=intent.trace_id,
        execution_intent_id=intent.execution_intent_id,
        execution_intent_hash=intent.execution_intent_hash,
        plan_id=plan.plan_id,
        plan_hash=plan.plan_hash,
        decision=decision,
        principal=principal,
        note=note,
        selected_recommendation_id=selected_recommendation_id,
    )
    # Invariants as code, even though the construction above already guarantees them.
    invariants.check_approval_matches_plan(approval, plan.plan_hash)
    invariants.check_approval_matches_intent(approval, intent)
    invariants.check_same_trace(approval, intent, plan)

    ledger.record_approval(session_id, approval)

    if decision == "approved":
        new_cp = checkpoint.with_platform_records(approval_id=approval.approval_id)
        new_cp = new_cp.transition(WorkflowState.READY_TO_EXECUTE) if checkpoint.state == WorkflowState.WAITING_FOR_APPROVAL else new_cp
        plan_status = "approved"
    elif decision == "rejected":
        new_cp = WorkflowCheckpoint.idle().with_platform_records(
            trace_id=checkpoint.trace_id, objective_id=checkpoint.objective_id, approval_id=approval.approval_id
        )
        plan_status = "rejected"
    else:  # changes_requested: stay staged; the next staging creates a new plan version
        new_cp = checkpoint.with_platform_records(approval_id=approval.approval_id)
        plan_status = "draft"

    try:
        ledger.record_plan(session_id, plan.model_copy(update={"status": plan_status}))
    except Exception as exc:  # status is informational; the approval record is the source of truth
        logger.debug("approval_service: could not update plan status: %s", exc)

    if save_checkpoint:
        history_manager.save_checkpoint(session_id, new_cp)
    logger.info(
        "[platform] approval trace=%s intent=%s decision=%s by=%s/%s",
        approval.trace_id, intent.execution_intent_id, decision, principal.identity_provider, principal.subject_id,
    )
    return ApprovalOutcome(approval=approval, checkpoint=new_cp, plan=plan)


def verify_approved_for_execution(
    session_id: str,
    checkpoint: WorkflowCheckpoint,
    *,
    ledger: Optional[LocalLedger] = None,
) -> HumanApproval:
    """Pre-execution re-check: the staged intent has an *approved* HumanApproval bound to its hash.

    Raises ``InvariantViolation`` (fail closed). Phase 3A moves this into
    ``ExecutionFabric.run()``; until then the broker path in ``/execute`` calls it.
    """
    ledger = ledger or get_ledger()
    if not checkpoint.pending_execution_intent_id:
        raise InvariantViolation("approval_required", "no execution intent is staged; nothing can have been approved")
    intent = ledger.load_intent(session_id, checkpoint.pending_execution_intent_id)
    if intent is None:
        raise InvariantViolation("approval_required", "staged execution intent not found")
    approval = ledger.load_approval(session_id, checkpoint.approval_id) if checkpoint.approval_id else None
    if approval is None:
        approvals = [a for a in ledger.approvals_for_intent(session_id, intent.execution_intent_id) if a.approved]
        approval = approvals[-1] if approvals else None
    invariants.check_approval_is_approved(approval)
    invariants.check_approval_matches_intent(approval, intent)  # type: ignore[arg-type]
    if checkpoint.pending_plan_hash and approval.plan_hash != checkpoint.pending_plan_hash:  # type: ignore[union-attr]
        raise InvariantViolation("approval_plan_hash", "approval was given for a different plan version")
    return approval  # type: ignore[return-value]


def approval_summary(outcome: ApprovalOutcome) -> Dict[str, Any]:
    a = outcome.approval
    return {
        "approval_id": a.approval_id,
        "trace_id": a.trace_id,
        "execution_intent_id": a.execution_intent_id,
        "execution_intent_hash": a.execution_intent_hash,
        "plan_id": a.plan_id,
        "plan_hash": a.plan_hash,
        "decision": a.decision,
        "principal": a.principal.model_dump(mode="json"),
        "timestamp": a.timestamp.isoformat(),
        "workflow_state": outcome.checkpoint.state.value,
    }
