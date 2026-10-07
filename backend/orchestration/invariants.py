"""
Execution-state invariants, as code.

Each invariant is a function that raises ``InvariantViolation``. The services
owning a transition call these; tests call them directly. Phase 0 defines
them; later phases wire them into checkpoint transitions, the intent builder,
the provider authorizer and the Execution Fabric (see plan rev 2 table).
"""

from __future__ import annotations

from typing import Iterable, Optional

from backend.contracts.execution_intent import ExecutionIntent
from backend.contracts.execution_recommendation import ExecutionRecommendation
from backend.contracts.execution_request import ExecutionRequest
from backend.contracts.human_approval import HumanApproval
from backend.contracts.scientific_objective import ScientificObjective
from backend.contracts.security_assessment import PolicyEnvelope, ScreeningResult

# Checkpoint states that must never lead to execution / executable intents.
# Kept as strings so this module does not import WorkflowCheckpoint (P2 adds the states there).
STATE_DENIED = "DENIED"
STATE_WAITING_FOR_SECURITY_REVIEW = "WAITING_FOR_SECURITY_REVIEW"
STATE_READY_TO_EXECUTE = "READY_TO_EXECUTE"
STATE_EXECUTING = "EXECUTING"
TERMINAL_BLOCKED_STATES = frozenset({STATE_DENIED})
NO_INTENT_STATES = frozenset({STATE_DENIED, STATE_WAITING_FOR_SECURITY_REVIEW})
EXECUTION_STATES = frozenset({STATE_READY_TO_EXECUTE, STATE_EXECUTING})


class InvariantViolation(RuntimeError):
    def __init__(self, invariant: str, detail: str):
        self.invariant = invariant
        self.detail = detail
        super().__init__(f"[{invariant}] {detail}")


def check_state_transition_allowed(current_state: str, new_state: str) -> None:
    """DENIED → cannot transition to execution."""
    if current_state in TERMINAL_BLOCKED_STATES and new_state in EXECUTION_STATES:
        raise InvariantViolation("denied_never_executes", f"{current_state} -> {new_state} is forbidden")


def check_intent_creation_allowed(current_state: str) -> None:
    """WAITING_FOR_SECURITY_REVIEW / DENIED → cannot create an executable ExecutionIntent."""
    if current_state in NO_INTENT_STATES:
        raise InvariantViolation("no_intent_in_blocked_state", f"cannot build ExecutionIntent while {current_state}")


def check_approval_matches_plan(approval: HumanApproval, current_plan_hash: str) -> None:
    """approval.plan_hash must match the current plan hash."""
    if approval.plan_hash != current_plan_hash:
        raise InvariantViolation("approval_plan_hash", "approval was given for a different plan version")


def check_approval_matches_intent(approval: HumanApproval, intent: ExecutionIntent) -> None:
    """approval.execution_intent_hash must match the submitted execution intent."""
    if approval.execution_intent_id != intent.execution_intent_id or not intent.matches(approval.execution_intent_hash):
        raise InvariantViolation("approval_intent_hash", "approval does not bind to this execution intent")


def check_request_matches_intent(request: ExecutionRequest, intent: ExecutionIntent) -> None:
    if request.execution_intent_id != intent.execution_intent_id or not intent.matches(request.execution_intent_hash):
        raise InvariantViolation("request_intent_hash", "execution request does not match the execution intent")
    if request.provider_id != intent.provider_id or request.capability_id != intent.capability_id:
        raise InvariantViolation("request_intent_provider", "request provider/capability differ from the intent")


def check_approval_is_approved(approval: Optional[HumanApproval]) -> None:
    """HumanApproval.decision != approved → ExecutionFabric.run() must fail closed."""
    if approval is None or not approval.approved:
        raise InvariantViolation("approval_required", "execution requires an approved HumanApproval")


def check_provider_in_candidate_set(intent: ExecutionIntent, recommendation: ExecutionRecommendation) -> None:
    """ExecutionIntent.provider_id must appear in the approved recommendation candidate set."""
    if intent.recommendation_id != recommendation.recommendation_id:
        raise InvariantViolation("intent_recommendation", "intent references a different recommendation")
    if intent.provider_id not in recommendation.candidate_set:
        raise InvariantViolation("provider_not_candidate", f"{intent.provider_id} not in recommendation candidate_set")


def check_selected_alternative(approval: HumanApproval, recommendation: ExecutionRecommendation) -> None:
    if approval.selected_recommendation_id and approval.selected_recommendation_id not in recommendation.candidate_set:
        raise InvariantViolation("selected_not_candidate", "approver selected a provider outside the candidate set")


def check_screening_satisfied(
    envelope: PolicyEnvelope, screening_results: Iterable[ScreeningResult], review_override: bool = False
) -> None:
    """screening_required → no execution without an accepted screening result or an authorized review override."""
    if not envelope.screening_required or review_override:
        return
    if not any(r.status == "accepted" for r in screening_results):
        raise InvariantViolation("screening_required", "no accepted screening result and no authorized review override")


def check_experimental_requires_approval(category: str, approval: Optional[HumanApproval]) -> None:
    if category == "experimental" and (approval is None or not approval.approved):
        raise InvariantViolation("experimental_requires_approval", "experimental execution requires an approved HumanApproval")


def check_idempotency_unique(idempotency_key: str, existing_keys: Iterable[str]) -> None:
    """Same idempotency_key → at most one provider submission. Caller re-attaches instead of resubmitting."""
    if idempotency_key in set(existing_keys):
        raise InvariantViolation("duplicate_submission", "a submission with this idempotency_key already exists")


def check_proposed_objective_not_activated(objective: ScientificObjective, user_action: bool) -> None:
    """Decision.proposed_objective stays 'proposed' until explicit user action."""
    if objective.parent_decision_id and objective.status != "proposed" and not user_action:
        raise InvariantViolation("proposed_objective_auto_activation", "a decision-proposed objective was activated without user action")


def check_same_trace(*records) -> None:
    traces = {getattr(r, "trace_id", None) for r in records if r is not None}
    if len(traces) > 1:
        raise InvariantViolation("trace_consistency", f"records span multiple trace_ids: {sorted(map(str, traces))}")
