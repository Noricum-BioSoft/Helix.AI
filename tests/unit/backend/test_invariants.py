from __future__ import annotations

import pytest

from backend.contracts.execution_intent import ExecutionIntent
from backend.contracts.execution_recommendation import CriteriaScores, ExecutionRecommendation
from backend.contracts.execution_request import ExecutionRequest, make_idempotency_key
from backend.contracts.human_approval import HumanApproval, Principal
from backend.contracts.ids import new_id, new_trace_id, stable_hash
from backend.contracts.scientific_objective import ScientificObjective
from backend.contracts.security_assessment import PolicyEnvelope, ScreeningResult
from backend.orchestration import invariants as inv

TRACE = new_trace_id()


def _intent(**over) -> ExecutionIntent:
    d = dict(execution_intent_id=new_id(), trace_id=TRACE, plan_id="p", plan_hash="ph", provider_id="local_compute",
             capability_id="local_compute:fastqc", input_manifest_hash=stable_hash([1]), execution_parameters_hash=stable_hash({}))
    d.update(over)
    return ExecutionIntent(**d)


def _approval(intent: ExecutionIntent, decision="approved", **over) -> HumanApproval:
    d = dict(approval_id=new_id(), trace_id=TRACE, execution_intent_id=intent.execution_intent_id,
             execution_intent_hash=intent.execution_intent_hash, plan_id=intent.plan_id, plan_hash=intent.plan_hash,
             decision=decision, principal=Principal(subject_id="u", identity_provider="dev_header", auth_method="header"))
    d.update(over)
    return HumanApproval(**d)


def test_denied_never_transitions_to_execution():
    with pytest.raises(inv.InvariantViolation, match="denied_never_executes"):
        inv.check_state_transition_allowed("DENIED", "READY_TO_EXECUTE")
    inv.check_state_transition_allowed("DENIED", "IDLE")
    inv.check_state_transition_allowed("WAITING_FOR_APPROVAL", "READY_TO_EXECUTE")


def test_no_intent_in_review_or_denied_state():
    for s in ("DENIED", "WAITING_FOR_SECURITY_REVIEW"):
        with pytest.raises(inv.InvariantViolation, match="no_intent_in_blocked_state"):
            inv.check_intent_creation_allowed(s)
    inv.check_intent_creation_allowed("PLANNING")


def test_approval_binds_to_plan_and_intent():
    intent = _intent()
    approval = _approval(intent)
    inv.check_approval_matches_plan(approval, "ph")
    inv.check_approval_matches_intent(approval, intent)
    with pytest.raises(inv.InvariantViolation, match="approval_plan_hash"):
        inv.check_approval_matches_plan(approval, "other")
    changed = _intent(execution_intent_id=intent.execution_intent_id, provider_id="carolina_cloud")
    with pytest.raises(inv.InvariantViolation, match="approval_intent_hash"):
        inv.check_approval_matches_intent(approval, changed)


def test_unapproved_or_missing_approval_fails_closed():
    intent = _intent()
    with pytest.raises(inv.InvariantViolation, match="approval_required"):
        inv.check_approval_is_approved(None)
    with pytest.raises(inv.InvariantViolation, match="approval_required"):
        inv.check_approval_is_approved(_approval(intent, decision="rejected"))
    inv.check_approval_is_approved(_approval(intent))


def test_request_must_match_intent():
    intent = _intent()
    ok = ExecutionRequest(execution_request_id=new_id(), trace_id=TRACE, execution_intent_id=intent.execution_intent_id,
                          execution_intent_hash=intent.execution_intent_hash, idempotency_key=make_idempotency_key(intent.execution_intent_hash),
                          provider_id=intent.provider_id, capability_id=intent.capability_id)
    inv.check_request_matches_intent(ok, intent)
    bad = ok.model_copy(update={"provider_id": "aws_emr"})
    with pytest.raises(inv.InvariantViolation, match="request_intent_provider"):
        inv.check_request_matches_intent(bad, intent)


def test_provider_must_be_in_candidate_set_and_selected_alternative_too():
    scores = CriteriaScores(**{k: 0.5 for k in CriteriaScores.model_fields})
    rec = ExecutionRecommendation(recommendation_id="r1", trace_id=TRACE, plan_id="p", provider_id="local_compute",
                                  capability_id="c", candidate_set=["local_compute", "nextflow"], criteria_scores=scores,
                                  score=0.5, confidence=0.9)
    inv.check_provider_in_candidate_set(_intent(recommendation_id="r1"), rec)
    with pytest.raises(inv.InvariantViolation, match="provider_not_candidate"):
        inv.check_provider_in_candidate_set(_intent(recommendation_id="r1", provider_id="aws_emr"), rec)
    with pytest.raises(inv.InvariantViolation, match="intent_recommendation"):
        inv.check_provider_in_candidate_set(_intent(recommendation_id="r2"), rec)
    intent = _intent(recommendation_id="r1")
    inv.check_selected_alternative(_approval(intent, selected_recommendation_id="nextflow"), rec)
    with pytest.raises(inv.InvariantViolation, match="selected_not_candidate"):
        inv.check_selected_alternative(_approval(intent, selected_recommendation_id="aws_emr"), rec)


def test_screening_required_needs_accepted_result_or_override():
    env = PolicyEnvelope(screening_required=True)
    with pytest.raises(inv.InvariantViolation, match="screening_required"):
        inv.check_screening_satisfied(env, [])
    with pytest.raises(inv.InvariantViolation):
        inv.check_screening_satisfied(env, [ScreeningResult(provider="mock", status="flagged")])
    inv.check_screening_satisfied(env, [ScreeningResult(provider="mock", status="accepted")])
    inv.check_screening_satisfied(env, [], review_override=True)
    inv.check_screening_satisfied(PolicyEnvelope(screening_required=False), [])


def test_experimental_requires_approval_and_idempotency_unique():
    with pytest.raises(inv.InvariantViolation, match="experimental_requires_approval"):
        inv.check_experimental_requires_approval("experimental", None)
    inv.check_experimental_requires_approval("computational", None)
    with pytest.raises(inv.InvariantViolation, match="duplicate_submission"):
        inv.check_idempotency_unique("k1", ["k0", "k1"])
    inv.check_idempotency_unique("k2", ["k0", "k1"])


def test_proposed_objective_and_trace_consistency():
    proposed = ScientificObjective(objective_id=new_id(), trace_id=TRACE, objective="o", question="q",
                                   status="active", parent_decision_id="d1")
    with pytest.raises(inv.InvariantViolation, match="proposed_objective_auto_activation"):
        inv.check_proposed_objective_not_activated(proposed, user_action=False)
    inv.check_proposed_objective_not_activated(proposed, user_action=True)
    a, b = _intent(), _intent(trace_id=new_trace_id())
    inv.check_same_trace(a, _approval(a), None)
    with pytest.raises(inv.InvariantViolation, match="trace_consistency"):
        inv.check_same_trace(a, b)
