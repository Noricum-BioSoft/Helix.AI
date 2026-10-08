"""
ORCH-001 — Scientific Orchestration Golden Path (acceptance test).

Created in Phase 0 as a skeleton: every step is skipped with the phase that
will implement it. The skip list may only shrink; by Phase 9 nothing is
skipped and the whole loop runs on deterministic fixtures in CI.

Flow: Objective → persisted Plan → Security Assessment (envelope) →
Capability resolution → Recommendation → Provider Authorization →
ExecutionIntent → Human Approval → ExecutionProvider → Result →
persist to DataWeaver → Observation → Interpretation → Decision.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.orch001

# Deterministic fixture plan used by the implemented steps.
PLAN = {
    "version": "v1",
    "steps": [
        {
            "id": "step1",
            "action_type": "run_analysis",
            "tool_name": "bulk_rnaseq_analysis",
            "arguments": {"counts_file": "s3://fixtures/counts.csv", "design": "~condition"},
            "description": "Differential expression, treated vs control",
        }
    ],
}
COMMAND = "Run differential expression on the counts matrix, treated vs control"


@pytest.fixture()
def loop(tmp_path, monkeypatch):
    """Shared state for the implemented steps: one session, one ledger, one trace."""
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")
    from backend.history_manager import history_manager

    history_manager.storage_dir = tmp_path / "sessions"
    history_manager.storage_dir.mkdir(parents=True, exist_ok=True)
    history_manager.sessions = {}
    history_manager._sessions_loaded = True

    from backend.orchestration.ledger import LocalLedger
    from backend.orchestration.plan_staging import stage_plan
    from backend.workflow_checkpoint import WorkflowCheckpoint

    ledger = LocalLedger(history_manager.storage_dir)
    staged = stage_plan("orch", COMMAND, PLAN, {"session_id": "orch"}, ledger=ledger)
    checkpoint = WorkflowCheckpoint.waiting_for_approval(pending_plan={"plan": PLAN, "command": COMMAND}).with_platform_records(
        **staged.checkpoint_fields()
    )
    return {"sid": "orch", "ledger": ledger, "staged": staged, "checkpoint": checkpoint}


# ── Implemented steps (P1) ───────────────────────────────────────────────────


def test_orch_001_objective_created(loop):
    from backend.contracts.ids import is_trace_id

    objective = loop["staged"].objective
    assert is_trace_id(objective.trace_id)
    assert objective.objective and objective.question
    assert objective.status == "active"
    assert loop["ledger"].load_objective(loop["sid"], objective.objective_id) == objective


def test_orch_001_plan_persisted_with_hash(loop):
    from backend.contracts.rationale import RationaleItem

    plan = loop["staged"].plan
    stored = loop["ledger"].load_plan(loop["sid"], plan.plan_id, plan.version)
    assert stored is not None and stored.plan_hash == plan.plan_hash and len(plan.plan_hash) == 64
    assert stored.objective_id == loop["staged"].objective.objective_id
    assert stored.trace_id == loop["staged"].objective.trace_id
    assert stored.plan_rationale and all(isinstance(r, RationaleItem) for r in stored.plan_rationale)
    assert stored.ir.steps[0].tool_name == "bulk_rnaseq_analysis"


def test_orch_001_execution_intent_built(loop):
    from backend.orchestration.intent_builder import build_intent

    intent = loop["staged"].intent
    assert intent.plan_hash == loop["staged"].plan.plan_hash
    assert intent.trace_id == loop["staged"].plan.trace_id
    assert intent.provider_id == "legacy:Local" and intent.capability_id == "local_compute:bulk_rnaseq_analysis"
    assert len(intent.execution_intent_hash) == 64
    # same plan on another provider is a different intent
    class _Infra:
        infrastructure = "EMR"

    assert build_intent(loop["staged"].plan, _Infra()).execution_intent_hash != intent.execution_intent_hash
    assert loop["ledger"].load_intent(loop["sid"], intent.execution_intent_id) == intent


def test_orch_001_approval_bound_to_intent_hash(loop):
    from backend.contracts.human_approval import Principal
    from backend.orchestration import invariants
    from backend.orchestration.approval_service import StaleApproval, decide_intent, verify_approved_for_execution

    sid, ledger, staged, cp = loop["sid"], loop["ledger"], loop["staged"], loop["checkpoint"]
    principal = Principal(subject_id="reviewer", identity_provider="dev_header", auth_method="header", roles=["approver"])
    with pytest.raises(StaleApproval):
        decide_intent(sid, staged.intent.execution_intent_id, "approved", principal, expected_intent_hash="0" * 64, checkpoint=cp, ledger=ledger, save_checkpoint=False)
    outcome = decide_intent(sid, staged.intent.execution_intent_id, "approved", principal, checkpoint=cp, ledger=ledger, save_checkpoint=False)
    approval = outcome.approval
    assert approval.execution_intent_hash == staged.intent.execution_intent_hash
    assert approval.plan_hash == staged.plan.plan_hash
    invariants.check_same_trace(staged.objective, staged.plan, staged.intent, approval)
    assert verify_approved_for_execution(sid, outcome.checkpoint, ledger=ledger) == approval
    # a different intent (e.g. provider changed) is not covered by this approval
    class _Infra:
        infrastructure = "EMR"

    from backend.orchestration.intent_builder import build_intent

    with pytest.raises(invariants.InvariantViolation):
        invariants.check_approval_matches_intent(approval, build_intent(staged.plan, _Infra()))


# ── Steps and the phase that un-skips them. Keep this list in sync with the plan.
STEPS = [
    ("security_assessment_with_envelope", "P2"),
    ("denied_state_blocks_intent", "P2"),
    ("capability_resolved_from_registry", "P3A"),
    ("provider_executes_local", "P3A"),
    ("provider_executes_mock_experimental_with_approval", "P3A"),
    ("retry_after_timeout_creates_no_duplicate", "P3A"),
    ("recommendation_with_candidate_set", "P4"),
    ("provider_authorization_against_envelope", "P4"),
    ("intent_provider_in_candidate_set", "P4"),
    ("one_real_computational_backend", "P3B/P3A-nextflow"),
    ("persisted_to_dataweaver", "P7"),
    ("ids_stable_across_stores", "P7"),
    ("trace_query_returns_full_loop", "P7"),
    ("observation_has_no_claims", "P8"),
    ("decision_references_evidence_ids", "P8"),
    ("proposed_objective_not_auto_activated", "P8"),
    ("provenance_from_decision_to_inputs", "P8"),
]


@pytest.mark.parametrize("step,phase", STEPS, ids=[s for s, _ in STEPS])
def test_orch_001_step(step: str, phase: str) -> None:
    pytest.skip(f"ORCH-001 step '{step}' is implemented in {phase}")


def test_orch_001_contracts_importable() -> None:
    """The only non-skipped assertion in Phase 0: the contract surface exists."""
    from backend.contracts import (  # noqa: F401
        evidence_assessment,
        execution_intent,
        execution_recommendation,
        execution_request,
        human_approval,
        provenance,
        provider_authorization,
        scientific_objective,
        scientific_plan,
        security_assessment,
    )
    from backend.execution.providers.base import ExecutionProvider  # noqa: F401
    from backend.orchestration import invariants  # noqa: F401
    from shared.capability_registry import CapabilityDescriptor  # noqa: F401
