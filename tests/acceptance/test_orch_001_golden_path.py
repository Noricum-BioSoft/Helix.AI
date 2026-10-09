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


# ── Security gate (P2) ───────────────────────────────────────────────────────


def test_orch_001_security_assessment_with_envelope(loop):
    """The staged plan carries a Secure Science assessment, and the intent is bound to it."""
    staged = loop["staged"]
    assessment = staged.assessment
    assert assessment is not None
    assert assessment.outcome == "ALLOW"
    assert assessment.plan_hash == staged.plan.plan_hash
    assert assessment.trace_id == staged.plan.trace_id
    assert assessment.policy_envelope is not None
    assert {p.policy_id for p in assessment.applied_policies} >= {"action_classification", "data_sensitivity"}
    # the intent records which assessment cleared it
    assert staged.intent is not None
    assert staged.intent.assessment_id == assessment.assessment_id
    assert staged.intent.assessment_hash == assessment.assessment_hash
    assert loop["ledger"].load_assessment(loop["sid"], assessment.assessment_id) == assessment


def test_orch_001_denied_state_blocks_intent(tmp_path, monkeypatch):
    """A DENY assessment yields no ExecutionIntent, and one cannot be forced."""
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")
    from backend.history_manager import history_manager

    history_manager.storage_dir = tmp_path / "sessions"
    history_manager.storage_dir.mkdir(parents=True, exist_ok=True)
    history_manager.sessions = {}
    history_manager._sessions_loaded = True

    from backend.orchestration.invariants import InvariantViolation
    from backend.orchestration.ledger import LocalLedger
    from backend.orchestration.plan_staging import build_intent_for, stage_plan
    from backend.security.checks.base import AssessmentContext

    ledger = LocalLedger(history_manager.storage_dir)
    context = AssessmentContext(
        uploaded_files=[
            {
                "name": "data.csv",
                "file_id": "f-data.csv",
                "policy_state": "cleared",
                "intake_policy": {"sensitivity_class": "internal", "scan_flags": ["suspicious_payload_pattern"]},
            }
        ]
    )
    staged = stage_plan("orch-deny", COMMAND, PLAN, {"session_id": "orch-deny"}, ledger=ledger, assessment_context=context)
    assert staged.security_outcome == "DENY"
    assert staged.blocked and staged.intent is None
    # the assessment is still recorded for audit
    assert ledger.load_assessment("orch-deny", staged.assessment.assessment_id) == staged.assessment
    # and nothing can turn a DENY into an intent
    with pytest.raises(InvariantViolation):
        build_intent_for(staged.plan, staged.assessment, None)


# ── Execution fabric (P3A): capability → provider → execution → result ───────


def _p3a_fabric(ledger, profile_name: str, *, runner=None):
    from backend.config.execution_profile import load_execution_profile
    from backend.execution.fabric import ExecutionFabric
    from backend.execution.providers.factory import build_providers
    from backend.execution.providers.local_compute import LocalComputeProvider
    from backend.execution.registry import CapabilityRegistry

    profile = load_execution_profile(profile_name, check_adapters=False)
    registry = CapabilityRegistry.from_config(profile)
    overrides = {}
    if runner is not None:
        overrides["local_compute"] = LocalComputeProvider(
            [d for d in registry.all() if d.provider == "local_compute"], {}, tool_runner=runner
        )
    providers = build_providers(profile, registry, **overrides)
    return registry, ExecutionFabric(registry, providers, ledger)


def _p3a_request(intent, approval):
    from backend.contracts.execution_request import ExecutionRequest, make_idempotency_key
    from backend.contracts.ids import new_id

    return ExecutionRequest(
        execution_request_id=new_id(), trace_id=intent.trace_id,
        execution_intent_id=intent.execution_intent_id, execution_intent_hash=intent.execution_intent_hash,
        approval_id=(approval.approval_id if approval else None),
        idempotency_key=make_idempotency_key(intent.execution_intent_hash),
        provider_id=intent.provider_id, capability_id=intent.capability_id,
    )


def _p3a_approval(intent, plan):
    from backend.contracts.human_approval import HumanApproval, Principal
    from backend.contracts.ids import new_id

    principal = Principal(subject_id="rev", identity_provider="dev_header", auth_method="header", roles=["approver"])
    return HumanApproval(
        approval_id=new_id(), trace_id=intent.trace_id, execution_intent_id=intent.execution_intent_id,
        execution_intent_hash=intent.execution_intent_hash, plan_id=plan.plan_id, plan_hash=plan.plan_hash,
        decision="approved", principal=principal,
    )


def test_orch_001_capability_resolved_from_registry(loop):
    _, fabric = _p3a_fabric(loop["ledger"], "local-only")
    cap_id = loop["staged"].intent.capability_id
    assert cap_id == "local_compute:bulk_rnaseq_analysis"
    assert fabric.registry.get(cap_id) is not None
    assert fabric.registry.resolve_alias("bulk_rnaseq_analysis") == cap_id


def test_orch_001_provider_executes_local(loop):
    staged = loop["staged"]
    approval = _p3a_approval(staged.intent, staged.plan)
    _, fabric = _p3a_fabric(loop["ledger"], "local-only", runner=lambda t, a: {"status": "success", "text": "ok"})
    run = fabric.run(loop["sid"], staged.intent, approval, _p3a_request(staged.intent, approval))
    assert run.status == "succeeded" and run.trace_id == staged.intent.trace_id


def test_orch_001_provider_executes_mock_experimental_with_approval(tmp_path, monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")
    from backend.contracts.ids import new_trace_id
    from backend.orchestration.intent_builder import build_intent
    from backend.orchestration.ledger import LocalLedger
    from backend.orchestration.objective_builder import build_objective_deterministic
    from backend.orchestration.plan_staging import build_scientific_plan

    ledger = LocalLedger(tmp_path)
    trace = new_trace_id()
    cmd = "Express the protein construct"
    objective = build_objective_deterministic(cmd, {"session_id": "lab"}, None, trace_id=trace)
    plan = build_scientific_plan(
        {"version": "v1", "steps": [{"id": "s1", "action_type": "run_analysis", "tool_name": "protein_expression", "arguments": {"construct": "pX-1"}}]},
        objective, cmd,
    )
    intent = build_intent(plan, None, capability_id="mock_experimental:protein_expression")
    approval = _p3a_approval(intent, plan)
    _, fabric = _p3a_fabric(ledger, "local-only")
    run = fabric.run("lab", intent, approval, _p3a_request(intent, approval))
    assert run.status == "succeeded"
    assert run.outputs and run.outputs[0].uri.startswith("mock://")


def test_orch_001_retry_after_timeout_creates_no_duplicate(loop):
    staged = loop["staged"]
    approval = _p3a_approval(staged.intent, staged.plan)
    calls = {"n": 0}

    def runner(tool, args):
        calls["n"] += 1
        return {"status": "success", "text": "ok"}

    _, fabric = _p3a_fabric(loop["ledger"], "local-only", runner=runner)
    req = _p3a_request(staged.intent, approval)
    run1 = fabric.run(loop["sid"], staged.intent, approval, req)
    run2 = fabric.run(loop["sid"], staged.intent, approval, req.model_copy(update={"execution_request_id": "req-2"}))
    assert run1.execution_run_id == run2.execution_run_id and calls["n"] == 1


# ── Steps and the phase that un-skips them. Keep this list in sync with the plan.
STEPS = [
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
