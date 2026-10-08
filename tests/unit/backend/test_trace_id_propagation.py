"""Phase 1 — one trace_id per loop, carried by every record and by the checkpoint across transitions."""

from __future__ import annotations

import pytest

from backend.contracts.human_approval import Principal
from backend.contracts.ids import is_trace_id
from backend.orchestration import invariants
from backend.orchestration.approval_service import (
    AlreadyDecided,
    StaleApproval,
    decide_intent,
    verify_approved_for_execution,
)
from backend.orchestration.ledger import LocalLedger
from backend.orchestration.plan_staging import stage_plan
from backend.workflow_checkpoint import WorkflowCheckpoint, WorkflowState

PLAN = {
    "version": "v1",
    "steps": [{"id": "step1", "tool_name": "bulk_rnaseq_analysis", "arguments": {"counts_file": "s3://b/counts.csv"}}],
}
PRINCIPAL = Principal(subject_id="alice", identity_provider="dev_header", auth_method="header")


@pytest.fixture(autouse=True)
def _isolate(tmp_path, monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")
    from backend.history_manager import history_manager

    history_manager.storage_dir = tmp_path / "sessions"
    history_manager.storage_dir.mkdir(parents=True, exist_ok=True)
    history_manager.sessions = {}
    history_manager._sessions_loaded = True


@pytest.fixture()
def ledger(tmp_path):
    return LocalLedger(tmp_path / "sessions")


def _staged_checkpoint(ledger, sid="s1"):
    staged = stage_plan(sid, "Run DE", PLAN, {"session_id": sid}, ledger=ledger)
    cp = WorkflowCheckpoint.waiting_for_approval(pending_plan={"plan": PLAN, "command": "Run DE"}).with_platform_records(
        **staged.checkpoint_fields()
    )
    return staged, cp


def test_all_records_of_a_loop_share_the_trace(ledger):
    staged, cp = _staged_checkpoint(ledger)
    outcome = decide_intent("s1", staged.intent.execution_intent_id, "approved", PRINCIPAL, checkpoint=cp, ledger=ledger, save_checkpoint=False)
    assert is_trace_id(staged.trace_id)
    invariants.check_same_trace(staged.objective, staged.plan, staged.intent, outcome.approval)
    by_trace = ledger.records_by_trace("s1", staged.trace_id)
    assert {k: len(v) for k, v in by_trace.items()} == {
        "objectives": 1, "plans": 1, "assessments": 1, "security_reviews": 0, "intents": 1, "approvals": 1,
    }
    assert outcome.checkpoint.trace_id == staged.trace_id
    assert outcome.checkpoint.state == WorkflowState.READY_TO_EXECUTE
    # the broker pre-check accepts this checkpoint
    approval = verify_approved_for_execution("s1", outcome.checkpoint, ledger=ledger)
    assert approval.approval_id == outcome.approval.approval_id


def test_checkpoint_platform_fields_round_trip_and_stay_absent_for_legacy():
    legacy = WorkflowCheckpoint.waiting_for_approval(pending_plan={"plan": PLAN})
    d = legacy.to_dict()
    assert not any(k in d for k in WorkflowCheckpoint.PLATFORM_FIELDS)  # byte-identical legacy checkpoints
    assert WorkflowCheckpoint.from_dict(d).trace_id is None

    cp = legacy.with_platform_records(trace_id="orch_" + "a" * 32, pending_plan_hash="h")
    rt = WorkflowCheckpoint.from_dict(cp.to_dict())
    assert rt.trace_id == "orch_" + "a" * 32 and rt.pending_plan_hash == "h"
    # transitions keep the trace
    assert rt.transition(WorkflowState.EXECUTING).trace_id == rt.trace_id
    with pytest.raises(ValueError):
        legacy.with_platform_records(not_a_field="x")


def test_trace_survives_a_revision(ledger):
    v1 = stage_plan("s1", "Run DE", PLAN, ledger=ledger)
    changed = {"version": "v1", "steps": [{"id": "step1", "tool_name": "bulk_rnaseq_analysis", "arguments": {"counts_file": "s3://b/other.csv"}}]}
    v2 = stage_plan("s1", "Run DE on the other file", changed, supersedes=v1.plan, ledger=ledger)
    assert v2.trace_id == v1.trace_id
    by_trace = ledger.records_by_trace("s1", v1.trace_id)
    assert len(by_trace["plans"]) == 2 and len(by_trace["intents"]) == 2 and len(by_trace["objectives"]) == 1


def test_decision_service_fails_closed(ledger):
    staged, cp = _staged_checkpoint(ledger)
    iid = staged.intent.execution_intent_id
    with pytest.raises(StaleApproval):
        decide_intent("s1", iid, "approved", PRINCIPAL, expected_intent_hash="0" * 64, checkpoint=cp, ledger=ledger, save_checkpoint=False)
    with pytest.raises(StaleApproval):
        decide_intent("s1", iid, "approved", PRINCIPAL, expected_plan_hash="0" * 64, checkpoint=cp, ledger=ledger, save_checkpoint=False)
    assert ledger.approvals_for_intent("s1", iid) == []
    # nothing approved yet → the broker pre-check refuses
    with pytest.raises(invariants.InvariantViolation):
        verify_approved_for_execution("s1", cp, ledger=ledger)
    # rejected → recorded, but execution still refused
    outcome = decide_intent("s1", iid, "rejected", PRINCIPAL, checkpoint=cp, ledger=ledger, save_checkpoint=False)
    assert not outcome.approval.approved and outcome.checkpoint.state == WorkflowState.IDLE
    with pytest.raises(invariants.InvariantViolation):
        verify_approved_for_execution("s1", cp.with_platform_records(approval_id=outcome.approval.approval_id), ledger=ledger)
    with pytest.raises(AlreadyDecided):
        decide_intent("s1", iid, "approved", PRINCIPAL, checkpoint=cp, ledger=ledger, save_checkpoint=False)
