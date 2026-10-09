"""
Phase 3A — ExecutionBroker → Fabric delegation seam (``try_fabric_execution``).

The broker calls this first. With the flag off, or when the session has no
approved READY_TO_EXECUTE intent, it returns ``None`` and the broker keeps its
legacy path (so existing flows are untouched). With the flag on and an approved
intent staged, it runs the tool through the fabric.
"""

from __future__ import annotations

import pytest

from backend.contracts.human_approval import HumanApproval, Principal
from backend.contracts.ids import new_id
from backend.execution.integration import try_fabric_execution
from backend.orchestration.intent_builder import build_intent
from backend.orchestration.ledger import LocalLedger
from backend.orchestration.objective_builder import build_objective_deterministic
from backend.orchestration.plan_staging import build_scientific_plan
from backend.workflow_checkpoint import WorkflowCheckpoint, WorkflowState


@pytest.fixture(autouse=True)
def _isolate(tmp_path, monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")
    monkeypatch.setattr("backend.main.dispatch_tool", _async_ok("tool ran"))
    from backend.history_manager import history_manager

    history_manager.storage_dir = tmp_path / "sessions"
    history_manager.storage_dir.mkdir(parents=True, exist_ok=True)
    history_manager.sessions = {}
    history_manager._sessions_loaded = True
    yield


def _async_ok(text: str):
    async def _dispatch(tool_name, arguments):
        return {"status": "success", "text": text, "downloadable_artifacts": []}

    return _dispatch


def _stage_ready_to_execute(session_id: str):
    from backend.history_manager import history_manager

    ledger = LocalLedger(history_manager.storage_dir)
    cmd = "Run bulk RNA-seq differential expression analysis"
    objective = build_objective_deterministic(cmd, {"session_id": session_id}, None)
    ledger.record_objective(session_id, objective)
    plan = build_scientific_plan(
        {"version": "v1", "steps": [{"id": "s1", "action_type": "run_analysis", "tool_name": "bulk_rnaseq_analysis", "arguments": {}}]},
        objective, cmd,
    )
    ledger.record_plan(session_id, plan)
    intent = build_intent(plan, None, capability_id="local_compute:bulk_rnaseq_analysis")
    ledger.record_intent(session_id, intent)
    principal = Principal(subject_id="alice", identity_provider="dev_header", auth_method="header", roles=["approver"])
    approval = HumanApproval(
        approval_id=new_id(), trace_id=intent.trace_id, execution_intent_id=intent.execution_intent_id,
        execution_intent_hash=intent.execution_intent_hash, plan_id=plan.plan_id, plan_hash=plan.plan_hash,
        decision="approved", principal=principal,
    )
    ledger.record_approval(session_id, approval)
    cp = WorkflowCheckpoint.waiting_for_approval(pending_plan={"plan": {}, "command": cmd}).with_platform_records(
        trace_id=intent.trace_id, objective_id=objective.objective_id,
        pending_plan_id=plan.plan_id, pending_plan_hash=plan.plan_hash,
        pending_execution_intent_id=intent.execution_intent_id, pending_execution_intent_hash=intent.execution_intent_hash,
        approval_id=approval.approval_id,
    ).transition(WorkflowState.READY_TO_EXECUTE)
    history_manager.save_checkpoint(session_id, cp)
    return intent


def test_flag_off_returns_none(monkeypatch):
    monkeypatch.delenv("HELIX_EXECUTION_FABRIC_V1", raising=False)
    _stage_ready_to_execute("sess-off")
    assert try_fabric_execution("bulk_rnaseq_analysis", {}, {"session_id": "sess-off"}) is None


def test_flag_on_no_intent_returns_none(monkeypatch):
    monkeypatch.setenv("HELIX_EXECUTION_FABRIC_V1", "1")
    assert try_fabric_execution("bulk_rnaseq_analysis", {}, {"session_id": "sess-none"}) is None


def test_flag_on_ready_intent_delegates(monkeypatch):
    monkeypatch.setenv("HELIX_EXECUTION_FABRIC_V1", "1")
    _stage_ready_to_execute("sess-on")
    result = try_fabric_execution("bulk_rnaseq_analysis", {}, {"session_id": "sess-on"})
    assert result is not None
    assert result["via"] == "execution_fabric"
    assert result["status"] == "success"


def test_flag_on_mismatched_tool_returns_none(monkeypatch):
    """A tool that does not match the staged intent's capability is not delegated."""
    monkeypatch.setenv("HELIX_EXECUTION_FABRIC_V1", "1")
    _stage_ready_to_execute("sess-mismatch")
    assert try_fabric_execution("single_cell_analysis", {}, {"session_id": "sess-mismatch"}) is None
