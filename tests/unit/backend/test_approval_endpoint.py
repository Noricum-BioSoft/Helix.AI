"""Phase 1.4/1.5 — approvals: HTTP endpoint and the NL path share one service and the same hash checks.

Route-level tests keep ``/execute`` orchestration real (mock-mode router from
conftest) and fake tool execution.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional

import pytest
from fastapi.testclient import TestClient

from backend.orchestration.ledger import LocalLedger

PLAN_CMD = "Run bulk RNA-seq differential expression analysis on s3://bucket/counts.csv"
HDR = {"X-Helix-User": "alice", "X-Helix-Roles": "scientist,approver"}


@pytest.fixture(autouse=True)
def _isolate(tmp_path, monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")
    monkeypatch.setenv("HELIX_DEMO_MODE", "0")
    monkeypatch.setenv("HELIX_SANDBOX_HOST_FALLBACK", "1")
    monkeypatch.setenv("HELIX_SCIENCE_GATE_V1", "1")
    monkeypatch.delenv("HELIX_AUTH_MODE", raising=False)
    monkeypatch.delenv("HELIX_ENV", raising=False)

    from backend.history_manager import history_manager

    history_manager.storage_dir = tmp_path / "sessions"
    history_manager.storage_dir.mkdir(parents=True, exist_ok=True)
    history_manager.sessions = {}
    history_manager._sessions_loaded = True

    import backend.main as _m

    _m._daily_prompt_counters.clear()

    # The unit conftest defaults the staging classifier to "execute directly";
    # these tests need the Plan → Approve → Execute gate, so stage everything.
    from backend.orchestration.approval_policy import StagingDecision

    monkeypatch.setattr(
        "backend.orchestration.approval_policy._classify_staging_intent",
        lambda *_a, **_k: StagingDecision(
            requires_approval=True, has_execute_intent=False, is_planning_request=True, method="test", reason="stage"
        ),
    )

    calls: Dict[str, int] = {"dispatch": 0}

    async def _fake_dispatch(tool: str, params: dict) -> dict:
        calls["dispatch"] += 1
        return {"status": "success", "text": f"mock {tool} complete", "run_id": "run-1"}

    monkeypatch.setattr("backend.main.dispatch_tool", _fake_dispatch)

    class _FakeBroker:
        async def execute_tool(self, request):
            calls["dispatch"] += 1
            return {"status": "success", "text": "plan executed", "type": "job", "run_id": "run-1"}

    monkeypatch.setattr("backend.main._get_execution_broker", lambda: _FakeBroker())
    yield calls


@pytest.fixture()
def client() -> TestClient:
    from backend.main import app

    return TestClient(app)


@pytest.fixture()
def ledger(tmp_path) -> LocalLedger:
    return LocalLedger(tmp_path / "sessions")


def _stage(client: TestClient, command: str = PLAN_CMD, session_id: Optional[str] = None) -> Dict[str, Any]:
    payload: Dict[str, Any] = {"command": command}
    if session_id:
        payload["session_id"] = session_id
    r = client.post("/execute", json=payload)
    assert r.status_code == 200, r.text
    out = r.json()
    assert out["status"] == "workflow_planned", out
    platform = out.get("platform") or (out.get("result") or {}).get("platform")
    assert platform, f"gate on but no platform block in response: {json.dumps(out)[:500]}"
    return {"session_id": out["session_id"], **platform}


# ── staging ──────────────────────────────────────────────────────────────────


def test_staging_persists_records_sharing_one_trace(client, ledger):
    staged = _stage(client)
    sid = staged["session_id"]
    assert staged["trace_id"].startswith("orch_")
    by_trace = ledger.records_by_trace(sid, staged["trace_id"])
    assert {k: len(v) for k, v in by_trace.items()} == {
        "objectives": 1, "plans": 1, "assessments": 1, "security_reviews": 0, "intents": 1, "approvals": 0,
    }
    assert by_trace["intents"][0].execution_intent_hash == staged["execution_intent_hash"]
    assert by_trace["plans"][0].plan_hash == staged["plan_hash"]

    from backend.history_manager import history_manager

    cp = history_manager.load_checkpoint(sid)
    assert cp.state.value == "WAITING_FOR_APPROVAL"
    assert cp.pending_execution_intent_id == staged["execution_intent_id"]
    assert cp.pending_plan_hash == staged["plan_hash"]
    assert cp.approval_id is None

    pending = client.get(f"/session/{sid}/intents/pending").json()
    assert pending["pending"]["intent"]["execution_intent_id"] == staged["execution_intent_id"]
    assert pending["pending"]["objective"]["trace_id"] == staged["trace_id"]


def test_gate_off_leaves_execute_unchanged(client, monkeypatch, ledger):
    monkeypatch.setenv("HELIX_SCIENCE_GATE_V1", "0")
    r = client.post("/execute", json={"command": PLAN_CMD}).json()
    assert r["status"] == "workflow_planned"
    assert "platform" not in r and "platform" not in (r.get("result") or {})
    sid = r["session_id"]
    assert not (Path(ledger.storage_dir) / sid / "platform").exists()
    assert client.get(f"/session/{sid}/intents/pending").status_code == 404
    approve = client.post("/execute", json={"command": "Approve.", "session_id": sid}).json()
    assert approve["status"] in {"success", "pipeline_submitted", "pipeline_executed", "job"}, approve


# ── HTTP endpoint ────────────────────────────────────────────────────────────


def test_http_approve_records_principal_and_transitions(client, ledger):
    staged = _stage(client)
    sid, iid = staged["session_id"], staged["execution_intent_id"]
    r = client.post(
        f"/session/{sid}/intents/{iid}/approve",
        json={"execution_intent_hash": staged["execution_intent_hash"], "plan_hash": staged["plan_hash"], "note": "LGTM"},
        headers=HDR,
    )
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["decision"] == "approved"
    assert body["workflow_state"] == "READY_TO_EXECUTE"
    assert body["principal"] == {
        "subject_id": "alice",
        "display_name": "alice",
        "identity_provider": "dev_header",
        "auth_method": "header",
        "roles": ["scientist", "approver"],
        "tenant_id": None,
        "project_id": None,
        "audit_signature": None,
    }
    approval = ledger.load_approval(sid, body["approval_id"])
    assert approval is not None and approval.approved
    assert approval.execution_intent_hash == staged["execution_intent_hash"]
    assert approval.plan_hash == staged["plan_hash"]
    assert approval.trace_id == staged["trace_id"]
    assert approval.note == "LGTM"
    # plan status follows the decision
    plan = ledger.load_plan(sid, staged["plan_id"], 1)
    assert plan.status == "approved"


def test_http_reject_and_request_changes(client, ledger):
    staged = _stage(client)
    sid, iid = staged["session_id"], staged["execution_intent_id"]
    r = client.post(f"/session/{sid}/intents/{iid}/reject", json={"note": "wrong file"}, headers=HDR)
    assert r.status_code == 200, r.text
    assert r.json()["decision"] == "rejected" and r.json()["workflow_state"] == "IDLE"
    assert ledger.load_plan(sid, staged["plan_id"], 1).status == "rejected"

    staged2 = _stage(client)
    sid2, iid2 = staged2["session_id"], staged2["execution_intent_id"]
    r2 = client.post(f"/session/{sid2}/intents/{iid2}/request-changes", json={"note": "use 8 threads"}, headers=HDR)
    assert r2.status_code == 200, r2.text
    assert r2.json()["decision"] == "changes_requested" and r2.json()["workflow_state"] == "WAITING_FOR_APPROVAL"


def test_http_intent_hash_mismatch_is_409_and_records_nothing(client, ledger):
    staged = _stage(client)
    sid, iid = staged["session_id"], staged["execution_intent_id"]
    r = client.post(f"/session/{sid}/intents/{iid}/approve", json={"execution_intent_hash": "f" * 64}, headers=HDR)
    assert r.status_code == 409, r.text
    assert r.json()["detail"]["code"] == "intent_hash_mismatch"
    assert ledger.approvals_for_intent(sid, iid) == []


def test_http_plan_hash_mismatch_is_409(client, ledger):
    staged = _stage(client)
    sid, iid = staged["session_id"], staged["execution_intent_id"]
    r = client.post(f"/session/{sid}/intents/{iid}/approve", json={"plan_hash": "0" * 64}, headers=HDR)
    assert r.status_code == 409, r.text
    assert r.json()["detail"]["code"] == "plan_hash_mismatch"
    assert ledger.approvals_for_intent(sid, iid) == []


def test_http_wrong_or_stale_intent_id_is_409(client):
    staged = _stage(client)
    sid = staged["session_id"]
    r = client.post(f"/session/{sid}/intents/not-the-staged-intent/approve", headers=HDR)
    assert r.status_code == 409
    assert r.json()["detail"]["code"] == "intent_not_staged"
    # Re-planning supersedes the intent: the old id can no longer be approved
    staged2 = _stage(client, "Run bulk RNA-seq differential expression analysis on s3://bucket/other_counts.csv", sid)
    assert staged2["execution_intent_id"] != staged["execution_intent_id"]
    r2 = client.post(f"/session/{sid}/intents/{staged['execution_intent_id']}/approve", headers=HDR)
    assert r2.status_code == 409


def test_http_requires_principal(client):
    staged = _stage(client)
    sid, iid = staged["session_id"], staged["execution_intent_id"]
    r = client.post(f"/session/{sid}/intents/{iid}/approve")
    assert r.status_code == 401


def test_http_double_decision_is_409(client):
    staged = _stage(client)
    sid, iid = staged["session_id"], staged["execution_intent_id"]
    assert client.post(f"/session/{sid}/intents/{iid}/approve", headers=HDR).status_code == 200
    r = client.post(f"/session/{sid}/intents/{iid}/reject", headers=HDR)
    assert r.status_code == 409 and r.json()["detail"]["code"] == "already_decided"


def test_http_unknown_session_is_404(client):
    assert client.post("/session/nope/intents/x/approve", headers=HDR).status_code == 404


# ── natural-language path ────────────────────────────────────────────────────


def test_nl_approval_yields_identical_record_and_executes(client, ledger, _isolate):
    staged = _stage(client)
    sid, iid = staged["session_id"], staged["execution_intent_id"]
    r = client.post("/execute", json={"command": "Approve.", "session_id": sid}).json()
    assert r["status"] in {"success", "pipeline_submitted", "pipeline_executed", "job"}, r
    assert _isolate["dispatch"] >= 1

    approvals = ledger.approvals_for_intent(sid, iid)
    assert len(approvals) == 1
    a = approvals[0]
    assert a.approved
    assert a.execution_intent_hash == staged["execution_intent_hash"]
    assert a.plan_hash == staged["plan_hash"]
    assert a.trace_id == staged["trace_id"]
    assert a.principal.identity_provider == "natural_language"
    assert a.principal.auth_method == "chat_message"
    assert a.principal.subject_id == f"session:{sid}"
    # same contract as the HTTP path: every field the HTTP record has, this one has
    assert set(a.model_dump()) == set(
        __import__("backend.contracts.human_approval", fromlist=["HumanApproval"]).HumanApproval.model_fields
    )


def test_nl_approval_refuses_plan_changed_after_staging(client, ledger, _isolate):
    staged = _stage(client)
    sid, iid = staged["session_id"], staged["execution_intent_id"]
    # Mutate the pending plan behind the user's back (what a buggy revision path could do)
    from backend.history_manager import history_manager

    cp = history_manager.load_checkpoint(sid)
    cp.pending_plan["plan"]["steps"][0]["arguments"]["counts_file"] = "s3://bucket/EVIL.fastq.gz"
    history_manager.save_checkpoint(sid, cp)
    pending = history_manager.sessions[sid].get("pending_plan")
    if isinstance(pending, dict):
        pending["plan"]["steps"][0]["arguments"]["counts_file"] = "s3://bucket/EVIL.fastq.gz"

    r = client.post("/execute", json={"command": "yes", "session_id": sid}).json()
    assert r["status"] == "approval_blocked", r
    assert r["raw_result"]["approval_error"] == "plan_changed"
    assert ledger.approvals_for_intent(sid, iid) == []
    assert _isolate["dispatch"] == 0
    assert history_manager.load_checkpoint(sid).state.value == "WAITING_FOR_APPROVAL"


def test_nl_approval_refuses_unrecorded_plan(client, ledger, _isolate):
    """A plan staged without an intent (e.g. gate flipped on mid-session) cannot be approved by chat."""
    staged = _stage(client)
    sid = staged["session_id"]
    from backend.history_manager import history_manager

    cp = history_manager.load_checkpoint(sid)
    history_manager.save_checkpoint(sid, cp.with_platform_records(pending_execution_intent_id=None))
    r = client.post("/execute", json={"command": "Approve.", "session_id": sid}).json()
    assert r["status"] == "approval_blocked", r
    assert _isolate["dispatch"] == 0


def test_http_approval_then_nl_execution_does_not_record_twice(client, ledger, _isolate):
    staged = _stage(client)
    sid, iid = staged["session_id"], staged["execution_intent_id"]
    assert client.post(f"/session/{sid}/intents/{iid}/approve", headers=HDR).status_code == 200
    r = client.post("/execute", json={"command": "Approve.", "session_id": sid}).json()
    assert r["status"] in {"success", "pipeline_submitted", "pipeline_executed", "job"}, r
    assert len(ledger.approvals_for_intent(sid, iid)) == 1
    assert ledger.approvals_for_intent(sid, iid)[0].principal.identity_provider == "dev_header"


def test_broker_precheck_blocks_rejected_intent(client, ledger, _isolate):
    staged = _stage(client)
    sid, iid = staged["session_id"], staged["execution_intent_id"]
    assert client.post(f"/session/{sid}/intents/{iid}/reject", headers=HDR).status_code == 200
    # Force the legacy pending plan to still exist and the checkpoint to look approvable
    from backend.history_manager import history_manager
    from backend.workflow_checkpoint import WorkflowCheckpoint, WorkflowState

    cp = history_manager.load_checkpoint(sid)
    rejected_id = cp.approval_id
    forged = WorkflowCheckpoint.waiting_for_approval(pending_plan=history_manager.sessions[sid]["pending_plan"]).with_platform_records(
        trace_id=staged["trace_id"],
        objective_id=staged["objective_id"],
        pending_plan_id=staged["plan_id"],
        pending_plan_hash=staged["plan_hash"],
        pending_execution_intent_id=iid,
        pending_execution_intent_hash=staged["execution_intent_hash"],
        approval_id=rejected_id,
    )
    history_manager.save_checkpoint(sid, forged)
    r = client.post("/execute", json={"command": "Approve.", "session_id": sid}).json()
    assert r["status"] in {"approval_blocked", "execution_blocked"}, r
    assert _isolate["dispatch"] == 0
