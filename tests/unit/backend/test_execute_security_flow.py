"""Phase 2 — Secure Science gate, end-to-end through ``/execute``.

Route-level flows: an unsafe upload DENYs a staged plan, a dual-use request
routes to manual review, a reviewer approve rebuilds the intent (review is not
approval — human approval still follows), a reviewer deny is terminal, PHI
tightens the envelope but still allows, and the gate off changes nothing.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import pytest
from fastapi.testclient import TestClient

from backend.orchestration.ledger import LocalLedger

CLEAN_CMD = "Run bulk RNA-seq differential expression analysis on s3://bucket/counts.csv"
OTHER_CMD = "Run bulk RNA-seq differential expression analysis on s3://bucket/other_counts.csv"
DUAL_USE_CMD = "Run bulk RNA-seq differential expression analysis for a gain-of-function study on s3://bucket/counts.csv"

SCIENTIST_HDR = {"X-Helix-User": "alice", "X-Helix-Roles": "scientist,approver"}
REVIEWER_HDR = {"X-Helix-User": "rev", "X-Helix-Roles": "security_reviewer"}


@pytest.fixture(autouse=True)
def _isolate(tmp_path, monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")
    monkeypatch.setenv("HELIX_DEMO_MODE", "0")
    monkeypatch.setenv("HELIX_SANDBOX_HOST_FALLBACK", "1")
    monkeypatch.setenv("HELIX_SCIENCE_GATE_V1", "1")
    monkeypatch.delenv("HELIX_AUTH_MODE", raising=False)
    monkeypatch.delenv("HELIX_ENV", raising=False)
    monkeypatch.delenv("HELIX_SCREENING_ADAPTER", raising=False)

    from backend.history_manager import history_manager

    history_manager.storage_dir = tmp_path / "sessions"
    history_manager.storage_dir.mkdir(parents=True, exist_ok=True)
    history_manager.sessions = {}
    history_manager._sessions_loaded = True

    import backend.main as _m

    _m._daily_prompt_counters.clear()

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


def _post(client: TestClient, command: str, session_id: Optional[str] = None) -> Dict[str, Any]:
    payload: Dict[str, Any] = {"command": command}
    if session_id:
        payload["session_id"] = session_id
    r = client.post("/execute", json=payload)
    assert r.status_code == 200, r.text
    return r.json()


def _platform(out: Dict[str, Any]) -> Dict[str, Any]:
    return out.get("platform") or (out.get("result") or {}).get("platform") or {}


def _new_session(client: TestClient) -> str:
    """Stage a clean plan to obtain a session id (no uploads yet)."""
    out = _post(client, CLEAN_CMD)
    assert out["status"] == "workflow_planned", out
    return out["session_id"]


def _inject_upload(session_id: str, sensitivity: str, flags: List[str] = (), policy_state: str = "cleared") -> None:
    from backend.history_manager import history_manager

    session = history_manager.get_session(session_id)
    assert session is not None
    meta = session.setdefault("metadata", {})
    meta.setdefault("uploaded_files", []).append(
        {
            "name": "data.csv",
            "file_id": "f-data.csv",
            "policy_state": policy_state,
            "intake_policy": {"sensitivity_class": sensitivity, "scan_flags": list(flags)},
        }
    )


# ── DENY ───────────────────────────────────────────────────────────────────


def test_unsafe_upload_denies_plan_and_blocks_execution(client, ledger, _isolate):
    sid = _new_session(client)
    _inject_upload(sid, "internal", flags=["suspicious_payload_pattern"])

    out = _post(client, OTHER_CMD, sid)
    assert out["status"] == "security_denied", out
    plat = _platform(out)
    assert plat["blocked"] is True and plat["security_outcome"] == "DENY"
    assert "execution_intent_id" not in plat and "review_url" not in plat

    from backend.history_manager import history_manager

    assert history_manager.load_checkpoint(sid).state.value == "DENIED"

    # execute_plan=True cannot bypass a denied plan
    r = client.post("/execute", json={"command": OTHER_CMD, "session_id": sid, "execute_plan": True})
    assert r.status_code == 200
    assert r.json()["status"] == "execution_blocked", r.json()
    assert _isolate["dispatch"] == 0


# ── REQUIRE_REVIEW → approve ─────────────────────────────────────────────────


def test_dual_use_routes_to_review_then_approve_rebuilds_intent(client, ledger, _isolate):
    out = _post(client, DUAL_USE_CMD)
    sid = out["session_id"]
    assert out["status"] == "security_review_required", out
    plat = _platform(out)
    assert plat["blocked"] is True and plat["security_outcome"] == "REQUIRE_REVIEW"
    assert plat["review_url"].endswith("/review")
    assessment_id = plat["assessment_id"]

    from backend.history_manager import history_manager

    assert history_manager.load_checkpoint(sid).state.value == "WAITING_FOR_SECURITY_REVIEW"

    # the assessment is durable and fetchable
    pending = client.get(f"/session/{sid}/assessments/pending").json()
    assert pending["assessment"]["assessment_id"] == assessment_id
    assert pending["workflow_state"] == "WAITING_FOR_SECURITY_REVIEW"

    # a non-reviewer cannot resolve it
    forbidden = client.post(
        f"/session/{sid}/assessments/{assessment_id}/review", json={"decision": "approve"}, headers=SCIENTIST_HDR
    )
    assert forbidden.status_code == 403, forbidden.text

    # a reviewer approve rebuilds the intent and returns to human approval
    r = client.post(
        f"/session/{sid}/assessments/{assessment_id}/review",
        json={"decision": "approve", "note": "defensive context, cleared"},
        headers=REVIEWER_HDR,
    )
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["decision"] == "approved" and body["resulting_outcome"] == "ALLOW_WITH_APPROVAL"
    assert body["execution_intent_id"]
    assert body["workflow_state"] == "WAITING_FOR_APPROVAL"

    # review is not approval — human approval still executes
    approve = _post(client, "Approve.", sid)
    assert approve["status"] in {"success", "pipeline_submitted", "pipeline_executed", "job"}, approve
    assert _isolate["dispatch"] >= 1


# ── REQUIRE_REVIEW → deny ────────────────────────────────────────────────────


def test_dual_use_review_deny_is_terminal(client, ledger, _isolate):
    out = _post(client, DUAL_USE_CMD)
    sid = out["session_id"]
    assessment_id = _platform(out)["assessment_id"]

    r = client.post(
        f"/session/{sid}/assessments/{assessment_id}/review",
        json={"decision": "deny", "note": "not cleared"},
        headers=REVIEWER_HDR,
    )
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["decision"] == "denied" and body["resulting_outcome"] == "DENY"
    assert body["execution_intent_id"] is None
    assert body["workflow_state"] == "DENIED"

    from backend.history_manager import history_manager

    assert history_manager.load_checkpoint(sid).state.value == "DENIED"
    # execution is barred even if the client forces execute_plan
    r2 = client.post("/execute", json={"command": DUAL_USE_CMD, "session_id": sid, "execute_plan": True})
    assert r2.status_code == 200 and r2.json()["status"] == "execution_blocked", r2.json()
    assert _isolate["dispatch"] == 0


def test_review_requires_principal(client):
    out = _post(client, DUAL_USE_CMD)
    sid = out["session_id"]
    assessment_id = _platform(out)["assessment_id"]
    r = client.post(f"/session/{sid}/assessments/{assessment_id}/review", json={"decision": "approve"})
    assert r.status_code == 401


# ── ALLOW_WITH_APPROVAL (PHI) ────────────────────────────────────────────────


def test_phi_upload_tightens_envelope_but_allows(client, ledger, _isolate):
    sid = _new_session(client)
    _inject_upload(sid, "restricted_human_data")

    out = _post(client, OTHER_CMD, sid)
    assert out["status"] == "workflow_planned", out
    plat = _platform(out)
    assert plat["blocked"] is False
    assert plat["security_outcome"] == "ALLOW_WITH_APPROVAL"
    assert plat["policy_envelope"]["external_execution_allowed"] is False
    assert plat["execution_intent_id"]

    approve = _post(client, "Approve.", sid)
    assert approve["status"] in {"success", "pipeline_submitted", "pipeline_executed", "job"}, approve
    assert _isolate["dispatch"] >= 1


# ── gate off ─────────────────────────────────────────────────────────────────


def test_gate_off_ignores_security(client, monkeypatch):
    monkeypatch.setenv("HELIX_SCIENCE_GATE_V1", "0")
    out = _post(client, DUAL_USE_CMD)
    assert out["status"] == "workflow_planned", out
    assert "platform" not in out and "platform" not in (out.get("result") or {})
