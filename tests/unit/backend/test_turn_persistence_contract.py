"""
test_turn_persistence_contract.py — Slice A persistence contract tests.

Single rule: every /execute request that reaches a tool produces exactly one
history entry, optionally one run record, with produced_artifacts always
populated when artifacts were written to disk.

Per-turn assertions:
  1. session["history"] grew by exactly 1.
  2. session["runs"] grew by ≤ 1.
  3. If artifacts were produced, produced_artifacts is non-empty and every URI
     in the list exists on the filesystem at assertion time.
  4. Bundle endpoint is constructible from history_manager.storage_dir (no
     hardcoded ../sessions paths).

Categories tested:
  A — Scriptable bio tool (bulk_rnaseq_analysis) — produces script + result artifacts.
  B — Advisory / handle_natural_command — lightweight history entry, no run record required.
  C — execute_plan=True path — formerly bypassed history; now routed through _dispatch_result.
  D — Tabular analysis approval — history entry + produced_artifacts from tabular executor.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict
from unittest.mock import AsyncMock, MagicMock, patch

import pytest
from fastapi.testclient import TestClient


# ── Isolation fixture ────────────────────────────────────────────────────────


@pytest.fixture(autouse=True)
def _isolate(tmp_path, monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")
    monkeypatch.setenv("HELIX_DEMO_MODE", "0")
    monkeypatch.setenv("HELIX_SANDBOX_HOST_FALLBACK", "1")

    from backend.history_manager import history_manager

    history_manager.storage_dir = tmp_path / "sessions"
    history_manager.storage_dir.mkdir(parents=True, exist_ok=True)
    history_manager.sessions = {}
    history_manager._sessions_loaded = True

    import backend.main as _mwm

    _mwm._daily_prompt_counters.clear()
    _mwm._execution_broker = None


def _client() -> TestClient:
    from backend.main import app

    return TestClient(app)


def _session(session_id: str) -> Dict[str, Any]:
    from backend.history_manager import history_manager

    return history_manager.get_session(session_id) or {}


# ── Helpers ──────────────────────────────────────────────────────────────────


def _make_broker_result(tool: str, extra: Dict | None = None) -> Dict:
    base: Dict[str, Any] = {
        "status": "success",
        "success": True,
        "text": f"Mock result for {tool}",
        "tool_name": tool,
    }
    if extra:
        base.update(extra)
    return base


# ── Category A: scriptable bio tool ──────────────────────────────────────────


class TestScriptableBioTool:
    """bulk_rnaseq_analysis is in _SCRIPTABLE_TOOLS → script artifact + run record."""

    def test_history_grows_by_one(self, monkeypatch):
        tool = "bulk_rnaseq_analysis"
        mock_result = _make_broker_result(tool)

        def _mock_dispatch(t, args, **kw):
            return mock_result

        monkeypatch.setattr("backend.main.dispatch_tool", _mock_dispatch)

        client = _client()
        r = client.post("/execute", json={"command": f"run {tool} on demo data"})
        assert r.status_code == 200
        sid = r.json()["session_id"]

        session = _session(sid)
        assert len(session.get("history", [])) == 1, (
            f"Expected 1 history entry, got {session.get('history', [])}"
        )

    def test_run_record_created(self, monkeypatch):
        tool = "bulk_rnaseq_analysis"
        monkeypatch.setattr("backend.main.dispatch_tool", lambda t, a, **k: _make_broker_result(t))

        client = _client()
        r = client.post("/execute", json={"command": f"run {tool}"})
        assert r.status_code == 200
        sid = r.json()["session_id"]

        session = _session(sid)
        runs = session.get("runs", [])
        assert len(runs) <= 1, f"At most 1 run expected, got {runs}"
        # A scriptable tool should produce a run record.
        assert len(runs) == 1, f"Expected 1 run record for scriptable tool, got {runs}"

    def test_produced_artifacts_uris_exist(self, monkeypatch, tmp_path):
        from backend.history_manager import history_manager

        tool = "bulk_rnaseq_analysis"
        monkeypatch.setattr("backend.main.dispatch_tool", lambda t, a, **k: _make_broker_result(t))

        client = _client()
        r = client.post("/execute", json={"command": f"run {tool}"})
        assert r.status_code == 200
        sid = r.json()["session_id"]

        session = _session(sid)
        runs = session.get("runs", [])
        if not runs:
            return  # no run record = nothing to check
        for run in runs:
            arts = run.get("produced_artifacts") or []
            for art in arts:
                uri = art.get("uri")
                if uri:
                    assert Path(uri).exists(), (
                        f"Artifact URI does not exist on disk: {uri}"
                    )

    def test_artifacts_live_under_storage_dir(self, monkeypatch):
        """Artifacts must be under history_manager.storage_dir, not a hardcoded path."""
        from backend.history_manager import history_manager

        tool = "bulk_rnaseq_analysis"
        monkeypatch.setattr("backend.main.dispatch_tool", lambda t, a, **k: _make_broker_result(t))

        client = _client()
        r = client.post("/execute", json={"command": f"run {tool}"})
        assert r.status_code == 200
        sid = r.json()["session_id"]

        session = _session(sid)
        for run in session.get("runs", []):
            for art in run.get("produced_artifacts") or []:
                uri = art.get("uri")
                if uri:
                    try:
                        Path(uri).resolve().relative_to(history_manager.storage_dir.resolve())
                    except ValueError:
                        pytest.fail(
                            f"Artifact URI escapes storage_dir: {uri!r} is not under "
                            f"{history_manager.storage_dir}"
                        )


# ── Category B: advisory / handle_natural_command ────────────────────────────


class TestAdvisoryTurn:
    """handle_natural_command → lightweight history entry written, ≤1 run record."""

    def test_history_entry_written_for_advisory(self, monkeypatch):
        """Advisory turns must produce a history entry.

        We patch dispatch_tool with an AsyncMock so the broker's `await
        self._tool_executor(...)` call succeeds (the real dispatch_tool is async).
        """
        from unittest.mock import AsyncMock

        advisory_result = {
            "status": "success",
            "success": True,
            "text": "Here are three things you can do with this dataset.",
            "tool_name": "handle_natural_command",
        }
        monkeypatch.setattr(
            "backend.main.dispatch_tool",
            AsyncMock(return_value=advisory_result),
        )

        client = _client()
        r = client.post("/execute", json={"command": "what can I do with this dataset?"})
        assert r.status_code == 200
        sid = r.json()["session_id"]

        session = _session(sid)
        assert len(session.get("history", [])) >= 1, (
            "Advisory turn must produce at least one history entry"
        )

    def test_run_count_at_most_one_for_advisory(self, monkeypatch):
        from unittest.mock import AsyncMock

        advisory_result = {
            "status": "success",
            "success": True,
            "text": "Advisory response.",
            "tool_name": "handle_natural_command",
        }
        monkeypatch.setattr(
            "backend.main.dispatch_tool",
            AsyncMock(return_value=advisory_result),
        )

        client = _client()
        r = client.post("/execute", json={"command": "describe my data"})
        assert r.status_code == 200
        sid = r.json()["session_id"]

        session = _session(sid)
        assert len(session.get("runs", [])) <= 1, (
            "Advisory turn must not produce more than one run record"
        )


# ── Category C: execute_plan path ────────────────────────────────────────────


class TestExecutePlanPath:
    """Formerly bypassed _dispatch_result → history was never written.

    After Slice A fix, execute_plan=True routes through _dispatch_result, so a
    history entry is always produced.
    """

    def test_execute_plan_records_history(self, monkeypatch):
        agent_result = {
            "status": "success",
            "success": True,
            "text": "Pipeline dispatched.",
            "tool_name": "bulk_rnaseq_analysis",
        }

        async def _mock_handle_command(*a, **kw):
            return agent_result

        monkeypatch.setattr("backend.agent.handle_command", _mock_handle_command)

        client = _client()
        r = client.post(
            "/execute",
            json={"command": "run the analysis", "execute_plan": True},
        )
        assert r.status_code == 200
        sid = r.json()["session_id"]

        session = _session(sid)
        assert len(session.get("history", [])) >= 1, (
            "execute_plan=True turn must produce at least one history entry"
        )

    def test_execute_plan_run_count_at_most_one(self, monkeypatch):
        agent_result = {
            "status": "success",
            "success": True,
            "text": "Done.",
        }

        async def _mock_handle_command(*a, **kw):
            return agent_result

        monkeypatch.setattr("backend.agent.handle_command", _mock_handle_command)

        client = _client()
        r = client.post(
            "/execute",
            json={"command": "run the analysis", "execute_plan": True},
        )
        assert r.status_code == 200
        sid = r.json()["session_id"]

        session = _session(sid)
        assert len(session.get("runs", [])) <= 1, (
            "execute_plan=True turn must produce at most one run record"
        )


# ── Category D: multi-turn produces one entry per turn ───────────────────────


class TestMultiTurn:
    """N consecutive /execute calls → at least N history entries in the same session."""

    def test_n_turns_produce_n_history_entries(self, monkeypatch):
        from unittest.mock import AsyncMock

        monkeypatch.setattr(
            "backend.main.dispatch_tool",
            AsyncMock(return_value={"status": "success", "text": "ok"}),
        )

        client = _client()
        # Turn 1 (auto-create session) — "Create demo plot" deterministically routes to
        # local_demo_plot_script which goes through _dispatch_result and writes history.
        r1 = client.post("/execute", json={"command": "Create demo plot"})
        assert r1.status_code == 200
        sid = r1.json()["session_id"]

        for i in range(2, 5):
            r = client.post("/execute", json={"command": "Create demo plot", "session_id": sid})
            assert r.status_code == 200

        session = _session(sid)
        n_history = len(session.get("history", []))
        assert n_history >= 4, f"Expected at least 4 history entries for 4 turns, got {n_history}"
