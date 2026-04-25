from __future__ import annotations

from unittest.mock import AsyncMock, MagicMock

import pytest
from fastapi.testclient import TestClient


@pytest.fixture(autouse=True)
def _isolate(tmp_path, monkeypatch):
    monkeypatch.setenv("HELIX_DEBUG_ROUTING", "1")
    monkeypatch.setenv("HELIX_DEMO_MODE", "0")
    monkeypatch.setenv("HELIX_SANDBOX_HOST_FALLBACK", "1")
    # Keyword routing for tests that stub the router; production uses LLM-first.
    monkeypatch.setenv("HELIX_LLM_ROUTER_FIRST", "0")

    from backend.history_manager import history_manager

    history_manager.storage_dir = tmp_path / "sessions"
    history_manager.storage_dir.mkdir(parents=True, exist_ok=True)
    history_manager.sessions = {}
    history_manager._sessions_loaded = True


def _client() -> TestClient:
    from backend.main import app

    return TestClient(app)


def test_phase2c_allowlist_short_circuits_agent(monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "0")

    mock_router = MagicMock()
    mock_router.route_command.return_value = ("toolbox_inventory", {})
    mock_router.route_command_with_shadow.return_value = ("toolbox_inventory", {})
    monkeypatch.setattr("backend.command_router.CommandRouter", lambda: mock_router)

    mock_dispatch = AsyncMock(return_value={"status": "success", "text": "ok"})
    monkeypatch.setattr("backend.main.dispatch_tool", mock_dispatch)

    mock_agent = AsyncMock(side_effect=AssertionError("Agent should not be called for allowlisted fast-path"))
    monkeypatch.setattr("backend.agent.handle_command", mock_agent)

    r = _client().post("/execute", json={"command": "list tools"})
    body = r.json()

    assert r.status_code == 200
    assert body.get("tool") == "toolbox_inventory"
    assert body.get("execution_path") == "phase2c_router"
    mock_dispatch.assert_awaited()


def test_phase2c_visualize_job_short_circuits_agent(monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "0")

    mock_router = MagicMock()
    _vis_ret = ("visualize_job_results", {"job_id": "99228f84-dea7-4424-8569-3f9c235a1547"})
    mock_router.route_command.return_value = _vis_ret
    mock_router.route_command_with_shadow.return_value = _vis_ret
    monkeypatch.setattr("backend.command_router.CommandRouter", lambda: mock_router)

    mock_dispatch = AsyncMock(
        return_value={"status": "success", "visualization_type": "results_viewer", "text": "ok"}
    )
    monkeypatch.setattr("backend.main.dispatch_tool", mock_dispatch)

    mock_agent = AsyncMock(side_effect=AssertionError("Agent should not be called for visualize_job_results fast-path"))
    monkeypatch.setattr("backend.agent.handle_command", mock_agent)

    r = _client().post("/execute", json={"command": "visualize results for that job"})
    body = r.json()

    assert r.status_code == 200
    assert body.get("tool") == "visualize_job_results"
    assert body.get("execution_path") == "phase2c_router"
    mock_dispatch.assert_awaited()


def test_non_allowlisted_router_result_uses_agent_path(monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "0")

    mock_router = MagicMock()
    mock_router.route_command.return_value = ("read_trimming", {"reads": "reads.fq"})
    monkeypatch.setattr("backend.command_router.CommandRouter", lambda: mock_router)

    mock_agent = AsyncMock(return_value={"status": "success", "text": "agent-handled"})
    monkeypatch.setattr("backend.agent.handle_command", mock_agent)

    r = _client().post("/execute", json={"command": "trim these reads"})
    body = r.json()

    assert r.status_code == 200
    assert body.get("tool") == "agent"
    assert body.get("execution_path") == "agent"
    mock_agent.assert_awaited()


def test_mock_mode_falls_back_to_router(monkeypatch):
    """In HELIX_MOCK_MODE=1 (agent disabled), the fallback router path runs.

    Keyword routing is enabled (HELIX_LLM_ROUTER_FIRST=0) so the mocked router
    can return a result without a live LLM.  Production uses LLM-first routing.
    """
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")

    mock_router = MagicMock()
    mock_router.route_command.return_value = ("read_trimming", {"reads": "reads.fq"})
    mock_router.route_command_with_shadow.return_value = ("read_trimming", {"reads": "reads.fq"})
    monkeypatch.setattr("backend.command_router.CommandRouter", lambda: mock_router)

    mock_broker = MagicMock()
    mock_broker.execute_tool = AsyncMock(return_value={"status": "success", "text": "fallback"})
    monkeypatch.setattr("backend.main._get_execution_broker", lambda: mock_broker)

    r = _client().post("/execute", json={"command": "trim these reads"})
    body = r.json()

    assert r.status_code == 200
    assert body.get("tool") == "read_trimming"
    assert body.get("execution_path") == "fallback_router"
    mock_broker.execute_tool.assert_awaited()
