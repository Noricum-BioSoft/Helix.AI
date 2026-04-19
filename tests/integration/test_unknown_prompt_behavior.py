from __future__ import annotations

import pytest
from fastapi.testclient import TestClient


@pytest.fixture(autouse=True)
def _isolate(tmp_path, monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")
    monkeypatch.setenv("HELIX_DEBUG_ROUTING", "1")
    monkeypatch.setenv("HELIX_DEMO_MODE", "0")
    monkeypatch.setenv("HELIX_SANDBOX_HOST_FALLBACK", "1")

    from backend.history_manager import history_manager

    history_manager.storage_dir = tmp_path / "sessions"
    history_manager.storage_dir.mkdir(parents=True, exist_ok=True)
    history_manager.sessions = {}
    history_manager._sessions_loaded = True


def _client() -> TestClient:
    from backend.main import app

    return TestClient(app)


def test_question_prompt_returns_safe_qa_message_without_tool_generation():
    r = _client().post("/execute", json={"command": "How do I run DESeq2 for factorial design?"})
    body = r.json()

    assert r.status_code == 200
    assert body.get("success") is True
    assert body.get("execution_path") == "qa_safe_fallback"
    assert "This looks like a question." in (body.get("text") or "")
    assert body.get("intent") == "qa"


def test_debug_routing_execution_path_is_exposed():
    r = _client().post("/execute", json={"command": "list tools"})
    body = r.json()

    assert r.status_code == 200
    assert body.get("execution_path") in {
        "phase2c_router",
        "fallback_router",
        "agent",
        "qa_safe_fallback",
    }
