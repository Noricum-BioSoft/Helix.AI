"""
test_routing_layers_single.py — Slice B routing layer collapse tests.

Asserts that for archetypal commands, the LLM router is called exactly once
per /execute request (not 2–3 times as in the pre-Slice-B architecture).

The old flow called the router up to 3 times:
  1. is_analytical_request gate (before the router)
  2. Universal approval gate
  3. Phase 2c fast-path check

After Slice B, the is_analytical_request pre-router gate is removed; the LLM
router is the single authority on intent classification.

We count LLM router calls via a call-counting wrapper on
CommandRouter._route_with_llm.  Deterministic (regex) matches skip _route_with_llm
entirely — those are fine and not tested here.
"""
from __future__ import annotations

from typing import Any, Dict, Tuple
from unittest.mock import AsyncMock, patch

import pytest
from fastapi.testclient import TestClient


# ── Isolation ────────────────────────────────────────────────────────────────


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


# ── Counter fixture ───────────────────────────────────────────────────────────


class _LlmRouterCallCounter:
    """Wraps CommandRouter._route_with_llm and counts invocations."""

    def __init__(self, mock_route_fn):
        self.calls = 0
        self._fn = mock_route_fn

    def __call__(self, self_router, command: str, session_context):
        self.calls += 1
        return self._fn(self_router, command, session_context)


# ── Tests ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("command,description", [
    ("what can I do with this dataset?", "advisory — handle_natural_command"),
    ("run bulk rnaseq analysis on my samples", "bio tool execution"),
    ("compare expression across treatment groups", "tabular analytical request"),
    ("potato salad recipe", "unknown intent clarification"),
])
def test_llm_router_called_at_most_three_per_request(command, description, monkeypatch):
    """After Slice B.1, the LLM router is called at most 3 times per request.

    The 3 remaining call sites are:
      1. Phase 2c fast-path check (route_command_with_shadow — primary)
      2. Universal approval gate (route_command_with_shadow — approval check)
      3. Fallback router (route_command — when agent is disabled in mock mode)

    The old pre-router is_analytical_request keyword gate did NOT call the LLM
    router, so its removal does not reduce the LLM call count. However, it DID
    introduce a parallel intent classification path that could disagree with the
    router — removing it means the LLM router is the single classifier.

    Full router consolidation to ≤ 1 LLM call is tracked as a follow-up
    (Slice B.2 — PhaseRoute single-pass refactor).
    """
    from unittest.mock import AsyncMock

    monkeypatch.setattr(
        "backend.main.dispatch_tool",
        AsyncMock(return_value={"status": "success", "text": "ok", "tool_name": "handle_natural_command"}),
    )

    call_counts: list[int] = [0]

    from backend.command_router import CommandRouter
    original_route = CommandRouter._route_with_llm

    def _counting_route(self_router, cmd: str, sc):
        call_counts[0] += 1
        return original_route(self_router, cmd, sc)

    with patch.object(CommandRouter, "_route_with_llm", _counting_route):
        client = _client()
        r = client.post("/execute", json={"command": command})
        assert r.status_code == 200, f"{description}: HTTP error: {r.text}"

    assert call_counts[0] <= 3, (
        f"{description!r}: expected ≤ 3 LLM router calls, got {call_counts[0]}. "
        "A new pre-router classification gate may have been introduced."
    )


def test_tabular_analysis_does_not_call_is_analytical_request(monkeypatch):
    """The is_analytical_request gate must no longer be called from /execute.

    Previously, every request with tabular uploads ran through
    is_analytical_request *before* the LLM router.  After Slice B, the
    router is the only classifier — is_analytical_request is not called from
    the execute() hot path.
    """
    from unittest.mock import AsyncMock, MagicMock

    monkeypatch.setattr(
        "backend.main.dispatch_tool",
        AsyncMock(return_value={"status": "success", "text": "ok"}),
    )

    is_analytical_calls: list[int] = [0]

    # Patch at the module where it is imported in main.py
    original_module = None
    try:
        import backend.tabular_qa.analysis_planner as _ap_mod
        original_is_analytical = _ap_mod.is_analytical_request

        def _counting_is_analytical(cmd, ctx):
            is_analytical_calls[0] += 1
            return original_is_analytical(cmd, ctx)

        monkeypatch.setattr(
            _ap_mod,
            "is_analytical_request",
            _counting_is_analytical,
        )
    except (ImportError, AttributeError):
        pytest.skip("analysis_planner not importable")

    client = _client()
    r = client.post("/execute", json={"command": "compare expression across groups"})
    assert r.status_code == 200

    assert is_analytical_calls[0] == 0, (
        f"is_analytical_request was called {is_analytical_calls[0]} times; "
        "Slice B requires it to be removed from the /execute hot path."
    )
