"""
test_user_flows.py — Slice D end-to-end flow tests, one per user-visible archetype.

Each test is a multi-turn /execute sequence with assertions on the user-visible
response shape.  A regression in any archetype here means a user-visible bug.

Archetypes:
  1. Tabular profile   — upload + "profile the workbook" → summary text
  2. Advisory          — "what can I do with this dataset?" → result.advisory set, text non-empty
  3. Unknown intent    — "potato" → unknown_intent clarification card
  4. Meta / Q&A        — after a run, "what's the difference…" → advisory, NOT a new plan
  5. Tabular planning  — analytical command with uploads → tabular_analysis_plan staged for approval

All tests run in HELIX_MOCK_MODE=1 (no LLM calls).  The LLM router and intent
classifier are auto-patched by conftest.py autouse fixtures.
"""
from __future__ import annotations

import json
from typing import Any, Dict
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


def _execute(client: TestClient, command: str, session_id: str | None = None, **kwargs) -> Dict[str, Any]:
    payload: Dict[str, Any] = {"command": command}
    if session_id:
        payload["session_id"] = session_id
    payload.update(kwargs)
    r = client.post("/execute", json=payload)
    assert r.status_code == 200, f"HTTP {r.status_code}: {r.text[:200]}"
    return r.json()


# ── Archetype 1: Tabular profile ─────────────────────────────────────────────


class TestTabularProfile:
    """upload + 'profile the workbook' → summary in result.text (non-empty)."""

    def test_profile_returns_text(self, monkeypatch):
        """Profiling command routes to tabular_qa and returns text summary."""
        profile_result = {
            "status": "success",
            "text": "The workbook has 150 rows and 12 columns. Key columns: patient_id, treatment, response.",
            "tool_name": "tabular_qa",
            "rows": 150,
            "cols": 12,
        }
        monkeypatch.setattr(
            "backend.main.dispatch_tool",
            AsyncMock(return_value=profile_result),
        )

        # Patch router to return tabular_qa for this command
        from backend.command_router import CommandRouter
        with patch.object(CommandRouter, "_route_with_llm", lambda self, c, sc: ("tabular_qa", {})):
            client = _client()
            resp = _execute(client, "profile the workbook")

        assert resp.get("success") is True, f"Expected success, got: {resp}"
        assert resp.get("text", ""), "result.text must be non-empty for a profile response"
        assert resp.get("status") != "error", f"Expected non-error status, got: {resp.get('status')}"

    def test_profile_records_history(self, monkeypatch):
        """Profile turn must be recorded in session history."""
        from backend.command_router import CommandRouter
        from backend.history_manager import history_manager

        monkeypatch.setattr(
            "backend.main.dispatch_tool",
            AsyncMock(return_value={"status": "success", "text": "150 rows, 12 columns.", "tool_name": "tabular_qa"}),
        )
        with patch.object(CommandRouter, "_route_with_llm", lambda self, c, sc: ("tabular_qa", {})):
            client = _client()
            resp = _execute(client, "profile the workbook")

        sid = resp["session_id"]
        session = history_manager.get_session(sid)
        assert len(session.get("history", [])) >= 1, "Profile turn must write a history entry"


# ── Archetype 2: Advisory ─────────────────────────────────────────────────────


class TestAdvisory:
    """'what can I do with this dataset?' → result.advisory set, result.text non-empty."""

    def _make_advisory_result(self) -> Dict:
        """Return a canonical HelixAdvisory JSON string as the 'text' field."""
        advisory = {
            "helix_type": "advisory",
            "title": "Your Analysis Options",
            "summary": "You have a gene expression dataset. Here are 3 workflows you could run.",
            "sections": [],
            "workflow_steps": [
                {"step": 1, "name": "DESeq2 differential expression", "description": "Identify differentially expressed genes."},
            ],
            "requirements": [],
            "questions_for_user": [],
            "next_steps": ["Upload count matrix", "Run differential expression analysis"],
        }
        return {
            "status": "success",
            "text": json.dumps(advisory),
            "tool_name": "handle_natural_command",
        }

    def test_advisory_sets_structured_field(self, monkeypatch):
        """When result.text is an advisory JSON, build_standard_response must set result.advisory."""
        from backend.command_router import CommandRouter
        from backend.intent_classifier import IntentDecision

        monkeypatch.setattr(
            "backend.main.dispatch_tool",
            AsyncMock(return_value=self._make_advisory_result()),
        )
        # Force "execute" intent so the fallback router dispatches to dispatch_tool
        # instead of short-circuiting to qa_safe_fallback (which the conftest mock
        # would set to "qa" for question-word commands like "what can I do…").
        monkeypatch.setattr(
            "backend.intent_classifier.classify_intent",
            lambda t, **kw: IntentDecision(intent="execute", reason="test_advisory"),
        )
        with patch.object(CommandRouter, "_route_with_llm", lambda s, c, sc: ("handle_natural_command", {})):
            client = _client()
            resp = _execute(client, "what can I do with this dataset?")

        # Core Slice C assertion: result.advisory is populated by the backend.
        advisory = resp.get("advisory")
        assert advisory is not None, (
            "result.advisory must be set when backend returns an advisory response. "
            "Slice C adds this field so the frontend doesn't need to JSON-parse result.text."
        )
        assert advisory.get("helix_type") == "advisory", f"Expected helix_type='advisory', got: {advisory}"
        assert advisory.get("title"), "Advisory must have a non-empty title"

    def test_advisory_text_is_non_empty_prose(self, monkeypatch):
        """result.text must be non-empty human-readable prose (not raw JSON)."""
        from backend.command_router import CommandRouter

        monkeypatch.setattr(
            "backend.main.dispatch_tool",
            AsyncMock(return_value=self._make_advisory_result()),
        )
        with patch.object(CommandRouter, "_route_with_llm", lambda s, c, sc: ("handle_natural_command", {})):
            client = _client()
            resp = _execute(client, "what can I do with this dataset?")

        text = resp.get("text", "")
        assert text, "result.text must be non-empty for advisory responses"
        # After Slice C, text is a markdown rendering of the advisory, not raw JSON.
        # The full JSON object should not be the sole content of text.
        is_just_json = text.strip().startswith("{") and text.strip().endswith("}")
        assert not is_just_json, (
            "result.text should be a markdown rendering, not the raw JSON advisory object. "
            "Slice C extracts the advisory to result.advisory and generates prose for result.text."
        )


# ── Archetype 3: Unknown intent ───────────────────────────────────────────────


class TestUnknownIntent:
    """'potato' → unknown_intent clarification card — tool='unknown_intent', text non-empty."""

    def test_unknown_returns_clarification(self, monkeypatch):
        from backend.command_router import CommandRouter

        with patch.object(CommandRouter, "_route_with_llm", lambda s, c, sc: (
            "unknown_intent",
            {"router_reasoning": {"suggested_steps": ["Run differential expression", "Run sequence alignment"]}},
        )):
            client = _client()
            resp = _execute(client, "potato salad")

        assert resp.get("tool") == "unknown_intent", (
            f"Expected tool='unknown_intent' for nonsense input, got: {resp.get('tool')!r}"
        )
        assert resp.get("text"), "Unknown intent must return a non-empty clarification text"

    def test_unknown_has_success_true(self, monkeypatch):
        """unknown_intent is a successful user-facing response (success=True)."""
        from backend.command_router import CommandRouter

        with patch.object(CommandRouter, "_route_with_llm", lambda s, c, sc: ("unknown_intent", {})):
            client = _client()
            resp = _execute(client, "xyz nonsense abcde")

        assert resp.get("success") is True, (
            "unknown_intent should set success=True — it's a deliberate clarification, not an error"
        )


# ── Archetype 4: Meta question after run ─────────────────────────────────────


class TestMetaQuestionAfterRun:
    """After a run: 'what's the difference between your next step and the previous one?'
    should produce advisory text, NOT a new analysis plan.
    """

    def test_meta_question_is_not_a_plan(self, monkeypatch):
        from backend.command_router import CommandRouter

        advisory_result = {
            "status": "success",
            "text": "The previous run used DESeq2 with alpha=0.05. The next step would filter for LFC > 1.",
            "tool_name": "handle_natural_command",
        }
        monkeypatch.setattr(
            "backend.main.dispatch_tool",
            AsyncMock(return_value=advisory_result),
        )

        # Route meta question to handle_natural_command (not tabular_analysis)
        with patch.object(CommandRouter, "_route_with_llm", lambda s, c, sc: ("handle_natural_command", {})):
            client = _client()

            # First turn: a bio run
            run_result = {"status": "success", "text": "Bulk RNA-seq complete.", "tool_name": "bulk_rnaseq_analysis"}
            monkeypatch.setattr("backend.main.dispatch_tool", AsyncMock(return_value=run_result))
            r1 = _execute(client, "run bulk rnaseq analysis on demo data")
            sid = r1["session_id"]

            # Second turn: meta question
            monkeypatch.setattr("backend.main.dispatch_tool", AsyncMock(return_value=advisory_result))
            with patch.object(CommandRouter, "_route_with_llm", lambda s, c, sc: ("handle_natural_command", {})):
                r2 = _execute(client, "what's the difference between your next step and the previous analysis?", session_id=sid)

        # Core assertion: meta question must NOT return a plan
        assert r2.get("tool") != "tabular_analysis_plan", (
            "Meta/comparative question must not produce a tabular_analysis_plan. "
            "Slice B removes the is_analytical_request pre-router gate that caused this regression."
        )
        assert r2.get("tool") != "__plan__", "Meta question must not produce a plan IR"
        assert r2.get("text"), "Meta question must return non-empty advisory text"

    def test_meta_question_is_advisory_not_workflow(self, monkeypatch):
        """workflow_state must be IDLE (not WAITING_FOR_APPROVAL) for meta questions."""
        from backend.command_router import CommandRouter

        advisory_result = {
            "status": "success",
            "text": "The previous analysis compared treatment vs control. The next step adds batch correction.",
            "tool_name": "handle_natural_command",
        }
        monkeypatch.setattr(
            "backend.main.dispatch_tool",
            AsyncMock(return_value=advisory_result),
        )

        with patch.object(CommandRouter, "_route_with_llm", lambda s, c, sc: ("handle_natural_command", {})):
            client = _client()
            resp = _execute(client, "what is the difference between those two analyses?")

        assert resp.get("workflow_state") not in ("waiting_for_approval", "planning"), (
            f"Meta question must not enter approval/planning workflow state, got: {resp.get('workflow_state')!r}"
        )


# ── Archetype 5: Tabular planning (approval gate) ─────────────────────────────


class TestTabularPlanning:
    """Analytical command with tabular context → tabular_analysis_plan staged for approval."""

    def _session_context_with_tabular_file(self) -> Dict:
        return {
            "uploaded_files": [
                {
                    "filename": "gene_expr.csv",
                    "schema_preview": {"family": "tabular"},
                }
            ]
        }

    def test_tabular_analytical_command_generates_plan(self, monkeypatch):
        """When router returns tabular_analysis AND plan_analysis succeeds, a plan is staged."""
        from backend.command_router import CommandRouter

        mock_plan = {
            "type": "tabular_analysis",
            "goal": "Compare expression across treatment groups",
            "steps": [
                {"id": "step1", "description": "Load data", "tool_name": "tabular_analysis"},
                {"id": "step2", "description": "Compare groups", "tool_name": "tabular_analysis"},
            ],
        }

        with patch.object(CommandRouter, "_route_with_llm", lambda s, c, sc: ("tabular_analysis", {})), \
             patch(
                 "backend.tabular_qa.analysis_planner.plan_analysis",
                 return_value={"status": "success", "plan": mock_plan},
             ):
            client = _client()
            resp = _execute(client, "compare expression across treatment groups in gene_expr.csv")

        # The response should be a plan awaiting approval
        assert resp.get("tool") in ("tabular_analysis_plan", "__plan__"), (
            f"Expected a plan tool, got: {resp.get('tool')!r}. "
            "Tabular analytical commands should generate a plan for approval."
        )
        assert resp.get("approval_required") is True or resp.get("execute_ready") is True, (
            "Plan response must signal that approval/execution is available"
        )
        assert resp.get("text"), "Plan response must have non-empty text"

    def test_plan_gate_does_not_fire_for_advisory(self, monkeypatch):
        """'what can I do with this dataset?' must NOT generate a tabular plan."""
        from backend.command_router import CommandRouter

        plan_analysis_calls: list[int] = [0]

        def _counting_plan_analysis(cmd, ctx):
            plan_analysis_calls[0] += 1
            return {"status": "success", "plan": {}}

        advisory_result = {
            "status": "success",
            "text": "Here are your options.",
            "tool_name": "handle_natural_command",
        }
        monkeypatch.setattr("backend.main.dispatch_tool", AsyncMock(return_value=advisory_result))

        with patch.object(CommandRouter, "_route_with_llm", lambda s, c, sc: ("handle_natural_command", {})), \
             patch("backend.tabular_qa.analysis_planner.plan_analysis", side_effect=_counting_plan_analysis):
            client = _client()
            resp = _execute(client, "what can I do with this dataset?")

        assert plan_analysis_calls[0] == 0, (
            f"plan_analysis was called {plan_analysis_calls[0]} times for an advisory command. "
            "After Slice B, plan_analysis is only called when the router returns tabular_analysis."
        )
        assert resp.get("tool") != "tabular_analysis_plan", (
            "Advisory commands must not generate a tabular_analysis_plan"
        )


# ── Archetype 6: Advisory quick-reply follow-up ───────────────────────────────


class TestAdvisoryFollowup:
    """Regression suite for the advisory → quick-reply interaction pattern.

    A quick-reply button submits its label text (e.g. "yes, inspect it")
    as a new /execute command in a session that has NO pending workflow plan.
    The backend must:
      - NOT classify the reply as a workflow approval
      - NOT return "There is no pending workflow to approve"
      - Route the reply as a normal analytical command
    """

    # Phrases the frontend can submit as quick-reply labels.
    # All start with "yes" to exercise the affirmative-prefix guard.
    QUICK_REPLY_PHRASES = [
        "yes, inspect it",
        "yes, profile it",
        "yes, proceed",
        "yes, run the analysis",
        "focus on differential analysis",
        "focus on visualization",
    ]

    def _advisory_dispatch(self) -> dict:
        return {
            "status": "success",
            "text": "Workbook profiling complete.",
            "tool_name": "tabular_qa",
        }

    @pytest.mark.parametrize("phrase", QUICK_REPLY_PHRASES)
    def test_quick_reply_does_not_hit_approval_gate(self, monkeypatch, phrase):
        """Submitting a quick-reply label with no pending plan must NOT return
        the 'no pending workflow' error — it must return a normal response."""
        from backend.command_router import CommandRouter
        from backend.intent_classifier import IntentDecision

        monkeypatch.setattr(
            "backend.main.dispatch_tool",
            AsyncMock(return_value=self._advisory_dispatch()),
        )
        monkeypatch.setattr(
            "backend.intent_classifier.classify_intent",
            lambda t, **kw: IntentDecision(intent="execute", reason="test_quickreply"),
        )
        with patch.object(CommandRouter, "_route_with_llm", lambda s, c, sc: ("tabular_qa", {})):
            client = _client()
            resp = _execute(client, phrase)

        # The approval-gate dead-end message must never appear.
        text = (resp.get("text") or "").lower()
        assert "no pending workflow" not in text, (
            f"Quick-reply '{phrase}' was misrouted to the approval gate. "
            "is_approval_command must return False when has_pending_plan=False."
        )
        assert resp.get("success") is True, f"Expected success response, got: {resp}"

    def test_quick_reply_followup_in_same_session(self, monkeypatch):
        """Two-turn test: advisory turn followed immediately by a quick-reply
        in the same session.  The follow-up must NOT hit the approval gate."""
        from backend.command_router import CommandRouter
        from backend.intent_classifier import IntentDecision
        import json

        # Turn 1: advisory
        advisory_payload = {
            "helix_type": "advisory",
            "title": "Analysis Options",
            "summary": "Here is what I can do.",
            "sections": [],
            "workflow_steps": [],
            "requirements": [],
            "questions_for_user": [
                {
                    "question": "Should I profile the workbook?",
                    "examples": ["yes, inspect it", "yes, profile it"],
                }
            ],
            "next_steps": [],
        }
        advisory_result = {
            "status": "success",
            "text": json.dumps(advisory_payload),
            "tool_name": "handle_natural_command",
        }
        action_result = {
            "status": "success",
            "text": "Profiling complete: 5 sheets, 200 rows.",
            "tool_name": "tabular_qa",
        }

        monkeypatch.setattr(
            "backend.intent_classifier.classify_intent",
            lambda t, **kw: IntentDecision(intent="execute", reason="test"),
        )

        call_count = [0]

        def _dispatch(*args, **kwargs):
            call_count[0] += 1
            if call_count[0] == 1:
                return advisory_result
            return action_result

        monkeypatch.setattr("backend.main.dispatch_tool", AsyncMock(side_effect=_dispatch))

        with patch.object(CommandRouter, "_route_with_llm", lambda s, c, sc: ("handle_natural_command", {})):
            client = _client()
            r1 = _execute(client, "what can I do with this dataset?")

        sid = r1["session_id"]

        with patch.object(CommandRouter, "_route_with_llm", lambda s, c, sc: ("tabular_qa", {})):
            r2 = _execute(client, "yes, inspect it", session_id=sid)

        text2 = (r2.get("text") or "").lower()
        assert "no pending workflow" not in text2, (
            "Quick-reply in turn 2 was misrouted to the approval gate. "
            "Session state is IDLE so is_approval_command must return False."
        )
        assert r2.get("success") is True, f"Turn 2 must succeed, got: {r2}"
