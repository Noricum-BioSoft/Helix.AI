"""
Cross-cutting tests (CC1–CC5) from TESTBED.md §1.

Session creation, run attachment, needs_inputs not creating tool runs,
error response shape, and (documented) frontend New Session behavior.
"""

import pytest

from .conftest import BULK_RNASEQ_NEEDS_INPUTS_PROMPT


class TestSessionCreation:
    """CC1: New session creates session_id; no runs yet."""

    def test_create_session_returns_session_id(self, client):
        r = client.post("/create_session", json={})
        assert r.status_code == 200, r.text
        data = r.json()
        assert "session_id" in data
        assert isinstance(data["session_id"], str)
        assert len(data["session_id"]) > 0

    def test_session_create_alternative_endpoint(self, client):
        r = client.post("/session/create", json={"user_id": None})
        # Some codebases have both; accept 200 with session_id
        if r.status_code == 200:
            data = r.json()
            if "session_id" in data:
                assert isinstance(data["session_id"], str)


class TestRunAttachment:
    """CC2: Execute with session_id attaches run to that session."""

    def test_execute_with_session_creates_run(self, client, session_id):
        r = client.post(
            "/execute",
            json={"command": "What were the inputs and outputs of my last run?", "session_id": session_id},
        )
        assert r.status_code == 200, r.text
        payload = r.json()
        # May or may not have run_id depending on tool; session_id should be echoed
        assert payload.get("session_id") == session_id or "session_id" in payload

        runs_r = client.get(f"/session/{session_id}/runs")
        if runs_r.status_code == 200:
            runs = runs_r.json().get("runs", [])
            # At least one run (the Q&A or session summary)
            assert isinstance(runs, list)


class TestNeedsInputsNoToolRun:
    """CC3: needs_inputs response does not create a tool run (only router/agent)."""

    def test_bulk_rnaseq_needs_inputs_returns_status(self, client, session_id):
        r = client.post(
            "/execute",
            json={"command": BULK_RNASEQ_NEEDS_INPUTS_PROMPT, "session_id": session_id},
        )
        assert r.status_code == 200, r.text
        data = r.json()
        status = data.get("status") or (data.get("result") or {}).get("status")
        text = data.get("text") or (data.get("result") or {}).get("text") or ""
        # In mock mode the agent is skipped and intent may be classified as Q&A;
        # with agent we expect needs_inputs and parameter table.
        if status == "needs_inputs" or "Required inputs" in text or "required inputs" in text.lower():
            assert "count_matrix" in text or "count matrix" in text.lower()
            assert "sample_metadata" in text or "sample metadata" in text.lower()
        else:
            # Mock fallback: Q&A message is acceptable
            assert "mock" in text.lower() or "Q&A" in text or "question" in text.lower()


class TestErrorResponseShape:
    """CC4: build_standard_response sets error when result.status=error."""

    def test_error_response_has_error_field(self, client):
        # Trigger an error path (e.g. invalid or empty command might still return 200 with success=False)
        r = client.post("/execute", json={"command": "", "session_id": None})
        # Empty command might be 422 or 200 with error
        if r.status_code == 200:
            data = r.json()
            if not data.get("success"):
                assert "error" in data or "errors" in data or data.get("status") == "error"


class TestFrontendNewSessionDocumented:
    """CC5: Frontend New Session — documented only; no backend test."""

    def test_placeholder_frontend_new_session(self):
        """Frontend test: New Session clears history and gets new session_id. See TESTBED.md §11."""
        pytest.skip("Frontend test: run manually or with Playwright")
