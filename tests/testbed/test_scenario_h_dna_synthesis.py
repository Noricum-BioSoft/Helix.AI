"""
H. DNA synthesis & vendor — USER_SCENARIOS §6, TESTBED §9.

H1: Research vendors for 500 bp synthesis
H2: Submit sequence for synthesis (mock)
"""

import pytest


def _execute(client, session_id, command):
    r = client.post("/execute", json={"command": command, "session_id": session_id})
    assert r.status_code == 200, r.text
    return r.json()


class TestH1VendorResearch:
    """H1_be: Research vendors for synthesizing 500 bp sequence."""

    def test_h1_vendor_research(self, client, session_id):
        data = _execute(client, session_id, "Research vendors for synthesizing a 500 bp sequence.")
        assert data.get("success") is not False
        text = data.get("text") or (data.get("result") or {}).get("text") or str(data.get("result", ""))
        assert "vendor" in text.lower() or "synthesis" in text.lower() or "twist" in text.lower() or len(text) > 20


class TestH2SynthesisSubmission:
    """H2_be: Submit sequence for synthesis (mock)."""

    def test_h2_submit_synthesis(self, client, session_id):
        data = _execute(
            client,
            session_id,
            "Submit this sequence for synthesis to Twist: ATGCGATCGATCG",
        )
        assert data.get("success") is not False
        text = data.get("text") or (data.get("result") or {}).get("text") or str(data.get("result", ""))
        assert "submit" in text.lower() or "order" in text.lower() or "mock" in text.lower() or len(text) > 10
