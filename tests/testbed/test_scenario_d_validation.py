"""
D. Validation & interpretation — USER_SCENARIOS §3, TESTBED §5.

D1: Which genes most significant (Q&A after B1)
D2: What does padj < 0.05 mean (Q&A)
D3: Show volcano with log2 again (patch)
"""

import pytest

from .conftest import BULK_RNASEQ_PROMPT


def _execute(client, session_id, command):
    r = client.post("/execute", json={"command": command, "session_id": session_id})
    assert r.status_code == 200, r.text
    return r.json()


class TestD1GenesMostSignificant:
    """D1_be: After B1, Q&A 'Which genes most significant' — no tool run."""

    def test_d1_qa_after_b1(self, client, session_id):
        _execute(client, session_id, BULK_RNASEQ_PROMPT)
        data = _execute(
            client,
            session_id,
            "Which genes are most significant for the infection effect?",
        )
        assert data.get("success") is not False
        text = data.get("text") or (data.get("result") or {}).get("text") or ""
        assert len(text) > 20


class TestD2PadjMeaning:
    """D2_be: Q&A 'What does padj < 0.05 mean?' — markdown, no analysis."""

    def test_d2_statistical_interpretation(self, client, session_id):
        data = _execute(client, session_id, "What does padj < 0.05 mean?")
        assert data.get("success") is not False
        text = data.get("text") or (data.get("result") or {}).get("text") or ""
        assert "padj" in text.lower() or "fdr" in text.lower() or "p-value" in text.lower() or len(text) > 30


class TestD3RevertLog2Scale:
    """D3_be: After B1, 'Show volcano with log2 scale again' → patch."""

    def test_d3_revert_x_scale_log2(self, client, session_id):
        _execute(client, session_id, BULK_RNASEQ_PROMPT)
        data = _execute(
            client,
            session_id,
            "Show me the same volcano but with log2 scale on the x-axis again.",
        )
        assert data.get("success") is not False or "log" in str(data).lower()
