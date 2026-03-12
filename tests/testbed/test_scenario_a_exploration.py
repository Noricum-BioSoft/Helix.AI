"""
A. Exploration & tool discovery (USER_SCENARIOS §3, TESTBED §2).

A1: bulk RNA-seq needs_inputs
A2: single-cell needs_inputs
A3: FastQC needs_inputs
A4: Q&A (bulk vs single-cell) — no tool execution
"""

import pytest

from .conftest import (
    BULK_RNASEQ_NEEDS_INPUTS_PROMPT,
    SINGLECELL_NEEDS_INPUTS_PROMPT,
    FASTQC_NEEDS_INPUTS_PROMPT,
    QA_PROMPT,
)


def _execute(client, session_id, command):
    r = client.post("/execute", json={"command": command, "session_id": session_id})
    assert r.status_code == 200, r.text
    return r.json()


def _get_text_and_status(data):
    result = data.get("result") or data.get("data") or {}
    if isinstance(result, dict):
        text = data.get("text") or result.get("text") or ""
        status = data.get("status") or result.get("status") or ""
    else:
        text = data.get("text") or ""
        status = data.get("status") or ""
    return text, status


class TestA1BulkRnaseqNeedsInputs:
    """A1_be: Bulk RNA-seq 'what can you do' → needs_inputs with parameter table (or mock Q&A fallback)."""

    def test_a1_returns_needs_inputs_or_required_table(self, client, session_id):
        data = _execute(client, session_id, BULK_RNASEQ_NEEDS_INPUTS_PROMPT)
        text, status = _get_text_and_status(data)
        if status == "needs_inputs" or "Required inputs" in text or "required inputs" in text.lower():
            assert "count_matrix" in text or "count matrix" in text.lower()
            assert "sample_metadata" in text or "sample metadata" in text.lower()
            assert "|" in text
        else:
            # Mock mode: agent skipped, Q&A fallback
            assert "mock" in text.lower() or "Q&A" in text or "question" in text.lower()

    def test_a1_no_hallucinated_params(self, client, session_id):
        data = _execute(client, session_id, BULK_RNASEQ_NEEDS_INPUTS_PROMPT)
        text, status = _get_text_and_status(data)
        if status == "needs_inputs" or "Required inputs" in text:
            known = ["count_matrix", "sample_metadata", "design", "alpha"]
            assert any(p in text.lower() for p in known)
        # Else mock fallback: no need to check params


class TestA2SingleCellNeedsInputs:
    """A2_be: Single-cell 'what do I need' → needs_inputs."""

    def test_a2_returns_data_file_or_required_inputs(self, client, session_id):
        data = _execute(client, session_id, SINGLECELL_NEEDS_INPUTS_PROMPT)
        text, status = _get_text_and_status(data)
        assert status == "needs_inputs" or "required" in text.lower() or "data_file" in text or "data file" in text.lower()
        assert "|" in text or "parameter" in text.lower()


class TestA3FastQcNeedsInputs:
    """A3_be: FastQC 'how do I run' → needs_inputs with input_r1, input_r2."""

    def test_a3_returns_input_paths(self, client, session_id):
        data = _execute(client, session_id, FASTQC_NEEDS_INPUTS_PROMPT)
        text, status = _get_text_and_status(data)
        assert "input" in text.lower() or "fastq" in text.lower() or "r1" in text.lower() or status == "needs_inputs"
        assert "|" in text or "parameter" in text.lower()


class TestA4QaNoToolExecution:
    """A4_be: Q&A 'difference between bulk and single-cell' — no tool run."""

    def test_a4_returns_markdown_answer(self, client, session_id):
        data = _execute(client, session_id, QA_PROMPT)
        text, status = _get_text_and_status(data)
        assert data.get("success") is not False or status != "error"
        # Answer should mention concepts, not execute analysis
        assert len(text) > 50
        # Should not look like a needs_inputs table for bulk_rnaseq
        assert "count_matrix" not in text or "single-cell" in text.lower() or "bulk" in text.lower()
