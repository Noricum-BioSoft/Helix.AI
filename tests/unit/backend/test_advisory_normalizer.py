"""Unit tests for backend.advisory_normalizer.

Verifies that normalize_advisory_text() converts all observed LLM output
shapes into the canonical HelixAdvisory JSON schema and that the result
always deserialises to a valid HelixAdvisory Pydantic model.
"""
from __future__ import annotations

import json

import pytest

from backend.advisory_normalizer import HelixAdvisory, normalize_advisory_text


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _parse(text: str) -> dict:
    return json.loads(text)


def _valid_advisory(d: dict) -> None:
    """Assert the dict is a valid HelixAdvisory and has required keys."""
    advisory = HelixAdvisory.model_validate(d)
    assert advisory.helix_type == "advisory"
    assert advisory.title
    assert advisory.summary


# ---------------------------------------------------------------------------
# Already-canonical input: must pass through unchanged
# ---------------------------------------------------------------------------

class TestPassthrough:
    def test_canonical_passes_through(self):
        canonical = json.dumps({
            "helix_type": "advisory",
            "title": "ChIP-seq Peak Calling",
            "summary": "Run MACS3 to identify H3K27ac binding sites.",
            "sections": [{"heading": "Overview", "content": "Use MACS3."}],
            "workflow_steps": [],
            "requirements": [],
            "questions_for_user": [],
            "next_steps": ["Upload your BAM files."],
        })
        result = normalize_advisory_text(canonical)
        d = _parse(result)
        assert d["helix_type"] == "advisory"
        assert d["title"] == "ChIP-seq Peak Calling"
        _valid_advisory(d)

    def test_fenced_canonical_strips_fences(self):
        canonical = json.dumps({
            "helix_type": "advisory",
            "title": "Test",
            "summary": "A summary.",
        })
        fenced = f"```json\n{canonical}\n```"
        result = normalize_advisory_text(fenced)
        d = _parse(result)
        assert d["helix_type"] == "advisory"
        _valid_advisory(d)


# ---------------------------------------------------------------------------
# Shape 1: legacy planning advisory  {classification, requirements, next_steps}
# ---------------------------------------------------------------------------

class TestShape1LegacyPlanning:
    def test_normalises_shape1(self):
        shape1 = json.dumps({
            "type": "planning_advisory",
            "title": "RNA-seq workflow",
            "summary": "Differential expression pipeline.",
            "classification": {
                "workflow_type": "bulk_rnaseq",
                "confidence": "high",
            },
            "requirements": [
                {"item": "FASTQ files", "details": "Raw reads"},
                {"item": "Reference genome"},
            ],
            "next_steps": ["Upload FASTQ", "Confirm parameters"],
        })
        result = normalize_advisory_text(shape1)
        d = _parse(result)
        _valid_advisory(d)
        assert len(d["requirements"]) == 2
        assert d["next_steps"] == ["Upload FASTQ", "Confirm parameters"]
        assert d["classification"] is not None


# ---------------------------------------------------------------------------
# Shape 2: legacy explanation  {type: "answer", sections}
# ---------------------------------------------------------------------------

class TestShape2LegacyExplanation:
    def test_normalises_shape2(self):
        shape2 = json.dumps({
            "type": "answer",
            "title": "What is DESeq2?",
            "summary": "DESeq2 is an R package for differential expression.",
            "sections": [
                {"heading": "Method", "content": "Uses NB distribution."},
                {"heading": "Output", "content": "Fold-change and p-values."},
            ],
            "next_steps": ["Install R", "Load counts matrix"],
        })
        result = normalize_advisory_text(shape2)
        d = _parse(result)
        _valid_advisory(d)
        assert len(d["sections"]) == 2
        assert d["sections"][0]["heading"] == "Method"


# ---------------------------------------------------------------------------
# Shape 3: gpt-4.1 nested  {plan: {...}, answer: {sections, next_steps}}
# ---------------------------------------------------------------------------

class TestShape3Gpt41Nested:
    def test_normalises_shape3(self):
        shape3 = json.dumps({
            "plan": {
                "title": "ChIP-seq Pipeline",
                "workflow_type": "chip_seq",
                "summary": "Peak calling and motif enrichment.",
            },
            "answer": {
                "sections": [
                    {"heading": "Peak calling", "content": "Use MACS3."},
                ],
                "next_steps": ["Upload BAM files."],
            },
        })
        result = normalize_advisory_text(shape3)
        d = _parse(result)
        _valid_advisory(d)
        assert "ChIP-seq" in d["title"] or "ChIP" in d["summary"]
        assert len(d["sections"]) >= 1
        assert d["next_steps"]

    def test_shape3_string_classification(self):
        """classification field may be a plain string, not a dict."""
        shape3 = json.dumps({
            "plan": {
                "title": "ATAC-seq",
                "workflow_type": "atac_seq",
                "summary": "Open chromatin.",
                "classification": "atac_seq",
            },
            "answer": {
                "sections": [],
                "next_steps": ["Upload files"],
            },
        })
        result = normalize_advisory_text(shape3)
        d = _parse(result)
        _valid_advisory(d)


# ---------------------------------------------------------------------------
# Non-advisory plain text: must be returned unchanged
# ---------------------------------------------------------------------------

class TestNonAdvisoryPassthrough:
    def test_plain_text_unchanged(self):
        text = "The analysis completed successfully."
        result = normalize_advisory_text(text)
        assert result == text

    def test_empty_string_unchanged(self):
        assert normalize_advisory_text("") == ""

    def test_non_advisory_json_unchanged(self):
        non_advisory = json.dumps({"status": "success", "rows": 42})
        result = normalize_advisory_text(non_advisory)
        assert result == non_advisory
