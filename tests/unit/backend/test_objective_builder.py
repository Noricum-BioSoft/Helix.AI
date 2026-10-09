"""Phase 1.1 — objective builder: deterministic in mock mode, IDs minted once, inputs from session."""

from __future__ import annotations

import pytest

from backend.contracts.ids import is_trace_id
from backend.orchestration import objective_builder
from backend.orchestration.objective_builder import build_objective, build_objective_deterministic


@pytest.fixture(autouse=True)
def _mock_mode(monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")
    monkeypatch.delenv("HELIX_OBJECTIVE_LLM", raising=False)
    monkeypatch.delenv("HELIX_EXECUTION_PROFILE", raising=False)


def test_objective_mints_trace_and_objective_ids():
    a = build_objective("Run FastQC on my reads", {"session_id": "s1"})
    b = build_objective("Run FastQC on my reads", {"session_id": "s1"})
    assert is_trace_id(a.trace_id) and is_trace_id(b.trace_id)
    assert a.trace_id != b.trace_id
    assert a.objective_id != b.objective_id
    assert a.status == "active"
    assert a.user_context["session_id"] == "s1"


def test_objective_reuses_given_trace_id():
    obj = build_objective("x", trace_id="orch_" + "0" * 32)
    assert obj.trace_id == "orch_" + "0" * 32


def test_objective_known_inputs_from_uploads_and_artifacts():
    ctx = {
        "uploaded_files": [
            {"filename": "reads.fastq.gz", "stored_path": "/tmp/reads.fastq.gz", "sha256": "abc"},
            "s3://bucket/counts.csv",
        ],
        "artifacts": {"art_1": {"uri": "/tmp/plot.png", "type": "plot"}},
    }
    obj = build_objective_deterministic("analyze", ctx)
    uris = {r.uri for r in obj.known_inputs}
    assert uris == {"/tmp/reads.fastq.gz", "s3://bucket/counts.csv", "/tmp/plot.png"}
    by_uri = {r.uri: r for r in obj.known_inputs}
    assert by_uri["/tmp/reads.fastq.gz"].content_hash == "abc"
    assert by_uri["/tmp/plot.png"].artifact_id == "art_1"


def test_objective_execution_constraints_follow_profile():
    obj = build_objective_deterministic("analyze")
    execution = obj.constraints.execution
    assert execution is not None
    assert execution.preferred_provider == "local_compute"
    assert "local_compute" in execution.allowed_providers
    # vendor neutrality: the default profile enables no cloud provider
    assert not any(p.startswith("aws") for p in execution.allowed_providers)


def test_objective_records_routed_tool_not_prose():
    obj = build_objective("Run FastQC", {}, {"tool": "fastqc_quality_analysis", "intent": "execute"})
    assert obj.user_context == {"routed_tool": "fastqc_quality_analysis", "intent": "execute"}


def test_llm_enrichment_disabled_in_mock_mode(monkeypatch):
    monkeypatch.setenv("HELIX_OBJECTIVE_LLM", "1")
    called = {"n": 0}

    def _boom(*_a, **_k):
        called["n"] += 1
        raise AssertionError("LLM must not be called in mock mode")

    monkeypatch.setattr(objective_builder, "_enrich_with_llm", _boom)
    build_objective("x")
    assert called["n"] == 0


def test_llm_enrichment_failure_falls_back(monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "0")
    monkeypatch.setenv("HELIX_OBJECTIVE_LLM", "1")

    def _boom(*_a, **_k):
        raise RuntimeError("no key")

    monkeypatch.setattr(objective_builder, "_enrich_with_llm", _boom)
    obj = build_objective("Compare expression between conditions")
    assert obj.question == "Compare expression between conditions"


def test_llm_enrichment_only_updates_allowed_fields(monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "0")
    monkeypatch.setenv("HELIX_OBJECTIVE_LLM", "1")

    class _Resp:
        content = (
            '{"question": "Which genes are DE?", "biological_system": "human PBMC", '
            '"desired_evidence": ["DE table"], "expected_outputs": ["volcano"], "success_criteria": ["padj<0.05"],'
            ' "objective_id": "hacked", "known_inputs": ["x"]}'
        )

    class _LLM:
        def invoke(self, _messages):
            return _Resp()

    monkeypatch.setattr("backend.orchestration.approval_classifier._get_llm", lambda: _LLM())
    obj = build_objective("Compare expression", {"uploaded_files": ["s3://b/counts.csv"]})
    assert obj.question == "Which genes are DE?"
    assert obj.biological_system == "human PBMC"
    assert obj.success_criteria == ["padj<0.05"]
    assert obj.objective_id != "hacked"
    assert [r.uri for r in obj.known_inputs] == ["s3://b/counts.csv"]
