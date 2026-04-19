import json
from pathlib import Path

from backend.artifact_resolver import build_alias_index, resolve_semantic_reference
from backend.main import _is_approval_command
from benchmarks.bio_benchmark_scorer import score_benchmark_run


def test_approval_command_accepts_prefixed_approve_phrase():
    assert _is_approval_command("Approve the correction and rerun.")


def test_artifact_resolver_builds_first_latest_deg_aliases():
    session = {
        "runs": [
            {"run_id": "run_a", "tool": "bulk_rnaseq_analysis"},
            {"run_id": "run_b", "tool": "bulk_rnaseq_analysis"},
        ],
        "artifacts": {},
    }
    idx = build_alias_index(session)
    assert idx["aliases"]["first_deg_results"]["run_id"] == "run_a"
    assert idx["aliases"]["latest_deg_results"]["run_id"] == "run_b"


def test_artifact_resolver_semantic_reference_resolution():
    session = {
        "runs": [{"run_id": "run_1", "tool": "bulk_rnaseq_analysis"}],
        "artifacts": {},
    }
    resolved = resolve_semantic_reference(session, "first DEG results")
    assert resolved["status"] == "resolved"
    assert resolved["target"]["run_id"] == "run_1"


def test_benchmark_scorer_emits_threshold_fields():
    spec = {
        "evaluation_dimensions": [
            {"id": "artifact_persistence", "severity": "high"},
        ],
        "scenarios": [
            {
                "id": "s1",
                "turns": [
                    {"id": "turn_01", "user": "Analyze data", "focus_dimensions": ["artifact_persistence"], "should_require_approval_before_execution": True},
                    {"id": "turn_02", "user": "Approve.", "focus_dimensions": ["artifact_persistence"]},
                ],
            }
        ],
    }
    run_data = {
        "session_id": "sid",
        "scenario_id": "s1",
        "rows": [
            {"turn_id": "turn_01", "tool": "__plan__", "status": "workflow_planned"},
            {"turn_id": "turn_02", "tool": "__plan__", "status": "success"},
        ],
    }
    scored = score_benchmark_run(run_data, spec)
    assert "raw_score" in scored and "max_score" in scored
    assert "overall_percentage" in scored and "threshold_met" in scored
    assert scored["raw_score"] <= scored["max_score"]
