"""Unit tests for the benchmark gate scorer."""

from __future__ import annotations

from benchmarks.scoring.use_case_gate_scorer import (
    GateThresholds,
    score_case_gates,
)


BASE_CASE = {
    "id": "test-case",
    "expected_route": "tabular_analysis",
    "forbidden_routes": ["sequence_alignment"],
    "forbidden_fallbacks": ["sequence_alignment"],
    "required_steps": ["ingest_tabular", "rank_and_select"],
    "required_outputs": ["table", "explanation"],
    "safety_requirements": ["explicit_assumptions"],
    "provenance_requirements": [
        "manifest_v1_present",
        "lineage_edges_present",
    ],
}


def _perfect_obs() -> dict:
    return {
        "observed_route": "tabular_analysis",
        "fallback_targets": [],
        "observed_steps": ["ingest_tabular", "rank_and_select"],
        "observed_output_types": ["table", "explanation"],
        "observed_safety_flags": ["explicit_assumptions"],
        "observed_provenance_signals": [
            "manifest_v1_present",
            "lineage_edges_present",
        ],
        "replay_status": "replay_match",
    }


def test_perfect_run_passes_gate() -> None:
    result = score_case_gates(BASE_CASE, _perfect_obs())
    assert result["gate_pass"] is True
    assert result["total_score"] == 1.0
    assert result["gate_fail_reasons"] == []


def test_forbidden_route_fails_gate() -> None:
    obs = _perfect_obs()
    obs["observed_route"] = "sequence_alignment"
    result = score_case_gates(BASE_CASE, obs)
    assert result["gate_pass"] is False
    assert "forbidden_or_wrong_route" in result["gate_fail_reasons"]


def test_forbidden_fallback_fails_gate() -> None:
    obs = _perfect_obs()
    obs["fallback_targets"] = ["sequence_alignment"]
    result = score_case_gates(BASE_CASE, obs)
    assert result["gate_pass"] is False
    assert "forbidden_fallback" in result["gate_fail_reasons"]


def test_replay_failure_fails_gate() -> None:
    obs = _perfect_obs()
    obs["replay_status"] = "replay_mismatch"
    result = score_case_gates(BASE_CASE, obs)
    assert result["gate_pass"] is False
    assert "replay_failed" in result["gate_fail_reasons"]


def test_incomplete_provenance_fails_gate() -> None:
    obs = _perfect_obs()
    obs["observed_provenance_signals"] = ["manifest_v1_present"]
    result = score_case_gates(BASE_CASE, obs)
    assert result["gate_pass"] is False
    assert "provenance_incomplete" in result["gate_fail_reasons"]


def test_min_score_threshold_gate() -> None:
    obs = _perfect_obs()
    obs["observed_steps"] = []
    obs["observed_output_types"] = []
    obs["observed_safety_flags"] = []
    result = score_case_gates(
        BASE_CASE,
        obs,
        thresholds=GateThresholds(min_total_score=0.95),
    )
    assert result["gate_pass"] is False
    assert "below_min_total_score" in result["gate_fail_reasons"]
