"""Unit gate tests for the three P0 tabular benchmark cases.

These tests run entirely without a live backend – they verify that the
gate scorer produces the correct PASS/FAIL decisions for synthetic
observations that reflect the tabular analysis MVP.

A test name prefix of ``test_pass_`` means the observations represent
ideal tabular-analysis execution and MUST score gate_pass=True.
A prefix of ``test_fail_`` means the observations represent a known-bad
execution (e.g. wrong route) and MUST score gate_pass=False.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import pytest
import yaml

from benchmarks.scoring.use_case_gate_scorer import GateThresholds, score_case_gates

REPO_ROOT = Path(__file__).resolve().parents[3]
CASES_DIR = REPO_ROOT / "benchmarks" / "cases"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load_case(relative_path: str) -> Dict[str, Any]:
    return yaml.safe_load((CASES_DIR / relative_path).read_text())


def _perfect_tabular_obs(steps: list[str], outputs: list[str] | None = None) -> Dict[str, Any]:
    return {
        "observed_route": "tabular_analysis",
        "fallback_targets": [],
        "observed_steps": steps,
        "observed_output_types": outputs or ["plot", "explanation"],
        "observed_safety_flags": ["uncertainty_note", "explicit_assumptions"],
        "observed_provenance_signals": [
            "manifest_v1_present",
            "lineage_edges_present",
            "fallback_events_present",
            "route_decision_trace_present",
        ],
        "replay_status": "replay_match",
    }


def _seq_alignment_obs() -> Dict[str, Any]:
    """Simulates a mis-routed run that fell back to sequence_alignment."""
    return {
        "observed_route": "sequence_alignment",
        "fallback_targets": ["sequence_alignment"],
        "observed_steps": [],
        "observed_output_types": [],
        "observed_safety_flags": [],
        "observed_provenance_signals": [],
        "replay_status": "no_replay",
    }


THRESHOLDS = GateThresholds(min_total_score=0.80)


# ---------------------------------------------------------------------------
# case-tabular-patient-cohort-correlation
# ---------------------------------------------------------------------------

CORRELATION_CASE = _load_case("tabular/tabular-patient-cohort-correlation.yaml")


def test_pass_correlation_case_perfect_run() -> None:
    obs = _perfect_tabular_obs(
        steps=["ingest_tabular", "compute_correlation", "generate_heatmap"],
        outputs=["plot", "explanation"],
    )
    result = score_case_gates(CORRELATION_CASE, obs, thresholds=THRESHOLDS)
    assert result["gate_pass"] is True, result["gate_fail_reasons"]
    assert result["total_score"] >= 0.80


def test_fail_correlation_case_wrong_route() -> None:
    obs = _seq_alignment_obs()
    result = score_case_gates(CORRELATION_CASE, obs, thresholds=THRESHOLDS)
    assert result["gate_pass"] is False
    # Scorer encodes this as the generic token "forbidden_or_wrong_route"
    assert "forbidden_or_wrong_route" in result["gate_fail_reasons"]


def test_fail_correlation_case_missing_required_steps() -> None:
    """Missing steps reduce score but don't trip a hard gate by themselves.
    The run fails *all checks* (passed_all_checks=False) even though gate_pass
    may still be True if the total weighted score stays above the threshold."""
    obs = _perfect_tabular_obs(
        steps=["ingest_tabular"],  # missing compute_correlation and generate_heatmap
        outputs=["plot", "explanation"],
    )
    result = score_case_gates(CORRELATION_CASE, obs, thresholds=THRESHOLDS)
    assert result["passed_all_checks"] is False
    assert result["checks"]["required_steps"]["ok"] is False


def test_fail_correlation_case_missing_plot_output() -> None:
    """Missing outputs reduce score but are a soft check; gate_pass depends on
    total weighted score.  We verify the specific check fails."""
    obs = _perfect_tabular_obs(
        steps=["ingest_tabular", "compute_correlation", "generate_heatmap"],
        outputs=["explanation"],  # plot missing
    )
    result = score_case_gates(CORRELATION_CASE, obs, thresholds=THRESHOLDS)
    assert result["checks"]["required_outputs"]["ok"] is False


# ---------------------------------------------------------------------------
# case-tabular-clinical-trial-group-comparison
# ---------------------------------------------------------------------------

GROUP_COMP_CASE = _load_case("tabular/tabular-clinical-trial-group-comparison.yaml")


def test_pass_group_comparison_perfect_run() -> None:
    obs = _perfect_tabular_obs(
        steps=["ingest_tabular", "group_statistics", "generate_plot"],
        outputs=["plot", "explanation"],
    )
    result = score_case_gates(GROUP_COMP_CASE, obs, thresholds=THRESHOLDS)
    assert result["gate_pass"] is True, result["gate_fail_reasons"]
    assert result["total_score"] >= 0.80


def test_fail_group_comparison_forbidden_fallback() -> None:
    obs = _perfect_tabular_obs(
        steps=["ingest_tabular", "group_statistics", "generate_plot"],
        outputs=["plot", "explanation"],
    )
    obs["fallback_targets"] = ["sequence_alignment"]
    result = score_case_gates(GROUP_COMP_CASE, obs, thresholds=THRESHOLDS)
    assert result["gate_pass"] is False
    # Scorer encodes this as the generic token "forbidden_fallback"
    assert "forbidden_fallback" in result["gate_fail_reasons"]


def test_fail_group_comparison_no_safety_flags() -> None:
    """No safety flags is a soft check; it fails the check but not necessarily
    the hard gate.  We verify the safety check failed."""
    obs = _perfect_tabular_obs(
        steps=["ingest_tabular", "group_statistics", "generate_plot"],
        outputs=["plot", "explanation"],
    )
    obs["observed_safety_flags"] = []  # stripped of all safety notes
    result = score_case_gates(GROUP_COMP_CASE, obs, thresholds=THRESHOLDS)
    assert result["checks"]["safety"]["ok"] is False


# ---------------------------------------------------------------------------
# case-routing-safety-proteomics-no-sequence-fallback
# ---------------------------------------------------------------------------

PROTEOMICS_CASE = _load_case("routing_safety/proteomics-no-sequence-fallback.yaml")


def test_pass_proteomics_routing_safety_perfect_run() -> None:
    obs = _perfect_tabular_obs(
        steps=["ingest_tabular", "run_pca", "generate_plot"],
        outputs=["plot", "explanation"],
    )
    result = score_case_gates(PROTEOMICS_CASE, obs, thresholds=THRESHOLDS)
    assert result["gate_pass"] is True, result["gate_fail_reasons"]


def test_fail_proteomics_routed_to_sequence_alignment() -> None:
    obs = _seq_alignment_obs()
    result = score_case_gates(PROTEOMICS_CASE, obs, thresholds=THRESHOLDS)
    assert result["gate_pass"] is False
    # Scorer uses the generic token "forbidden_or_wrong_route"
    assert "forbidden_or_wrong_route" in result["gate_fail_reasons"]


def test_fail_proteomics_missing_provenance() -> None:
    obs = _perfect_tabular_obs(
        steps=["ingest_tabular", "run_pca", "generate_plot"],
        outputs=["plot", "explanation"],
    )
    obs["observed_provenance_signals"] = []  # all provenance stripped
    result = score_case_gates(PROTEOMICS_CASE, obs, thresholds=THRESHOLDS)
    assert result["gate_pass"] is False


# ---------------------------------------------------------------------------
# Suite-level gate: run_suite over all three P0 cases
# ---------------------------------------------------------------------------


def test_run_suite_passes_with_all_ideal_observations(tmp_path: Path) -> None:
    """run_suite should return release_gate_status='pass' when every P0
    case gets a perfect observation set."""
    import json

    from benchmarks.scoring.run_suite import run_suite

    observations = {
        "case-tabular-patient-cohort-correlation": _perfect_tabular_obs(
            steps=["ingest_tabular", "compute_correlation", "generate_heatmap"],
            outputs=["plot", "explanation"],
        ),
        "case-tabular-clinical-trial-group-comparison": _perfect_tabular_obs(
            steps=["ingest_tabular", "group_statistics", "generate_plot"],
            outputs=["plot", "explanation"],
        ),
        "case-routing-safety-proteomics-no-sequence-fallback": _perfect_tabular_obs(
            steps=["ingest_tabular", "run_pca", "generate_plot"],
            outputs=["plot", "explanation"],
        ),
    }

    obs_path = tmp_path / "observations.json"
    obs_path.write_text(json.dumps(observations))

    report_path = tmp_path / "gate_report.json"

    # Only score the three P0 cases to avoid failures from unrelated YAML
    # cases whose fixtures / route expectations differ.
    p0_cases_dir = tmp_path / "p0_cases"
    p0_cases_dir.mkdir()
    for src in [
        CASES_DIR / "tabular" / "tabular-patient-cohort-correlation.yaml",
        CASES_DIR / "tabular" / "tabular-clinical-trial-group-comparison.yaml",
        CASES_DIR / "routing_safety" / "proteomics-no-sequence-fallback.yaml",
    ]:
        (p0_cases_dir / src.name).write_text(src.read_text())

    report = run_suite(
        p0_cases_dir,
        obs_path,
        report_path=report_path,
        thresholds=THRESHOLDS,
    )

    assert report["release_gate_status"] == "pass", json.dumps(report, indent=2)
    assert report["gate_fail_count"] == 0
    assert report["total_cases"] == 3


def test_run_suite_fails_when_critical_case_is_misrouted(tmp_path: Path) -> None:
    """run_suite should return release_gate_status='fail' when a release_gate
    critical case is routed to sequence_alignment."""
    import json

    from benchmarks.scoring.run_suite import run_suite

    observations = {
        "case-tabular-patient-cohort-correlation": _seq_alignment_obs(),
        "case-tabular-clinical-trial-group-comparison": _perfect_tabular_obs(
            steps=["ingest_tabular", "group_statistics", "generate_plot"],
        ),
        "case-routing-safety-proteomics-no-sequence-fallback": _perfect_tabular_obs(
            steps=["ingest_tabular", "run_pca", "generate_plot"],
        ),
    }

    obs_path = tmp_path / "observations.json"
    obs_path.write_text(json.dumps(observations))

    report_path = tmp_path / "gate_report.json"

    p0_cases_dir = tmp_path / "p0_cases"
    p0_cases_dir.mkdir()
    for src in [
        CASES_DIR / "tabular" / "tabular-patient-cohort-correlation.yaml",
        CASES_DIR / "tabular" / "tabular-clinical-trial-group-comparison.yaml",
        CASES_DIR / "routing_safety" / "proteomics-no-sequence-fallback.yaml",
    ]:
        (p0_cases_dir / src.name).write_text(src.read_text())

    report = run_suite(
        p0_cases_dir,
        obs_path,
        report_path=report_path,
        thresholds=THRESHOLDS,
    )

    assert report["release_gate_status"] == "fail"
    assert report["gate_fail_count"] >= 1
    assert "case-tabular-patient-cohort-correlation" in report["critical_failures"]
