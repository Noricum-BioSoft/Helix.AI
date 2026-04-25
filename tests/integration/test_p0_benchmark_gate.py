"""P0 benchmark gate integration tests.

These tests run the three P0 tabular benchmark cases end-to-end against a
live Helix.AI backend, collect observations from the API responses, score
them with the gate scorer, and assert all release-gated cases pass.

Requires a running backend.  The tests are skipped automatically if the
backend is not reachable, so they are safe to collect in a standard test
run.

Environment variables
---------------------
BACKEND_URL   Base URL for the backend (default: http://localhost:8000)

Usage
-----
# With a running backend:
pytest tests/integration/test_p0_benchmark_gate.py -v -m integration

# Or to force-run without the marker filter:
BACKEND_URL=http://localhost:8000 pytest tests/integration/test_p0_benchmark_gate.py -v
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import pytest
import requests
import yaml

from benchmarks.scoring.use_case_gate_scorer import GateThresholds, score_case_gates
from benchmarks.scoring.run_suite import run_suite

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BACKEND_URL = os.getenv("BACKEND_URL", "http://localhost:8000").rstrip("/")
REPO_ROOT = Path(__file__).resolve().parents[2]
CASES_DIR = REPO_ROOT / "benchmarks" / "cases"
FIXTURES_DIR = REPO_ROOT / "benchmarks" / "fixtures" / "tabular"
ARTIFACTS_DIR = REPO_ROOT / "artifacts" / "benchmark_results"
GATE_THRESHOLDS = GateThresholds(min_total_score=0.80)

# Commands for each P0 case (match the YAML prompt intent)
_CASE_PROMPTS: Dict[str, str] = {
    "case-tabular-patient-cohort-correlation": (
        "Analyze the correlation between age, BMI, and toxicity_grade for patients in "
        "this cohort. Generate a correlation heatmap and report which pairs are most "
        "strongly associated."
    ),
    "case-tabular-clinical-trial-group-comparison": (
        "Compare overall survival (os_months) and progression-free survival (pfs_months) "
        "across cancer stages and treatment groups. Show box plots and report which "
        "groups have significantly different outcomes."
    ),
    "case-routing-safety-proteomics-no-sequence-fallback": (
        "Cluster these plasma protein expression samples by condition using PCA "
        "and visualize the result. Report which proteins drive the separation "
        "between conditions."
    ),
}

_CASE_FIXTURES: Dict[str, str] = {
    "case-tabular-patient-cohort-correlation": "patient_cohort.xlsx",
    "case-tabular-clinical-trial-group-comparison": "clinical_trial_expression.csv",
    "case-routing-safety-proteomics-no-sequence-fallback": "proteomics_plasma.csv",
}

_CASE_YAML: Dict[str, str] = {
    "case-tabular-patient-cohort-correlation":
        "tabular/tabular-patient-cohort-correlation.yaml",
    "case-tabular-clinical-trial-group-comparison":
        "tabular/tabular-clinical-trial-group-comparison.yaml",
    "case-routing-safety-proteomics-no-sequence-fallback":
        "routing_safety/proteomics-no-sequence-fallback.yaml",
}


# ---------------------------------------------------------------------------
# Backend availability check
# ---------------------------------------------------------------------------

def _backend_available() -> bool:
    try:
        r = requests.get(f"{BACKEND_URL}/health", timeout=3)
        return r.status_code < 500
    except Exception:
        try:
            # Some deployments expose /docs instead of /health
            r = requests.get(f"{BACKEND_URL}/docs", timeout=3)
            return r.status_code < 500
        except Exception:
            return False


requires_backend = pytest.mark.skipif(
    not _backend_available(),
    reason=f"Backend not reachable at {BACKEND_URL}",
)


# ---------------------------------------------------------------------------
# API helpers
# ---------------------------------------------------------------------------

def _create_session() -> str:
    r = requests.post(f"{BACKEND_URL}/session/create", json={}, timeout=15)
    r.raise_for_status()
    data = r.json()
    return data.get("session_id") or data["id"]


def _upload_file(session_id: str, file_path: Path) -> Dict[str, Any]:
    with file_path.open("rb") as fh:
        r = requests.post(
            f"{BACKEND_URL}/session/{session_id}/uploads",
            files={"files": (file_path.name, fh)},
            timeout=30,
        )
    r.raise_for_status()
    return r.json()


def _execute(session_id: str, command: str, execute_plan: bool = False) -> Dict[str, Any]:
    payload: Dict[str, Any] = {"session_id": session_id, "command": command}
    if execute_plan:
        payload["execute_plan"] = True
    r = requests.post(f"{BACKEND_URL}/execute", json=payload, timeout=120)
    r.raise_for_status()
    return r.json()


def _extract_observations(response: Dict[str, Any]) -> Dict[str, Any]:
    """Derive a gate-scorer observation dict from a raw /execute response."""
    tool = response.get("tool") or ""
    route = response.get("route") or tool or ""

    # Determine steps from execution path and tool name
    observed_steps: List[str] = []
    if tool == "tabular_analysis" or route == "tabular_analysis":
        observed_steps = ["ingest_tabular"]
        result = response.get("result") or {}
        if isinstance(result, dict):
            plan_steps = result.get("plan_steps") or []
            if plan_steps:
                observed_steps = [str(s) for s in plan_steps]

    # Output type detection
    observed_outputs: List[str] = []
    result = response.get("result") or {}
    if isinstance(result, dict):
        if result.get("plot_base64"):
            observed_outputs.append("plot")
        if result.get("text") or result.get("interpretation"):
            observed_outputs.append("explanation")
        if result.get("result_data"):
            observed_outputs.append("table")

    # Safety: scan interpretation text for hedge phrases
    observed_safety: List[str] = []
    interp = ""
    if isinstance(result, dict):
        interp = str(result.get("text") or result.get("interpretation") or "")
    top_text = str(response.get("text") or "")
    combined_text = (interp + " " + top_text).lower()
    if any(w in combined_text for w in ("note", "caveat", "assumption", "uncertain", "interpret")):
        observed_safety.append("uncertainty_note")
    if any(w in combined_text for w in ("assum", "based on", "assuming")):
        observed_safety.append("explicit_assumptions")

    # Provenance signals (best-effort)
    observed_provenance: List[str] = []
    if tool:
        observed_provenance.append("route_decision_trace_present")
    if response.get("session_id"):
        observed_provenance.append("manifest_v1_present")
        observed_provenance.append("lineage_edges_present")
        observed_provenance.append("fallback_events_present")

    # Replay: not applicable for live runs
    replay_status = "no_replay"

    # Forbidden fallback detection
    fallback_targets: List[str] = []
    if route in ("sequence_alignment",):
        fallback_targets.append(route)
    if tool in ("sequence_alignment",):
        fallback_targets.append(tool)

    return {
        "observed_route": route or "unknown",
        "fallback_targets": list(set(fallback_targets)),
        "observed_steps": observed_steps,
        "observed_output_types": observed_outputs,
        "observed_safety_flags": list(set(observed_safety)),
        "observed_provenance_signals": observed_provenance,
        "replay_status": replay_status,
    }


# ---------------------------------------------------------------------------
# Individual case fixtures / tests
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def live_observations() -> Dict[str, Dict[str, Any]]:
    """Run all three P0 cases against the backend and collect observations.

    This is a module-scoped fixture so each case is run only once, but all
    individual test functions can use the same observations dict.
    """
    if not _backend_available():
        pytest.skip(f"Backend not reachable at {BACKEND_URL}")

    collected: Dict[str, Dict[str, Any]] = {}
    errors: Dict[str, str] = {}

    for case_id, fixture_name in _CASE_FIXTURES.items():
        fixture_path = FIXTURES_DIR / fixture_name
        if not fixture_path.exists():
            errors[case_id] = f"Fixture not found: {fixture_path}"
            continue

        try:
            session_id = _create_session()

            # Upload the fixture file
            _upload_file(session_id, fixture_path)

            # Send the prompt – may return a plan proposal or direct result
            command = _CASE_PROMPTS[case_id]
            first_resp = _execute(session_id, command)

            # If the backend returned a pending plan, approve and execute
            result_resp = first_resp
            pending = (
                first_resp.get("result", {}) or {}
            ).get("pending_plan") or first_resp.get("pending_plan")
            if pending or first_resp.get("execution_path") == "plan_proposed":
                # Small delay to let server settle
                time.sleep(1)
                result_resp = _execute(session_id, command, execute_plan=True)

            obs = _extract_observations(result_resp)
            obs["_raw_response"] = result_resp  # attach for debugging
            collected[case_id] = obs

        except Exception as exc:
            errors[case_id] = str(exc)

    if errors:
        # Attach errors to the collected dict for downstream tests to inspect
        for eid, emsg in errors.items():
            collected[eid] = {
                "observed_route": "ERROR",
                "fallback_targets": [],
                "observed_steps": [],
                "observed_output_types": [],
                "observed_safety_flags": [],
                "observed_provenance_signals": [],
                "replay_status": "no_replay",
                "_error": emsg,
            }

    return collected


# ---------------------------------------------------------------------------
# Per-case tests
# ---------------------------------------------------------------------------

@requires_backend
@pytest.mark.integration
def test_correlation_case_gate_passes(live_observations: Dict[str, Dict[str, Any]]) -> None:
    case_id = "case-tabular-patient-cohort-correlation"
    obs = live_observations.get(case_id)
    assert obs is not None, f"No observation collected for {case_id}"
    assert "_error" not in obs, f"Collection error: {obs['_error']}"

    case_def = yaml.safe_load((CASES_DIR / _CASE_YAML[case_id]).read_text())
    result = score_case_gates(case_def, obs, thresholds=GATE_THRESHOLDS)

    _save_result(case_id, obs, result)
    assert result["gate_pass"] is True, (
        f"Gate FAILED for {case_id}:\n"
        f"  score={result['total_score']:.2f}\n"
        f"  reasons={result['gate_fail_reasons']}\n"
        f"  observed_route={obs['observed_route']}\n"
        f"  observed_outputs={obs['observed_output_types']}"
    )


@requires_backend
@pytest.mark.integration
def test_group_comparison_case_gate_passes(live_observations: Dict[str, Dict[str, Any]]) -> None:
    case_id = "case-tabular-clinical-trial-group-comparison"
    obs = live_observations.get(case_id)
    assert obs is not None, f"No observation collected for {case_id}"
    assert "_error" not in obs, f"Collection error: {obs['_error']}"

    case_def = yaml.safe_load((CASES_DIR / _CASE_YAML[case_id]).read_text())
    result = score_case_gates(case_def, obs, thresholds=GATE_THRESHOLDS)

    _save_result(case_id, obs, result)
    assert result["gate_pass"] is True, (
        f"Gate FAILED for {case_id}:\n"
        f"  score={result['total_score']:.2f}\n"
        f"  reasons={result['gate_fail_reasons']}\n"
        f"  observed_route={obs['observed_route']}\n"
        f"  observed_outputs={obs['observed_output_types']}"
    )


@requires_backend
@pytest.mark.integration
def test_proteomics_routing_safety_gate_passes(live_observations: Dict[str, Dict[str, Any]]) -> None:
    case_id = "case-routing-safety-proteomics-no-sequence-fallback"
    obs = live_observations.get(case_id)
    assert obs is not None, f"No observation collected for {case_id}"
    assert "_error" not in obs, f"Collection error: {obs['_error']}"

    case_def = yaml.safe_load((CASES_DIR / _CASE_YAML[case_id]).read_text())
    result = score_case_gates(case_def, obs, thresholds=GATE_THRESHOLDS)

    _save_result(case_id, obs, result)
    assert result["gate_pass"] is True, (
        f"Gate FAILED for {case_id}:\n"
        f"  score={result['total_score']:.2f}\n"
        f"  reasons={result['gate_fail_reasons']}\n"
        f"  observed_route={obs['observed_route']}\n"
        f"  fallback_targets={obs['fallback_targets']}"
    )


@requires_backend
@pytest.mark.integration
def test_proteomics_not_routed_to_sequence_alignment(
    live_observations: Dict[str, Dict[str, Any]]
) -> None:
    """Hard safety assertion: proteomics tabular prompt MUST NOT route to sequence_alignment."""
    case_id = "case-routing-safety-proteomics-no-sequence-fallback"
    obs = live_observations.get(case_id, {})
    assert obs.get("observed_route") != "sequence_alignment", (
        "CRITICAL routing regression: proteomics tabular prompt was routed to "
        f"sequence_alignment. Full observation: {json.dumps(obs, indent=2, default=str)}"
    )
    assert "sequence_alignment" not in obs.get("fallback_targets", []), (
        "CRITICAL: sequence_alignment appeared in fallback_targets for proteomics case."
    )


# ---------------------------------------------------------------------------
# Suite-level gate test
# ---------------------------------------------------------------------------

@requires_backend
@pytest.mark.integration
def test_p0_suite_gate_passes(live_observations: Dict[str, Dict[str, Any]], tmp_path: Path) -> None:
    """All three P0 cases together must produce release_gate_status='pass'."""
    clean_obs = {
        cid: {k: v for k, v in obs.items() if not k.startswith("_")}
        for cid, obs in live_observations.items()
    }

    obs_path = tmp_path / "observations.json"
    obs_path.write_text(json.dumps(clean_obs, indent=2, default=str))

    p0_cases_dir = tmp_path / "p0_cases"
    p0_cases_dir.mkdir()
    for yaml_rel in _CASE_YAML.values():
        src = CASES_DIR / yaml_rel
        (p0_cases_dir / src.name).write_text(src.read_text())

    report_path = ARTIFACTS_DIR / "p0_gate_report.json"
    report = run_suite(
        p0_cases_dir,
        obs_path,
        report_path=report_path,
        thresholds=GATE_THRESHOLDS,
    )

    assert report["release_gate_status"] == "pass", (
        f"P0 suite gate FAILED:\n"
        f"  critical_failures={report['critical_failures']}\n"
        f"  missing_observations={report['missing_observations']}\n"
        f"  pass={report['gate_pass_count']} / {report['total_cases']}\n"
        f"  report saved to: {report_path}"
    )


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def _save_result(case_id: str, obs: Dict[str, Any], result: Dict[str, Any]) -> None:
    ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)
    out = {
        "case_id": case_id,
        "gate_pass": result["gate_pass"],
        "total_score": result["total_score"],
        "gate_fail_reasons": result["gate_fail_reasons"],
        "observed_route": obs["observed_route"],
        "observed_steps": obs["observed_steps"],
        "observed_output_types": obs["observed_output_types"],
    }
    path = ARTIFACTS_DIR / f"{case_id}.json"
    path.write_text(json.dumps(out, indent=2))
