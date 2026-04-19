"""Aggregate gate scoring over a suite of scenario YAML files.

Inputs:
- a directory of scenario files (or list of files),
- a JSON file containing run observations keyed by `case_id`.

Outputs:
- a consolidated JSON report under `artifacts/benchmark_results/gate_report.json`.
- a non-zero exit code if any release-gated scenario fails.

Run observations schema (per case_id):
{
  "observed_route": str,
  "fallback_targets": [str, ...],
  "observed_steps": [str, ...],
  "observed_output_types": [str, ...],
  "observed_safety_flags": [str, ...],
  "observed_provenance_signals": [str, ...],
  "replay_status": str
}
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import yaml

from benchmarks.scoring.use_case_gate_scorer import (
    GateThresholds,
    score_case_gates,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CASES_DIR = REPO_ROOT / "benchmarks" / "cases"
DEFAULT_REPORT_DIR = REPO_ROOT / "artifacts" / "benchmark_results"


def _iter_scenario_files(root: Path) -> List[Path]:
    root = root.resolve()
    if root.is_file():
        return [root]
    return sorted(p.resolve() for p in root.rglob("*.yaml") if p.is_file())


def _safe_relpath(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _load_yaml(path: Path) -> Dict[str, Any]:
    return yaml.safe_load(path.read_text()) or {}


def _is_release_gated(case_def: Dict[str, Any]) -> bool:
    gate = case_def.get("release_gate") or {}
    return bool(gate.get("must_pass")) or str(gate.get("severity", "")).lower() == "critical"


def run_suite(
    cases_dir: Path,
    observations_path: Path,
    *,
    report_path: Path | None = None,
    thresholds: GateThresholds | None = None,
) -> Dict[str, Any]:
    scenario_files = _iter_scenario_files(cases_dir)
    observations: Dict[str, Dict[str, Any]] = json.loads(observations_path.read_text())

    results: List[Dict[str, Any]] = []
    missing_observations: List[str] = []
    critical_fail: List[str] = []

    for path in scenario_files:
        case_def = _load_yaml(path)
        case_id = case_def.get("id")
        if not case_id:
            continue
        obs = observations.get(case_id)
        if obs is None:
            missing_observations.append(case_id)
            continue
        result = score_case_gates(case_def, obs, thresholds=thresholds)
        result["scenario_path"] = _safe_relpath(path)
        result["release_gated"] = _is_release_gated(case_def)
        results.append(result)
        if result["release_gated"] and not result["gate_pass"]:
            critical_fail.append(case_id)

    total = len(results)
    passed = sum(1 for r in results if r["gate_pass"])
    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "suite_root": _safe_relpath(cases_dir.resolve()),
        "total_cases": total,
        "gate_pass_count": passed,
        "gate_fail_count": total - passed,
        "missing_observations": missing_observations,
        "critical_failures": critical_fail,
        "release_gate_status": "fail" if (critical_fail or missing_observations) else ("pass" if total else "empty"),
        "results": results,
    }

    target = report_path or (DEFAULT_REPORT_DIR / "gate_report.json")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2))
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases-dir", type=Path, default=DEFAULT_CASES_DIR)
    parser.add_argument("--observations", type=Path, required=True)
    parser.add_argument("--report", type=Path, default=None)
    parser.add_argument("--min-total-score", type=float, default=0.80)
    args = parser.parse_args()

    thresholds = GateThresholds(min_total_score=args.min_total_score)
    report = run_suite(args.cases_dir, args.observations, report_path=args.report, thresholds=thresholds)

    print(f"total={report['total_cases']} pass={report['gate_pass_count']} fail={report['gate_fail_count']}")
    if report["missing_observations"]:
        print("missing_observations=" + ",".join(report["missing_observations"]))
    if report["critical_failures"]:
        print("critical_failures=" + ",".join(report["critical_failures"]))

    return 0 if report["release_gate_status"] == "pass" else 1


if __name__ == "__main__":
    sys.exit(main())
