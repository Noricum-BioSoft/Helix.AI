from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict


REPO_ROOT = Path(__file__).resolve().parents[2]
ARTIFACTS_DIR = REPO_ROOT / "artifacts"
WORKLOG_PATH = ARTIFACTS_DIR / "worklog.md"
TEST_RESULTS_DIR = ARTIFACTS_DIR / "test_results"
FAILURE_SUMMARIES_DIR = ARTIFACTS_DIR / "failure_summaries"
BENCHMARK_RESULTS_DIR = ARTIFACTS_DIR / "benchmark_results"
RELEASE_READINESS_PATH = ARTIFACTS_DIR / "release_readiness.json"


def _require_path_exists(path: Path, label: str) -> None:
    if not path.exists():
        raise ValueError(f"Missing required artifact: {label} ({path})")


def _require_release_contract(data: Dict[str, Any]) -> None:
    required = {
        "timestamp",
        "working_tree",
        "suites_executed",
        "pass_fail_counts",
        "benchmark_metrics",
        "threshold_comparison",
        "blocker_summary",
        "final_readiness_value",
        "release_ready",
    }
    missing = sorted(required.difference(data.keys()))
    if missing:
        raise ValueError(f"release_readiness.json missing required keys: {missing}")

    threshold = data.get("threshold_comparison") or {}
    if not isinstance(threshold, dict):
        raise ValueError("threshold_comparison must be an object")
    if "meets_all_thresholds" not in threshold:
        raise ValueError("threshold_comparison.meets_all_thresholds is required")
    if "details" not in threshold:
        raise ValueError("threshold_comparison.details is required")

    blocker_summary = str(data.get("blocker_summary") or "").strip()
    if not blocker_summary:
        raise ValueError("blocker_summary must be non-empty")

    readiness_value = str(data.get("final_readiness_value") or "").strip().lower()
    if readiness_value not in {"ready", "not_ready", "blocked"}:
        raise ValueError("final_readiness_value must be one of: ready, not_ready, blocked")


def validate() -> int:
    _require_path_exists(WORKLOG_PATH, "artifacts/worklog.md")
    _require_path_exists(TEST_RESULTS_DIR, "artifacts/test_results/")
    _require_path_exists(FAILURE_SUMMARIES_DIR, "artifacts/failure_summaries/")
    _require_path_exists(BENCHMARK_RESULTS_DIR, "artifacts/benchmark_results/")
    _require_path_exists(RELEASE_READINESS_PATH, "artifacts/release_readiness.json")

    try:
        data = json.loads(RELEASE_READINESS_PATH.read_text())
    except Exception as exc:  # noqa: BLE001
        raise ValueError(f"release_readiness.json is invalid JSON: {exc}") from exc

    _require_release_contract(data)
    print("Artifacts readiness validation passed.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(validate())
    except Exception as exc:  # noqa: BLE001
        print(f"Artifacts readiness validation failed: {exc}")
        raise SystemExit(1)
