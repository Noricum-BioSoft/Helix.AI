from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
ARTIFACTS_DIR = REPO_ROOT / "artifacts"
WORKLOG_PATH = ARTIFACTS_DIR / "worklog.md"
TEST_RESULTS_DIR = ARTIFACTS_DIR / "test_results"
FAILURE_SUMMARIES_DIR = ARTIFACTS_DIR / "failure_summaries"
BENCHMARK_RESULTS_DIR = ARTIFACTS_DIR / "benchmark_results"
RELEASE_READINESS_PATH = ARTIFACTS_DIR / "release_readiness.json"


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _ensure_file(path: Path, content: str) -> None:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)


def _default_release_readiness() -> dict:
    return {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "working_tree": "unknown",
        "suites_executed": [],
        "pass_fail_counts": {
            "passed": 0,
            "failed": 0,
            "skipped": 0,
        },
        "benchmark_metrics": {},
        "threshold_comparison": {
            "meets_all_thresholds": False,
            "details": [],
        },
        "blockers": [
            "Release readiness not yet evaluated for this working tree."
        ],
        "blocker_summary": "Release readiness not yet evaluated.",
        "final_readiness_value": "not_ready",
        "release_ready": False,
    }


def _merge_release_defaults(existing: dict) -> dict:
    merged = _default_release_readiness()
    merged.update(existing or {})
    threshold = merged.get("threshold_comparison")
    if not isinstance(threshold, dict):
        threshold = {}
    threshold_defaults = _default_release_readiness()["threshold_comparison"]
    threshold_defaults.update(threshold)
    merged["threshold_comparison"] = threshold_defaults
    return merged


def bootstrap() -> None:
    _ensure_dir(ARTIFACTS_DIR)
    _ensure_dir(TEST_RESULTS_DIR)
    _ensure_dir(FAILURE_SUMMARIES_DIR)
    _ensure_dir(BENCHMARK_RESULTS_DIR)

    _ensure_file(
        WORKLOG_PATH,
        "# Worklog\n\n- Initialized artifact tracking.\n",
    )
    _ensure_file(
        FAILURE_SUMMARIES_DIR / "latest.md",
        "# Latest Failure Summary\n\nNo failure summary recorded yet.\n",
    )
    RELEASE_READINESS_PATH.parent.mkdir(parents=True, exist_ok=True)
    if RELEASE_READINESS_PATH.exists():
        try:
            existing = json.loads(RELEASE_READINESS_PATH.read_text())
        except Exception:  # noqa: BLE001
            existing = {}
        merged = _merge_release_defaults(existing if isinstance(existing, dict) else {})
        RELEASE_READINESS_PATH.write_text(json.dumps(merged, indent=2))
    else:
        RELEASE_READINESS_PATH.write_text(
            json.dumps(_default_release_readiness(), indent=2)
        )


def main() -> int:
    bootstrap()
    print("Artifacts bootstrap completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
