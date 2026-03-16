from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, List

import requests
import yaml

from backend.orchestration.benchmark_controller import (
    BenchmarkLoopConfig,
    run_benchmark_iteration_loop,
)
from benchmarks.bio_benchmark_scorer import score_benchmark_run


def _run_once(api_base: str, benchmark_yaml: Path) -> Dict[str, Any]:
    spec = yaml.safe_load(benchmark_yaml.read_text()) or {}
    scenario = (spec.get("scenarios") or [{}])[0]
    turns = scenario.get("turns") or []
    sid = requests.post(f"{api_base}/session/create", json={}, timeout=20).json()["session_id"]
    rows: List[Dict[str, Any]] = []
    for t in turns:
        prompt = t.get("user", "")
        try:
            r = requests.post(
                f"{api_base}/execute",
                json={"command": prompt, "session_id": sid},
                timeout=180,
            )
            data = r.json()
        except requests.exceptions.Timeout:
            data = {
                "tool": "timeout",
                "status": "timeout",
                "execute_ready": False,
                "text": "Request timed out after 180s.",
                "error": "ReadTimeout",
            }
        except Exception as exc:  # noqa: BLE001
            data = {
                "tool": "error",
                "status": "error",
                "execute_ready": False,
                "text": f"Benchmark runner error: {exc}",
                "error": str(exc),
            }
        rows.append(
            {
                "turn_id": t.get("id"),
                "user": prompt,
                "tool": data.get("tool"),
                "status": data.get("status"),
                "execute_ready": data.get("execute_ready"),
                "text": (data.get("text") or "")[:500],
                "error": data.get("error"),
            }
        )
    return {"session_id": sid, "scenario_id": scenario.get("id"), "rows": rows}


def main() -> None:
    parser = argparse.ArgumentParser(description="Run benchmark improvement loop.")
    parser.add_argument(
        "--api-base",
        default="http://localhost:8000",  # default matches the standard dev-server port
    )
    parser.add_argument("--benchmark-yaml", default="tests/bioinformatics_benchmark.yaml")
    parser.add_argument("--out-dir", default="benchmarks/runs")
    parser.add_argument("--max-iters", type=int, default=3)
    parser.add_argument("--stall-iters", type=int, default=2)
    parser.add_argument("--stall-delta-pct", type=float, default=0.5)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    bench_yaml = Path(args.benchmark_yaml)
    spec = yaml.safe_load(bench_yaml.read_text()) or {}

    cfg = BenchmarkLoopConfig(
        max_iters=args.max_iters,
        threshold_pct=80.0,
        stall_delta_pct=args.stall_delta_pct,
        stall_iters=args.stall_iters,
    )
    loop_result = run_benchmark_iteration_loop(
        run_once=lambda: _run_once(args.api_base, bench_yaml),
        score_once=lambda raw: score_benchmark_run(raw, spec),
        out_dir=out_dir,
        config=cfg,
    )
    history = loop_result.get("history") or []
    for h in history:
        print(f"Iteration {h['iteration']}: Overall Percentage: {h['percentage']}%")
    print(f"Stop reason: {loop_result.get('stop_reason')}")


if __name__ == "__main__":
    main()

