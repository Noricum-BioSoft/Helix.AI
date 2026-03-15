from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import requests
import yaml

from benchmarks.bio_benchmark_scorer import score_benchmark_run


def _run_once(api_base: str, benchmark_yaml: Path) -> Dict[str, Any]:
    spec = yaml.safe_load(benchmark_yaml.read_text()) or {}
    scenario = (spec.get("scenarios") or [{}])[0]
    turns = scenario.get("turns") or []
    sid = requests.post(f"{api_base}/session/create", json={}, timeout=20).json()["session_id"]
    rows: List[Dict[str, Any]] = []
    for t in turns:
        prompt = t.get("user", "")
        r = requests.post(
            f"{api_base}/execute",
            json={"command": prompt, "session_id": sid},
            timeout=180,
        )
        data = r.json()
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


def _write_markdown_report(scored: Dict[str, Any], out_file: Path) -> None:
    lines = [
        f"# Iteration Report ({out_file.stem})",
        "",
        f"- Generated: {datetime.now(timezone.utc).isoformat()}",
        f"- Session: {scored.get('session_id')}",
        f"- Scenario: {scored.get('scenario_id')}",
        f"- Raw Score: {scored.get('raw_score')} / {scored.get('max_score')}",
        f"- Overall Percentage: {scored.get('overall_percentage')}%",
        f"- Threshold Met: {'Yes' if scored.get('threshold_met') else 'No'}",
        f"- Critical failures: {len(scored.get('critical_failures') or [])}",
        "",
        "## Critical Failures",
    ]
    cfs = scored.get("critical_failures") or []
    if not cfs:
        lines.append("- None")
    else:
        for cf in cfs:
            lines.append(
                f"- `{cf.get('turn_id')}` tool=`{cf.get('tool')}` status=`{cf.get('status')}`: {cf.get('reason')}"
            )
    out_file.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run benchmark improvement loop.")
    parser.add_argument("--api-base", default="http://localhost:8001")
    parser.add_argument("--benchmark-yaml", default="tests/bioinformatics_benchmark.yaml")
    parser.add_argument("--out-dir", default="benchmarks/runs")
    parser.add_argument("--max-iters", type=int, default=3)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    bench_yaml = Path(args.benchmark_yaml)
    spec = yaml.safe_load(bench_yaml.read_text()) or {}

    for i in range(1, args.max_iters + 1):
        raw = _run_once(args.api_base, bench_yaml)
        raw_file = out_dir / f"iteration_{i}_raw.json"
        scored_file = out_dir / f"iteration_{i}_scored.json"
        report_file = out_dir / f"iteration_{i}_report.md"

        raw_file.write_text(json.dumps(raw, indent=2))
        scored = score_benchmark_run(raw, spec)
        scored_file.write_text(json.dumps(scored, indent=2))
        _write_markdown_report(scored, report_file)

        print(f"Iteration {i}: Raw Score: {scored['raw_score']} / {scored['max_score']}")
        print(f"Iteration {i}: Overall Percentage: {scored['overall_percentage']}%")
        print(f"Iteration {i}: Threshold Met: {'Yes' if scored['threshold_met'] else 'No'}")

        if scored.get("threshold_met"):
            break


if __name__ == "__main__":
    main()

