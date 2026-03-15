from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional


@dataclass
class BenchmarkLoopConfig:
    max_iters: int = 6
    threshold_pct: float = 80.0
    stall_delta_pct: float = 0.5
    stall_iters: int = 2


def derive_lessons_and_targets(scored: Dict[str, Any]) -> Dict[str, List[str]]:
    critical = scored.get("critical_failures") or []
    lessons: List[str] = []
    targets: List[str] = []
    if not critical:
        lessons.append("No critical failures in this iteration.")
        targets.append("Preserve behavior with regression tests.")
        return {"lessons_learned": lessons, "remediation_targets": targets}

    by_reason: Dict[str, int] = {}
    for cf in critical:
        reason = str(cf.get("reason") or "unknown")
        by_reason[reason] = by_reason.get(reason, 0) + 1
    for reason, count in sorted(by_reason.items(), key=lambda x: x[1], reverse=True):
        lessons.append(f"{count} critical failure(s): {reason}")
        if "approval" in reason:
            targets.append("Harden approval staging/execution consistency across execution paths.")
        elif "inputs" in reason or "resolver" in reason:
            targets.append("Improve semantic artifact binding and historical reference resolution.")
        elif "planner" in reason:
            targets.append("Improve action-plan generation and execution-intent handling.")
        else:
            targets.append("Harden execution path robustness and fallback behavior.")
    return {"lessons_learned": lessons, "remediation_targets": sorted(set(targets))}


def write_iteration_report(scored: Dict[str, Any], out_file: Path, extra: Optional[Dict[str, Any]] = None) -> None:
    extra = extra or {}
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
    ]
    if extra.get("stop_reason"):
        lines.append(f"- Stop reason: {extra['stop_reason']}")
    lines += ["", "## Lessons Learned"]
    for item in extra.get("lessons_learned", []) or ["-"]:
        lines.append(f"- {item}")
    lines += ["", "## Remediation Targets"]
    for item in extra.get("remediation_targets", []) or ["-"]:
        lines.append(f"- {item}")
    lines += ["", "## Critical Failures"]
    critical = scored.get("critical_failures") or []
    if not critical:
        lines.append("- None")
    else:
        for cf in critical:
            lines.append(
                f"- `{cf.get('turn_id')}` tool=`{cf.get('tool')}` status=`{cf.get('status')}`: {cf.get('reason')}"
            )
    out_file.write_text("\n".join(lines) + "\n")


def run_benchmark_iteration_loop(
    *,
    run_once: Callable[[], Dict[str, Any]],
    score_once: Callable[[Dict[str, Any]], Dict[str, Any]],
    out_dir: Path,
    config: Optional[BenchmarkLoopConfig] = None,
) -> Dict[str, Any]:
    cfg = config or BenchmarkLoopConfig()
    out_dir.mkdir(parents=True, exist_ok=True)
    history: List[Dict[str, Any]] = []
    stall_counter = 0
    prev_pct: Optional[float] = None
    stop_reason = "max_iters_reached"
    last_scored: Dict[str, Any] = {}

    for i in range(1, cfg.max_iters + 1):
        raw = run_once()
        scored = score_once(raw)
        lessons_bundle = derive_lessons_and_targets(scored)

        raw_file = out_dir / f"iteration_{i}_raw.json"
        scored_file = out_dir / f"iteration_{i}_scored.json"
        report_file = out_dir / f"iteration_{i}_report.md"
        raw_file.write_text(json.dumps(raw, indent=2))
        scored_file.write_text(json.dumps(scored, indent=2))

        pct = float(scored.get("overall_percentage") or 0.0)
        crit = scored.get("critical_failures") or []
        threshold_met = pct >= cfg.threshold_pct and len(crit) == 0

        if threshold_met:
            stop_reason = "threshold_met_and_no_critical_failures"
            write_iteration_report(
                scored,
                report_file,
                {**lessons_bundle, "stop_reason": stop_reason},
            )
            history.append({"iteration": i, "percentage": pct, "critical_failures": len(crit)})
            last_scored = scored
            (out_dir / "latest_iteration.json").write_text(
                json.dumps(
                    {
                        "iteration": i,
                        "overall_percentage": pct,
                        "critical_failures": len(crit),
                        "threshold_met": True,
                        "stop_reason": stop_reason,
                        "raw_file": raw_file.name,
                        "scored_file": scored_file.name,
                        "report_file": report_file.name,
                    },
                    indent=2,
                )
            )
            break

        if prev_pct is not None:
            if (pct - prev_pct) <= cfg.stall_delta_pct:
                stall_counter += 1
            else:
                stall_counter = 0
        prev_pct = pct

        if stall_counter >= cfg.stall_iters:
            stop_reason = "progress_stalled"
            write_iteration_report(
                scored,
                report_file,
                {**lessons_bundle, "stop_reason": stop_reason},
            )
            history.append({"iteration": i, "percentage": pct, "critical_failures": len(crit)})
            last_scored = scored
            (out_dir / "latest_iteration.json").write_text(
                json.dumps(
                    {
                        "iteration": i,
                        "overall_percentage": pct,
                        "critical_failures": len(crit),
                        "threshold_met": False,
                        "stop_reason": stop_reason,
                        "raw_file": raw_file.name,
                        "scored_file": scored_file.name,
                        "report_file": report_file.name,
                    },
                    indent=2,
                )
            )
            break

        write_iteration_report(scored, report_file, lessons_bundle)
        history.append({"iteration": i, "percentage": pct, "critical_failures": len(crit)})
        last_scored = scored
        (out_dir / "latest_iteration.json").write_text(
            json.dumps(
                {
                    "iteration": i,
                    "overall_percentage": pct,
                    "critical_failures": len(crit),
                    "threshold_met": False,
                    "stop_reason": "in_progress",
                    "raw_file": raw_file.name,
                    "scored_file": scored_file.name,
                    "report_file": report_file.name,
                },
                indent=2,
            )
        )

    loop_summary = {"stop_reason": stop_reason, "history": history, "last_scored": last_scored}
    (out_dir / "loop_summary.json").write_text(json.dumps(loop_summary, indent=2))
    return loop_summary

