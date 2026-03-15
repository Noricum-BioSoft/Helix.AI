from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml


@dataclass
class ScoreConfig:
    threshold_pct: float = 80.0


def _load_yaml(path: Path) -> Dict[str, Any]:
    return yaml.safe_load(path.read_text()) or {}


def _index_turn_defs(benchmark_spec: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    turns: Dict[str, Dict[str, Any]] = {}
    for scenario in benchmark_spec.get("scenarios", []) or []:
        for turn in scenario.get("turns", []) or []:
            tid = turn.get("id")
            if tid:
                turns[tid] = turn
    return turns


def _score_turn(turn_def: Dict[str, Any], turn_obs: Dict[str, Any]) -> Dict[str, Any]:
    """
    Heuristic scorer for benchmark turn logs.
    Scale: 0 (fail), 1 (partial), 2 (pass).
    """
    status = str(turn_obs.get("status") or "").lower()
    tool = turn_obs.get("tool")
    prompt = (turn_def.get("user") or "").strip().lower()
    approval_required = bool(turn_def.get("should_require_approval_before_execution"))
    is_approve_turn = prompt.startswith("approve")

    score = 0
    rationale: List[str] = []
    root_causes: List[str] = []

    # Baseline transport/result sanity.
    if status in {"success", "workflow_planned", "needs_inputs"}:
        score = 1
        rationale.append("response returned without transport failure")
    else:
        rationale.append("hard error or malformed response")
        root_causes.append("executor_gap")

    # Approval-aware scoring.
    if approval_required and not is_approve_turn:
        if status == "workflow_planned":
            score = 2
            rationale.append("approval respected via staged plan")
        elif status == "success":
            score = 0
            rationale.append("executed despite approval requirement")
            root_causes.append("approval_loop_gap")
        elif status == "needs_inputs":
            score = 0
            rationale.append("blocked instead of staging approval-ready plan")
            root_causes.append("planner_gap")
    elif is_approve_turn:
        if status == "success":
            score = 2
            rationale.append("approve turn executed successfully")
        else:
            score = 0
            rationale.append("approve turn did not execute pending plan")
            root_causes.append("approval_loop_gap")
    else:
        # Execution-expected turns.
        if status == "success":
            score = 2
            rationale.append("execution turn succeeded")
        elif status == "workflow_planned":
            score = 0
            rationale.append("execution turn returned planning only")
            root_causes.append("planner_gap")
        elif status == "needs_inputs":
            score = 0
            rationale.append("execution turn blocked on unresolved inputs")
            root_causes.append("reference_resolver_gap")

    return {
        "turn_id": turn_def.get("id"),
        "prompt": turn_def.get("user"),
        "tool": tool,
        "status": status or None,
        "score": score,
        "max_score": 2,
        "focus_dimensions": turn_def.get("focus_dimensions", []) or [],
        "root_cause_categories": sorted(set(root_causes)),
        "rationale": "; ".join(rationale),
        "evidence": {
            "tool": tool,
            "status": status or None,
            "execute_ready": turn_obs.get("execute_ready"),
            "error": turn_obs.get("error") or turn_obs.get("raw_error"),
            "text_snippet": (turn_obs.get("text") or "")[:240],
            "duration_s": turn_obs.get("dur"),
        },
    }


def score_benchmark_run(
    run_data: Dict[str, Any],
    benchmark_spec: Dict[str, Any],
    *,
    config: Optional[ScoreConfig] = None,
) -> Dict[str, Any]:
    cfg = config or ScoreConfig()
    turn_defs = _index_turn_defs(benchmark_spec)
    turns_obs = run_data.get("rows", []) or []
    scenario_id = run_data.get("scenario_id") or (
        (benchmark_spec.get("scenarios", [{}])[0] or {}).get("id")
    )

    per_turn: List[Dict[str, Any]] = []
    for obs in turns_obs:
        tid = obs.get("turn_id")
        td = turn_defs.get(tid)
        if not td:
            continue
        per_turn.append(_score_turn(td, obs))

    raw = sum(t["score"] for t in per_turn)
    max_score = sum(t["max_score"] for t in per_turn) or 1
    pct = round((raw / max_score) * 100.0, 2)

    # Per-dimension.
    dim_scores: Dict[str, List[int]] = {}
    for t in per_turn:
        for d in t.get("focus_dimensions", []) or []:
            dim_scores.setdefault(d, []).append(int(t["score"]))
    per_dimension = {
        d: round(sum(vals) / len(vals), 3)
        for d, vals in dim_scores.items()
        if vals
    }

    # Severity-aware critical failures.
    severity_map = {
        d.get("id"): d.get("severity", "medium")
        for d in benchmark_spec.get("evaluation_dimensions", []) or []
        if isinstance(d, dict)
    }
    critical_failures: List[Dict[str, Any]] = []
    for t in per_turn:
        if t["score"] == 0:
            failed_high_dims = [
                d for d in t.get("focus_dimensions", []) if severity_map.get(d) == "high"
            ]
            critical_failures.append(
                {
                    "turn_id": t["turn_id"],
                    "tool": t["tool"],
                    "status": t["status"],
                    "high_severity_dimensions": failed_high_dims,
                    "reason": t["rationale"],
                }
            )

    threshold_met = pct >= cfg.threshold_pct and len(critical_failures) == 0
    scenario_result = (
        "Pass"
        if threshold_met
        else ("Partial" if pct >= 50.0 else "Fail")
    )

    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "scenario_id": scenario_id,
        "session_id": run_data.get("session_id"),
        "per_turn_scores": per_turn,
        "per_dimension_scores": per_dimension,
        "per_scenario_scores": [
            {
                "scenario_id": scenario_id,
                "raw_score": raw,
                "max_score": max_score,
                "overall_percentage": pct,
                "result": scenario_result,
            }
        ],
        "raw_score": raw,
        "max_score": max_score,
        "overall_percentage": pct,
        "critical_failures": critical_failures,
        "threshold": cfg.threshold_pct,
        "threshold_met": threshold_met,
    }


def score_run_file(run_file: Path, benchmark_yaml: Path) -> Dict[str, Any]:
    run_data = json.loads(run_file.read_text())
    spec = _load_yaml(benchmark_yaml)
    return score_benchmark_run(run_data, spec)

