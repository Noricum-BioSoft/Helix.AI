from __future__ import annotations

import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import yaml


@dataclass
class ScoreConfig:
    threshold_pct: float = 80.0


TOOL_CLASS_BY_NAME = {
    "__plan__": "plan",
    "bio_rerun": "rerun",
    "patch_and_rerun": "rerun",
    "bio_diff_runs": "comparison",
    "lookup_go_term": "lookup",
    "query_uniprot": "lookup",
    "bulk_rnaseq_analysis": "analysis",
    "single_cell_analysis": "analysis",
    "handle_natural_command": "agent",
}

UNRESOLVED_PATTERNS = [
    "could not resolve",
    "could not find prior",
    "run at least two persisted analyses first",
    "no historical artifact found",
    "no comparable run found",
    "please run the base analysis first",
    "could not resolve comparable run ids",
]

GENERIC_FALLBACK_PATTERNS = [
    "please run the base analysis first",
    "run additional analyses in this session",
    "could not resolve",
    "enable the agent",
    "rephrase as",
]

NON_ACTIONABLE_PATTERNS = [
    "re-run complete.",
    "analysis complete.",
    "done.",
]


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


def _infer_expectation_category(turn_def: Dict[str, Any]) -> str:
    prompt = str(turn_def.get("user") or "").strip().lower()
    behaviors = [str(b).lower() for b in (turn_def.get("expected_behaviors") or [])]
    if prompt.startswith("approve"):
        return "approval_execution_expected"
    if turn_def.get("should_require_approval_before_execution"):
        if any("detect_under_specified_request" in b for b in behaviors):
            return "clarification_expected"
        return "planning_expected"
    if any("resolve_exact_historical_state" in b for b in behaviors):
        return "historical_recreation_expected"
    if any("compare" in b for b in behaviors):
        return "comparison_expected"
    if any("go_enrichment" in b or "pathway_enrichment" in b for b in behaviors):
        return "downstream_enrichment_expected"
    if any("rerun" in b for b in behaviors):
        return "rerun_expected"
    if any("plot" in b or "visualization" in b or "heatmap" in b or "pca" in b for b in behaviors):
        return "visualization_update_expected"
    if any("clarif" in b for b in behaviors):
        return "clarification_expected"
    return "execution_expected"


def _tool_class(tool_name: Any) -> str:
    return TOOL_CLASS_BY_NAME.get(str(tool_name or ""), "unknown")


def _is_empty_or_weak_text(text: str) -> bool:
    t = (text or "").strip()
    return not t


def _detect_failure_patterns(
    *,
    turn_def: Dict[str, Any],
    turn_obs: Dict[str, Any],
    expectation_category: str,
) -> Set[str]:
    patterns: Set[str] = set()
    status = str(turn_obs.get("status") or "").lower()
    tool = str(turn_obs.get("tool") or "")
    tool_class = _tool_class(tool)
    text = str(turn_obs.get("text") or "")
    text_l = text.lower().strip()
    hard_error = bool(turn_obs.get("error") or turn_obs.get("raw_error"))

    if hard_error or status in {"error", "failed", "workflow_failed", "timeout"}:
        patterns.add("hard_execution_failure")

    if _is_empty_or_weak_text(text):
        patterns.add("empty_output")
    elif len(text_l) < 24:
        patterns.add("non_actionable_output")

    if any(p in text_l for p in NON_ACTIONABLE_PATTERNS):
        patterns.add("non_actionable_output")

    if "synthetic demo data used" in text_l:
        patterns.add("synthetic_success_only")

    if any(p in text_l for p in UNRESOLVED_PATTERNS):
        patterns.add("unresolved_historical_reference")
        patterns.add("missing_artifact_resolution")

    if any(p in text_l for p in GENERIC_FALLBACK_PATTERNS):
        patterns.add("fallback_instead_of_execution")

    if expectation_category in {
        "execution_expected",
        "rerun_expected",
        "comparison_expected",
        "historical_recreation_expected",
        "downstream_enrichment_expected",
        "visualization_update_expected",
        "approval_execution_expected",
    } and status == "workflow_planned":
        patterns.add("plan_instead_of_execution")

    if expectation_category in {"comparison_expected", "historical_recreation_expected"}:
        if tool_class not in {"comparison", "rerun", "agent"}:
            patterns.add("wrong_tool_class")
        if "unresolved_historical_reference" in patterns:
            patterns.add("comparison_not_performed")
        if text_l and all(
            token not in text_l for token in ("compare", "comparison", "changed", "difference", "run a", "run b")
        ):
            patterns.add("output_not_specific")

    if expectation_category == "downstream_enrichment_expected":
        if tool_class == "lookup":
            patterns.add("enrichment_substituted_with_lookup")
            patterns.add("wrong_tool_class")
        if "enrichment" not in text_l and "go" not in text_l:
            patterns.add("output_not_specific")

    if expectation_category in {"rerun_expected", "visualization_update_expected"}:
        if tool_class not in {"rerun", "analysis", "agent"}:
            patterns.add("wrong_tool_class")
        if "fallback_instead_of_execution" in patterns:
            patterns.add("rerun_not_performed")

    if expectation_category == "approval_execution_expected":
        if status != "success":
            patterns.add("fallback_instead_of_execution")
        if "executed" not in text_l and "submitted" not in text_l:
            patterns.add("output_not_specific")

    if expectation_category == "planning_expected":
        if status != "workflow_planned":
            patterns.add("wrong_tool_class")
        if "pipeline plan" not in text_l and "steps" not in text_l:
            patterns.add("non_actionable_output")

    if expectation_category == "clarification_expected":
        if status not in {"workflow_planned", "success", "needs_inputs"}:
            patterns.add("non_actionable_output")
        if not any(t in text_l for t in ("plan", "steps", "assumption", "input", "clarif", "workflow")):
            patterns.add("output_not_specific")

    return patterns


def _assess_tool_alignment(expectation_category: str, tool_name: Any, failure_patterns: Set[str]) -> str:
    tool_class = _tool_class(tool_name)
    if "wrong_tool_class" in failure_patterns or "enrichment_substituted_with_lookup" in failure_patterns:
        return "misaligned"
    if expectation_category == "planning_expected":
        return "aligned" if tool_class == "plan" else "partial"
    if expectation_category == "comparison_expected":
        return "aligned" if tool_class in {"comparison", "agent"} else "partial"
    if expectation_category in {"rerun_expected", "visualization_update_expected"}:
        return "aligned" if tool_class in {"rerun", "analysis", "agent"} else "partial"
    return "aligned" if tool_class != "unknown" else "partial"


def _assess_execution_vs_planning(expectation_category: str, status: str, failure_patterns: Set[str]) -> str:
    if "plan_instead_of_execution" in failure_patterns:
        return "incorrect_planning_instead_of_execution"
    if expectation_category == "planning_expected" and status == "workflow_planned":
        return "correct_planning"
    if expectation_category == "approval_execution_expected" and status == "success":
        return "approved_execution_happened"
    if expectation_category in {
        "execution_expected",
        "rerun_expected",
        "comparison_expected",
        "downstream_enrichment_expected",
        "visualization_update_expected",
        "historical_recreation_expected",
    } and status == "success":
        return "execution_happened"
    return "uncertain"


def _assess_artifact_resolution(expectation_category: str, failure_patterns: Set[str]) -> str:
    if expectation_category in {"comparison_expected", "historical_recreation_expected", "rerun_expected"}:
        if "missing_artifact_resolution" in failure_patterns or "unresolved_historical_reference" in failure_patterns:
            return "failed"
        return "resolved_or_not_required"
    return "not_primary"


def _assess_output_usefulness(text: str, failure_patterns: Set[str]) -> str:
    if "empty_output" in failure_patterns:
        return "empty"
    if "non_actionable_output" in failure_patterns or "output_not_specific" in failure_patterns:
        return "weak"
    if len((text or "").strip()) < 40:
        return "minimal"
    return "specific"


def _score_from_patterns(
    *,
    expectation_category: str,
    status: str,
    failure_patterns: Set[str],
) -> int:
    if "hard_execution_failure" in failure_patterns:
        return 0
    if expectation_category == "planning_expected":
        if status == "workflow_planned" and not {"non_actionable_output", "empty_output"} & failure_patterns:
            return 2
        return 1 if status in {"workflow_planned", "needs_inputs"} else 0
    if expectation_category == "clarification_expected":
        if status in {"workflow_planned", "success", "needs_inputs"} and not (
            failure_patterns & {"empty_output", "non_actionable_output", "output_not_specific"}
        ):
            return 2
        if status in {"workflow_planned", "success", "needs_inputs"}:
            return 1
        return 0

    severe_patterns = {
        "plan_instead_of_execution",
        "unresolved_historical_reference",
        "missing_artifact_resolution",
        "comparison_not_performed",
        "rerun_not_performed",
        "enrichment_substituted_with_lookup",
        "empty_output",
        "fallback_instead_of_execution",
    }
    medium_patterns = {
        "wrong_tool_class",
        "non_actionable_output",
        "output_not_specific",
        "synthetic_success_only",
    }
    if failure_patterns & severe_patterns:
        return 0
    if failure_patterns & medium_patterns:
        return 1
    if status == "success":
        return 2
    return 1


def _patterns_to_root_causes(patterns: Set[str]) -> List[str]:
    mapping = {
        "unresolved_historical_reference": "reference_resolver_gap",
        "missing_artifact_resolution": "asset_registry_gap",
        "fallback_instead_of_execution": "executor_gap",
        "wrong_tool_class": "planner_gap",
        "plan_instead_of_execution": "planner_gap",
        "comparison_not_performed": "lineage_gap",
        "rerun_not_performed": "executor_gap",
        "enrichment_substituted_with_lookup": "planner_gap",
        "empty_output": "error_handling_gap",
        "non_actionable_output": "error_handling_gap",
        "output_not_specific": "lineage_gap",
        "synthetic_success_only": "state_session_gap",
        "hard_execution_failure": "executor_gap",
    }
    return sorted({mapping[p] for p in patterns if p in mapping})


def _score_turn(turn_def: Dict[str, Any], turn_obs: Dict[str, Any]) -> Dict[str, Any]:
    """
    Heuristic scorer for benchmark turn logs.
    Scale: 0 (fail), 1 (partial), 2 (pass).
    """
    status = str(turn_obs.get("status") or "").lower()
    tool = turn_obs.get("tool")
    text = str(turn_obs.get("text") or "")
    expectation_category = _infer_expectation_category(turn_def)
    failure_patterns = _detect_failure_patterns(
        turn_def=turn_def,
        turn_obs=turn_obs,
        expectation_category=expectation_category,
    )
    score = _score_from_patterns(
        expectation_category=expectation_category,
        status=status,
        failure_patterns=failure_patterns,
    )
    root_causes = _patterns_to_root_causes(failure_patterns)
    tool_alignment = _assess_tool_alignment(expectation_category, tool, failure_patterns)
    execution_assessment = _assess_execution_vs_planning(expectation_category, status, failure_patterns)
    artifact_assessment = _assess_artifact_resolution(expectation_category, failure_patterns)
    output_assessment = _assess_output_usefulness(text, failure_patterns)

    rationale_parts = [
        f"expectation={expectation_category}",
        f"tool_alignment={tool_alignment}",
        f"execution_vs_planning={execution_assessment}",
        f"artifact_resolution={artifact_assessment}",
        f"output_usefulness={output_assessment}",
    ]
    if failure_patterns:
        rationale_parts.append("deductions=" + ",".join(sorted(failure_patterns)))
    else:
        rationale_parts.append("deductions=none")

    return {
        "turn_id": turn_def.get("id"),
        "prompt": turn_def.get("user"),
        "tool": tool,
        "status": status or None,
        "score": score,
        "max_score": 2,
        "expectation_category": expectation_category,
        "focus_dimensions": turn_def.get("focus_dimensions", []) or [],
        "root_cause_categories": sorted(set(root_causes)),
        "failure_patterns": sorted(failure_patterns),
        "assessments": {
            "tool_action_alignment": tool_alignment,
            "artifact_state_resolution": artifact_assessment,
            "execution_vs_planning": execution_assessment,
            "output_usefulness": output_assessment,
        },
        "rationale": "; ".join(rationale_parts),
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

    # Defensive invariant checks — catch silent mis-scoring bugs early.
    assert 0 <= raw <= max_score, (
        f"Scorer invariant violated: raw_score {raw} is outside [0, {max_score}]"
    )
    assert 0.0 <= pct <= 100.0, (
        f"Scorer invariant violated: overall_percentage {pct} is outside [0, 100]"
    )
    for _t in per_turn:
        assert _t["score"] <= _t["max_score"], (
            f"Turn {_t['turn_id']}: score {_t['score']} > max_score {_t['max_score']}"
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

