"""
Reviewer: produces a concise natural-language summary of what happened in a run.

No LLM required — template-driven from structured run_data.
"""
from __future__ import annotations

from typing import Any, Dict, List


def review(run_data: Dict[str, Any]) -> str:
    """
    Generate a human-readable review summary for the run.

    Covers: what ran, what changed, key results, decision, what to do next.
    """
    lines: List[str] = []

    run_id = run_data.get("run_id", "unknown")
    objective = run_data.get("objective", "")
    hypothesis = run_data.get("hypothesis", "")
    changes = run_data.get("changes", "")
    steps_run = run_data.get("steps_run") or []
    metrics = run_data.get("metrics") or {}
    decision = run_data.get("decision", "")
    next_steps = run_data.get("next_steps") or []
    comparison = run_data.get("comparison") or {}
    validation = run_data.get("validation_results") or {}
    eda = run_data.get("eda_summary") or {}
    config = run_data.get("config") or {}
    train_info = run_data.get("train_info") or {}

    lines.append(f"## Run Summary: {run_id}\n")
    lines.append(f"**Objective:** {objective}")
    if hypothesis:
        lines.append(f"**Hypothesis:** {hypothesis}")
    if changes:
        lines.append(f"**What changed:** {changes}")
    lines.append(f"**Steps completed:** {', '.join(steps_run) or 'none'}\n")

    # Data
    shape = eda.get("shape", [])
    if shape:
        lines.append(f"**Dataset shape after cleaning:** {shape[0]} rows × {shape[1]} columns")

    val_errors = validation.get("n_errors", 0)
    val_warnings = validation.get("n_warnings", 0)
    if val_errors or val_warnings:
        lines.append(f"**Validation:** {val_errors} error(s), {val_warnings} warning(s)")
    else:
        lines.append("**Validation:** All checks passed ✅")

    # Model & metrics
    task_type = train_info.get("task_type", "")
    model_class = train_info.get("model_class", "")
    if task_type and model_class:
        n_features = train_info.get("n_features", "?")
        lines.append(
            f"**Baseline model:** {model_class} ({task_type}) on {n_features} numeric features"
        )

    if metrics:
        metric_parts = [f"{k}={v:.4f}" for k, v in metrics.items() if isinstance(v, (int, float))]
        lines.append(f"**Metrics:** {', '.join(metric_parts)}")

    # Comparison
    if comparison and comparison.get("prior_best"):
        is_new = comparison.get("is_new_best", False)
        pm = comparison.get("primary_metric", "")
        delta = (comparison.get("deltas") or {}).get(pm)
        if is_new:
            lines.append(f"**Result:** ✅ New best run! ({pm} Δ={delta:+.4f})")
        else:
            lines.append(f"**Result:** ℹ️ No improvement vs prior best ({pm} Δ={delta:+.4f})")
    elif decision == "first_run":
        lines.append("**Result:** First completed run — baseline established.")

    # Decision
    decision_labels = {
        "new_best": "Continue — new best achieved",
        "first_run": "Continue — baseline established",
        "no_change": "Investigate — no improvement; try a different approach",
        "regressed": "Investigate — metrics regressed; check for data issues",
    }
    lines.append(f"\n**Decision:** {decision_labels.get(decision, decision)}")

    # Next steps
    if next_steps:
        lines.append("\n**Recommended next experiments:**")
        for i, step in enumerate(next_steps, 1):
            if isinstance(step, dict):
                score = step.get("score")
                score_str = f" (score: {score:.2f})" if score is not None else ""
                lines.append(f"  {i}. **{step.get('name')}**{score_str} — {step.get('description', '')}")
            else:
                lines.append(f"  {i}. {step}")

    return "\n".join(lines)
