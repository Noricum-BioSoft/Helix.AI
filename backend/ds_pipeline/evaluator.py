"""
Evaluator: wraps the evaluate pipeline step and determines a run decision.

Decision logic:
  - "new_best"   — current run beats the prior best on the primary metric
  - "no_change"  — metrics within 1% of prior best
  - "regressed"  — metrics are worse than prior best
  - "first_run"  — no prior best exists
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from backend.ds_pipeline.pipelines.evaluate import compare_with_prior_best


def evaluate_run(
    metrics: Dict[str, Any],
    experiment_log: List[Dict[str, Any]],
    primary_metric: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Compare current metrics against the prior best from experiment_log.

    Returns a dict with 'comparison', 'decision', and 'prior_best_run_id'.
    """
    prior_best_metrics, prior_best_run_id = _find_prior_best(
        experiment_log, primary_metric
    )

    comparison = compare_with_prior_best(metrics, prior_best_metrics, primary_metric)
    pm = comparison.get("primary_metric")

    if prior_best_metrics is None:
        decision = "first_run"
    elif comparison.get("is_new_best"):
        decision = "new_best"
    else:
        delta = (comparison.get("deltas") or {}).get(pm, 0) if pm else 0
        try:
            rel_change = abs(float(delta)) / max(abs(float(prior_best_metrics.get(pm, 1))), 1e-9)
        except (TypeError, ValueError):
            rel_change = 0.0
        decision = "no_change" if rel_change < 0.01 else "regressed"

    return {
        "comparison": comparison,
        "decision": decision,
        "prior_best_run_id": prior_best_run_id,
    }


def _find_prior_best(
    experiment_log: List[Dict[str, Any]],
    primary_metric: Optional[str],
) -> tuple[Optional[Dict[str, Any]], Optional[str]]:
    if not experiment_log:
        return None, None

    best_row = None
    best_val: Optional[float] = None
    lower_is_better = primary_metric in {"mse", "rmse", "mae"} if primary_metric else False

    for row in experiment_log:
        metric_key = f"metric_{primary_metric}" if primary_metric else None
        if metric_key and metric_key in row:
            try:
                val = float(row[metric_key])
                if best_val is None:
                    best_val = val
                    best_row = row
                elif lower_is_better and val < best_val:
                    best_val = val
                    best_row = row
                elif not lower_is_better and val > best_val:
                    best_val = val
                    best_row = row
            except (TypeError, ValueError):
                continue

    if best_row is None:
        best_row = experiment_log[-1]

    prior_metrics: Dict[str, Any] = {}
    for k, v in best_row.items():
        if k.startswith("metric_"):
            metric_name = k[len("metric_"):]
            try:
                prior_metrics[metric_name] = float(v)
            except (TypeError, ValueError):
                pass

    return (prior_metrics or None), best_row.get("run_id")
