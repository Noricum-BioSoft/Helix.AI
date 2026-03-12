"""
Evaluation pipeline step: compute metrics and compare against the prior best run.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional


def evaluate_model(
    model: Any,
    train_info: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Compute evaluation metrics from the fitted model and train_info produced by train.py.

    Returns a metrics dict suitable for storing in run_data['metrics'].
    """
    task_type = train_info.get("task_type", "regression")
    X_test = train_info.get("X_test")
    y_test = train_info.get("y_test")

    if X_test is None or y_test is None:
        return {"error": "Missing X_test / y_test in train_info."}

    y_pred = model.predict(X_test)
    metrics: Dict[str, Any] = {}

    if task_type == "classification":
        metrics = _classification_metrics(y_test, y_pred, model, X_test, train_info)
    else:
        metrics = _regression_metrics(y_test, y_pred)

    return metrics


def _classification_metrics(
    y_test: Any,
    y_pred: Any,
    model: Any,
    X_test: Any,
    train_info: Dict[str, Any],
) -> Dict[str, Any]:
    from sklearn.metrics import (
        accuracy_score, f1_score, precision_score, recall_score, roc_auc_score,
    )
    import numpy as np

    classes = train_info.get("classes") or []
    binary = len(classes) == 2
    avg = "binary" if binary else "weighted"

    metrics: Dict[str, Any] = {
        "accuracy": round(float(accuracy_score(y_test, y_pred)), 4),
        "f1": round(float(f1_score(y_test, y_pred, average=avg, zero_division=0)), 4),
        "precision": round(float(precision_score(y_test, y_pred, average=avg, zero_division=0)), 4),
        "recall": round(float(recall_score(y_test, y_pred, average=avg, zero_division=0)), 4),
    }

    if binary and hasattr(model, "predict_proba"):
        try:
            y_prob = model.predict_proba(X_test)[:, 1]
            metrics["roc_auc"] = round(float(roc_auc_score(y_test, y_prob)), 4)
        except Exception:
            pass

    return metrics


def _regression_metrics(y_test: Any, y_pred: Any) -> Dict[str, Any]:
    from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
    import numpy as np

    mse = float(mean_squared_error(y_test, y_pred))
    return {
        "mse": round(mse, 6),
        "rmse": round(float(np.sqrt(mse)), 6),
        "mae": round(float(mean_absolute_error(y_test, y_pred)), 6),
        "r2": round(float(r2_score(y_test, y_pred)), 4),
    }


def compare_with_prior_best(
    current_metrics: Dict[str, Any],
    prior_best_metrics: Optional[Dict[str, Any]],
    primary_metric: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Compare current metrics against the prior best.

    Returns a comparison dict with deltas and whether this run is a new best.
    """
    if not prior_best_metrics:
        return {
            "is_new_best": True,
            "prior_best": None,
            "deltas": {},
            "primary_metric": primary_metric,
        }

    deltas: Dict[str, Any] = {}
    for key in current_metrics:
        if key in prior_best_metrics:
            try:
                delta = float(current_metrics[key]) - float(prior_best_metrics[key])
                deltas[key] = round(delta, 6)
            except (TypeError, ValueError):
                pass

    # Determine primary metric automatically if not specified
    if not primary_metric:
        primary_metric = _pick_primary_metric(current_metrics)

    is_new_best = False
    if primary_metric and primary_metric in current_metrics and primary_metric in prior_best_metrics:
        try:
            # Higher is better for accuracy/f1/r2/roc_auc; lower is better for mse/rmse/mae
            lower_is_better = primary_metric in {"mse", "rmse", "mae"}
            curr_val = float(current_metrics[primary_metric])
            prev_val = float(prior_best_metrics[primary_metric])
            is_new_best = (curr_val < prev_val) if lower_is_better else (curr_val > prev_val)
        except (TypeError, ValueError):
            pass

    return {
        "is_new_best": is_new_best,
        "prior_best": prior_best_metrics,
        "deltas": deltas,
        "primary_metric": primary_metric,
    }


def _pick_primary_metric(metrics: Dict[str, Any]) -> Optional[str]:
    preferred_order = ["roc_auc", "f1", "accuracy", "r2", "rmse", "mse", "mae"]
    for m in preferred_order:
        if m in metrics:
            return m
    keys = [k for k in metrics if isinstance(metrics.get(k), (int, float))]
    return keys[0] if keys else None
