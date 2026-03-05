"""
Experiment planner: proposes the top 3 next experiments using a heuristic scoring model.

Score = (expected_impact * confidence) / cost

No LLM required — fully deterministic heuristics based on current run results.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

# Candidate experiment templates ranked by typical ROI
_CANDIDATE_POOL: List[Dict[str, Any]] = [
    {
        "name": "feature_engineering",
        "description": "Engineer new features: polynomial interactions, log transforms, ratio features.",
        "base_impact": 0.8,
        "base_confidence": 0.7,
        "base_cost": 0.4,
        "triggers": [],
    },
    {
        "name": "hyperparameter_tuning",
        "description": "Grid-search / cross-validate key hyperparameters of the baseline model.",
        "base_impact": 0.6,
        "base_confidence": 0.8,
        "base_cost": 0.5,
        "triggers": [],
    },
    {
        "name": "handle_class_imbalance",
        "description": "Address class imbalance via oversampling (SMOTE) or class weights.",
        "base_impact": 0.9,
        "base_confidence": 0.9,
        "base_cost": 0.3,
        "triggers": ["imbalance"],
    },
    {
        "name": "imputation_strategy",
        "description": "Try a different imputation strategy (KNN or iterative imputer) for high-missingness columns.",
        "base_impact": 0.7,
        "base_confidence": 0.7,
        "base_cost": 0.3,
        "triggers": ["high_missingness"],
    },
    {
        "name": "categorical_encoding",
        "description": "Encode categorical features (one-hot, target encoding, ordinal) for inclusion in the model.",
        "base_impact": 0.7,
        "base_confidence": 0.8,
        "base_cost": 0.3,
        "triggers": ["has_categoricals"],
    },
    {
        "name": "outlier_removal",
        "description": "Detect and remove or Winsorise outliers to reduce noise.",
        "base_impact": 0.5,
        "base_confidence": 0.6,
        "base_cost": 0.3,
        "triggers": [],
    },
    {
        "name": "ensemble_model",
        "description": "Replace the linear baseline with a gradient-boosting ensemble (e.g. RandomForest, GBM).",
        "base_impact": 0.9,
        "base_confidence": 0.7,
        "base_cost": 0.7,
        "triggers": [],
    },
    {
        "name": "temporal_split",
        "description": "Switch from random to temporal train/test split to prevent future data leakage.",
        "base_impact": 0.8,
        "base_confidence": 0.95,
        "base_cost": 0.2,
        "triggers": ["time_col"],
    },
    {
        "name": "feature_selection",
        "description": "Remove low-importance or collinear features to reduce overfitting.",
        "base_impact": 0.5,
        "base_confidence": 0.6,
        "base_cost": 0.4,
        "triggers": [],
    },
    {
        "name": "collect_more_data",
        "description": "Collect or generate more labelled samples to improve generalisation.",
        "base_impact": 1.0,
        "base_confidence": 0.5,
        "base_cost": 1.0,
        "triggers": ["small_dataset"],
    },
]


def propose_next_steps(
    run_data: Dict[str, Any],
    experiment_log: Optional[List[Dict[str, Any]]] = None,
    top_k: int = 3,
) -> List[Dict[str, Any]]:
    """
    Return the top-k next experiment proposals ordered by score.

    Score = (expected_impact * confidence) / cost
    """
    active_triggers = _detect_triggers(run_data)
    already_tried = _already_tried_names(experiment_log or [])

    scored: List[Dict[str, Any]] = []
    for cand in _CANDIDATE_POOL:
        name = cand["name"]
        if name in already_tried:
            continue

        impact = cand["base_impact"]
        confidence = cand["base_confidence"]
        cost = cand["base_cost"]

        # Boost triggered candidates
        if cand["triggers"] and any(t in active_triggers for t in cand["triggers"]):
            confidence = min(confidence + 0.15, 1.0)
            impact = min(impact + 0.1, 1.0)

        score = (impact * confidence) / max(cost, 0.01)
        scored.append({
            "name": name,
            "description": cand["description"],
            "expected_impact": round(impact, 3),
            "confidence": round(confidence, 3),
            "cost": round(cost, 3),
            "score": round(score, 3),
            "triggered_by": [t for t in cand["triggers"] if t in active_triggers],
        })

    scored.sort(key=lambda x: x["score"], reverse=True)
    return scored[:top_k]


def _detect_triggers(run_data: Dict[str, Any]) -> List[str]:
    triggers: List[str] = []

    eda = run_data.get("eda_summary") or {}
    validation = run_data.get("validation_results") or {}
    config = run_data.get("config") or {}

    # Small dataset
    shape = eda.get("shape", [])
    if shape and len(shape) >= 1 and int(shape[0]) < 500:
        triggers.append("small_dataset")

    # High missingness flagged in validation
    val_results = validation.get("results", [])
    if any(r.get("check") == "missingness" for r in val_results):
        triggers.append("high_missingness")

    # Has categorical columns
    if eda.get("n_categorical", 0) > 0:
        triggers.append("has_categoricals")

    # Time column configured
    if config.get("time_col"):
        triggers.append("time_col")

    # Class imbalance: classification task with uneven target distribution
    task_type = (run_data.get("train_info") or {}).get("task_type", "")
    if task_type == "classification":
        target_dist = eda.get("target_distribution", {})
        vc = target_dist.get("value_counts", {})
        if vc:
            counts = list(vc.values())
            if max(counts) / max(sum(counts), 1) > 0.8:
                triggers.append("imbalance")

    return triggers


def _already_tried_names(experiment_log: List[Dict[str, Any]]) -> List[str]:
    tried = []
    for row in experiment_log:
        name = row.get("objective") or row.get("changes") or ""
        tried.append(name.lower())
    return tried
