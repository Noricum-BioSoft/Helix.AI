"""
Baseline model training pipeline step.

Auto-detects task type (classification vs regression) from the target column
and trains a minimal sklearn model:
  - Classification → LogisticRegression (max_iter=1000, solver="lbfgs")
  - Regression     → Ridge (alpha=1.0)

Only numeric columns are used as features; object columns are dropped.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Tuple


def _detect_task(series: Any) -> str:
    """Return 'classification' or 'regression' based on target column heuristics."""
    if series.dtype == object:
        return "classification"
    n_unique = series.nunique()
    if n_unique <= 10 and n_unique / max(len(series), 1) < 0.05:
        return "classification"
    return "regression"


def train_baseline(
    df: Any,
    *,
    target_col: str,
    task_type: str = "auto",
    test_size: float = 0.2,
    random_seed: int = 42,
    model_kwargs: Optional[Dict[str, Any]] = None,
) -> Tuple[Any, Any, Any, Dict[str, Any]]:
    """
    Fit a baseline model.

    Returns
    -------
    model      : fitted sklearn estimator
    X_train    : training features DataFrame
    X_test     : test features DataFrame
    train_info : dict with task_type, feature_cols, train/test sizes, model params
    """
    import pandas as pd
    from sklearn.model_selection import train_test_split
    from sklearn.linear_model import LogisticRegression, Ridge
    from sklearn.preprocessing import LabelEncoder

    if target_col not in df.columns:
        raise ValueError(f"Target column '{target_col}' not found in dataset.")

    # Feature selection: numeric columns only, excluding target
    feature_df = df.select_dtypes(include="number").drop(
        columns=[target_col], errors="ignore"
    )
    feature_cols = feature_df.columns.tolist()

    if not feature_cols:
        raise ValueError("No numeric feature columns available for training.")

    X = feature_df.values
    y_raw = df[target_col]

    # Auto-detect or validate task type
    if task_type == "auto":
        task_type = _detect_task(y_raw)

    if task_type == "classification":
        le = LabelEncoder()
        y = le.fit_transform(y_raw.astype(str))
        classes = le.classes_.tolist()
    else:
        y = y_raw.values.astype(float)
        le = None
        classes = None

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_seed
    )

    kwargs = model_kwargs or {}
    if task_type == "classification":
        model = LogisticRegression(max_iter=1000, solver="lbfgs", **kwargs)
    else:
        model = Ridge(alpha=kwargs.get("alpha", 1.0))

    model.fit(X_train, y_train)

    train_info: Dict[str, Any] = {
        "task_type": task_type,
        "feature_cols": feature_cols,
        "n_features": len(feature_cols),
        "n_train": len(X_train),
        "n_test": len(X_test),
        "model_class": type(model).__name__,
        "model_params": model.get_params(),
        "classes": classes,
        "label_encoder": le,
        "y_train": y_train,
        "y_test": y_test,
        "X_train": X_train,
        "X_test": X_test,
    }
    return model, X_train, X_test, train_info
