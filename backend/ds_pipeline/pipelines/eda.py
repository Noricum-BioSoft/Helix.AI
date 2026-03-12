"""
EDA pipeline step: compute descriptive statistics and highlight key patterns.

Uses DuckDB for aggregations (scales to larger files) and pandas for summaries.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional


def run_eda(
    conn: Any,
    df: Any,
    *,
    target_col: Optional[str] = None,
    max_corr_pairs: int = 10,
) -> Dict[str, Any]:
    """
    Compute EDA summary statistics.

    Returns a dict with:
      - shape: (n_rows, n_cols)
      - dtypes: column → dtype string
      - describe: pandas describe() output as dict
      - missingness: column → pct missing
      - target_distribution: value counts / histogram bins for the target
      - top_correlations: top correlated feature pairs with target
      - numeric_cols / categorical_cols: lists of column names
    """
    import pandas as pd

    numeric_cols = df.select_dtypes(include="number").columns.tolist()
    categorical_cols = df.select_dtypes(exclude="number").columns.tolist()

    # Descriptive stats (DuckDB SUMMARIZE is fast for large tables)
    try:
        describe_dict = conn.execute("SUMMARIZE data").df().to_dict(orient="records")
    except Exception:
        describe_dict = df.describe(include="all").to_dict()

    # Missingness
    missingness: Dict[str, float] = {
        col: round(float(df[col].isnull().mean()), 4) for col in df.columns
    }

    # Target distribution
    target_dist: Dict[str, Any] = {}
    if target_col and target_col in df.columns:
        col_data = df[target_col].dropna()
        if col_data.dtype == object or col_data.nunique() <= 20:
            target_dist = {
                "type": "categorical",
                "value_counts": col_data.value_counts().to_dict(),
            }
        else:
            hist, edges = _histogram(col_data)
            target_dist = {
                "type": "numeric",
                "min": float(col_data.min()),
                "max": float(col_data.max()),
                "mean": float(col_data.mean()),
                "median": float(col_data.median()),
                "std": float(col_data.std()),
                "histogram": {"counts": hist, "edges": edges},
            }

    # Top correlations with target (numeric features only)
    top_correlations: List[Dict[str, Any]] = []
    if target_col and target_col in numeric_cols:
        feature_cols = [c for c in numeric_cols if c != target_col]
        if feature_cols:
            corr = df[feature_cols + [target_col]].corr()[target_col].drop(target_col)
            top = corr.abs().nlargest(min(max_corr_pairs, len(corr)))
            top_correlations = [
                {"feature": col, "correlation": round(float(corr[col]), 4)}
                for col in top.index
            ]

    return {
        "shape": list(df.shape),
        "dtypes": {col: str(df[col].dtype) for col in df.columns},
        "describe": describe_dict,
        "missingness": missingness,
        "target_distribution": target_dist,
        "top_correlations": top_correlations,
        "numeric_cols": numeric_cols,
        "categorical_cols": categorical_cols,
        "n_numeric": len(numeric_cols),
        "n_categorical": len(categorical_cols),
    }


def _histogram(series: Any, bins: int = 20) -> tuple[list, list]:
    """Return histogram counts and bin edges as plain Python lists."""
    try:
        import numpy as np
        counts, edges = np.histogram(series.dropna(), bins=bins)
        return counts.tolist(), [round(float(e), 6) for e in edges.tolist()]
    except Exception:
        return [], []
