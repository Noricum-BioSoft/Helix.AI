"""
Cleaning pipeline step: removes duplicates and imputes missing values.

Conservative defaults that preserve data integrity and are fully reversible.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Tuple


def clean(
    df: Any,
    conn: Any,
    *,
    drop_duplicates: bool = True,
    numeric_fill_strategy: str = "median",  # "median" | "mean" | "zero" | "none"
    categorical_fill_strategy: str = "mode",  # "mode" | "missing" | "none"
    drop_high_missingness: float = 0.9,
) -> Tuple[Any, Any, Dict[str, Any]]:
    """
    Clean the DataFrame in place.

    Returns
    -------
    conn : updated DuckDB connection (data table re-registered)
    df   : cleaned pandas DataFrame
    report : dict summarising what was changed
    """
    import pandas as pd

    original_shape = df.shape
    report: Dict[str, Any] = {
        "original_rows": original_shape[0],
        "original_cols": original_shape[1],
        "dropped_columns": [],
        "dropped_duplicate_rows": 0,
        "imputed_columns": {},
    }

    # 1. Drop columns that exceed the missingness hard stop
    cols_to_drop = [
        col for col in df.columns
        if df[col].isnull().mean() >= drop_high_missingness
    ]
    if cols_to_drop:
        df = df.drop(columns=cols_to_drop)
        report["dropped_columns"] = cols_to_drop

    # 2. Remove duplicate rows
    if drop_duplicates:
        n_before = len(df)
        df = df.drop_duplicates()
        report["dropped_duplicate_rows"] = n_before - len(df)

    # 3. Impute missing values
    for col in df.columns:
        n_missing = int(df[col].isnull().sum())
        if n_missing == 0:
            continue

        if pd.api.types.is_numeric_dtype(df[col]):
            if numeric_fill_strategy == "median":
                fill_val = df[col].median()
            elif numeric_fill_strategy == "mean":
                fill_val = df[col].mean()
            elif numeric_fill_strategy == "zero":
                fill_val = 0
            else:
                continue
            df[col] = df[col].fillna(fill_val)
            report["imputed_columns"][col] = {
                "strategy": numeric_fill_strategy,
                "n_imputed": n_missing,
                "fill_value": float(fill_val) if fill_val is not None else None,
            }
        else:
            if categorical_fill_strategy == "mode":
                mode_vals = df[col].mode()
                fill_val = mode_vals.iloc[0] if len(mode_vals) > 0 else "MISSING"
            elif categorical_fill_strategy == "missing":
                fill_val = "MISSING"
            else:
                continue
            df[col] = df[col].fillna(fill_val)
            report["imputed_columns"][col] = {
                "strategy": categorical_fill_strategy,
                "n_imputed": n_missing,
                "fill_value": str(fill_val),
            }

    report["final_rows"] = len(df)
    report["final_cols"] = len(df.columns)

    # Re-register cleaned DataFrame in DuckDB
    conn.unregister("data") if hasattr(conn, "unregister") else None
    conn.register("data", df)

    return conn, df, report
