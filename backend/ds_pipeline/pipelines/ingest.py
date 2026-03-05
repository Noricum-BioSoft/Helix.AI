"""
Data ingestion: load a CSV file into a DuckDB in-memory connection and a pandas DataFrame.

Returns lightweight metadata about the dataset (shape, column types, sample rows)
without holding large data in memory beyond what's needed.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Tuple


def ingest_csv(path: str | Path) -> Tuple[Any, Any, Dict[str, Any]]:
    """
    Load a CSV into DuckDB and pandas.

    Returns
    -------
    conn : duckdb.DuckDBPyConnection
        In-memory DuckDB connection with the data registered as the 'data' table.
    df : pandas.DataFrame
        Full dataset as a DataFrame.
    metadata : dict
        Shape, column names/types, nullability, and 5-row sample.
    """
    import duckdb
    import pandas as pd

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Data file not found: {path}")

    df = pd.read_csv(path)
    conn = duckdb.connect(database=":memory:")
    conn.register("data", df)

    col_info = []
    for col in df.columns:
        col_info.append({
            "name": col,
            "dtype": str(df[col].dtype),
            "n_missing": int(df[col].isnull().sum()),
            "pct_missing": round(float(df[col].isnull().mean()), 4),
            "n_unique": int(df[col].nunique()),
        })

    metadata: Dict[str, Any] = {
        "path": str(path),
        "n_rows": len(df),
        "n_cols": len(df.columns),
        "columns": col_info,
        "sample": df.head(5).to_dict(orient="records"),
        "size_bytes": path.stat().st_size,
    }
    return conn, df, metadata
