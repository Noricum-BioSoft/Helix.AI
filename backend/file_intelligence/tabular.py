"""Upload-time profiler for tabular files (CSV / TSV / Excel)."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional


def profile_tabular(path: str | Path, *, sheet: Optional[str] = None) -> Dict[str, Any]:
    """
    Return a rich FileProfile for a tabular file.

    Calls ``ingest_tabular`` for loading, then computes per-column statistics
    (min/max/mean for numerics, value_counts for low-cardinality categoricals).
    """
    import pandas as pd
    from backend.ds_pipeline.pipelines.ingest import ingest_tabular

    path = Path(path)

    try:
        _conn, df, meta = ingest_tabular(path, sheet=sheet)
    except Exception as exc:
        return {
            "format": path.suffix.lstrip("."),
            "family": "tabular",
            "n_records": None,
            "summary": {},
            "schema": {},
            "sample": [],
            "available_sheets": None,
            "raw_metadata": {},
            "profiler_error": str(exc),
        }

    col_profiles = []
    for col_info in meta["columns"]:
        col = col_info["name"]
        series = df[col]
        cp: Dict[str, Any] = {
            "name": col,
            "dtype": col_info["dtype"],
            "n_missing": col_info["n_missing"],
            "pct_missing": col_info["pct_missing"],
            "n_unique": col_info["n_unique"],
        }
        if pd.api.types.is_numeric_dtype(series):
            described = series.describe()
            cp.update({
                "min": _safe_scalar(described.get("min")),
                "max": _safe_scalar(described.get("max")),
                "mean": _safe_scalar(described.get("mean")),
                "std": _safe_scalar(described.get("std")),
                "median": _safe_scalar(series.median()),
            })
        elif col_info["n_unique"] <= 20:
            cp["top_values"] = series.value_counts().head(10).to_dict()
        col_profiles.append(cp)

    return {
        "format": meta["source_format"],
        "family": "tabular",
        "n_records": meta["n_rows"],
        "summary": {
            "n_rows": meta["n_rows"],
            "n_cols": meta["n_cols"],
            "size_bytes": meta["size_bytes"],
            "source_sheet": meta.get("source_sheet"),
        },
        "schema": {"columns": col_profiles},
        "sample": meta["sample"],
        "available_sheets": meta.get("available_sheets"),
        "raw_metadata": {},
        "profiler_error": None,
    }


def _safe_scalar(val: Any) -> Any:
    """Convert numpy scalar to Python native type for JSON serialisation."""
    try:
        import numpy as np
        if isinstance(val, (np.integer,)):
            return int(val)
        if isinstance(val, (np.floating,)):
            return None if np.isnan(val) else float(val)
    except Exception:
        pass
    return val
