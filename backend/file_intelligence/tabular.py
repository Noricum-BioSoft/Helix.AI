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

    result = {
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
    return _sanitize_json(result)


def _safe_scalar(val: Any) -> Any:
    """Convert a scalar to a JSON-safe Python type.

    Handles numpy integers/floats, Python float NaN/Inf, and passes through
    everything else unchanged.  NaN and Inf are mapped to None (JSON null).
    """
    import math

    if val is None:
        return None
    # Handle Python-native float NaN / Inf (e.g. from pandas operations on
    # all-null columns where the result is float('nan') not np.nan).
    if isinstance(val, float):
        return None if (math.isnan(val) or math.isinf(val)) else val
    try:
        import numpy as np
        if isinstance(val, np.integer):
            return int(val)
        if isinstance(val, np.floating):
            return None if (np.isnan(val) or np.isinf(val)) else float(val)
        if isinstance(val, np.bool_):
            return bool(val)
    except Exception:
        pass
    return val


def _sanitize_json(obj: Any) -> Any:
    """Recursively replace NaN/Inf scalars in any nested dict/list/tuple.

    This is the last-mile guard applied to the entire profiler output before
    it is handed to FastAPI's JSON serialiser.  Python's json.dumps raises a
    ValueError for float('nan'), which propagates *past* the CORS middleware
    and causes the browser to receive a 500 with no CORS headers.
    """
    if isinstance(obj, dict):
        return {_sanitize_json(k): _sanitize_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize_json(v) for v in obj]
    return _safe_scalar(obj) if isinstance(obj, (float, int)) else obj
