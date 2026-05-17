"""Resolve Excel sheet names for tabular analysis from plans and user commands."""
from __future__ import annotations

import re
from typing import Any, Dict, List, Optional

# sheet 'ts_final' | sheet ts_final | sheet_name='ts_final' | read_excel(..., sheet_name="ts_final")
_SHEET_PATTERNS = (
    re.compile(r"sheet\s+['\"]?([\w][\w.\-]*)['\"]?", re.IGNORECASE),
    re.compile(r"sheet_name\s*=\s*['\"]([\w][\w.\-]*)['\"]", re.IGNORECASE),
    re.compile(r"read_excel\([^)]*sheet_name\s*=\s*['\"]([\w][\w.\-]*)['\"]", re.IGNORECASE),
    re.compile(r"from\s+(?:the\s+)?['\"]?([\w][\w.\-]*)['\"]?\s+sheet\b", re.IGNORECASE),
)


def extract_sheet_name(text: Optional[str]) -> Optional[str]:
    """Return the first Excel sheet name mentioned in *text*, if any."""
    if not text or not str(text).strip():
        return None
    blob = str(text)
    for pat in _SHEET_PATTERNS:
        m = pat.search(blob)
        if m:
            return m.group(1)
    return None


def extract_sheet_from_plan(plan: Dict[str, Any]) -> Optional[str]:
    """Collect sheet hints from plan fields and step operations."""
    if not isinstance(plan, dict):
        return None
    explicit = plan.get("sheet")
    if isinstance(explicit, str) and explicit.strip():
        return explicit.strip()
    chunks: List[str] = [str(plan.get("goal") or ""), str(plan.get("title") or "")]
    for step in plan.get("steps") or []:
        if isinstance(step, dict):
            chunks.append(str(step.get("operation") or ""))
            chunks.append(str(step.get("description") or ""))
            chunks.append(str(step.get("name") or ""))
    return extract_sheet_name("\n".join(chunks))


def resolve_execution_sheet(
    *,
    plan: Optional[Dict[str, Any]] = None,
    command: Optional[str] = None,
    profile: Optional[Dict[str, Any]] = None,
    available_sheets: Optional[List[str]] = None,
) -> Optional[str]:
    """
    Choose the worksheet to load for execution.

    Priority: explicit plan.sheet → plan step text → user command → upload profile default.
    When *available_sheets* is provided, the resolved name must exist in the workbook.
    """
    candidates: List[Optional[str]] = []
    if plan:
        candidates.append(plan.get("sheet") if isinstance(plan.get("sheet"), str) else None)
        candidates.append(extract_sheet_from_plan(plan))
    candidates.append(extract_sheet_name(command))
    if profile:
        summary = profile.get("summary") or {}
        candidates.append(summary.get("source_sheet"))

    chosen: Optional[str] = None
    for c in candidates:
        if c and str(c).strip():
            chosen = str(c).strip()
            break

    if chosen and available_sheets:
        if chosen not in available_sheets:
            # Case-insensitive fallback
            lower_map = {s.lower(): s for s in available_sheets}
            chosen = lower_map.get(chosen.lower(), chosen)
    return chosen


def refresh_schema_preview_for_sheet(local_path: str, sheet: str) -> Dict[str, Any]:
    """Re-profile a tabular file for a specific Excel sheet."""
    try:
        from backend.file_intelligence.tabular import profile_tabular

        return profile_tabular(local_path, sheet=sheet)
    except Exception:
        return _minimal_schema_preview(local_path, sheet)


def _minimal_schema_preview(local_path: str, sheet: str) -> Dict[str, Any]:
    """Lightweight schema preview when full profiler (DuckDB) is unavailable."""
    import pandas as pd
    from backend.ds_pipeline.pipelines.ingest import list_sheets

    df = pd.read_excel(local_path, sheet_name=sheet)
    columns = [
        {
            "name": str(c),
            "dtype": str(df[c].dtype),
            "n_missing": int(df[c].isna().sum()),
            "pct_missing": round(float(df[c].isna().mean()), 4),
            "n_unique": int(df[c].nunique(dropna=True)),
        }
        for c in df.columns
    ]
    return {
        "format": "excel",
        "family": "tabular",
        "n_records": len(df),
        "summary": {
            "n_rows": len(df),
            "n_cols": len(df.columns),
            "source_sheet": sheet,
        },
        "schema": {"columns": columns},
        "sample": df.head(3).to_dict(orient="records"),
        "available_sheets": list_sheets(local_path),
        "raw_metadata": {},
        "profiler_error": None,
    }


def unwrap_tabular_result(serialized: Any) -> Any:
    """Unwrap executor serialisation envelope {type, value}."""
    if isinstance(serialized, dict) and serialized.get("type") == "dict" and "value" in serialized:
        return serialized["value"]
    return serialized


def tabular_result_indicates_failure(serialized: Any) -> Optional[str]:
    """
    Return a user-facing error message when sandbox ``result`` signals logical failure.
    """
    inner = unwrap_tabular_result(serialized)
    if not isinstance(inner, dict):
        return None
    status = str(inner.get("status") or "").lower()
    if status == "missing_required_columns":
        missing = inner.get("missing_columns") or []
        avail = inner.get("available_columns") or []
        return (
            "Required columns are missing for this analysis: "
            f"{missing}. Loaded columns include: {avail[:12]}"
            + ("…" if len(avail) > 12 else "")
        )
    if status in {"error", "failed"}:
        return str(inner.get("message") or inner.get("error") or "Tabular analysis failed.")
    return None
