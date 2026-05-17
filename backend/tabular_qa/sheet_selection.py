"""Resolve Excel sheet names for tabular analysis from plans and user commands."""
from __future__ import annotations

import re
from typing import Any, Dict, List, Optional

# English words that follow "sheet" in prose but are not worksheet names.
_SHEET_STOPWORDS = frozenset(
    {
        "from",
        "the",
        "a",
        "an",
        "in",
        "on",
        "at",
        "to",
        "for",
        "of",
        "and",
        "or",
        "with",
        "using",
        "use",
        "load",
        "named",
        "called",
        "is",
        "was",
        "are",
        "be",
        "this",
        "that",
        "your",
        "my",
        "our",
        "uploaded",
        "workbook",
        "file",
        "data",
        "table",
        "named",
    }
)

# Order matters: specific patterns before generic `sheet <word>`.
_SHEET_PATTERNS = (
    re.compile(r"sheet_name\s*=\s*['\"]([\w][\w.\-]*)['\"]", re.IGNORECASE),
    re.compile(r"read_excel\([^)]*sheet_name\s*=\s*['\"]([\w][\w.\-]*)['\"]", re.IGNORECASE),
    re.compile(r"sheet\s+from\s+['\"]?([\w][\w.\-]*)['\"]?", re.IGNORECASE),
    re.compile(r"from\s+(?:the\s+)?['\"]?([\w][\w.\-]*)['\"]?\s+sheet\b", re.IGNORECASE),
    re.compile(r"on\s+(?:the\s+)?['\"]?([\w][\w.\-]*)['\"]?\s+sheet\b", re.IGNORECASE),
    re.compile(r"sheet\s+['\"]([\w][\w.\-]*)['\"]", re.IGNORECASE),
    re.compile(r"sheet\s+['\"]?([\w][\w.\-]*)['\"]?", re.IGNORECASE),
)


def _normalize_sheet_token(name: str) -> str:
    return str(name).strip()


def _is_plausible_sheet_name(name: str, available_sheets: Optional[List[str]] = None) -> bool:
    token = _normalize_sheet_token(name)
    if not token:
        return False
    if token.lower() in _SHEET_STOPWORDS:
        return False
    if available_sheets:
        if token in available_sheets:
            return True
        lower_map = {s.lower(): s for s in available_sheets}
        return token.lower() in lower_map
    # Without a workbook list, require at least one underscore or digit (ts_final, Sheet1).
    return bool(re.search(r"[_\d]", token)) or token.lower().endswith("_final")


def _match_sheet_names(text: str, available_sheets: Optional[List[str]] = None) -> List[str]:
    found: List[str] = []
    seen: set[str] = set()
    for pat in _SHEET_PATTERNS:
        for m in pat.finditer(text):
            raw = _normalize_sheet_token(m.group(1))
            if not raw or raw.lower() in seen:
                continue
            seen.add(raw.lower())
            if _is_plausible_sheet_name(raw, available_sheets):
                if available_sheets and raw not in available_sheets:
                    lower_map = {s.lower(): s for s in available_sheets}
                    raw = lower_map.get(raw.lower(), raw)
                found.append(raw)
    return found


def extract_sheet_name(
    text: Optional[str],
    *,
    available_sheets: Optional[List[str]] = None,
) -> Optional[str]:
    """Return the best Excel sheet name mentioned in *text*, if any."""
    if not text or not str(text).strip():
        return None
    matches = _match_sheet_names(str(text), available_sheets)
    return matches[0] if matches else None


def extract_sheet_from_plan(
    plan: Dict[str, Any],
    *,
    available_sheets: Optional[List[str]] = None,
) -> Optional[str]:
    """Collect sheet hints from plan fields and step operations."""
    if not isinstance(plan, dict):
        return None
    explicit = plan.get("sheet")
    if isinstance(explicit, str) and explicit.strip():
        name = _normalize_sheet_token(explicit)
        if _is_plausible_sheet_name(name, available_sheets):
            if available_sheets and name not in available_sheets:
                lower_map = {s.lower(): s for s in available_sheets}
                return lower_map.get(name.lower(), name)
            return name
    chunks: List[str] = [str(plan.get("goal") or ""), str(plan.get("title") or "")]
    for step in plan.get("steps") or []:
        if isinstance(step, dict):
            chunks.append(str(step.get("operation") or ""))
            chunks.append(str(step.get("description") or ""))
            chunks.append(str(step.get("name") or ""))
    return extract_sheet_name("\n".join(chunks), available_sheets=available_sheets)


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
    When *available_sheets* is provided, only names that exist in the workbook are returned.
    """
    if available_sheets is None and profile:
        available_sheets = profile.get("available_sheets")

    candidates: List[Optional[str]] = []
    if plan:
        if isinstance(plan.get("sheet"), str):
            candidates.append(plan.get("sheet"))
        candidates.append(extract_sheet_from_plan(plan, available_sheets=available_sheets))
        # Scan full plan JSON for any valid sheet mention (planner may embed ts_final in steps).
        try:
            import json

            for name in _match_sheet_names(json.dumps(plan), available_sheets):
                candidates.append(name)
        except (TypeError, ValueError):
            pass
    if command:
        candidates.append(extract_sheet_name(command, available_sheets=available_sheets))
        for name in _match_sheet_names(command, available_sheets):
            candidates.append(name)
    if profile:
        summary = profile.get("summary") or {}
        candidates.append(summary.get("source_sheet"))

    for c in candidates:
        if not c:
            continue
        name = _normalize_sheet_token(c)
        if not _is_plausible_sheet_name(name, available_sheets):
            continue
        if available_sheets:
            if name in available_sheets:
                return name
            lower_map = {s.lower(): s for s in available_sheets}
            if name.lower() in lower_map:
                return lower_map[name.lower()]
            continue
        return name

    # No valid named sheet — fall back to first profile sheet only if nothing else matched.
    if profile:
        summary = profile.get("summary") or {}
        default = summary.get("source_sheet")
        if default and (not available_sheets or default in available_sheets):
            return default
    return None


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
