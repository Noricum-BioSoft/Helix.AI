"""Deterministic tumor-vs-normal target ranking for immunopeptidomics-style tables."""
from __future__ import annotations

import re
from typing import Any, Dict, List, Optional

_REQUIRED = ("median_tumor", "max_median_gtex")
_OPTIONAL_OUTPUT = ("gene_symbol", "pep", "cancer")

_RANK_CUES = (
    "median_tumor",
    "max_median_gtex",
    "tumor_normal",
    "tumor vs normal",
    "tumor-to-normal",
    "tumor to normal",
)


def wants_tumor_normal_ratio_ranking(
    plan: Optional[Dict[str, Any]] = None,
    command: Optional[str] = None,
) -> bool:
    """True when the user/plan asks for median_tumor vs max_median_gtex ranking."""
    parts: List[str] = []
    if command:
        parts.append(command)
    if plan:
        parts.append(str(plan.get("goal") or ""))
        parts.append(str(plan.get("title") or ""))
        for step in plan.get("steps") or []:
            if isinstance(step, dict):
                parts.append(str(step.get("operation") or ""))
                parts.append(str(step.get("description") or ""))
    blob = "\n".join(parts).lower()
    if not any(cue.replace("_", " ") in blob or cue in blob for cue in _RANK_CUES):
        return False
    return "median_tumor" in blob and "max_median_gtex" in blob


def _extract_top_k(text: str, default: int = 50) -> int:
    m = re.search(r"\btop\s+(\d+)\b", text, re.IGNORECASE)
    if m:
        try:
            return max(1, min(500, int(m.group(1))))
        except ValueError:
            pass
    return default


def run_tumor_normal_ratio_ranking(
    df: Any,
    *,
    top_k: int = 50,
    command: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Compute tumor_normal_ratio = median_tumor / max_median_gtex (safe divide) and rank descending.

    Returns a JSON-serialisable dict stored as sandbox ``result``.
    """
    import numpy as np
    import pandas as pd

    missing = [c for c in _REQUIRED if c not in df.columns]
    if missing:
        return {
            "status": "missing_required_columns",
            "missing_columns": missing,
            "available_columns": [str(c) for c in df.columns],
        }

    work = df.copy()
    work["median_tumor"] = pd.to_numeric(work["median_tumor"], errors="coerce")
    work["max_median_gtex"] = pd.to_numeric(work["max_median_gtex"], errors="coerce")

    denom = work["max_median_gtex"]
    work["tumor_normal_ratio"] = np.where(denom > 0, work["median_tumor"] / denom, np.nan)

    ranked = work.sort_values("tumor_normal_ratio", ascending=False, na_position="last")
    if command:
        top_k = _extract_top_k(command, top_k)

    out_cols = [c for c in _OPTIONAL_OUTPUT if c in ranked.columns]
    out_cols.append("tumor_normal_ratio")
    top = ranked[out_cols].head(top_k)

    records: List[Dict[str, Any]] = []
    for rec in top.to_dict(orient="records"):
        row: Dict[str, Any] = {}
        for k, v in rec.items():
            if v is None or (isinstance(v, float) and np.isnan(v)):
                row[k] = None
            elif k == "tumor_normal_ratio":
                row[k] = float(v)
            else:
                row[k] = str(v) if not isinstance(v, (int, float, bool)) else v
        records.append(row)

    return {
        "status": "ok",
        "n_rows_input": int(len(df)),
        "n_valid_ratios": int(pd.notna(work["tumor_normal_ratio"]).sum()),
        "n_zero_or_missing_normal_medians": int(
            (denom <= 0).fillna(False).sum() + denom.isna().sum()
        ),
        "top_k": top_k,
        "top_by_tumor_normal_ratio": records,
    }


def format_ranking_interpretation(result: Dict[str, Any], *, sheet: Optional[str] = None) -> str:
    """Short prose summary for deterministic ranking (no LLM)."""
    sheet_note = f" from sheet '{sheet}'" if sheet else ""
    top = result.get("top_by_tumor_normal_ratio") or []
    if not top:
        return (
            f"No ranked targets could be produced{sheet_note}. "
            "Check that median_tumor and max_median_gtex are populated."
        )
    lead = top[0]
    gene = lead.get("gene_symbol") or lead.get("pep") or "top row"
    ratio = lead.get("tumor_normal_ratio")
    intro = (
        f"Ranked {result.get('n_valid_ratios', 0)} rows by tumor-to-normal expression ratio"
        f" (median_tumor / max_median_gtex){sheet_note}."
    )
    if ratio is not None:
        intro += f" Highest ratio: {gene} (ratio={float(ratio):.4g})."
    tail = (
        f" Showing top {len(top)} candidates. Higher ratios suggest stronger tumor enrichment "
        "relative to the highest GTEx normal-tissue median. Near-zero max_median_gtex can inflate ratios."
    )
    return intro + tail
