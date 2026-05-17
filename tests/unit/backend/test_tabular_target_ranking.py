"""Tests for sheet-aware tabular execution and tumor/normal target ranking."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict
from unittest.mock import MagicMock

import pytest

from backend.tabular_qa.analysis_executor import execute_analysis_plan
from backend.tabular_qa.sheet_selection import (
    extract_sheet_name,
    resolve_execution_sheet,
    tabular_result_indicates_failure,
)
from backend.tabular_qa.target_ranking import (
    run_tumor_normal_ratio_ranking,
    wants_tumor_normal_ratio_ranking,
)


def _write_media5_style_xlsx(dir_: Path) -> Path:
    import openpyxl

    wb = openpyxl.Workbook()
    ws0 = wb.active
    ws0.title = "gene_morpheus_mem"
    ws0.append(["gene", "cancer", "normal"])
    ws0.append(["GENE_A", 1.0, 0.1])
    ts = wb.create_sheet("ts_final")
    ts.append(["gene_symbol", "pep", "cancer", "median_tumor", "max_median_gtex"])
    ts.append(["G1", "PEP1", "BRCA", 10.0, 1.0])
    ts.append(["G2", "PEP2", "SKCM", 5.0, 5.0])
    ts.append(["G3", "PEP3", "BRCA", 20.0, 0.5])
    p = dir_ / "media5.xlsx"
    wb.save(str(p))
    return p


def test_extract_sheet_name_from_command():
    cmd = "Using sheet ts_final compute median_tumor / max_median_gtex and rank descending"
    assert extract_sheet_name(cmd) == "ts_final"


def test_extract_sheet_name_sheet_from_phrase():
    """Regression: 'sheet from ts_final' must not capture stopword 'from'."""
    avail = ["gene_morpheus_mem", "ts_final"]
    cmd = "Load data from sheet ts_final and rank median_tumor / max_median_gtex"
    assert extract_sheet_name(cmd, available_sheets=avail) == "ts_final"
    assert resolve_execution_sheet(
        command="Using the uploaded workbook, load sheet from ts_final",
        available_sheets=avail,
    ) == "ts_final"


def test_resolve_execution_sheet_prefers_plan_over_profile():
    plan = {"sheet": "ts_final", "goal": "rank targets"}
    profile = {"summary": {"source_sheet": "gene_morpheus_mem"}}
    assert resolve_execution_sheet(plan=plan, profile=profile) == "ts_final"


def test_run_tumor_normal_ratio_ranking_orders_descending():
    import pandas as pd

    df = pd.DataFrame(
        {
            "gene_symbol": ["G1", "G2", "G3"],
            "pep": ["A", "B", "C"],
            "median_tumor": [10.0, 5.0, 20.0],
            "max_median_gtex": [1.0, 5.0, 0.5],
        }
    )
    out = run_tumor_normal_ratio_ranking(df, top_k=3)
    assert out["status"] == "ok"
    top = out["top_by_tumor_normal_ratio"]
    assert top[0]["gene_symbol"] == "G3"
    assert top[0]["tumor_normal_ratio"] == pytest.approx(40.0)


def test_tabular_result_indicates_failure_on_missing_columns():
    msg = tabular_result_indicates_failure(
        {
            "type": "dict",
            "value": {
                "status": "missing_required_columns",
                "missing_columns": ["median_tumor"],
                "available_columns": ["cancer"],
            },
        }
    )
    assert msg is not None
    assert "median_tumor" in msg


def test_execute_analysis_plan_uses_ts_final_deterministic_path(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Executor must load ts_final (not first sheet) and return ranked targets without LLM."""
    xlsx = _write_media5_style_xlsx(tmp_path)
    profile = {
        "family": "tabular",
        "summary": {"source_sheet": "gene_morpheus_mem", "n_rows": 1, "n_cols": 3},
        "available_sheets": ["gene_morpheus_mem", "ts_final"],
        "schema": {"columns": [{"name": "gene", "dtype": "object"}]},
    }
    plan: Dict[str, Any] = {
        "type": "tabular_analysis",
        "title": "Tumor Normal Ratio Ranking",
        "goal": "Rank by median_tumor / max_median_gtex on sheet ts_final",
        "sheet": "ts_final",
        "steps": [],
        "_files": [
            {
                "name": "media5.xlsx",
                "local_path": str(xlsx),
                "schema_preview": profile,
            }
        ],
    }
    ctx = {"uploaded_files": plan["_files"]}
    cmd = (
        "On sheet ts_final compute median_tumor / max_median_gtex, "
        "rank descending, top 10 cancer targets"
    )

    def _no_llm():
        raise AssertionError("LLM should not be called for deterministic ranking")

    monkeypatch.setattr(
        "backend.tabular_qa.analysis_executor._get_llm",
        _no_llm,
    )

    out = execute_analysis_plan(plan, "sess-rank", ctx, original_command=cmd)
    assert out["status"] == "success"
    assert out.get("execution_sheet") == "ts_final"
    assert wants_tumor_normal_ratio_ranking(plan, cmd)
    payload = out["result_data"]["value"]
    assert payload["status"] == "ok"
    assert payload["top_by_tumor_normal_ratio"][0]["gene_symbol"] == "G3"
