"""Unit tests for the universal tabular ingest module."""
from __future__ import annotations

import io
import tempfile
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Helpers to create tiny in-memory fixture files on disk
# ---------------------------------------------------------------------------

def _write_csv(dir_: Path, name: str = "data.csv") -> Path:
    p = dir_ / name
    p.write_text("gene,tumor,normal\nGENE_A,10.0,2.5\nGENE_B,5.0,5.0\nGENE_C,8.0,1.0\n")
    return p


def _write_tsv(dir_: Path, name: str = "data.tsv") -> Path:
    p = dir_ / name
    p.write_text("gene\ttumor\tnormal\nGENE_A\t10.0\t2.5\nGENE_B\t5.0\t5.0\n")
    return p


def _write_xlsx(dir_: Path, name: str = "data.xlsx") -> Path:
    import openpyxl

    wb = openpyxl.Workbook()
    ws1 = wb.active
    ws1.title = "Sheet1"
    ws1.append(["gene", "tumor", "normal"])
    ws1.append(["GENE_A", 10.0, 2.5])
    ws1.append(["GENE_B", 5.0, 5.0])
    ws2 = wb.create_sheet("ts_final")
    ws2.append(["gene", "median_tumor", "max_median_gtex"])
    ws2.append(["TARGET_1", 8.5, 1.2])
    ws2.append(["TARGET_2", 6.0, 3.0])
    p = dir_ / name
    wb.save(str(p))
    return p


# ---------------------------------------------------------------------------
# ingest_tabular — CSV
# ---------------------------------------------------------------------------

class TestIngestTabularCsv:
    def test_returns_duckdb_conn_df_metadata(self, tmp_path):
        from backend.ds_pipeline.pipelines.ingest import ingest_tabular

        path = _write_csv(tmp_path)
        conn, df, meta = ingest_tabular(path)

        assert df.shape == (3, 3)
        assert list(df.columns) == ["gene", "tumor", "normal"]
        assert meta["n_rows"] == 3
        assert meta["n_cols"] == 3
        assert meta["source_format"] == "csv"
        assert meta["source_sheet"] is None

    def test_duckdb_query_works(self, tmp_path):
        from backend.ds_pipeline.pipelines.ingest import ingest_tabular

        path = _write_csv(tmp_path)
        conn, _, _ = ingest_tabular(path)
        result = conn.execute("SELECT gene FROM data ORDER BY tumor DESC").fetchall()
        assert result[0][0] == "GENE_A"

    def test_file_not_found_raises(self, tmp_path):
        from backend.ds_pipeline.pipelines.ingest import ingest_tabular

        with pytest.raises(FileNotFoundError):
            ingest_tabular(tmp_path / "missing.csv")


# ---------------------------------------------------------------------------
# ingest_tabular — TSV
# ---------------------------------------------------------------------------

class TestIngestTabularTsv:
    def test_tsv_loads_correctly(self, tmp_path):
        from backend.ds_pipeline.pipelines.ingest import ingest_tabular

        path = _write_tsv(tmp_path)
        conn, df, meta = ingest_tabular(path)

        assert df.shape == (2, 3)
        assert meta["source_format"] == "tsv"

    def test_explicit_delimiter_override(self, tmp_path):
        from backend.ds_pipeline.pipelines.ingest import ingest_tabular

        p = tmp_path / "pipe.txt"
        p.write_text("gene|tumor\nGENE_A|10.0\n")
        _, df, _ = ingest_tabular(p, delimiter="|")
        assert list(df.columns) == ["gene", "tumor"]


# ---------------------------------------------------------------------------
# ingest_tabular — Excel
# ---------------------------------------------------------------------------

class TestIngestTabularExcel:
    def test_loads_first_sheet_by_default(self, tmp_path):
        from backend.ds_pipeline.pipelines.ingest import ingest_tabular

        path = _write_xlsx(tmp_path)
        _, df, meta = ingest_tabular(path)

        assert meta["source_format"] == "excel"
        assert meta["source_sheet"] == "Sheet1"
        assert list(df.columns) == ["gene", "tumor", "normal"]

    def test_loads_named_sheet(self, tmp_path):
        from backend.ds_pipeline.pipelines.ingest import ingest_tabular

        path = _write_xlsx(tmp_path)
        _, df, meta = ingest_tabular(path, sheet="ts_final")

        assert meta["source_sheet"] == "ts_final"
        assert "median_tumor" in df.columns

    def test_invalid_sheet_raises(self, tmp_path):
        from backend.ds_pipeline.pipelines.ingest import ingest_tabular

        path = _write_xlsx(tmp_path)
        with pytest.raises(ValueError, match="not found"):
            ingest_tabular(path, sheet="NoSuchSheet")

    def test_available_sheets_in_metadata(self, tmp_path):
        from backend.ds_pipeline.pipelines.ingest import ingest_tabular

        path = _write_xlsx(tmp_path)
        _, _, meta = ingest_tabular(path)
        assert set(meta["available_sheets"]) == {"Sheet1", "ts_final"}


# ---------------------------------------------------------------------------
# list_sheets
# ---------------------------------------------------------------------------

class TestListSheets:
    def test_excel_returns_sheet_names(self, tmp_path):
        from backend.ds_pipeline.pipelines.ingest import list_sheets

        path = _write_xlsx(tmp_path)
        sheets = list_sheets(path)
        assert "ts_final" in sheets

    def test_csv_returns_default(self, tmp_path):
        from backend.ds_pipeline.pipelines.ingest import list_sheets

        path = _write_csv(tmp_path)
        assert list_sheets(path) == ["default"]


# ---------------------------------------------------------------------------
# ingest_csv shim (backwards compatibility)
# ---------------------------------------------------------------------------

class TestIngestCsvShim:
    def test_shim_delegates_to_ingest_tabular(self, tmp_path):
        from backend.ds_pipeline.pipelines.ingest import ingest_csv, ingest_tabular

        path = _write_csv(tmp_path)
        _, df_shim, _ = ingest_csv(path)
        _, df_direct, _ = ingest_tabular(path)
        assert df_shim.equals(df_direct)


# ---------------------------------------------------------------------------
# Command router — tabular_analysis route
# ---------------------------------------------------------------------------

class TestCommandRouterTabularAnalysis:
    def _router(self):
        from backend.command_router import CommandRouter
        return CommandRouter()

    def _ctx(self, files=None):
        return {"session_id": "test-session", "uploaded_files": files or []}

    def test_rank_command_with_xlsx_file_routes_tabular(self):
        router = self._router()
        cmd = "rank genes by median_tumor/max_median_gtex from ts_final sheet in data.xlsx"
        tool, params = router.route_command(cmd, self._ctx())
        assert tool == "tabular_analysis"
        assert "xlsx" in params.get("data_path", "").lower() or params.get("data_path") != ""

    def test_sort_command_with_tsv_routes_tabular(self):
        router = self._router()
        cmd = "sort data.tsv by tumor expression descending"
        tool, params = router.route_command(cmd, self._ctx())
        assert tool == "tabular_analysis"

    def test_ratio_command_with_session_xlsx_routes_tabular(self):
        router = self._router()
        cmd = "calculate ratio of tumor vs normal expression"
        ctx = self._ctx(files=[{"filename": "experiment.xlsx", "path": "/uploads/experiment.xlsx"}])
        tool, params = router.route_command(cmd, ctx)
        assert tool == "tabular_analysis"

    def test_sheet_name_extracted(self):
        router = self._router()
        cmd = "rank genes from sheet 'ts_final' in results.xlsx"
        tool, params = router.route_command(cmd, self._ctx())
        assert tool == "tabular_analysis"
        assert params.get("sheet") == "ts_final"

    def test_plain_sequence_command_does_not_route_tabular(self):
        router = self._router()
        cmd = "align these protein sequences and show me the alignment"
        tool, _ = router.route_command(cmd, self._ctx())
        assert tool != "tabular_analysis"

    def test_tabular_route_does_not_fire_without_file_or_op(self):
        router = self._router()
        cmd = "analyze my data"
        tool, _ = router.route_command(cmd, self._ctx())
        assert tool != "tabular_analysis"
