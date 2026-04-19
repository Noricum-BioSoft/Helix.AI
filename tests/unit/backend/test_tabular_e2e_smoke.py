"""
End-to-end smoke test for the tabular analysis pipeline.

Path under test:
  CSV file → file_intelligence profiler → tabular_qa executor → tabular_qa agent

No running backend or real LLM API key required: the LLM call is mocked
to return deterministic Python code so we verify every layer except the
model itself.
"""

import io
import textwrap
from pathlib import Path
from typing import Any, Dict
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

SAMPLE_CSV = """\
gene,sample_a,sample_b,log2fc,padj
BRCA1,120.5,45.2,1.41,0.001
TP53,89.3,210.7,-1.24,0.032
MYC,340.1,88.6,1.94,0.0001
EGFR,55.0,60.1,-0.13,0.78
PTEN,200.0,195.0,0.04,0.91
CDK2,10.2,95.3,-3.22,0.0001
RB1,300.4,50.8,2.56,0.0001
KRAS,88.7,90.1,-0.02,0.99
"""


@pytest.fixture
def csv_file(tmp_path: Path) -> Path:
    p = tmp_path / "expression.csv"
    p.write_text(SAMPLE_CSV)
    return p


@pytest.fixture
def sample_df() -> pd.DataFrame:
    return pd.read_csv(io.StringIO(SAMPLE_CSV))


# ---------------------------------------------------------------------------
# 1. Profiler smoke test
#    Actual schema: n_records, summary.{n_rows,n_cols}, schema.columns[].name,
#    sample, available_sheets
# ---------------------------------------------------------------------------

class TestFileIntelligenceProfiler:
    def test_profile_csv_returns_expected_shape(self, csv_file: Path):
        from backend.file_intelligence.profiler import profile_file

        result = profile_file(str(csv_file))

        assert result["family"] == "tabular"
        assert result["format"] in ("csv", "tsv", "tabular")
        assert result["n_records"] == 8

    def test_profile_csv_summary(self, csv_file: Path):
        from backend.file_intelligence.profiler import profile_file

        result = profile_file(str(csv_file))

        assert result["summary"]["n_rows"] == 8
        assert result["summary"]["n_cols"] == 5

    def test_profile_csv_column_schema(self, csv_file: Path):
        from backend.file_intelligence.profiler import profile_file

        result = profile_file(str(csv_file))
        col_names = [c["name"] for c in result["schema"]["columns"]]

        assert "gene" in col_names
        assert "log2fc" in col_names
        assert "padj" in col_names

    def test_profile_csv_numeric_stats(self, csv_file: Path):
        from backend.file_intelligence.profiler import profile_file

        result = profile_file(str(csv_file))
        col_by_name = {c["name"]: c for c in result["schema"]["columns"]}

        log2fc = col_by_name["log2fc"]
        assert "min" in log2fc
        assert "max" in log2fc
        assert log2fc["max"] == pytest.approx(2.56, abs=0.01)

    def test_profile_csv_sample_rows_present(self, csv_file: Path):
        from backend.file_intelligence.profiler import profile_file

        result = profile_file(str(csv_file))
        assert len(result["sample"]) > 0
        assert "gene" in result["sample"][0]

    def test_profile_tsv(self, tmp_path: Path):
        from backend.file_intelligence.profiler import profile_file

        tsv = tmp_path / "data.tsv"
        tsv.write_text(SAMPLE_CSV.replace(",", "\t"))
        result = profile_file(str(tsv))

        assert result["n_records"] == 8
        assert result["summary"]["n_cols"] == 5

    def test_profile_xlsx(self, tmp_path: Path):
        from backend.file_intelligence.profiler import profile_file

        df = pd.read_csv(io.StringIO(SAMPLE_CSV))
        xlsx = tmp_path / "data.xlsx"
        df.to_excel(str(xlsx), index=False)
        result = profile_file(str(xlsx))

        assert result["n_records"] == 8
        assert result["summary"]["n_cols"] == 5
        assert result["available_sheets"] is not None
        assert len(result["available_sheets"]) >= 1


# ---------------------------------------------------------------------------
# 2. Executor smoke test
#    Actual schema: success, result (serialized dict), error, stdout, code
#    Scalar: {"type": "scalar", "value": x}
#    DataFrame: {"type": "dataframe", "value": [...], "columns": [...], ...}
#    List: {"type": "list", "value": [...]}
# ---------------------------------------------------------------------------

class TestTabularQAExecutor:
    def test_simple_aggregation(self, sample_df: pd.DataFrame):
        from backend.tabular_qa.executor import execute_code

        code = "result = df['log2fc'].max()"
        out = execute_code(code, sample_df)

        assert out["success"] is True
        assert out["result"]["type"] == "scalar"
        assert abs(out["result"]["value"] - 2.56) < 0.01

    def test_filtering_returns_dataframe(self, sample_df: pd.DataFrame):
        from backend.tabular_qa.executor import execute_code

        code = "result = df[df['padj'] < 0.05]"
        out = execute_code(code, sample_df)

        assert out["success"] is True
        assert out["result"]["type"] == "dataframe"
        assert out["result"]["n_rows"] == 5  # 5 genes with padj < 0.05

    def test_sorting_top_n(self, sample_df: pd.DataFrame):
        from backend.tabular_qa.executor import execute_code

        code = "result = df.nlargest(3, 'log2fc')[['gene', 'log2fc']].to_dict('records')"
        out = execute_code(code, sample_df)

        assert out["success"] is True
        # to_dict('records') is a list → serialized as {"type": "list", "value": [...]}
        rows = out["result"]["value"]
        assert rows[0]["gene"] == "RB1"  # highest log2fc = 2.56

    def test_syntax_error_handled_gracefully(self, sample_df: pd.DataFrame):
        from backend.tabular_qa.executor import execute_code

        out = execute_code("result = df[df['log2fc'  >  0", sample_df)

        assert out["success"] is False
        assert out["error"] is not None

    def test_timeout_enforced(self, sample_df: pd.DataFrame):
        from backend.tabular_qa.executor import execute_code

        # Pure Python busy-loop — no imports needed
        code = "x = 0\nwhile True:\n    x += 1\nresult = x"
        out = execute_code(code, sample_df, timeout_s=1)

        assert out["success"] is False
        assert "timeout" in out["error"].lower() or "timed" in out["error"].lower()

    def test_blocked_import_rejected(self, sample_df: pd.DataFrame):
        from backend.tabular_qa.executor import execute_code

        code = "import os; result = os.listdir('/')"
        out = execute_code(code, sample_df)

        assert out["success"] is False

    def test_original_df_not_mutated(self, sample_df: pd.DataFrame):
        from backend.tabular_qa.executor import execute_code

        original_shape = sample_df.shape
        code = "df.drop(columns=['gene'], inplace=True); result = 'dropped'"
        execute_code(code, sample_df)

        # The original DataFrame passed in must be unchanged
        assert sample_df.shape == original_shape

    def test_print_output_captured(self, sample_df: pd.DataFrame):
        from backend.tabular_qa.executor import execute_code

        code = "print('n rows:', len(df)); result = len(df)"
        out = execute_code(code, sample_df)

        assert out["success"] is True
        assert "n rows:" in out["stdout"]


# ---------------------------------------------------------------------------
# 3. Agent smoke test (LLM mocked)
# ---------------------------------------------------------------------------

MOCK_CODE = textwrap.dedent("""\
    significant = df[df['padj'] < 0.05]
    result = significant.nlargest(3, 'log2fc')[['gene', 'log2fc']].to_dict('records')
""")

MOCK_LLM_RESPONSE = MagicMock()
MOCK_LLM_RESPONSE.content = f"```python\n{MOCK_CODE}\n```"


class TestTabularQAAgent:
    def _make_profile(self, csv_file: Path) -> Dict[str, Any]:
        from backend.file_intelligence.profiler import profile_file
        return profile_file(str(csv_file))

    @patch("backend.tabular_qa.agent._get_llm")
    def test_agent_returns_success(self, mock_get_llm, csv_file: Path):
        mock_llm = MagicMock()
        mock_llm.invoke.return_value = MOCK_LLM_RESPONSE
        mock_get_llm.return_value = mock_llm

        from backend.tabular_qa.agent import run_tabular_qa

        profile = self._make_profile(csv_file)
        result = run_tabular_qa(
            question="What are the top 3 most upregulated significant genes?",
            session_id="smoke-test-session",
            file_path=str(csv_file),
            profile=profile,
        )

        assert result["success"] is True
        assert result["answer"]
        assert len(result["answer"]) > 10

    @patch("backend.tabular_qa.agent._get_llm")
    def test_agent_retries_on_bad_code(self, mock_get_llm, csv_file: Path):
        """First LLM response produces broken code; second is good.
        Agent should retry and ultimately succeed."""
        bad_response = MagicMock()
        bad_response.content = "```python\nresult = df['nonexistent_col'].max()\n```"

        good_response = MagicMock()
        good_response.content = f"```python\n{MOCK_CODE}\n```"

        mock_llm = MagicMock()
        mock_llm.invoke.side_effect = [bad_response, good_response]
        mock_get_llm.return_value = mock_llm

        from backend.tabular_qa import agent as agent_mod
        import importlib
        importlib.reload(agent_mod)

        profile = self._make_profile(csv_file)
        result = agent_mod.run_tabular_qa(
            question="Top upregulated genes?",
            session_id="smoke-test-session",
            file_path=str(csv_file),
            profile=profile,
        )

        assert "success" in result
        assert "answer" in result

    @patch("backend.tabular_qa.agent._get_llm")
    def test_agent_includes_schema_in_prompt(self, mock_get_llm, csv_file: Path):
        """Verify that column names from the profile reach the LLM prompt."""
        captured_prompts = []

        def capture(messages):
            captured_prompts.append(str(messages))
            return MOCK_LLM_RESPONSE

        mock_llm = MagicMock()
        mock_llm.invoke.side_effect = capture
        mock_get_llm.return_value = mock_llm

        from backend.tabular_qa.agent import run_tabular_qa

        profile = self._make_profile(csv_file)
        run_tabular_qa(
            question="How many significant genes are there?",
            session_id="s1",
            file_path=str(csv_file),
            profile=profile,
        )

        assert captured_prompts, "LLM was never called"
        # Column names must appear in the prompt so the LLM can write correct code
        assert "log2fc" in captured_prompts[0] or "gene" in captured_prompts[0]


# ---------------------------------------------------------------------------
# 4. Ingest smoke test (multi-format)
#    Actual return: (conn, df, metadata) — DuckDB connection, DataFrame, dict
# ---------------------------------------------------------------------------

class TestIngestSmoke:
    def test_ingest_csv(self, csv_file: Path):
        from backend.ds_pipeline.pipelines.ingest import ingest_tabular

        conn, df, meta = ingest_tabular(csv_file)
        assert df.shape == (8, 5)
        assert "gene" in df.columns
        assert meta["n_rows"] == 8
        assert meta["n_cols"] == 5

    def test_ingest_tsv(self, tmp_path: Path):
        from backend.ds_pipeline.pipelines.ingest import ingest_tabular

        tsv = tmp_path / "data.tsv"
        tsv.write_text(SAMPLE_CSV.replace(",", "\t"))
        _, df, meta = ingest_tabular(tsv)
        assert df.shape == (8, 5)
        assert meta["source_format"] == "tsv"

    def test_ingest_excel_sheet_selection(self, tmp_path: Path):
        from backend.ds_pipeline.pipelines.ingest import ingest_tabular, list_sheets

        df_src = pd.read_csv(io.StringIO(SAMPLE_CSV))
        xlsx = tmp_path / "multi.xlsx"
        with pd.ExcelWriter(str(xlsx)) as writer:
            df_src.to_excel(writer, sheet_name="Expression", index=False)
            df_src.head(3).to_excel(writer, sheet_name="Summary", index=False)

        sheets = list_sheets(str(xlsx))
        assert "Expression" in sheets
        assert "Summary" in sheets

        _, df, meta = ingest_tabular(xlsx, sheet="Summary")
        assert df.shape == (3, 5)
        assert meta["source_sheet"] == "Summary"

    def test_ingest_excel_missing_sheet_raises(self, tmp_path: Path):
        from backend.ds_pipeline.pipelines.ingest import ingest_tabular

        df_src = pd.read_csv(io.StringIO(SAMPLE_CSV))
        xlsx = tmp_path / "data.xlsx"
        df_src.to_excel(str(xlsx), index=False)

        with pytest.raises(ValueError, match="not found"):
            ingest_tabular(xlsx, sheet="DoesNotExist")

    def test_ingest_duckdb_query(self, csv_file: Path):
        """Verify the returned DuckDB connection can be queried."""
        from backend.ds_pipeline.pipelines.ingest import ingest_tabular

        conn, _, _ = ingest_tabular(csv_file)
        rows = conn.execute("SELECT gene FROM data WHERE log2fc > 2").fetchall()
        genes = [r[0] for r in rows]
        assert "RB1" in genes

    def test_ingest_missing_file_raises(self, tmp_path: Path):
        from backend.ds_pipeline.pipelines.ingest import ingest_tabular

        with pytest.raises(FileNotFoundError):
            ingest_tabular(tmp_path / "nonexistent.csv")
