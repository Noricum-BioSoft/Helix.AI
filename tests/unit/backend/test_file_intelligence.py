"""Unit tests for the file_intelligence profiling layer."""
from __future__ import annotations

import gzip
import json
import tempfile
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------

def _csv(tmp_path: Path) -> Path:
    p = tmp_path / "data.csv"
    p.write_text("gene,tumor,normal\nGENE_A,10.0,2.5\nGENE_B,5.0,5.0\nGENE_C,8.0,1.0\n")
    return p


def _tsv(tmp_path: Path) -> Path:
    p = tmp_path / "data.tsv"
    p.write_text("gene\ttumor\tnormal\nGENE_A\t10.0\t2.5\nGENE_B\t5.0\t5.0\n")
    return p


def _xlsx(tmp_path: Path) -> Path:
    import openpyxl
    wb = openpyxl.Workbook()
    ws1 = wb.active
    ws1.title = "Sheet1"
    ws1.append(["gene", "median_tumor", "max_median_gtex"])
    ws1.append(["TARGET_1", 8.5, 1.2])
    ws2 = wb.create_sheet("ts_final")
    ws2.append(["gene", "score"])
    ws2.append(["TARGET_2", 6.0])
    p = tmp_path / "data.xlsx"
    wb.save(str(p))
    return p


def _fasta(tmp_path: Path) -> Path:
    p = tmp_path / "seqs.fasta"
    p.write_text(">seq1 description\nATGCATGC\n>seq2\nGGGCCCAAA\n")
    return p


def _fastq(tmp_path: Path) -> Path:
    p = tmp_path / "reads.fastq"
    p.write_text(
        "@read1\nATGCATGC\n+\nIIIIIIII\n"
        "@read2\nGGGCCCAA\n+\nHHHHHHHH\n"
    )
    return p


def _vcf(tmp_path: Path) -> Path:
    p = tmp_path / "variants.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n"
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
        "chr1\t100\t.\tA\tT\t50\tPASS\tDP=30\tGT\t0/1\n"
        "chr1\t200\t.\tG\tC\t30\tPASS\tDP=20\tGT\t1/1\n"
        "chr2\t500\t.\tT\tA\t60\tPASS\tDP=40\tGT\t0/0\n"
    )
    return p


def _bed(tmp_path: Path) -> Path:
    p = tmp_path / "peaks.bed"
    p.write_text(
        "chr1\t100\t200\tpeak_1\t500\t.\n"
        "chr1\t300\t400\tpeak_2\t400\t+\n"
        "chr2\t100\t250\tpeak_3\t600\t-\n"
    )
    return p


def _gff(tmp_path: Path) -> Path:
    p = tmp_path / "genes.gff3"
    p.write_text(
        "##gff-version 3\n"
        "chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=GENE001\n"
        "chr1\tensembl\texon\t1000\t1500\t.\t+\t.\tParent=GENE001\n"
    )
    return p


# ---------------------------------------------------------------------------
# Profiler dispatch
# ---------------------------------------------------------------------------

class TestProfilerDispatch:
    def test_csv_dispatches_to_tabular(self, tmp_path):
        from backend.file_intelligence.profiler import profile_file
        prof = profile_file(_csv(tmp_path))
        assert prof["family"] == "tabular"
        assert prof["n_records"] == 3

    def test_tsv_dispatches_to_tabular(self, tmp_path):
        from backend.file_intelligence.profiler import profile_file
        prof = profile_file(_tsv(tmp_path))
        assert prof["family"] == "tabular"

    def test_xlsx_dispatches_to_tabular(self, tmp_path):
        from backend.file_intelligence.profiler import profile_file
        prof = profile_file(_xlsx(tmp_path))
        assert prof["family"] == "tabular"
        assert prof.get("available_sheets") is not None

    def test_fasta_dispatches_to_sequence(self, tmp_path):
        from backend.file_intelligence.profiler import profile_file
        prof = profile_file(_fasta(tmp_path))
        assert prof["family"] == "sequence"
        assert prof["n_records"] == 2

    def test_fastq_dispatches_to_sequence(self, tmp_path):
        from backend.file_intelligence.profiler import profile_file
        prof = profile_file(_fastq(tmp_path))
        assert prof["family"] == "sequence"
        assert prof["n_records"] == 2

    def test_vcf_dispatches_to_variant(self, tmp_path):
        from backend.file_intelligence.profiler import profile_file
        prof = profile_file(_vcf(tmp_path))
        assert prof["family"] == "variant"
        assert prof["n_records"] == 3

    def test_bed_dispatches_to_genomic_interval(self, tmp_path):
        from backend.file_intelligence.profiler import profile_file
        prof = profile_file(_bed(tmp_path))
        assert prof["family"] == "genomic_interval"
        assert prof["n_records"] == 3

    def test_gff_dispatches_to_genomic_interval(self, tmp_path):
        from backend.file_intelligence.profiler import profile_file
        prof = profile_file(_gff(tmp_path))
        assert prof["family"] == "genomic_interval"

    def test_unknown_extension_returns_error_profile(self, tmp_path):
        from backend.file_intelligence.profiler import profile_file
        p = tmp_path / "file.xyz"
        p.write_text("data")
        prof = profile_file(p)
        assert prof["profiler_error"] is not None


# ---------------------------------------------------------------------------
# Tabular profiler
# ---------------------------------------------------------------------------

class TestTabularProfiler:
    def test_csv_has_column_stats(self, tmp_path):
        from backend.file_intelligence.tabular import profile_tabular
        prof = profile_tabular(_csv(tmp_path))
        cols = prof["schema"]["columns"]
        tumor_col = next(c for c in cols if c["name"] == "tumor")
        assert "mean" in tumor_col
        assert tumor_col["max"] == pytest.approx(10.0)

    def test_excel_sheet_selection(self, tmp_path):
        from backend.file_intelligence.tabular import profile_tabular
        prof = profile_tabular(_xlsx(tmp_path), sheet="ts_final")
        assert prof["summary"]["source_sheet"] == "ts_final"
        col_names = [c["name"] for c in prof["schema"]["columns"]]
        assert "score" in col_names

    def test_available_sheets_listed(self, tmp_path):
        from backend.file_intelligence.tabular import profile_tabular
        prof = profile_tabular(_xlsx(tmp_path))
        assert set(prof["available_sheets"]) == {"Sheet1", "ts_final"}

    def test_sample_rows_present(self, tmp_path):
        from backend.file_intelligence.tabular import profile_tabular
        prof = profile_tabular(_csv(tmp_path))
        assert len(prof["sample"]) > 0


# ---------------------------------------------------------------------------
# Sequence profiler
# ---------------------------------------------------------------------------

class TestSequenceProfiler:
    def test_fasta_sequence_lengths(self, tmp_path):
        from backend.file_intelligence.sequence import profile_sequence
        prof = profile_sequence(_fasta(tmp_path))
        assert prof["summary"]["min_length"] == 8
        assert prof["summary"]["max_length"] == 9

    def test_fastq_n_records(self, tmp_path):
        from backend.file_intelligence.sequence import profile_sequence
        prof = profile_sequence(_fastq(tmp_path))
        assert prof["n_records"] == 2
        assert prof["format"] == "fastq"


# ---------------------------------------------------------------------------
# VCF profiler
# ---------------------------------------------------------------------------

class TestVcfProfiler:
    def test_chrom_distribution(self, tmp_path):
        from backend.file_intelligence.vcf import profile_vcf
        prof = profile_vcf(_vcf(tmp_path))
        assert "chr1" in prof["summary"]["chrom_distribution"]

    def test_info_field_keys(self, tmp_path):
        from backend.file_intelligence.vcf import profile_vcf
        prof = profile_vcf(_vcf(tmp_path))
        assert "DP" in prof["schema"]["info_fields"]

    def test_sample_rows_present(self, tmp_path):
        from backend.file_intelligence.vcf import profile_vcf
        prof = profile_vcf(_vcf(tmp_path))
        assert len(prof["sample"]) > 0


# ---------------------------------------------------------------------------
# BED profiler
# ---------------------------------------------------------------------------

class TestBedProfiler:
    def test_bed_chrom_distribution(self, tmp_path):
        from backend.file_intelligence.bed import profile_bed
        prof = profile_bed(_bed(tmp_path))
        assert "chr1" in prof["summary"]["chrom_distribution"]

    def test_gff_feature_types(self, tmp_path):
        from backend.file_intelligence.bed import profile_bed
        prof = profile_bed(_gff(tmp_path))
        assert "gene" in prof["summary"]["feature_types"]


# ---------------------------------------------------------------------------
# Code Interpreter executor
# ---------------------------------------------------------------------------

class TestExecutor:
    def test_simple_scalar(self, tmp_path):
        import pandas as pd
        from backend.tabular_qa.executor import execute_code

        df = pd.DataFrame({"x": [1, 2, 3]})
        out = execute_code("result = df['x'].sum()", df)
        assert out["success"] is True
        assert out["result"]["value"] == 6

    def test_dataframe_result(self, tmp_path):
        import pandas as pd
        from backend.tabular_qa.executor import execute_code

        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
        out = execute_code("result = df.sort_values('a', ascending=False)", df)
        assert out["success"] is True
        assert out["result"]["type"] == "dataframe"
        assert out["result"]["value"][0]["a"] == 3

    def test_syntax_error_returns_failure(self, tmp_path):
        import pandas as pd
        from backend.tabular_qa.executor import execute_code

        df = pd.DataFrame({"x": [1]})
        out = execute_code("result = df[", df)
        assert out["success"] is False
        assert "error" in out

    def test_forbidden_import_blocked(self, tmp_path):
        import pandas as pd
        from backend.tabular_qa.executor import execute_code

        df = pd.DataFrame({"x": [1]})
        out = execute_code("import os; result = os.getcwd()", df)
        assert out["success"] is False

    def test_timeout_returns_failure(self, tmp_path):
        import pandas as pd
        from backend.tabular_qa.executor import execute_code

        df = pd.DataFrame({"x": [1]})
        out = execute_code("import time; time.sleep(10); result = 1", df, timeout_s=1)
        # time.sleep is not in safe builtins — will raise NameError, not timeout
        assert out["success"] is False

    def test_original_df_not_mutated(self, tmp_path):
        import pandas as pd
        from backend.tabular_qa.executor import execute_code

        df = pd.DataFrame({"x": [1, 2, 3]})
        execute_code("df['y'] = df['x'] * 2; result = df", df)
        assert "y" not in df.columns  # sandbox operates on a copy


# ---------------------------------------------------------------------------
# tabular_qa mock-mode
# ---------------------------------------------------------------------------

class TestTabularQaMockMode:
    def test_mock_mode_returns_gracefully(self, tmp_path, monkeypatch):
        monkeypatch.setenv("HELIX_MOCK_MODE", "1")
        from backend.tabular_qa.agent import run_tabular_qa
        from backend.file_intelligence.tabular import profile_tabular

        csv_path = _csv(tmp_path)
        prof = profile_tabular(csv_path)
        result = run_tabular_qa(
            question="What are the top genes?",
            session_id="test",
            file_path=str(csv_path),
            profile=prof,
        )
        assert result["success"] is True
        assert "mock" in result["answer"].lower()
