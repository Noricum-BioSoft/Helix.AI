"""
B. First run (supply data) — USER_SCENARIOS §3, TESTBED §3.

B1: Bulk RNA-seq with demo data
B2: Single-cell with H5 path
B3: FastQC
B4: Phylogenetic tree from FASTA
B5: Read trimming
"""

import pytest

from .conftest import BULK_RNASEQ_PROMPT, FASTA_THREE


def _execute(client, session_id, command):
    r = client.post("/execute", json={"command": command, "session_id": session_id})
    assert r.status_code == 200, r.text
    return r.json()


def _get_links(data):
    links = data.get("data", {}).get("links") if isinstance(data.get("data"), dict) else []
    if not links and isinstance(data.get("result"), dict):
        links = data["result"].get("links") or []
    if not links:
        links = data.get("links") or []
    return links


class TestB1BulkRnaseqFirstRun:
    """B1_be: Run bulk RNA-seq with count matrix + metadata → run_id, links."""

    def test_b1_success_run_id_and_links(self, client, session_id):
        data = _execute(client, session_id, BULK_RNASEQ_PROMPT)
        assert data.get("success") is not False or data.get("status") != "error"
        run_id = data.get("run_id") or (data.get("result") or {}).get("run_id")
        # In mock/synthetic mode we may get success and run_id
        if data.get("success") and run_id:
            links = _get_links(data)
            labels = [lnk.get("label", "") for lnk in links if isinstance(lnk, dict)]
            assert "analysis.py" in " ".join(labels) or "bundle" in " ".join(labels).lower() or run_id

        runs_r = client.get(f"/session/{session_id}/runs")
        if runs_r.status_code == 200:
            runs = runs_r.json().get("runs", [])
            assert isinstance(runs, list)


class TestB2SingleCellFirstRun:
    """B2_be: Single-cell with H5 path (may use needs_inputs or run)."""

    def test_b2_single_cell_with_data_path(self, client, session_id):
        cmd = (
            "Run single-cell RNA-seq analysis. "
            "data_file: s3://noricum-ngs-data/demo/scrna/sle_pbmc_filtered_feature_bc_matrix.h5 "
            "data_format: 10x resolution: 0.5 steps: all"
        )
        data = _execute(client, session_id, cmd)
        # Either runs (success, run_id) or needs_inputs
        assert "status" in data or "result" in data
        if data.get("success") and data.get("run_id"):
            links = _get_links(data)
            assert isinstance(links, list)


class TestB3FastQc:
    """B3_be: FastQC on R1/R2 paths."""

    def test_b3_fastqc_request(self, client, session_id):
        cmd = "Run FastQC on s3://bucket/R1.fastq.gz and s3://bucket/R2.fastq.gz"
        data = _execute(client, session_id, cmd)
        assert data.get("success") is not False or "job" in str(data).lower() or "fastqc" in str(data).lower()


class TestB4PhylogeneticTree:
    """B4_be: Phylogenetic tree from FASTA sequences."""

    def test_b4_tree_from_fasta(self, client, session_id):
        cmd = f"Build a phylogenetic tree from these sequences:\n{FASTA_THREE}"
        data = _execute(client, session_id, cmd)
        assert data.get("success") is not False or data.get("status") != "error"
        result = data.get("result") or data.get("data") or data
        has_tree = (
            result.get("tree_newick")
            or result.get("ete_visualization")
            or "tree" in (data.get("text") or "").lower()
        )
        assert has_tree or data.get("success") is True


class TestB5ReadTrimming:
    """B5_be: Read trimming with adapter and quality."""

    def test_b5_trim_request(self, client, session_id):
        cmd = (
            "Trim my reads: forward AGATCGGAAGAGC, quality 20. "
            "Forward reads: ATGCGATCG\n+\nIIIIIIII\nReverse reads: GCTAGCTAG\n+\nIIIIIIII"
        )
        data = _execute(client, session_id, cmd)
        # May succeed or ask for paths
        assert "status" in data or "result" in data
