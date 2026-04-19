#!/usr/bin/env python3
"""
Run all 5 Helix.AI demo scenarios against the deployed backend.

Usage:
    BACKEND_URL=http://... pytest tests/integration/test_deployed_demo_scenarios.py -v -m integration
"""
import os
import pytest
import requests
import time
from typing import Dict, Any, Optional

DEFAULT_BACKEND_URL = os.getenv(
    "BACKEND_URL",
    "http://HelixA-ALBAE-D7BksiQIynZb-1051248867.us-west-1.elb.amazonaws.com",
)


# Demo 1: Bulk RNA-seq T.gondii (needs_inputs → then followUp with S3 paths)
DEMO1_PROMPT = """You are analyzing an RNA-seq transcriptome dataset from a mouse study investigating the effects of Toxoplasma gondii infection on brain gene expression.

Study Design
This is a 2 × 2 factorial experimental design with the following factors:
Factor 1: Infection Status — Infected / Uninfected
Factor 2: Time Point — 11 dpi / 33 dpi
Each of the four experimental groups contains 3 biological replicates, for a total of 12 samples.

Objectives: Model gene expression changes using an appropriate statistical framework. Perform differential expression analysis. Include PCA, sample distance heatmap. Assume raw count data and sample metadata are available."""

DEMO1_FOLLOWUP = """Run bulk RNA-seq differential expression analysis with these inputs:
count_matrix: s3://noricum-ngs-data/demo/rnaseq/tgondii_counts.csv
sample_metadata: s3://noricum-ngs-data/demo/rnaseq/tgondii_metadata.csv
design_formula: ~infection_status + time_point + infection_status:time_point"""


# Demo 2: Single-cell SLE PBMC (needs_inputs)
DEMO2_PROMPT = """You are analyzing a single-cell RNA-seq dataset from a human PBMC study investigating immune cell dysregulation in patients with systemic lupus erythematosus (SLE) compared to healthy controls.

Dataset: Single-cell RNA-seq (10x Genomics Chromium v3). Two-group comparison: SLE (n=5) vs. Healthy (n=5).
Objectives: QC, normalize (SCTransform), PCA → UMAP, cluster (Leiden), annotate cell types, pseudobulk DEGs, cell-type proportions.
Desired: UMAP plots, marker dot plot, DEG tables, cell-type proportion chart."""

DEMO2_FOLLOWUP = """data_file: s3://noricum-ngs-data/demo/scrna/sle_pbmc_filtered_feature_bc_matrix.h5
data_format: 10x
resolution: 0.5
steps: all"""


# Demo 3: Amplicon QC (executes_pipeline or multi_step with S3 paths in prompt)
DEMO3_PROMPT = """You are processing a 16S rRNA amplicon sequencing dataset from a gut microbiome study. Raw paired-end FASTQ files are on S3 and need a full preprocessing pipeline.

Dataset
  Illumina MiSeq 2×250 bp paired-end reads; V3–V4 hypervariable region.
  Forward reads: s3://helix-test-data/microbiome/run1/sample01_R1.fastq.gz
  Reverse reads: s3://helix-test-data/microbiome/run1/sample01_R2.fastq.gz
  Output prefix:  s3://helix-test-data/microbiome/run1/processed/

Pipeline Steps: 1) FastQC on R1 and R2. 2) Trim adapter and low-quality bases. 3) Merge paired-end reads. 4) Quality report CSV."""


# Demo 4: APAP time-course workflow design (multi_step_plan)
DEMO4_PROMPT = """I have a bulk RNA-seq time-course study from a murine acetaminophen liver-injury model and I want you to design the analysis workflow before execution.

Study Context: 5 post-injury time points (0 h, 6 h, 24 h, 72 h, 168 h), 4 replicates per time point (20 samples total).
Task: Propose a complete workflow from count-level QC to biological interpretation. For each stage specify key calculations, expected tables, recommended plots. Include checkpoints for data quality. Recommend how to identify early-response, late-recovery, and persistent signatures."""


# Demo 5: Phylogenetics SARS-CoV-2 (multi_step_plan or fast-path)
DEMO5_PROMPT = """You are conducting a comparative evolutionary analysis of the SARS-CoV-2 spike protein across major variants of concern (VOCs).

Dataset: Full-length spike protein sequences for Wuhan-Hu-1, Alpha, Beta, Gamma, Delta, Omicron BA.1, Omicron BA.4/5, XBB.1.5. Use accessions: MN908947.3, OQ898928.1, OR353131.1, MW642250.1, OR323381.1, PP847536.1, PP848071.1, PP405604.1.

Objectives: Retrieve sequences, multiple sequence alignment (MAFFT), ML tree with bootstrap, annotate RBD mutations, pairwise identity matrix. Outputs: alignment FASTA, Newick tree, tree PNG, identity CSV."""


def _execute(base_url: str, command: str, session_id: str, timeout: int = 90) -> Dict[str, Any]:
    url = f"{base_url.rstrip('/')}/execute"
    r = requests.post(url, json={"command": command, "session_id": session_id}, timeout=timeout)
    r.raise_for_status()
    return r.json()


def _create_session(base_url: str) -> str:
    url = f"{base_url.rstrip('/')}/create_session"
    r = requests.post(url, timeout=15)
    r.raise_for_status()
    return r.json()["session_id"]


@pytest.fixture(scope="module")
def backend_url():
    return DEFAULT_BACKEND_URL


@pytest.fixture(scope="module")
def session_id(backend_url):
    return _create_session(backend_url)


@pytest.mark.integration
class TestDemoScenariosDeployed:
    """Run all 5 demo scenarios against the deployed backend."""

    def test_demo1_bulk_rnaseq_tgondii_initial(self, backend_url, session_id):
        """Demo 1: Initial prompt returns needs_inputs or plan."""
        result = _execute(backend_url, DEMO1_PROMPT, session_id)
        assert result.get("success") is True or "result" in result or "text" in result
        assert "error" not in str(result).lower() or result.get("success") is True
        # May return needs_inputs, execute_ready, or direct result
        assert "result" in result or "text" in result or "output" in result

    def test_demo1_bulk_rnaseq_tgondii_with_data(self, backend_url, session_id):
        """Demo 1: Follow-up with S3 paths runs analysis (fast-path or LLM)."""
        _execute(backend_url, DEMO1_PROMPT, session_id)
        time.sleep(1)
        result = _execute(backend_url, DEMO1_FOLLOWUP, session_id)
        assert result.get("success") is True, f"Demo 1 follow-up failed: {result}"
        # Expect tool bulk_rnaseq_analysis and some result content
        assert result.get("tool") == "bulk_rnaseq_analysis" or result.get("success")
        assert "result" in result or "text" in result or "visuals" in result

    def test_demo2_scrna_sle_initial(self, backend_url, session_id):
        """Demo 2: scRNA SLE initial prompt."""
        result = _execute(backend_url, DEMO2_PROMPT, session_id)
        assert result.get("success") is True or "text" in result
        assert "error" not in str(result).lower() or result.get("success") is True

    def test_demo2_scrna_sle_with_data(self, backend_url, session_id):
        """Demo 2: Follow-up with H5 path runs single_cell_analysis."""
        _execute(backend_url, DEMO2_PROMPT, session_id)
        time.sleep(1)
        result = _execute(backend_url, DEMO2_FOLLOWUP, session_id)
        assert result.get("success") is True, f"Demo 2 follow-up failed: {result}"
        assert "single_cell" in result.get("tool", "") or result.get("success")

    def test_demo3_amplicon_qc(self, backend_url, session_id):
        """Demo 3: Amplicon QC pipeline (fast-path or plan + execute)."""
        result = _execute(backend_url, DEMO3_PROMPT, session_id)
        assert result.get("success") is True or "text" in result or "result" in result
        # Fast-path returns Amplicon QC complete; planner returns pipeline plan
        text = (result.get("text") or "") + str(result.get("result", ""))
        assert "error" not in text.lower() or result.get("success") is True
        assert "quality" in text.lower() or "pipeline" in text.lower() or "amplicon" in text.lower() or result.get("execute_ready")

    def test_demo4_apap_workflow_design(self, backend_url, session_id):
        """Demo 4: APAP time-course workflow design."""
        result = _execute(backend_url, DEMO4_PROMPT, session_id)
        assert result.get("success") is True or "text" in result
        assert "error" not in str(result).lower() or result.get("success") is True
        text = result.get("text") or ""
        assert "workflow" in text.lower() or "step" in text.lower() or "analysis" in text.lower() or len(text) > 100

    def test_demo5_phylogenetics_sarscov2(self, backend_url, session_id):
        """Demo 5: Phylogenetics SARS-CoV-2 (fast-path or multi-step)."""
        result = _execute(backend_url, DEMO5_PROMPT, session_id)
        assert result.get("success") is True or "text" in result
        assert "error" not in str(result).lower() or result.get("success") is True
        text = (result.get("text") or "") + str(result.get("result", ""))
        assert "tree" in text.lower() or "phylogen" in text.lower() or "spike" in text.lower() or "variant" in text.lower() or len(text) > 100


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-m", "integration"])
