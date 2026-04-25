"""
Smoke tests — LLM command router.

Verifies that the live LLM correctly maps natural-language commands to the
expected tool names.  These tests run against the real OpenAI API and are
skipped automatically when no key is present (see conftest.py).

Each test case is a (command, expected_tool) pair.  The assertion checks the
exact tool name returned by CommandRouter._route_with_llm so we can catch
regressions in the prompt / model without running a full end-to-end workflow.
"""
from __future__ import annotations

import os

import pytest

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def router():
    """Return a CommandRouter instance with live LLM enabled."""
    # Ensure mock mode is off (conftest._disable_mock_mode handles env, but
    # CommandRouter reads it at call time, so this double-check is harmless).
    os.environ["HELIX_MOCK_MODE"] = "0"
    from backend.command_router import CommandRouter
    return CommandRouter()


# ---------------------------------------------------------------------------
# Routing cases
# ---------------------------------------------------------------------------

ROUTING_CASES = [
    # ── Genomics / transcriptomics ───────────────────────────────────────────
    (
        "Run differential expression analysis on my RNA-seq counts using DESeq2",
        "bulk_rnaseq_analysis",
    ),
    (
        "I have scRNA-seq data from 10x Genomics. Run clustering and UMAP.",
        "single_cell_analysis",
    ),
    (
        "Run FastQC quality control on my paired-end FASTQ files",
        "fastqc_quality_analysis",
    ),
    (
        "Trim adapters and low-quality bases from my reads",
        "read_trimming",
    ),
    (
        "Align these two DNA sequences and show me the alignment",
        "sequence_alignment",
    ),
    (
        "Fetch the sequence for accession NM_001301717 from NCBI",
        "fetch_ncbi_sequence",
    ),
    (
        "Build a phylogenetic tree from my FASTA sequences",
        "phylogenetic_tree",
    ),
    # ── Tabular / data science ───────────────────────────────────────────────
    (
        "Sort my Excel file by the 'expression' column and show the top 20 rows",
        "tabular_analysis",
    ),
    (
        "Can you answer a question about my CSV: what is the average logFC value? Just tell me the number.",
        "tabular_qa",
    ),
    (
        "Run an exploratory data analysis on my uploaded dataset",
        "ds_run_analysis",
    ),
    # ── Iterative / session ops ──────────────────────────────────────────────
    (
        "Re-run the analysis with alpha=0.01 instead of 0.05",
        "bio_rerun",
    ),
    (
        "Edit the existing volcano plot script to highlight all genes where logFC > 2 in red",
        "patch_and_rerun",
    ),
    (
        "What were the inputs and outputs of the first run in this session?",
        "session_run_io_summary",
    ),
    # ── Utility ─────────────────────────────────────────────────────────────
    (
        "What tools do you have available?",
        "toolbox_inventory",
    ),
]


@pytest.mark.parametrize("command,expected_tool", ROUTING_CASES)
def test_router_maps_command_to_correct_tool(router, command, expected_tool):
    """LLM router returns the expected tool for a well-known command."""
    tool_name, _ = router._route_with_llm(command, {})
    assert tool_name == expected_tool, (
        f"Command: {command!r}\n"
        f"Expected tool: {expected_tool!r}\n"
        f"Got tool:      {tool_name!r}"
    )
