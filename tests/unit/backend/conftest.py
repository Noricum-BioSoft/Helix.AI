"""
Shared pytest fixtures for backend unit tests.

These fixtures provide minimal, reusable session/run/artifact state objects
that mirror the structures used in the live session history manager.  Using
them keeps individual test files lean and makes it easy to extend regression
coverage without repeating boilerplate setup.
"""

from __future__ import annotations

import os
from unittest.mock import patch

import pytest


# ---------------------------------------------------------------------------
# LLM mock-mode compatibility layer
#
# Unit tests run with HELIX_MOCK_MODE=1 (no live LLM).  The LLM-based staging
# and intent classifiers raise in mock mode, which would cause unrelated tests
# to fail with "LLM disabled" errors.  These autouse fixtures provide safe
# defaults so dispatch/routing tests focus on what they're actually testing.
#
# Tests that explicitly verify staging or intent classification behaviour
# mock the relevant functions themselves (via @patch), which overrides these
# defaults at the per-test level.
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _mock_staging_classifier(request, monkeypatch):
    """Auto-patch LLM-dependent classification in HELIX_MOCK_MODE.

    Prevents "LLM disabled" errors in dispatch/routing tests that don't need
    real LLM behaviour:

    - ``_classify_staging_intent`` (in approval_policy) → "do not stage"
    - ``approval_policy.is_approval_command`` → False (for staging checks)
    - ``approval_classifier._get_llm`` → mock that returns ``{"approval": false}``
      so ambiguous phrases that slip past early-rejection don't raise.

    Per-test overrides:
    - Tests that explicitly use ``@patch("...approval_classifier._get_llm", ...)``
      will override the _get_llm mock (test @patch runs inside this fixture's
      context and takes precedence).
    - Tests that explicitly use ``@patch("...approval_policy._classify_staging_intent", ...)``
      will likewise override the staging mock.
    """
    if os.getenv("HELIX_MOCK_MODE", "0") != "1":
        yield
        return

    from unittest.mock import MagicMock
    from backend.orchestration.approval_policy import StagingDecision

    _default_staging = StagingDecision(
        requires_approval=False,
        has_execute_intent=True,
        is_planning_request=False,
        method="mock",
        reason="helix_mock_mode",
    )

    # Mock LLM that always returns "not approval" for the approval classifier.
    _false_approval_llm = MagicMock()
    _false_approval_llm.invoke.return_value.content = '{"approval": false}'

    with patch(
        "backend.orchestration.approval_policy._classify_staging_intent",
        return_value=_default_staging,
    ), patch(
        "backend.orchestration.approval_policy.is_approval_command",
        return_value=False,
    ), patch(
        "backend.orchestration.approval_classifier._get_llm",
        return_value=_false_approval_llm,
    ):
        yield


@pytest.fixture(autouse=True)
def _mock_intent_classifier(request, monkeypatch):
    """Auto-patch classify_intent when running in HELIX_MOCK_MODE.

    Returns "execute" intent by default.  Tests that need specific intent
    behaviour must mock ``classify_intent`` themselves.
    """
    if os.getenv("HELIX_MOCK_MODE", "0") != "1":
        yield
        return

    from backend.intent_classifier import IntentDecision

    _default = IntentDecision(intent="execute", reason="helix_mock_mode")

    with patch(
        "backend.intent_classifier.classify_intent",
        return_value=_default,
    ):
        yield


@pytest.fixture(autouse=True)
def _mock_command_router():
    """Auto-patch CommandRouter._route_with_llm in HELIX_MOCK_MODE.

    Provides deterministic LLM-free routing for unit tests.  Tests that
    directly verify routing decisions should patch ``_route_with_llm``
    themselves via ``@patch`` or ``monkeypatch``.
    """
    if os.getenv("HELIX_MOCK_MODE", "0") != "1":
        yield
        return

    import re as _re
    from backend.command_router import CommandRouter

    def _mock_route(self, command: str, session_context):
        c = (command or "").lower()
        sid = (session_context or {}).get("session_id", "")
        # Vague / open-ended analysis
        if any(x in c for x in ["what is going on", "tell me what", "and tell me"]):
            return "handle_natural_command", {"command": command, "session_id": sid}
        # Demo plot (test/local tool)
        if any(x in c for x in ["demo plot", "demo scatter", "local demo plot", "create demo plot", "generate demo plot"]):
            return "local_demo_plot_script", {}
        # Plot / viz updates
        if any(x in c for x in ["highlight", "on the volcano", "on the plot", "update the plot", "color the"]):
            return "patch_and_rerun", {}
        # Plot axis scale
        if any(x in c for x in ["log scale", "linear scale", "switch to log", "switch to linear", "change the plots from"]):
            return "patch_and_rerun", {"change_request": command}
        # Iterative rerun / parameter change
        if any(x in c for x in ["re-run", "rerun", "run again", "redo", "repeat with", "change parameter", "update parameter", "re-analyse", "reanalyze"]):
            return "bio_rerun", {}
        # Session run I/O summary
        if any(x in c for x in ["inputs/outputs of", "inputs and outputs of", "what were the inputs", "what were the outputs"]):
            return "session_run_io_summary", {}
        # Diff / compare runs
        if any(x in c for x in ["recreate", "historical", "differences between", "compare runs", "bio_diff"]) or (
            _re.search(r"\bcompare\b", c) and _re.search(r"\bruns?\b", c)
        ):
            return "bio_diff_runs", {}
        # Data science / EDA pipeline
        if any(x in c for x in ["analyze my data", "analyse my data", "run eda", "data science run", "analyze dataset", "analyse dataset"]):
            return "ds_run_analysis", {}
        # Tabular analysis (specific ops with file extension present)
        if (
            any(x in c for x in ["rank ", "sort ", "ratio ", "filter ", "aggregate ", "sum ", "mean ", "count ", "pivot "])
            and any(x in c for x in [".xlsx", ".xls", ".csv", ".tsv"])
        ):
            return "tabular_analysis", {}
        # Tabular analysis via session file
        if (
            any(x in c for x in ["rank ", "sort ", "ratio ", "filter ", "aggregate "])
            and (session_context or {}).get("uploaded_files")
        ):
            return "tabular_analysis", {}
        # NCBI sequence fetch
        if any(x in c for x in ["fetch sequence", "get sequence", "ncbi", "accession", "genbank", "refseq", "download sequence"]):
            return "fetch_ncbi_sequence", {}
        # Sequence alignment (commands containing actual DNA/RNA sequences)
        if _re.search(r"\b[ACGTNU]{6,}\b", command.upper()) and any(x in c for x in ["align", "alignment", "check mismatch", "compare"]):
            return "sequence_alignment", {}
        # Single-cell
        if any(x in c for x in ["scrna", "single cell", "single-cell", "pbmc", ".h5ad", "anndata", "umap", "seurat"]):
            return "single_cell_analysis", {}
        # FastQC / amplicon
        if any(x in c for x in ["fastqc", "amplicon qc", "quality control", "multiqc"]):
            return "fastqc_quality_analysis", {}
        # Bulk RNA-seq
        if any(x in c for x in ["bulk rna", "rnaseq", "rna-seq", "deseq", "transcriptomics", "differential expression", "time-course"]):
            return "bulk_rnaseq_analysis", {}
        # Phylogenetics
        if any(x in c for x in ["phylo", "evolutionary tree", "phylogenetic"]):
            return "phylogenetic_tree", {}
        # Read ops
        if any(x in c for x in ["trim reads", "trim", "trimming", "adapter"]):
            return "read_trimming", {}
        if any(x in c for x in ["merge reads", "merge paired"]):
            return "read_merging", {}
        # Toolbox
        if any(x in c for x in ["list tools", "toolbox", "what tools", "capabilities", "what can you do"]):
            return "toolbox_inventory", {}
        # Visualise existing job
        if any(x in c for x in ["visualize the results", "show results of", "visualize results"]):
            return "visualize_job_results", {}
        # Default
        return "handle_natural_command", {"command": command, "session_id": sid}

    with patch.object(CommandRouter, "_route_with_llm", _mock_route):
        yield


# ─────────────────────────────────────────────────────────────────────────────
# Session state fixtures
# ─────────────────────────────────────────────────────────────────────────────


@pytest.fixture
def empty_session() -> dict:
    """Bare session with no runs and no artifacts."""
    return {"runs": [], "artifacts": []}


@pytest.fixture
def single_run_session() -> dict:
    """Session with one completed bulk RNA-seq run."""
    return {
        "runs": [
            {
                "run_id": "run_v1",
                "tool": "bulk_rnaseq_analysis",
                "command": "Run DEG analysis on counts.csv",
                "status": "success",
                "tool_args": {
                    "count_matrix": "counts.csv",
                    "sample_metadata": "meta.csv",
                    "design_formula": "~condition",
                    "alpha": 0.05,
                },
            }
        ],
        "artifacts": [
            {
                "artifact_id": "art_deg_v1",
                "title": "DEG_table_v1",
                "type": "deg_table",
                "state_tags": ["initial"],
                "source_run_id": "run_v1",
            }
        ],
    }


@pytest.fixture
def two_run_session() -> dict:
    """Session with two completed runs — used for current-vs-first comparisons."""
    return {
        "runs": [
            {
                "run_id": "run_v1",
                "tool": "bulk_rnaseq_analysis",
                "command": "Run initial DEG analysis",
                "status": "success",
            },
            {
                "run_id": "run_v2",
                "tool": "bulk_rnaseq_analysis",
                "command": "Rerun after metadata correction",
                "status": "success",
            },
        ],
        "artifacts": [
            {
                "artifact_id": "art_v1",
                "title": "DEG_table_v1",
                "type": "deg_table",
                "state_tags": ["initial"],
                "source_run_id": "run_v1",
            },
            {
                "artifact_id": "art_v2",
                "title": "DEG_table_v2",
                "type": "deg_table",
                "state_tags": ["corrected", "metadata"],
                "source_run_id": "run_v2",
            },
        ],
    }


@pytest.fixture
def batch_exclusion_session() -> dict:
    """Session used for turn_16 'use before batch exclusion' patterns."""
    return {
        "runs": [],
        "artifacts": [
            {
                "artifact_id": "a1",
                "title": "cleaned_dataset_v1",
                "state_tags": ["cleaned", "metadata_corrected"],
                "source_run_id": "run_clean",
            },
            {
                "artifact_id": "a2",
                "title": "cleaned_dataset_after_batch_exclusion_v2",
                "state_tags": ["cleaned", "batch_exclusion"],
                "source_run_id": "run_excluded",
            },
        ],
    }


@pytest.fixture
def historical_figures_session() -> dict:
    """Session used for turn_19 'corrected metadata before fold-change fix' patterns."""
    return {
        "runs": [],
        "artifacts": [
            {
                "artifact_id": "fig_v2",
                "title": "figure_set_corrected_metadata_v2",
                "state_tags": ["corrected", "metadata"],
                "source_run_id": "run_v2",
            },
            {
                "artifact_id": "fig_v3",
                "title": "figure_set_corrected_metadata_after_fold_change_bug_fix_v3",
                "state_tags": ["corrected", "metadata", "fold_change_bug_fix"],
                "source_run_id": "run_v3",
            },
        ],
    }


# ─────────────────────────────────────────────────────────────────────────────
# Benchmark scoring helpers
# ─────────────────────────────────────────────────────────────────────────────


@pytest.fixture
def bench_spec_factory():
    """
    Factory that builds a minimal benchmark spec around a single turn.

    Usage::

        def test_foo(bench_spec_factory):
            spec = bench_spec_factory(
                user="Approve.",
                behaviors=["execute_approved_workflow"],
                should_require_approval=False,
            )
    """

    def _make(
        user: str,
        behaviors: list[str],
        *,
        should_require_approval: bool = False,
        turn_id: str = "t1",
    ) -> dict:
        return {
            "evaluation_dimensions": [
                {"id": "artifact_persistence", "severity": "high"},
                {"id": "workflow_selection", "severity": "high"},
                {"id": "scientific_intent_resolution", "severity": "high"},
            ],
            "scenarios": [
                {
                    "id": "s1",
                    "turns": [
                        {
                            "id": turn_id,
                            "user": user,
                            "expected_behaviors": behaviors,
                            "should_require_approval_before_execution": should_require_approval,
                            "focus_dimensions": [
                                "workflow_selection",
                                "scientific_intent_resolution",
                            ],
                        }
                    ],
                }
            ],
        }

    return _make


@pytest.fixture
def score_turn(bench_spec_factory):
    """
    Fixture that scores a single turn row and returns the per-turn score dict.

    Usage::

        def test_foo(score_turn):
            result = score_turn(
                user="Approve.",
                behaviors=["execute_approved_workflow"],
                row={"turn_id": "t1", "tool": "__plan__", "status": "success", "text": "Done."},
            )
            assert result["score"] == 2
    """
    from benchmarks.bio_benchmark_scorer import score_benchmark_run

    def _score(
        user: str,
        behaviors: list[str],
        row: dict,
        *,
        should_require_approval: bool = False,
    ) -> dict:
        spec = bench_spec_factory(
            user,
            behaviors,
            should_require_approval=should_require_approval,
        )
        run_data = {"session_id": "sid", "scenario_id": "s1", "rows": [row]}
        scored = score_benchmark_run(run_data, spec)
        return scored["per_turn_scores"][0]

    return _score
