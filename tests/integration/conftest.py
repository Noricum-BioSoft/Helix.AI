"""Shared pytest fixtures for integration tests.

All LLM classifiers (command router, staging classifier, intent classifier,
approval classifier) are mocked so integration tests run deterministically in
HELIX_MOCK_MODE=1 without live API calls.

Mocking philosophy
------------------
- ``CommandRouter._route_with_llm``: replaced with a deterministic dispatch
  table keyed on command content.  This is acceptable test infrastructure —
  it mirrors what the real LLM would return for well-known commands.
- ``_classify_staging_intent``: returns a simple staging decision based on
  whether the command is a vague planning request or a concrete execution.
- ``intent_classifier.classify_intent``: returns "qa" for question-like
  commands, "execute" for everything else.
- ``approval_classifier._get_llm``: returns a mock LLM that always classifies
  as non-approval.  Note: the keyword fast-path in ``classify_approval``
  handles real approval phrases like "Approve." *before* the LLM is called,
  so the mock is only reached for ambiguous inputs.
"""
from __future__ import annotations

import os
import re
from typing import Any, Dict, Optional, Tuple
from unittest.mock import MagicMock, patch

import pytest


# ---------------------------------------------------------------------------
# Deterministic command → tool dispatch (mirrors the real LLM behaviour)
# ---------------------------------------------------------------------------

def _mock_route(self, command: str, session_context: Dict[str, Any]) -> Tuple[str, Dict[str, Any]]:
    """Return (tool_name, base_params) for *command* without LLM."""
    c = command.lower()
    sid = session_context.get("session_id", "")

    # Workflow design / planning requests → planning handler (must precede specific tools)
    if any(x in c for x in [
        "design the workflow",
        "design workflow",
        "propose the workflow",
        "plan the workflow",
        "plan this workflow",
        "design the analysis",
        "what is going on",
        "tell me what",
        "and tell me",
    ]):
        return "handle_natural_command", {"command": command, "session_id": sid}

    # Plot / visualization updates (highlight, color, update plot, etc.)
    if any(x in c for x in ["highlight", "on the volcano", "on the plot", "update the plot", "color the"]):
        return "patch_and_rerun", {}

    # Historical recreation / diff / compare-runs
    if any(
        x in c
        for x in [
            "recreate the figure",
            "historical recreation",
            "differences between my runs",
            "show differences between",
            "compare.*runs",
            "bio_diff_runs",
        ]
    ) or (re.search(r"\bcompare\b", c) and re.search(r"\bruns?\b", c)):
        return "bio_diff_runs", {}

    # Rerun / iterate
    if any(x in c for x in ["outlier", "exclude them and rerun", "rerun with", "re-run", "reanalyze"]):
        return "bio_rerun", {}

    # Single-cell RNA-seq
    if any(x in c for x in ["scrna", "scRNA-seq", "single cell", "single-cell", "pbmc", ".h5"]):
        return "single_cell_analysis", {}

    # FastQC / amplicon QC (check before generic fastq / S3 to pick specific tool)
    if any(x in c for x in ["fastqc", "amplicon qc", "quality control", "multiqc"]):
        return "fastqc_quality_analysis", {}

    # Bulk RNA-seq
    if any(
        x in c
        for x in [
            "bulk rna-seq",
            "bulk rna seq",
            "rnaseq",
            "rna-seq",
            "deseq",
            "transcriptomics",
            "differential expression",
            "time-course",
        ]
    ):
        return "bulk_rnaseq_analysis", {}

    # Phylogenetics
    if any(x in c for x in ["phylo", "evolutionary tree", "phylogenetic"]):
        return "phylogenetic_tree", {}

    # Tool inventory
    if any(x in c for x in ["list tools", "toolbox", "what tools", "capabilities", "what can you"]):
        return "toolbox_inventory", {}

    # Visualise existing job
    if any(x in c for x in ["visualize the results", "show results of", "visualize results"]):
        return "visualize_job_results", {}

    # Default: general command / planning handler
    return "handle_natural_command", {"command": command, "session_id": sid}


# ---------------------------------------------------------------------------
# Deterministic staging decision
# ---------------------------------------------------------------------------

def _mock_staging(command: str) -> "StagingDecision":  # type: ignore[name-defined]
    from backend.orchestration.approval_policy import StagingDecision

    c = (command or "").lower()

    # Vague "tell me what is going on"-style requests → planning
    if (
        "and tell me what is going on" in c
        or "design the workflow" in c
        or (
            "analyze" in c
            and "dataset" in c
            and "s3://" not in c
            and not c.strip().startswith("run")
        )
    ):
        return StagingDecision(
            requires_approval=False,
            has_execute_intent=False,
            is_planning_request=True,
            method="mock",
            reason="mock_planning_request",
        )

    return StagingDecision(
        requires_approval=False,
        has_execute_intent=True,
        is_planning_request=False,
        method="mock",
        reason="mock_execute",
    )


# ---------------------------------------------------------------------------
# Autouse fixture
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _mock_llm_classifiers():
    """Patch all LLM classifiers when HELIX_MOCK_MODE=1 (the default in tests)."""
    if os.getenv("HELIX_MOCK_MODE", "0") != "1":
        yield
        return

    from backend.command_router import CommandRouter
    from backend.intent_classifier import IntentDecision

    # Approval classifier LLM mock (keyword fast-path still fires first for
    # real approval phrases like "Approve." so this is only reached for
    # ambiguous short inputs).
    _non_approval_llm = MagicMock()
    _non_approval_llm.invoke.return_value.content = '{"approval": false}'

    # Intent mock: questions → qa, everything else → execute
    _QUESTION_WORDS = re.compile(
        r"^(?:what|how|why|where|when|who|which|whose|is|are|am|was|were|"
        r"can|could|may|might|will|would|shall|should|do|does|did|has|have|had)\b",
        re.IGNORECASE,
    )

    def _classify_intent(text, *, session_context=None, workflow_state=None):
        t = (text or "").strip()
        if t.endswith("?") or _QUESTION_WORDS.match(t):
            return IntentDecision(intent="qa", reason="mock_question")
        return IntentDecision(intent="execute", reason="mock_execute")

    with (
        patch.object(CommandRouter, "_route_with_llm", _mock_route),
        patch(
            "backend.orchestration.approval_policy._classify_staging_intent",
            side_effect=_mock_staging,
        ),
        patch(
            "backend.intent_classifier.classify_intent",
            side_effect=_classify_intent,
        ),
        patch(
            "backend.orchestration.approval_classifier._get_llm",
            return_value=_non_approval_llm,
        ),
    ):
        yield
