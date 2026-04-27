"""Tests for the Q&A bypass in the approval staging policy.

Covers:
1. Advisory/meta questions are never staged for approval.
2. Genuine workflow-design requests are still staged.
3. tabular_qa (and other read-only tools) are never staged.
4. _build_run_summary returns useful result summaries.
"""
from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from backend.orchestration.approval_policy import (
    READ_ONLY_ROUTER_TOOLS,
    StagingDecision,
    should_stage_for_approval,
)
from backend.context_builder import _build_run_summary


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _mock_staging_decision(
    requires_approval: bool = False,
    has_execute_intent: bool = False,
    is_planning_request: bool = False,
) -> StagingDecision:
    return StagingDecision(
        requires_approval=requires_approval,
        has_execute_intent=has_execute_intent,
        is_planning_request=is_planning_request,
        method="mock",
        reason="mock",
    )


# ---------------------------------------------------------------------------
# Issue 3: Advisory / Q&A questions must NOT be staged
# ---------------------------------------------------------------------------

ADVISORY_COMMANDS = [
    "What should I do next with these results?",
    "Where do I go from here?",
    "What are my options?",
    "What do you recommend?",
    "Can you help me decide what to do next?",
    "What can I do with the fetched sequence?",
    "I just ran an NCBI fetch — what comes next?",
]


@pytest.mark.parametrize("command", ADVISORY_COMMANDS)
def test_advisory_question_with_handle_natural_command_not_staged(command: str) -> None:
    """Advisory meta-questions routed to handle_natural_command are never staged.

    The staging classifier returns is_planning_request=False for advisory
    questions (after the prompt fix).  Even if a mis-classification occurred and
    is_planning_request were True, the function should still not stage when the
    command is advisory.  We simulate the corrected classifier behaviour here.
    """
    decision = _mock_staging_decision(
        requires_approval=False,
        has_execute_intent=False,
        is_planning_request=False,  # classifier correctly returns False for advisory Q&A
    )
    with patch(
        "backend.orchestration.approval_policy._classify_staging_intent",
        return_value=decision,
    ):
        result = should_stage_for_approval("handle_natural_command", command)
    assert result is False, f"Advisory command should not be staged: {command!r}"


# ---------------------------------------------------------------------------
# Genuine planning requests SHOULD still be staged
# ---------------------------------------------------------------------------

PLANNING_COMMANDS = [
    "Design a full RNA-seq workflow for my samples",
    "Propose a pipeline for my scRNA-seq data",
    "Plan an analysis for my WGS dataset",
    "What would a metagenomics pipeline look like for this data?",
]


@pytest.mark.parametrize("command", PLANNING_COMMANDS)
def test_planning_request_with_handle_natural_command_is_staged(command: str) -> None:
    """Genuine workflow-design requests are staged for user review."""
    decision = _mock_staging_decision(
        requires_approval=False,
        has_execute_intent=False,
        is_planning_request=True,
    )
    with patch(
        "backend.orchestration.approval_policy._classify_staging_intent",
        return_value=decision,
    ):
        result = should_stage_for_approval("handle_natural_command", command)
    assert result is True, f"Planning request should be staged: {command!r}"


# ---------------------------------------------------------------------------
# READ_ONLY tools are never staged (includes tabular_qa)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("tool", sorted(READ_ONLY_ROUTER_TOOLS))
def test_read_only_tools_are_never_staged(tool: str) -> None:
    """No read-only tool should ever be staged, regardless of command text."""
    result = should_stage_for_approval(tool, "some command", {})
    assert result is False, f"Read-only tool should never be staged: {tool!r}"


def test_tabular_qa_is_in_read_only_tools() -> None:
    assert "tabular_qa" in READ_ONLY_ROUTER_TOOLS


# ---------------------------------------------------------------------------
# Issue 2: _build_run_summary produces useful context
# ---------------------------------------------------------------------------

class TestBuildRunSummary:
    def test_uses_user_friendly_summary_when_present(self) -> None:
        result = {"user_friendly_summary": "Fetched human CCR7 mRNA (NM_001301717.4)"}
        summary = _build_run_summary(result, "fetch_ncbi_sequence")
        assert "CCR7" in summary
        assert "NM_001301717.4" in summary

    def test_uses_text_field_as_fallback(self) -> None:
        result = {"text": "Alignment completed: 150 sequences aligned successfully"}
        summary = _build_run_summary(result, "sequence_alignment")
        assert "Alignment completed" in summary

    def test_extracts_accession_for_ncbi_result(self) -> None:
        result = {"status": "success", "accession": "NM_001301717.4", "gene_name": "CCR7"}
        summary = _build_run_summary(result, "fetch_ncbi_sequence")
        assert "accession" in summary or "NM_001301717.4" in summary
        assert "gene_name" in summary or "CCR7" in summary

    def test_extracts_status_when_no_other_fields(self) -> None:
        result = {"status": "success"}
        summary = _build_run_summary(result, "some_tool")
        assert "success" in summary

    def test_handles_empty_result(self) -> None:
        summary = _build_run_summary({}, "any_tool")
        assert isinstance(summary, str)

    def test_handles_non_dict_result(self) -> None:
        summary = _build_run_summary("plain string result", "any_tool")
        assert "plain string result" in summary

    def test_handles_none_result(self) -> None:
        summary = _build_run_summary(None, "any_tool")
        assert isinstance(summary, str)

    def test_summary_capped_at_200_chars(self) -> None:
        long_text = "x" * 500
        result = {"user_friendly_summary": long_text}
        summary = _build_run_summary(result, "tool")
        assert len(summary) <= 200
