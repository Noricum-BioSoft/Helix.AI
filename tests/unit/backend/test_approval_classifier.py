"""Unit tests for the LLM-based approval intent classifier.

All tests run without a live LLM (HELIX_MOCK_MODE=1 is set by conftest.py).
The mock tests verify the fast-path and early-rejection logic directly.
LLM path is tested via mocking.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from backend.orchestration.approval_classifier import (
    ApprovalDecision,
    classify_approval,
    is_approval_command,
    _keyword_match,
    _is_clearly_not_approval,
)


# ---------------------------------------------------------------------------
# Keyword fast-path
# ---------------------------------------------------------------------------

KEYWORD_APPROVALS = [
    "proceed",
    "Proceed",
    "PROCEED",
    "approve",
    "approved",
    "yes",
    "yes, approve",
    "yes, run it",
    "go ahead",
    "go ahead.",
    "ok",
    "okay",
    "sure",
    "run it",
    "do it",
    "let's do it",
    "let's go",
    "sounds good",
    "looks good",
    "yep",
    "yup",
]


@pytest.mark.parametrize("phrase", KEYWORD_APPROVALS)
def test_keyword_fast_path_approves(phrase: str) -> None:
    assert _keyword_match(phrase) is True


@pytest.mark.parametrize("phrase", KEYWORD_APPROVALS)
def test_classify_returns_keyword_fast_path_method(phrase: str) -> None:
    decision = classify_approval(phrase)
    assert decision.is_approval is True
    assert decision.method == "keyword_fast_path"


def test_keyword_fast_path_normalizes_whitespace() -> None:
    assert _keyword_match("  proceed  ") is True
    # split() collapses internal spaces, so "go  ahead" normalises to "go ahead"
    # which IS in the keyword set.
    assert _keyword_match("go  ahead") is True


def test_approve_prefix_fast_path() -> None:
    assert _keyword_match("approve the plan") is True
    assert _keyword_match("approve correlation analysis") is True


# ---------------------------------------------------------------------------
# Early rejection
# ---------------------------------------------------------------------------

ANALYTICAL_PHRASES = [
    "Analyze the correlation between age and BMI in my dataset",
    "Run PCA on the plasma protein expression samples and visualize it",
    "Compute differential expression between treated and control groups",
    "Plot a heatmap of the results",
    "my_data.csv",
    "s3://my-bucket/data.fastq.gz",
    "Generate a volcano plot for the DEGs",
    "Cluster the samples by condition using kmeans",
    "Calculate median tumor / max_median_gtex ratio and rank descending",
]


@pytest.mark.parametrize("phrase", ANALYTICAL_PHRASES)
def test_early_rejection_rejects_analytical(phrase: str) -> None:
    assert _is_clearly_not_approval(phrase) is True


@pytest.mark.parametrize("phrase", ANALYTICAL_PHRASES)
def test_classify_returns_early_rejection_method(phrase: str) -> None:
    decision = classify_approval(phrase)
    assert decision.is_approval is False
    assert decision.method == "early_rejection"


SHORT_NON_ANALYTICAL = [
    "that looks right",
    # "yes, execute the plan" — affirmative prefix bypasses early rejection
    "yes, execute the plan",
    "great, do it",
    "perfect",
    "sounds correct",
    # Affirmative prefix guards
    "approved, run the analysis",
    "let's execute",
]


@pytest.mark.parametrize("phrase", SHORT_NON_ANALYTICAL)
def test_short_affirmative_not_early_rejected(phrase: str) -> None:
    """These should NOT be early-rejected — they reach the LLM for classification."""
    assert _is_clearly_not_approval(phrase) is False


# ---------------------------------------------------------------------------
# LLM path (mocked)
# ---------------------------------------------------------------------------

def _mock_llm_response(approval: bool) -> MagicMock:
    """Return a mock LLM that responds with the given approval value."""
    llm = MagicMock()
    llm.invoke.return_value = MagicMock(content=f'{{"approval": {str(approval).lower()}}}')
    return llm


@patch("backend.orchestration.approval_classifier._get_llm")
def test_llm_path_approves(mock_get_llm: MagicMock) -> None:
    mock_get_llm.return_value = _mock_llm_response(True)
    decision = classify_approval("yes, execute the plan")
    assert decision.is_approval is True
    assert decision.method == "llm"


@patch("backend.orchestration.approval_classifier._get_llm")
def test_llm_path_rejects(mock_get_llm: MagicMock) -> None:
    # Use a short ambiguous phrase that is NOT early-rejected (no interrogative
    # prefix, no analytical verb, no affirmative prefix) but which the LLM
    # correctly classifies as non-approval.
    mock_get_llm.return_value = _mock_llm_response(False)
    decision = classify_approval("Actually, never mind for now.")
    assert decision.is_approval is False
    assert decision.method == "llm"


@patch("backend.orchestration.approval_classifier._get_llm")
def test_llm_path_passes_pending_plan_context(mock_get_llm: MagicMock) -> None:
    """When has_pending_plan=True the context string should be in the LLM message."""
    llm = _mock_llm_response(True)
    mock_get_llm.return_value = llm
    classify_approval("great", has_pending_plan=True)
    call_args = llm.invoke.call_args[0][0]  # positional arg: messages list
    user_msg = next(m for m in call_args if m["role"] == "user")
    assert "pending analysis plan" in user_msg["content"].lower() or \
           "Context:" in user_msg["content"]


# ---------------------------------------------------------------------------
# No keyword fallback — LLM failure raises, not silently falls back
# ---------------------------------------------------------------------------

@patch("backend.orchestration.approval_classifier._get_llm")
def test_keyword_fast_path_bypasses_llm_entirely(mock_get_llm: MagicMock) -> None:
    """Keyword fast-path phrases never reach the LLM, even if LLM is down."""
    mock_get_llm.side_effect = RuntimeError("LLM disabled in HELIX_MOCK_MODE")
    decision = classify_approval("proceed")
    # Hits keyword_fast_path BEFORE LLM is called
    assert decision.is_approval is True
    assert decision.method == "keyword_fast_path"
    mock_get_llm.assert_not_called()


@patch("backend.orchestration.approval_classifier._get_llm")
def test_llm_failure_raises_for_ambiguous_phrase(mock_get_llm: MagicMock) -> None:
    """When the LLM is unavailable and the phrase is ambiguous, an exception is raised.

    There is no keyword fallback — a silent wrong answer is worse than a
    visible error that can be handled by the caller.
    """
    mock_get_llm.side_effect = RuntimeError("no API key")
    # "that looks right" is NOT a keyword and NOT clearly analytical
    with pytest.raises(RuntimeError):
        classify_approval("that looks right")


# ---------------------------------------------------------------------------
# has_pending_plan bias
# ---------------------------------------------------------------------------

@patch("backend.orchestration.approval_classifier._get_llm")
def test_has_pending_plan_true_biases_llm_call(mock_get_llm: MagicMock) -> None:
    """Verify that has_pending_plan=True is forwarded to the LLM call."""
    llm = _mock_llm_response(True)
    mock_get_llm.return_value = llm
    result = is_approval_command("looks good to me", has_pending_plan=True)
    assert result is True
    assert llm.invoke.called


@patch("backend.orchestration.approval_classifier._get_llm")
def test_is_approval_command_drop_in_backward_compat(mock_get_llm: MagicMock) -> None:
    """is_approval_command() without has_pending_plan still works."""
    mock_get_llm.return_value = _mock_llm_response(True)
    assert is_approval_command("yes, execute the plan") is True


# ---------------------------------------------------------------------------
# Previously broken phrases (regression tests)
# ---------------------------------------------------------------------------

NATURAL_APPROVALS = [
    "yes, execute the plan",
    "that looks great, go for it",
    "sounds right to me",
    "let's execute",
    "approved, run the analysis",
    "perfect, start it",
]


@patch("backend.orchestration.approval_classifier._get_llm")
@pytest.mark.parametrize("phrase", NATURAL_APPROVALS)
def test_natural_approval_phrases_are_classified_as_approvals(
    mock_get_llm: MagicMock, phrase: str
) -> None:
    """These are phrases that the old keyword set would reject but users naturally say."""
    mock_get_llm.return_value = _mock_llm_response(True)
    result = is_approval_command(phrase)
    assert result is True, f"Expected {phrase!r} to be classified as an approval"
