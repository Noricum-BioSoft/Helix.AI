"""
Smoke tests — LLM approval classifier.

Verifies that the live LLM correctly distinguishes approval phrases from
non-approval phrases.  The classifier must never silently fall back to a
keyword list — if it misclassifies, that is a real regression.
"""
from __future__ import annotations

import pytest

# ---------------------------------------------------------------------------
# Cases: (phrase, should_approve)
# ---------------------------------------------------------------------------

APPROVAL_CASES = [
    # Clear approvals
    ("Yes, go ahead.", True),
    ("Approve", True),
    ("proceed", True),
    ("looks good, execute it", True),
    ("run it", True),
    ("yes please", True),
    ("confirmed", True),
    ("that looks correct, proceed", True),
    # Clear rejections / questions / modifications
    ("No, cancel.", False),
    ("Wait, change the alpha to 0.01 first", False),
    ("What does step 2 do exactly?", False),
    ("Can you explain why you chose DESeq2?", False),
    ("Stop", False),
    ("Let me think about it", False),
    ("Actually, use the other file instead", False),
    # Edge cases that should NOT be approvals
    ("Run the analysis on my data", False),
    ("Execute a bulk RNA-seq on counts.csv", False),
]


@pytest.fixture(scope="module")
def classifier():
    from backend.orchestration.approval_classifier import classify_approval
    return classify_approval


@pytest.mark.parametrize("phrase,expected_approval", APPROVAL_CASES)
def test_approval_classifier(classifier, phrase, expected_approval):
    """LLM classifier correctly approves or rejects each phrase."""
    result = classifier(phrase, has_pending_plan=True)
    assert result.is_approval == expected_approval, (
        f"Phrase: {phrase!r}\n"
        f"Expected is_approval={expected_approval}, got is_approval={result.is_approval}\n"
        f"Reason: {result.reason}"
    )
