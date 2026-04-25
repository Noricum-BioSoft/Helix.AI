"""Unit tests for the LLM-based intent classifier.

All tests mock the LLM — there are no heuristic fallbacks to rely on.
Tests verify that:
1. Session-aware shortcuts work without any LLM call.
2. The LLM is called for unambiguous classification tasks.
3. LLM responses are correctly mapped to execute/qa decisions.
4. LLM failure raises an exception (no silent fallback).
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from backend.intent_classifier import classify_intent, IntentDecision


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _mock_llm(intent_labels: list[str]) -> MagicMock:
    """Return a mock LLM that responds with the given intent labels."""
    import json
    llm = MagicMock()
    llm.invoke.return_value = MagicMock(
        content=json.dumps({"prompt": "test", "intent": intent_labels})
    )
    return llm


# ---------------------------------------------------------------------------
# Session-aware state-machine shortcuts (no LLM needed)
# ---------------------------------------------------------------------------

def test_waiting_for_approval_state_short_circuits_to_execute():
    """WAITING_FOR_APPROVAL state → always execute, no LLM call."""
    with patch("backend.intent_classifier._get_llm") as mock_get_llm:
        mock_get_llm.side_effect = AssertionError("LLM should not be called for state shortcut")
        decision = classify_intent("ok", workflow_state="WAITING_FOR_APPROVAL")
    assert decision.intent == "execute"
    assert "checkpoint_state" in decision.reason


def test_waiting_for_clarification_state_short_circuits_to_execute():
    with patch("backend.intent_classifier._get_llm") as mock_get_llm:
        mock_get_llm.side_effect = AssertionError("LLM should not be called for state shortcut")
        decision = classify_intent("yes", workflow_state="WAITING_FOR_CLARIFICATION")
    assert decision.intent == "execute"


def test_waiting_for_inputs_state_short_circuits_to_execute():
    with patch("backend.intent_classifier._get_llm") as mock_get_llm:
        mock_get_llm.side_effect = AssertionError("LLM should not be called for state shortcut")
        decision = classify_intent("here it is", workflow_state="WAITING_FOR_INPUTS")
    assert decision.intent == "execute"


def test_rerun_cue_in_session_with_history_short_circuits_to_execute():
    """Prior runs + rerun cue → execute without LLM."""
    session = {"runs": [{"run_id": "r1", "tool": "bulk_rnaseq_analysis"}]}
    with patch("backend.intent_classifier._get_llm") as mock_get_llm:
        mock_get_llm.side_effect = AssertionError("LLM should not be called for rerun shortcut")
        decision = classify_intent("run it again", session_context=session)
    assert decision.intent == "execute"
    assert "rerun" in decision.reason


# ---------------------------------------------------------------------------
# LLM-driven classification
# ---------------------------------------------------------------------------

@patch("backend.intent_classifier._get_llm")
def test_llm_action_label_maps_to_execute(mock_get_llm):
    mock_get_llm.return_value = _mock_llm(["action"])
    decision = classify_intent("Run FastQC on my reads")
    assert decision.intent == "execute"
    assert "llm_classified_action" in decision.reason
    mock_get_llm.assert_called()


@patch("backend.intent_classifier._get_llm")
def test_llm_question_label_maps_to_qa(mock_get_llm):
    mock_get_llm.return_value = _mock_llm(["question"])
    decision = classify_intent("What is a FASTQ file?")
    assert decision.intent == "qa"
    assert "llm_classified_question" in decision.reason


@patch("backend.intent_classifier._get_llm")
def test_llm_action_overrides_question_label(mock_get_llm):
    """When both action and question labels present, action wins."""
    mock_get_llm.return_value = _mock_llm(["question", "action", "data"])
    decision = classify_intent("Can you run FastQC on this data?")
    assert decision.intent == "execute"


@patch("backend.intent_classifier._get_llm")
def test_llm_s3_uri_command_classified_as_execute(mock_get_llm):
    mock_get_llm.return_value = _mock_llm(["action", "data"])
    decision = classify_intent("analyze s3://my-bucket/sample_R1.fastq and s3://my-bucket/sample_R2.fastq")
    assert decision.intent == "execute"


@patch("backend.intent_classifier._get_llm")
def test_llm_is_called_with_user_text(mock_get_llm):
    """Verify the user command is passed to the LLM invoke call."""
    llm = _mock_llm(["action"])
    mock_get_llm.return_value = llm
    classify_intent("align these reads")
    call_args = llm.invoke.call_args[0][0]  # messages list
    user_msg = next(m for m in call_args if m["role"] == "user")
    assert "align these reads" in user_msg["content"]


# ---------------------------------------------------------------------------
# No heuristic fallback — LLM failure raises
# ---------------------------------------------------------------------------

@patch("backend.intent_classifier._get_llm")
def test_llm_failure_raises_exception(mock_get_llm):
    """If the LLM is unavailable, an exception propagates — no silent fallback."""
    mock_get_llm.side_effect = RuntimeError("LLM disabled in HELIX_MOCK_MODE")
    with pytest.raises(RuntimeError, match="LLM disabled"):
        classify_intent("Run FastQC on my reads")
