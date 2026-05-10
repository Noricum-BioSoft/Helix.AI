"""
Unit tests for the tabular analysis planning gate (is_analytical_request).

Key invariant: meta/conversational/comparative messages must NOT be classified
as analytical execution requests, regardless of message length.  Only explicit
analytical cue words should trigger planning.
"""

import pytest
from backend.tabular_qa.analysis_planner import is_analytical_request


# --------------------------------------------------------------------------
# Shared session context fixture with one tabular upload
# --------------------------------------------------------------------------

TABULAR_SESSION = {
    "uploaded_files": [
        {
            "name": "data.xlsx",
            "schema_preview": {"family": "tabular"},
        }
    ]
}

NO_UPLOAD_SESSION: dict = {"uploaded_files": []}


# --------------------------------------------------------------------------
# Should NOT trigger the planner gate (fall through to LLM router)
# --------------------------------------------------------------------------

@pytest.mark.parametrize("command", [
    # The exact failure case from session 2af19038
    "what's the difference between your next suggested step and the previous analysis?",
    # General meta / clarification turns
    "what do you mean by that?",
    "can you explain what the next step involves?",
    "how does that differ from what you just did?",
    "why did you suggest that particular approach?",
    "which of those options would you recommend and why?",
    "what would happen if I did something different?",
    "I'm not sure I understand the goal of that step",
    "can you tell me more about the suggested next steps?",
    # Follow-up / suggestion requests
    "what should I do next?",
    "what are my options from here?",
    "give me suggestions for the next step",
    "what else can we learn from this dataset?",
    # Simple factual questions (already excluded by simple_patterns but reconfirm)
    "how many rows are there?",
    "what are the column names?",
    # Long but non-analytical
    "could you please tell me a bit more about what the previous analysis actually found in terms of the evidence categories?",
])
def test_meta_and_conversational_turns_do_not_trigger_planner(command: str) -> None:
    assert not is_analytical_request(command, TABULAR_SESSION), (
        f"Expected meta/conversational command NOT to trigger planning gate, but it did: {command!r}"
    )


# --------------------------------------------------------------------------
# SHOULD trigger the planner gate (clear analytical cue words present)
# --------------------------------------------------------------------------

@pytest.mark.parametrize("command", [
    "compare expression levels across the treatment and control groups",
    "cluster samples by their gene expression profiles",
    "identify the top ranked genes by log2 fold change",
    "run a correlation analysis between all numeric columns",
    "visualize the distribution of evidence categories as a bar chart",
    "group by cancer type and compute the mean evidence score",
    "analyze differential expression between responders and non-responders",
    "can you identify biomarkers associated with treatment response?",
    "build a heatmap of the top 50 differentially expressed genes",
    "stratify patients by age group and compare survival outcomes",
])
def test_explicit_analytical_cues_trigger_planner(command: str) -> None:
    assert is_analytical_request(command, TABULAR_SESSION), (
        f"Expected analytical command to trigger planning gate, but it did not: {command!r}"
    )


# --------------------------------------------------------------------------
# No tabular upload → always returns False regardless of content
# --------------------------------------------------------------------------

def test_no_tabular_upload_always_false() -> None:
    analytical = "cluster samples by gene expression and build a heatmap"
    assert not is_analytical_request(analytical, NO_UPLOAD_SESSION)
    assert not is_analytical_request(analytical, {})
