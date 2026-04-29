"""
Smoke tests — LLM intent classifier.

Verifies that the live LLM correctly classifies user commands as either
'execute' (analytical action) or 'qa' (question / conversational).
"""
from __future__ import annotations

import pytest

INTENT_CASES = [
    # Execute intents
    ("Run differential expression on my RNA-seq data", "execute"),
    ("Analyze my uploaded CSV file", "execute"),
    ("Build a phylogenetic tree from sequences.fasta", "execute"),
    ("Trim adapters from my FASTQ reads", "execute"),
    ("Re-run with alpha=0.01", "execute"),
    # QA intents
    ("What is DESeq2?", "qa"),
    ("Why did you choose this normalization method?", "qa"),
    ("What does UMAP stand for?", "qa"),
    ("How many samples do I have?", "qa"),
    ("What tools are available?", "qa"),
]


@pytest.fixture(scope="module")
def classify():
    from backend.intent_classifier import classify_intent
    return classify_intent


@pytest.mark.parametrize("command,expected_intent", INTENT_CASES)
def test_intent_classifier(classify, command, expected_intent):
    """LLM intent classifier returns the correct intent category."""
    result = classify(command)
    assert result.intent == expected_intent, (
        f"Command: {command!r}\n"
        f"Expected intent={expected_intent!r}, got intent={result.intent!r}\n"
        f"Reason: {result.reason}"
    )
