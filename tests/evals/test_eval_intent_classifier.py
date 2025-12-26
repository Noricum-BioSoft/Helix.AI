import pytest

from backend.intent_classifier import classify_intent
from tests.evals.eval_utils import load_cases


@pytest.mark.parametrize(
    "case",
    load_cases("tests/evals/cases/intent_classification.jsonl"),
    ids=lambda c: c.id,
)
def test_intent_classifier_eval(case):
    text = case.input.get("text", "")
    decision = classify_intent(text)
    assert decision.intent == case.expect["intent"]





