import pytest

from backend.command_router import CommandRouter
from tests.evals.eval_utils import assert_subset, load_cases


@pytest.mark.parametrize(
    "case",
    load_cases("tests/evals/cases/router_tool_mapping.jsonl"),
    ids=lambda c: c.id,
)
def test_router_tool_mapping_eval(case):
    router = CommandRouter()
    command = case.input.get("command", "")
    session_context = case.input.get("session_context", {}) or {}

    # Plan evals (multi-step workflows)
    if "plan_tool_names" in case.expect:
        plan = router.route_plan(command, session_context)
        tool_names = [s.get("tool_name") for s in (plan or {}).get("steps", [])]
        assert tool_names == case.expect["plan_tool_names"]
        return

    tool_name, params = router.route_command(command, session_context)
    assert tool_name == case.expect["tool_name"]

    expected_subset = case.expect.get("parameters_subset") or {}
    if expected_subset:
        assert_subset(expected_subset, params or {})





