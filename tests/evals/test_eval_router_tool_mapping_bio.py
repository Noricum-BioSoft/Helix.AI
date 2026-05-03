import pytest

from backend.command_router import CommandRouter
from tests.evals.eval_utils import assert_subset, load_cases


@pytest.mark.parametrize(
    "case",
    load_cases("tests/evals/cases/bioinformatics_router_tool_mapping.jsonl"),
    ids=lambda c: c.id,
)
def test_router_tool_mapping_bio_eval(case):
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

    # Tests can accept either a single expected tool or any-of-N for cases
    # where multiple routings are semantically reasonable (e.g. a phylogenetic
    # request with raw accessions could go to phylogenetic_tree OR
    # fetch_ncbi_sequence-first multi-step).
    expected_tool = case.expect.get("tool_name")
    expected_oneof = case.expect.get("tool_name_oneof")
    if expected_oneof:
        assert tool_name in expected_oneof, (
            f"Expected tool in {expected_oneof}, got {tool_name!r}"
        )
    else:
        assert tool_name == expected_tool, (
            f"Expected tool {expected_tool!r}, got {tool_name!r}"
        )

    expected_subset = case.expect.get("parameters_subset") or {}
    if expected_subset:
        assert_subset(expected_subset, params or {})

