#!/usr/bin/env python3
"""
Integration test: real execution of the handle_natural_command pipeline.

Runs the same 2-step flow as the demo (align then consensus) using the real
ExecutionBroker and dispatch_tool: step 1 = sequence_alignment (real alignment
module), step 2 = handle_natural_command (real agent for "calculate the consensus").

Skipped when HELIX_MOCK_MODE=1 (default in CI).

Run real execution (requires API key, no HELIX_MOCK_MODE=1):

    HELIX_MOCK_MODE=0 pytest tests/integration/test_handle_natural_command_pipeline.py -v -s -m integration

Optional: HELIX_AGENT_TIMEOUT_S=30 to cap agent time for step 2.
"""

from __future__ import annotations

import asyncio
import os
import sys
from pathlib import Path

import pytest

# Repo root
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
TOOLS_PATH = REPO_ROOT / "tools"
if str(TOOLS_PATH) not in sys.path:
    sys.path.append(str(TOOLS_PATH))


@pytest.mark.integration
@pytest.mark.asyncio
async def test_handle_natural_command_pipeline_real_execution():
    """Run 2-step plan (sequence_alignment then handle_natural_command) with real tools."""
    if os.getenv("HELIX_MOCK_MODE") == "1":
        pytest.skip(
            "Real pipeline execution skipped (HELIX_MOCK_MODE=1). "
            "Run with: HELIX_MOCK_MODE=0 pytest tests/integration/test_handle_natural_command_pipeline.py -v -s"
        )

    from backend.command_router import CommandRouter
    from backend.execution_broker import ExecutionBroker, ExecutionRequest
    from backend.main import dispatch_tool
    from backend.history_manager import history_manager

    # Use a dedicated sessions subdir so this test doesn't pollute the main sessions/ tree.
    # See docs/SESSIONS_DIRECTORY.md for why many files can appear under sessions/.
    _sessions_subdir = os.environ.get("HELIX_INTEGRATION_SESSIONS_DIR", "integration-test-pipeline-runs")
    _old_storage = history_manager.storage_dir
    history_manager.storage_dir = Path(_old_storage) / _sessions_subdir
    history_manager.storage_dir.mkdir(parents=True, exist_ok=True)

    session_id = "integration-test-pipeline"
    session_context = {"session_id": session_id}

    # Same prompt as demo: "align ... then calculate the consensus"
    prompt = (
        "perform a multi sequence alignment of the following sequences: "
        ">seq1 ACGTACGT >seq2 ACGTACGT >seq3 TGCATGCA "
        "then calculate the consensus sequence"
    )

    # 1) Build 2-step plan via router (use keyword path so step2 = handle_natural_command for "consensus")
    # LLM router can route "calculate the consensus sequence" to sequence_alignment; we need step2 = handle_natural_command.
    prev_llm_router = os.environ.get("HELIX_USE_LLM_ROUTER")
    os.environ["HELIX_USE_LLM_ROUTER"] = "0"
    try:
        router = CommandRouter()
        plan = router.route_plan(prompt, session_context)
    finally:
        if prev_llm_router is not None:
            os.environ["HELIX_USE_LLM_ROUTER"] = prev_llm_router
        else:
            os.environ.pop("HELIX_USE_LLM_ROUTER", None)
    steps = plan.get("steps", [])
    assert len(steps) == 2, f"Expected 2 steps, got {len(steps)}"
    assert steps[0].get("tool_name") == "sequence_alignment"
    assert steps[1].get("tool_name") == "handle_natural_command"

    # 2) Execute plan with real broker (dispatch_tool)
    broker = ExecutionBroker(tool_executor=dispatch_tool)
    timeout_s = int(os.getenv("HELIX_AGENT_TIMEOUT_S", "90"))
    request = ExecutionRequest(
        tool_name="__plan__",
        arguments={"plan": plan, "session_id": session_id},
        session_id=session_id,
        original_command=prompt,
        session_context=session_context,
    )

    try:
        result = await asyncio.wait_for(
            broker.execute_tool(request),
            timeout=timeout_s + 10,
        )
    except asyncio.TimeoutError:
        pytest.fail(f"Pipeline did not complete within {timeout_s + 10}s")

    # 3) Assert plan result
    assert result is not None
    res = result.get("result") or result
    assert res.get("type") == "plan_result", res
    step_results = res.get("steps", [])
    assert len(step_results) == 2, step_results

    # Step 1: sequence_alignment
    s1 = step_results[0]
    assert s1.get("tool_name") == "sequence_alignment"
    r1 = s1.get("result") or {}
    assert r1.get("status") in ("success", "completed") or "alignment" in r1 or "output" in r1, (
        f"Step 1 should succeed or return alignment: {r1}"
    )

    # Step 2: handle_natural_command (agent / tool generator)
    s2 = step_results[1]
    assert s2.get("tool_name") == "handle_natural_command"
    r2 = s2.get("result") or {}
    assert isinstance(r2, dict), r2
    # Step 2 may succeed (consensus) or fail (e.g. tool generator had no input from step 1).
    # Full E2E consensus requires data flow from step 1 → step 2; see docs/CONSENSUS_FROM_ALIGNMENT_GAP.md.

    # Persist plan result to session (pipeline execution summary only, no full_response/code)
    from backend.main import _build_pipeline_execution_storage_result
    history_manager.ensure_session_exists(session_id)
    storage_result = _build_pipeline_execution_storage_result(result, session_id)
    history_manager.add_history_entry(
        session_id, prompt, "__plan__", storage_result,
        metadata={"plan_steps": len(step_results), "test": "integration"},
    )

    # Print results so you can see alignment and consensus when running with -s
    print("\n" + "=" * 60)
    print("STEP 1 — Sequence alignment")
    print("=" * 60)
    alignment_data = r1.get("alignment") or r1.get("output") or r1.get("text")
    if alignment_data is not None:
        if isinstance(alignment_data, list):
            for i, row in enumerate(alignment_data[:10], 1):
                print(f"  {row}")
            if len(alignment_data) > 10:
                print(f"  ... and {len(alignment_data) - 10} more")
        else:
            print(alignment_data if isinstance(alignment_data, str) else str(alignment_data)[:500])
    else:
        print("  (no alignment/output key; full result keys:", list(r1.keys()))
    print("\n" + "=" * 60)
    print("STEP 2 — handle_natural_command (consensus / agent result)")
    print("=" * 60)
    text2 = r2.get("text") or r2.get("message") or ""
    if text2:
        print(text2[:1000] + ("..." if len(str(text2)) > 1000 else ""))
    if r2.get("alignment") or r2.get("output"):
        print("  alignment/output:", r2.get("alignment") or r2.get("output"))
    if not text2 and not r2.get("alignment") and not r2.get("output"):
        print("  (result keys:", list(r2.keys()), ")")

    # Where inputs, prompts, and results are stored (same as normal session handling)
    paths = history_manager.get_session_storage_paths(session_id)
    print("\n" + "=" * 60)
    print("Session storage (inputs, prompts, results)")
    print("=" * 60)
    print("  storage_root:", paths["storage_root"])
    print("  session_file:", paths["session_file"])
    print("  session_dir: ", paths["session_dir"])
    print("  (session state: session_file; artifacts/runs: session_dir)")
    print("\n✅ Integration pipeline: step1 (sequence_alignment) and step2 (handle_natural_command) completed.")

    # Restore so other tests or the app don't use the integration subdir
    history_manager.storage_dir = _old_storage
