#!/usr/bin/env python3
"""
Demo: proof that handle_natural_command and multi-step flow work as documented.

Uses the MSA + consensus prompt. Run from repo root:
  python -m tests.demo_scenarios.demo_handle_natural_command_flow
  or
  pytest tests/demo_scenarios/demo_handle_natural_command_flow.py -v -s

Step 6 executes the 2-step plan through the ExecutionBroker (mock tool executor)
so we prove the pipeline runs end-to-end, not just routing.
"""
from __future__ import annotations

import asyncio
import os
import sys

# Repo root on path
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

# Same prompt as in conversation (MSA + consensus with FASTA)
MSA_CONSENSUS_PROMPT = """perform a multi sequence alignment and calculate the consensus sequence of the following sequences: 
>seq1 
GGACGGTTACCGGTCAACCGGAAGTTCGCTGTACTACATATATATCTGGATACTCGGTCGATGAATTCGATAGCGATCGTCGTGGTCAAAACTTCTTT 
>seq2 
TCGTGTGGTGGATTTGACTGAAGGCCATTCTAGTGTCCTTTAAGCCGTGACCTTCGCAACCGGCCATCCACTCACGTTCGATTAGTCGCCCGCCCAGCAGGTA 
>seq3 
ACAGCAGAGGTGTTATTCCCTATGGTGCCCTTTTTATAGGTCATCTTGACGGCTCTAGACATAATCTGACAGCGATGCTTACAGGACTAGTCACGATACGGTGGT"""


def _looks_like_workflow(cmd: str) -> bool:
    """Replicate main_with_mcp logic so we can show why single-sentence prompt does not get a plan."""
    c = (cmd or "").lower()
    if any(tok in c for tok in ["apply code patch", "apply patch:", "replace script with", "replace the script with"]) or "```" in c:
        return False
    return any(tok in c for tok in [" and then ", " then ", "->", "→", "\n", ";"]) and len(c) > 20


def run_demo():
    from backend.command_router import CommandRouter

    session_context = {"session_id": "demo-session"}
    router = CommandRouter()

    print("=" * 70)
    print("PROOF: handle_natural_command and multi-step flow")
    print("=" * 70)
    print()

    # --- 1. Router on the exact prompt (keyword path: no "then", so no workflow split) ---
    print("1) ROUTER on the exact prompt (single sentence, no 'then')")
    print("-" * 60)
    # Force keyword path for deterministic proof (no LLM)
    os.environ["HELIX_USE_LLM_ROUTER"] = "0"
    tool_name, params = router.route_command(MSA_CONSENSUS_PROMPT, session_context)
    print(f"   route_command(...) -> tool_name = {repr(tool_name)}")
    print(f"   params keys: {list(params.keys())}")
    if params.get("router_reasoning"):
        print(f"   router_reasoning: {params['router_reasoning']}")
    print()
    assert tool_name == "handle_natural_command", f"Expected handle_natural_command, got {tool_name}"

    # --- 2. Workflow detection ---
    print("2) WORKFLOW DETECTION (_looks_like_workflow)")
    print("-" * 60)
    looks = _looks_like_workflow(MSA_CONSENSUS_PROMPT)
    print(f"   _looks_like_workflow(prompt) = {looks}")
    if looks:
        print("   (Prompt contains newlines in FASTA, so workflow delimiter triggers. See step 4 for clean 2-step.)")
    else:
        print("   (No ' then ' / newline / ; -> no plan from router; agent gets full command.)")
    print()

    # --- 3. Execute path: when router returns handle_natural_command ---
    print("3) EXECUTE PATH WHEN ROUTER RETURNS handle_natural_command")
    print("-" * 60)
    print("   Phase2c: router returned handle_natural_command -> skip direct tool call.")
    if looks:
        print("   _looks_like_workflow is True -> route_plan() would be used (splits only on 'then'/';'/->).")
    else:
        print("   _looks_like_workflow is False -> handle_command(full_command) (agent plans internally).")
    print("   Either way, no single toolbox tool is forced; agent/plan handles multi-step.")
    print()

    # --- 4. Exact prompt: route_plan no longer splits on FASTA newlines (prompt normalization) ---
    print("4) PLAN FROM EXACT PROMPT (no split on newlines; 1 step)")
    print("-" * 60)
    plan_exact = router.route_plan(MSA_CONSENSUS_PROMPT, session_context)
    steps_exact = plan_exact.get("steps", [])
    print(f"   route_plan(exact_prompt) -> {len(steps_exact)} step(s) (only 'then'/';' split; FASTA preserved)")
    for i, step in enumerate(steps_exact[:5], 1):
        print(f"   Step {i}: tool_name={step.get('tool_name')}, desc={str(step.get('description', ''))[:55]}...")
    if len(steps_exact) > 5:
        print(f"   ... and {len(steps_exact) - 5} more steps")
    assert len(steps_exact) == 1, "Exact prompt should be 1 step after prompt-normalization (no newline split)"
    print()

    # --- 5. Same intent WITH " then " -> clean 2-step plan (toolbox + handle_natural_command) ---
    prompt_with_then = (
        "perform a multi sequence alignment of the following sequences: >seq1 ACGT >seq2 TGCA "
        "then calculate the consensus sequence"
    )
    print("5) SAME INTENT WITH ' then ' -> CLEAN 2-STEP PLAN")
    print("-" * 60)
    looks_then = _looks_like_workflow(prompt_with_then)
    print(f"   _looks_like_workflow(prompt_with_then) = {looks_then}")
    plan = router.route_plan(prompt_with_then, session_context)
    print(f"   route_plan(...) -> plan version = {plan.get('version')}, steps = {len(plan.get('steps', []))}")
    for i, step in enumerate(plan.get("steps", []), 1):
        print(f"   Step {i}: id={step.get('id')}, tool_name={step.get('tool_name')}")
        print(f"           description = {str(step.get('description', ''))[:60]}...")
    print()
    assert looks_then is True
    assert len(plan["steps"]) == 2
    # First part -> sequence_alignment; second part "calculate the consensus sequence" -> handle_natural_command
    step1_tool = plan["steps"][0].get("tool_name")
    step2_tool = plan["steps"][1].get("tool_name")
    assert step1_tool == "sequence_alignment", f"Step 1 (align) should be sequence_alignment, got {step1_tool}"
    assert step2_tool == "handle_natural_command", f"Step 2 (consensus) should be handle_natural_command, got {step2_tool}"
    print("   -> Step 1 = toolbox (sequence_alignment), Step 2 = handle_natural_command (agent/code-gen).")
    print()

    # --- 6. Execute the 2-step plan through the broker (mock executor) ---
    print("6) EXECUTE THE 2-STEP PLAN (ExecutionBroker + mock tool executor)")
    print("-" * 60)
    ran = asyncio.run(_execute_plan_demo(plan))
    print(f"   Broker executed plan: {len(ran)} steps ran.")
    print(f"   Step results: {list(ran.keys())}")
    print()
    assert ran.get("step1", {}).get("tool_name") == "sequence_alignment", ran
    assert ran.get("step2", {}).get("tool_name") == "handle_natural_command", ran
    print("   -> Pipeline execution verified (mock executor; no real LLM/tools).")
    print()

    print("=" * 70)
    print("All assertions passed. The system behaves as documented.")
    print("=" * 70)


async def _execute_plan_demo(plan: dict) -> dict:
    """Run the 2-step plan through ExecutionBroker with a mock tool executor."""
    from backend.execution_broker import ExecutionBroker, ExecutionRequest

    async def mock_tool_executor(tool_name: str, args: dict):
        if tool_name == "sequence_alignment":
            return {"status": "success", "aligned_sequences": ">aln\nACGT\nTGCA\n"}
        if tool_name == "handle_natural_command":
            return {"status": "success", "text": "Consensus computed (mock)."}
        return {"status": "success"}

    broker = ExecutionBroker(tool_executor=mock_tool_executor)
    out = await broker.execute_tool(
        ExecutionRequest(
            tool_name="__plan__",
            arguments={"plan": plan, "session_id": "demo-session"},
            session_id="demo-session",
            original_command="align then consensus",
            session_context={},
        )
    )
    steps_by_id = {}
    if out.get("result", {}).get("type") == "plan_result":
        for s in out["result"].get("steps", []):
            steps_by_id[s["id"]] = {
                "tool_name": s.get("tool_name"),
                "status": s.get("result", {}).get("status"),
            }
    return steps_by_id


def test_demo_handle_natural_command_flow():
    """Pytest entry point: run_demo() must complete without assertion errors."""
    run_demo()


if __name__ == "__main__":
    run_demo()
