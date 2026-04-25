import pytest


def test_command_router_route_plan_emits_plan_ir(monkeypatch):
    # Enable keyword routing so route_plan can decompose the multi-step command.
    # Production uses HELIX_LLM_ROUTER_FIRST=1; this test opts into keyword mode.
    monkeypatch.setenv("HELIX_LLM_ROUTER_FIRST", "0")
    from backend.command_router import CommandRouter

    router = CommandRouter()
    plan = router.route_plan("mutate sequence ATGC then align sequences", session_context={})
    assert plan["version"] == "v1"
    assert len(plan["steps"]) == 2
    assert plan["steps"][0]["id"] == "step1"
    assert "tool_name" in plan["steps"][0]


@pytest.mark.asyncio
async def test_broker_executes_plan_sync_with_refs(monkeypatch):
    from backend.execution_broker import ExecutionBroker, ExecutionRequest

    async def exec_tool(tool_name: str, args: dict):
        if tool_name == "mutate_sequence":
            return {"status": "success", "value": "ATGC_MUT"}
        if tool_name == "sequence_alignment":
            # should receive resolved arg
            return {"status": "success", "aligned": args.get("sequences")}
        return {"status": "success"}

    broker = ExecutionBroker(tool_executor=exec_tool)

    plan = {
        "version": "v1",
        "steps": [
            {"id": "step1", "tool_name": "mutate_sequence", "arguments": {"sequence": "ATGC", "num_variants": 1}},
            {
                "id": "step2",
                "tool_name": "sequence_alignment",
                "arguments": {"sequences": {"$ref": "steps.step1.result.value"}},
            },
        ],
    }

    out = await broker.execute_tool(
        ExecutionRequest(
            tool_name="__plan__",
            arguments={"plan": plan, "session_id": "s"},
            session_id="s",
            original_command="x",
            session_context={},
        )
    )

    assert out["mode"] == "sync"
    assert out["result"]["type"] == "plan_result"
    assert out["result"]["steps"][-1]["result"]["aligned"] == "ATGC_MUT"


@pytest.mark.asyncio
async def test_broker_submits_plan_async(monkeypatch):
    from backend.execution_broker import ExecutionBroker, ExecutionRequest

    monkeypatch.setenv("HELIX_ASYNC_BYTES_THRESHOLD", "1")

    async def exec_tool(tool_name: str, args: dict):
        return {"status": "success"}

    broker = ExecutionBroker(tool_executor=exec_tool)

    async def fake_submit_plan(req, plan):
        return {"type": "job", "status": "submitted", "job_id": "job-plan"}

    monkeypatch.setattr(broker, "_submit_plan_emr_job", fake_submit_plan)

    plan = {"version": "v1", "steps": [{"id": "step1", "tool_name": "mutate_sequence", "arguments": {"path": "/tmp/x"}}]}

    # Mock infrastructure decision agent import to raise exception so it falls back to threshold-based routing
    original_import = __import__
    def mock_import(name, *args, **kwargs):
        if name == "backend.infrastructure_decision_agent":
            raise ImportError("Infrastructure agent disabled for test")
        return original_import(name, *args, **kwargs)
    
    monkeypatch.setattr("builtins.__import__", mock_import)

    # Force size > threshold so routing becomes async (threshold is 10 by default, or from env)
    monkeypatch.setenv("HELIX_ASYNC_BYTES_THRESHOLD", "1")  # Set low threshold
    monkeypatch.setattr(broker, "_try_get_local_size", lambda p: 2)  # Return size > threshold

    out = await broker.execute_tool(
        ExecutionRequest(
            tool_name="__plan__",
            arguments={"plan": plan, "session_id": "s"},
            session_id="s",
            original_command="x",
            session_context={},
        )
    )

    assert out["mode"] == "async"
    assert out["result"]["job_id"] == "job-plan"






