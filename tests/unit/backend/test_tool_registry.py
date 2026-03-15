import asyncio

from backend.orchestration.tool_registry import dispatch_via_registry


def test_tool_registry_handles_needs_inputs_gate():
    async def _run():
        handled, result = await dispatch_via_registry(
            "bulk_rnaseq_analysis",
            {"needs_inputs": True},
            needs_inputs_builder=lambda tool, args: {"status": "needs_inputs", "tool": tool, "args": args},
        )
        assert handled is True
        assert result["status"] == "needs_inputs"
        assert result["tool"] == "bulk_rnaseq_analysis"

    asyncio.run(_run())

