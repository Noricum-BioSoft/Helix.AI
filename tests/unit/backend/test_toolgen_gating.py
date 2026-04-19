import pytest


@pytest.mark.asyncio
async def test_unknown_tool_does_not_trigger_toolgen_for_qa(monkeypatch):
    """
    Phase 4: ensure tool-generator-agent is only used when intent is execute.
    """
    from backend.main import dispatch_tool

    called = {"toolgen": False}

    async def fake_toolgen(*args, **kwargs):
        called["toolgen"] = True
        return {"status": "success"}

    # Patch the tool generator symbol that dispatch_tool imports.
    import backend.tool_generator_agent as tga

    monkeypatch.setattr(tga, "generate_and_execute_tool", fake_toolgen)

    out = await dispatch_tool("some_unknown_tool", {"command": "What is bioinformatics?"})
    assert out.get("tool_generated") is False
    assert called["toolgen"] is False


