import pytest


@pytest.mark.asyncio
async def test_unknown_tool_does_not_trigger_toolgen_for_qa(monkeypatch):
    """Phase 4: ensure tool-generator-agent is only used when intent is execute.

    The intent classifier is mocked to return Q&A intent so the test doesn't
    depend on a live LLM.
    """
    from backend.main import dispatch_tool
    from backend.intent_classifier import IntentDecision

    called = {"toolgen": False}

    async def fake_toolgen(*args, **kwargs):
        called["toolgen"] = True
        return {"status": "success"}

    # Mock intent classifier to return Q&A (no LLM needed in test).
    monkeypatch.setattr(
        "backend.intent_classifier.classify_intent",
        lambda *a, **kw: IntentDecision(intent="qa", reason="mocked_qa"),
    )

    # Patch the tool generator symbol that dispatch_tool imports.
    import backend.tool_generator_agent as tga

    monkeypatch.setattr(tga, "generate_and_execute_tool", fake_toolgen)

    out = await dispatch_tool("some_unknown_tool", {"command": "What is bioinformatics?"})
    assert out.get("tool_generated") is False
    assert called["toolgen"] is False


