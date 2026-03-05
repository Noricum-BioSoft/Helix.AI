import os
import pytest
import asyncio
from backend.agent import handle_command

pytestmark = pytest.mark.integration


@pytest.mark.asyncio
@pytest.mark.parametrize("command, expected_text, expected_plot_title", [
    ("Align the sequences", "Sequences aligned successfully.", "Alignment Visualization"),
    ("Mutate a sequence and generate 96 variants", "Generated 96 mutated variants.", "Mutation Variants"),
    ("Build a phylogenetic tree", "Unknown command.", None)
])
async def test_handle_command(command:str, expected_text:str, expected_plot_title:str):

    if os.getenv("HELIX_MOCK_MODE") == "1" or not os.getenv("OPENAI_API_KEY"):
        pytest.skip("Requires live OpenAI access (integration-only).")

    result = await handle_command(command)

    assert "text" in result
    assert result["text"] == expected_text

    if expected_plot_title:
        assert "plot" in result
        assert "layout" in result["plot"]
        assert result["plot"]["layout"]["title"] == expected_plot_title
    else:
        assert "plot" not in result


def test_run_agent():
    if os.getenv("HELIX_MOCK_MODE") == "1" or not os.getenv("OPENAI_API_KEY"):
        pytest.skip("Requires live OpenAI access (integration-only).")

    # Import relevant functionality
    from langchain.chat_models import init_chat_model
    from langgraph.checkpoint.memory import MemorySaver
    from langgraph.prebuilt import create_react_agent
    from backend.agent import sequence_alignment

    # Create the agent
    memory = MemorySaver()
    model = init_chat_model("openai:gpt-4o")
    tools = [sequence_alignment]
    agent_executor = create_react_agent(model, tools, checkpointer=memory)

    # Use the agent
    config = {"configurable": {"thread_id": "abc123"}}

    input_message = {
        "role": "user",
        "content": "Align the following sequences: ACGT AACGT",
    }
    for step in agent_executor.stream(
            {"messages": [input_message]}, config, stream_mode="values"
    ):
        step["messages"][-1].pretty_print()