# /backend/agent.py

import asyncio

from typing import Dict

from langgraph.prebuilt import create_react_agent
from langgraph.checkpoint.memory import MemorySaver

from langchain.chat_models import init_chat_model
from langchain_core.tools import tool

from pydantic import BaseModel, Field

# load the environment
from dotenv import load_dotenv
load_dotenv()

from langchain_core.prompts import PromptTemplate
from langchain.globals import set_verbose, set_debug

set_verbose(True)
set_debug(True)

PROMPT_TEMPLATE = """Execute the following command. You have access to the following tools:

{tools}

Use the following format:

Command: the input command you must execute
Thought: you should always think about what to do
Action: the action to take, should be one of [{tool_names}]
Action Input: the input to the action
Observation: the result of the action

Thought: I now know the final answer
Final Answer: the final answer to the original input question

Begin!

Command: {input}
Thought: {agent_scratchpad}"""


class OutputFormatter(BaseModel):
    """Always use this tool to structure your response to the user."""

    input: str = Field(description="The input of the tool/function call")
    output: str = Field(description="The output of the tool/function call")
    plot: str = Field(description="An optional plot visualizing the output")


@tool
def sequence_alignment(sequences: str) -> Dict:
    """performs a sequence alignment on a given set of sequences."""

    seq_align_dict = [
        {"name": "seq1", "sequence": "ACTG--TTGAC"},
        {"name": "seq2", "sequence": "ACTGCA-T--C"},
        {"name": "seq3", "sequence": "ACTGCAATGAC"},
    ]

    from pymsaviz import MsaViz

    with open("msa_file.fa", "w") as fp:
        fp.write(">seq01\nACTG--TTGAC")
        fp.write(">seq02\nACTGCA-T--C")
        fp.write(">seq03\nACTGCAATGAC")

    with open("msa_file.fa", "r") as fp:
        msa_data = fp. readlines()

    return {
        "text": "Sequences aligned successfully.",
        "input": sequences,
        "output": msa_data,
        "plot": {
            "data": [{"x": [1, 2, 3], "y": [3, 3, 3], "type": "bar"}],
            "layout": {"title": "Alignment Visualization"},
        },
    }


@tool
def mutate_sequence(sequence: str, num_variants: int = 96) -> Dict:
    """mutates a given sequence and returns 96 variants per default."""

    seq_mut_dict = [
        {"name": "seq01", "sequence": "ACTG"},
        {"name": "seq02", "sequence": "ACGA"},
        {"name": "seq03", "sequence": "ACGT"},
        {"name": "seq04", "sequence": "ACGC"},
        {"name": "seq05", "sequence": "ACTG"},
        {"name": "seq06", "sequence": "TCGA"},
        {"name": "seq07", "sequence": "CCGT"},
        {"name": "seq08", "sequence": "GCGC"},
    ]

    return {
        "text": "Sequence mutated successfully.",
        "input": {
            "sequence": sequence,
            "variants": num_variants,
        },
        "output": seq_mut_dict,
        "plot": {
            "data": [{"x": [1, 2, 3], "y": [1, 1, 1], "type": "bar"}],
            "layout": {"title": "Mutation Visualization"},
        },
    }


# llm = init_chat_model("deepseek:deepseek-chat", temperature=0)

from langchain_deepseek import ChatDeepSeek

llm = ChatDeepSeek(
    model="deepseek-chat",
    temperature=0,
    max_tokens=None,
    timeout=None,
    max_retries=2,
)

memory = MemorySaver()
model = init_chat_model("openai:gpt-4o")
config = {"configurable": {"thread_id": "abc123"}}
agent = create_react_agent(
    model=model,
    tools=[sequence_alignment, mutate_sequence],
    checkpointer=memory,
    # prompt=PromptTemplate.from_template(PROMPT_TEMPLATE),
)


async def handle_command(command: str):
    print(f"[handle_command] command: {command}")

    input_message = {
        "role": "user",
        "content": command,
    }

    result = agent.invoke({"messages": [input_message]}, config)
    print(f"[handle_command] result: {result}")

    return result
    #
    # from pymsaviz import MsaViz, get_msa_testdata
    #
    # msa_file = get_msa_testdata("HIGD2A.fa")
    # mv = MsaViz(msa_file, wrap_length=60, show_count=True)
    #
    # fig = mv.plotfig()
    #
    # # mv.savefig("api_example01.png")


    return {
        "text": result.tool_input,
        "plot": fig
    }


if __name__ == "__main__":
    seq_align_results = asyncio.run(handle_command("align the given sequences"))
    print(seq_align_results)

    mut_results = asyncio.run(handle_command("mutate the given sequence."))
    print(mut_results)

