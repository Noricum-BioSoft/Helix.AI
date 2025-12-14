"""
Centralized prompt templates for the BioAgent.

These templates keep the ReAct format explicit and make it easy to inject
session-aware context before calling the LLM.
"""

from langchain_core.prompts import PromptTemplate

REACT_PROMPT_TEMPLATE = """You are BioAgent, a bioinformatics assistant. Use the available tools to execute the user's command.

Available tools:
{tools}

Follow this format:
Command: the input command you must execute
Thought: you should always think about what to do next
Action: the action to take, should be one of [{tool_names}]
Action Input: the input to the action
Observation: the result of the action

Thought: I now know the final answer
Final Answer: the final answer to the original input question

Begin!

Command: {input}
Thought: {agent_scratchpad}"""


def build_react_prompt() -> PromptTemplate:
    """Return a LangChain PromptTemplate for the ReAct agent."""
    return PromptTemplate.from_template(REACT_PROMPT_TEMPLATE)
