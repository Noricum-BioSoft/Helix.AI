# /backend/agent.py

import asyncio

from typing import Dict, List

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
    """Performs a sequence alignment on a given set of sequences."""
    
    # Import the alignment function
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
    from alignment import run_alignment_tool
    
    result = run_alignment_tool(sequences)
    
    return {
        "text": result.get("text", "Sequences aligned successfully."),
        "input": sequences,
        "output": result.get("alignment", []),
        "statistics": result.get("statistics", {}),
        "plot": {
            "data": [{"x": [1, 2, 3], "y": [3, 3, 3], "type": "bar"}],
            "layout": {"title": "Alignment Visualization"},
        },
    }


@tool
def mutate_sequence(sequence: str, num_variants: int = 96) -> Dict:
    """Mutates a given sequence and returns the specified number of variants. Use this when you need to CREATE new sequence variants from an input sequence."""
    
    # Import the mutation function
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
    from mutations import run_mutation_raw
    
    result = run_mutation_raw(sequence, num_variants)
    
    return {
        "text": result.get("text", "Sequence mutated successfully."),
        "input": {
            "sequence": sequence,
            "variants": num_variants,
        },
        "output": result.get("statistics", {}),
        "plot": result.get("plot", {
            "data": [{"x": [1, 2, 3], "y": [1, 1, 1], "type": "bar"}],
            "layout": {"title": "Mutation Visualization"},
        }),
    }


@tool
def dna_vendor_research(command: str, sequence_length: int = None, quantity: str = "standard") -> Dict:
    """Research DNA synthesis vendors and testing options for experimental validation. Use this when the user wants to ORDER sequences, find VENDORS, or research TESTING options. Keywords: order, vendor, synthesis, test, assay, expression, function, binding."""
    
    # Import the research function
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
    from dna_vendor_research import run_dna_vendor_research_raw
    
    result = run_dna_vendor_research_raw(command, sequence_length, quantity)
    
    return {
        "text": f"DNA vendor research completed: {result.get('message', 'Research done')}",
        "input": {
            "command": command,
            "sequence_length": sequence_length,
            "quantity": quantity
        },
        "output": result,
        "plot": {
            "data": [{"x": ["Vendors", "Testing Options"], "y": [result.get('total_vendors', 0), result.get('total_testing_options', 0)], "type": "bar"}],
            "layout": {"title": "DNA Vendor Research Results"},
        },
    }


@tool
def phylogenetic_tree(aligned_sequences: str) -> Dict:
    """Create phylogenetic tree visualization from aligned sequences."""
    
    print(f"🔧 Agent phylogenetic_tree called with aligned_sequences: '{aligned_sequences}'")
    
    # Import the phylogenetic tree function
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
    from phylogenetic_tree import run_phylogenetic_tree
    
    result = run_phylogenetic_tree(aligned_sequences)
    
    print(f"🔧 Agent phylogenetic_tree result: {result}")
    
    return {
        "text": result.get("text", "Phylogenetic tree created successfully."),
        "input": aligned_sequences,
        "output": result.get("statistics", {}),
        "plot": result.get("plot", {
            "data": [{"x": [1, 2, 3], "y": [1, 1, 1], "type": "bar"}],
            "layout": {"title": "Phylogenetic Tree"},
        }),
    }


@tool
def sequence_selection(aligned_sequences: str, selection_type: str = "random", num_sequences: int = 1) -> Dict:
    """Select sequences from aligned sequences based on various criteria (random, best_conservation, lowest_gaps, highest_gc, longest, shortest)."""
    
    print(f"🔧 Agent sequence_selection called with:")
    print(f"🔧 aligned_sequences: '{aligned_sequences}'")
    print(f"🔧 selection_type: '{selection_type}'")
    print(f"🔧 num_sequences: {num_sequences}")
    
    # Import the sequence selection function
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
    from sequence_selection import run_sequence_selection_raw
    
    result = run_sequence_selection_raw(aligned_sequences, selection_type, num_sequences)
    
    print(f"🔧 Agent sequence_selection result: {result}")
    
    return {
        "text": result.get("text", "Sequence selection completed successfully."),
        "input": {
            "aligned_sequences": aligned_sequences,
            "selection_type": selection_type,
            "num_sequences": num_sequences
        },
        "output": result.get("selected_sequences", []),
        "selection_criteria": result.get("selection_criteria", {}),
    }


@tool
def synthesis_submission(sequences: str, vendor_preference: str = None, quantity: str = "standard", delivery_time: str = "standard") -> Dict:
    """Submit sequences for DNA synthesis and get pricing quote. Quantity options: standard, large, custom. Delivery options: rush, standard, economy."""
    
    # Import the synthesis submission function
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
    from synthesis_submission import run_synthesis_submission_raw
    
    result = run_synthesis_submission_raw(sequences, vendor_preference, quantity, delivery_time)
    
    return {
        "text": result.get("text", "Synthesis submission completed successfully."),
        "input": {
            "sequences": sequences,
            "vendor_preference": vendor_preference,
            "quantity": quantity,
            "delivery_time": delivery_time
        },
        "output": result.get("quote", {}),
        "validation_results": result.get("validation_results", {}),
    }


@tool
def create_session() -> Dict:
    """Create a new session for tracking user interactions. Use this when the user wants to start a new session or create a new experiment."""
    
    import uuid
    
    session_id = str(uuid.uuid4())
    
    return {
        "text": f"Session created successfully with ID: {session_id}",
        "input": {},
        "output": {"session_id": session_id},
        "plot": {}
    }


@tool
def plasmid_visualization(vector_name: str = None, cloning_sites: str = "", insert_sequence: str = "", 
                          full_plasmid_sequence: str = None, insert_position: int = None) -> Dict:
    """Generate plasmid visualization data.
    
    Supports two modes:
    1. Full plasmid: When full_plasmid_sequence is provided, visualize the complete plasmid
    2. Insert mode: When vector_name and insert_sequence are provided, create plasmid with insert
    """
    
    # Import the plasmid visualization function
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
    from plasmid_visualizer import run_plasmid_visualization_raw
    
    result = run_plasmid_visualization_raw(
        vector_name=vector_name,
        cloning_sites=cloning_sites,
        insert_sequence=insert_sequence,
        full_plasmid_sequence=full_plasmid_sequence,
        insert_position=insert_position
    )
    
    return {
        "text": result.get("text", "Plasmid visualization created successfully."),
        "input": {
            "vector_name": vector_name,
            "cloning_sites": cloning_sites,
            "insert_sequence": insert_sequence,
            "full_plasmid_sequence": full_plasmid_sequence,
            "insert_position": insert_position
        },
        "output": result.get("plasmid_data", {}),
        "visualization_type": result.get("visualization_type", "circular_plasmid"),
    }

@tool
def plasmid_for_representatives(representatives: List[str], aligned_sequences: str, vector_name: str = "pUC19", cloning_sites: str = "EcoRI, BamHI, HindIII") -> Dict:
    """Create plasmid visualizations for representative sequences from clustering analysis."""
    
    # Import the plasmid visualization function
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
    from plasmid_visualizer import create_plasmid_for_representatives
    
    result = create_plasmid_for_representatives(representatives, aligned_sequences, vector_name, cloning_sites)
    
    return {
        "text": result.get("text", "Plasmid visualizations for representatives created successfully."),
        "plasmid_results": result.get("plasmid_results", []),
        "total_representatives": result.get("total_representatives", 0),
        "vector_name": result.get("vector_name", vector_name),
        "cloning_sites": result.get("cloning_sites", cloning_sites),
    }


# llm = init_chat_model("deepseek:deepseek-chat", temperature=0)

import os
from langchain_deepseek import ChatDeepSeek

# Initialize LLM with proper error handling
try:
    if os.getenv("DEEPSEEK_API_KEY"):
        llm = ChatDeepSeek(
            model="deepseek-chat",
            temperature=0,
            max_tokens=None,
            timeout=None,
            max_retries=2,
        )
        print("✅ DeepSeek LLM initialized with API key")
    elif os.getenv("OPENAI_API_KEY"):
        # Fallback to OpenAI if DeepSeek API key is not available
        llm = init_chat_model("openai:gpt-4o", temperature=0)
        print("⚠️  DEEPSEEK_API_KEY not found, using OpenAI as fallback")
    else:
        raise ValueError("No API keys found. Please set either DEEPSEEK_API_KEY or OPENAI_API_KEY in your environment variables.")
except Exception as e:
    print(f"❌ Error initializing LLM: {e}")
    print("📋 Please check ENVIRONMENT_SETUP.md for configuration instructions")
    raise

memory = MemorySaver()
# Use the same LLM instance for consistency
model = llm
config = {"configurable": {"thread_id": "abc123"}}

agent = create_react_agent(
    model=model,
    tools=[sequence_alignment, mutate_sequence, dna_vendor_research, phylogenetic_tree, sequence_selection, synthesis_submission, plasmid_visualization, plasmid_for_representatives],
    checkpointer=memory,
)


async def handle_command(command: str, session_id: str = "default", session_context: Dict = None):
    print(f"[handle_command] command: {command}")
    print(f"[handle_command] session_id: {session_id}")

    # NOTE: Memory clearing is left disabled so the LangGraph MemorySaver can
    # retain session-scoped context. If token growth becomes an issue, consider
    # adding a configurable flag to reset selectively rather than wiping on
    # every request.

    vendor_response = _maybe_handle_vendor_request(command)
    if vendor_response is not None:
        return vendor_response

    # For other commands, use the agent with session-specific config
    input_message = {
        "role": "user",
        "content": command,
    }

    # Use session-specific config
    session_config = {"configurable": {"thread_id": session_id}}
    result = agent.invoke({"messages": [input_message]}, session_config)
    print(f"[handle_command] result: {result}")

    return result


if __name__ == "__main__":
    seq_align_results = asyncio.run(handle_command("align the given sequences"))
    print(seq_align_results)

    mut_results = asyncio.run(handle_command("mutate the given sequence."))
    print(mut_results)


def _maybe_handle_vendor_request(command: str):
    """Detect vendor-focused requests and handle them directly."""
    vendor_keywords = ['order', 'vendor', 'synthesis', 'test', 'assay', 'expression', 'function', 'binding']
    if not any(keyword in command.lower() for keyword in vendor_keywords):
        return None

    print("[handle_command] Detected vendor research command, using dna_vendor_research tool")

    sequence_length = None
    if '96' in command or 'variants' in command:
        sequence_length = 1000

    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
    from dna_vendor_research import run_dna_vendor_research_raw

    result = run_dna_vendor_research_raw(command, sequence_length, 'large')

    return {
        "text": f"DNA vendor research completed: {result.get('message', 'Research done')}",
        "input": command,
        "output": result,
        "plot": {
            "data": [
                {
                    "x": ["Vendors", "Testing Options"],
                    "y": [result.get('total_vendors', 0), result.get('total_testing_options', 0)],
                    "type": "bar",
                }
            ],
            "layout": {"title": "DNA Vendor Research Results"},
        },
    }

