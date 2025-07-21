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
def plasmid_visualization(vector_name: str, cloning_sites: str, insert_sequence: str) -> Dict:
    """Generate plasmid visualization data from vector name, cloning sites, and insert sequence."""
    
    # Import the plasmid visualization function
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
    from plasmid_visualizer import run_plasmid_visualization_raw
    
    result = run_plasmid_visualization_raw(vector_name, cloning_sites, insert_sequence)
    
    return {
        "text": result.get("text", "Plasmid visualization created successfully."),
        "input": {
            "vector_name": vector_name,
            "cloning_sites": cloning_sites,
            "insert_sequence": insert_sequence
        },
        "output": result.get("plasmid_data", {}),
        "visualization_type": result.get("visualization_type", "circular_plasmid"),
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
model = init_chat_model("openai:gpt-4o", max_tokens=4000)
config = {"configurable": {"thread_id": "abc123"}}

# Create a custom prompt that helps the agent choose the right tool
CUSTOM_PROMPT = """You are a bioinformatics assistant with access to the following tools:

1. sequence_alignment - Use for aligning DNA/RNA sequences
2. mutate_sequence - Use for CREATING new sequence variants from an input sequence
3. dna_vendor_research - Use for ORDERING sequences, finding VENDORS, or researching TESTING options
4. phylogenetic_tree - Use for creating phylogenetic trees from aligned sequences
5. sequence_selection - Use for selecting sequences from alignments based on criteria
6. synthesis_submission - Use for submitting sequences for DNA synthesis

IMPORTANT: 
- When the user mentions "order", "vendor", "synthesis", "test", "assay", "expression", "function", "binding" - use the dna_vendor_research tool.
- When the user wants to CREATE new variants from a sequence, use mutate_sequence.
- When the user wants to COMPARE existing sequences, use sequence_alignment.
- When the user wants to create a phylogenetic tree from aligned sequences, use phylogenetic_tree.
- When the user wants to select sequences from an alignment, use sequence_selection.
- When the user wants to submit sequences for synthesis, use synthesis_submission.

SEQUENCE SELECTION vs ALIGNMENT:
- Use sequence_selection when the user says: "pick", "select", "choose", "from the", "randomly pick", "select from", "choose from"
- Use sequence_alignment when the user says: "align", "alignment", "compare sequences", "multiple sequence alignment"
- If the user says "pick 1 sequence" or "select 3 sequences" - this is sequence_selection, NOT alignment

The user's request: {input}

Think step by step about what tool to use, then respond accordingly."""

agent = create_react_agent(
    model=model,
    tools=[sequence_alignment, mutate_sequence, dna_vendor_research, phylogenetic_tree, sequence_selection, synthesis_submission, plasmid_visualization],
    checkpointer=memory,
)


async def handle_command(command: str, session_id: str = "default", session_context: Dict = None):
    print(f"[handle_command] command: {command}")
    print(f"[handle_command] session_id: {session_id}")

    # Clear agent memory to prevent context accumulation
    try:
        # Clear the memory for this session to prevent token bloat
        memory.clear()
        print(f"[handle_command] Cleared agent memory for session {session_id}")
    except Exception as e:
        print(f"[handle_command] Could not clear memory: {e}")

    # Check if this is a vendor research command
    vendor_keywords = ['order', 'vendor', 'synthesis', 'test', 'assay', 'expression', 'function', 'binding']
    if any(keyword in command.lower() for keyword in vendor_keywords):
        print("[handle_command] Detected vendor research command, using dna_vendor_research tool")
        
        # Extract sequence length from command if mentioned
        sequence_length = None
        if '96' in command or 'variants' in command:
            sequence_length = 1000  # Default length for variants
        
        # Use the vendor research tool directly
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
                "data": [{"x": ["Vendors", "Testing Options"], "y": [result.get('total_vendors', 0), result.get('total_testing_options', 0)], "type": "bar"}],
                "layout": {"title": "DNA Vendor Research Results"},
            },
        }

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

