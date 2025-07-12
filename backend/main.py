from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from agent import handle_command
from command_router import CommandRouter
import uuid
from typing import Dict, Optional

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Global session storage
sessions: Dict[str, Dict] = {}

class CommandRequest(BaseModel):
    command: str
    session_id: Optional[str] = None

class SessionResponse(BaseModel):
    session_id: str

# Initialize command router
command_router = CommandRouter()

async def handle_command_direct(command: str, session_id: str = "default"):
    """Handle commands directly using keyword matching instead of ChatGPT."""
    print(f"[handle_command_direct] command: {command}")
    print(f"[handle_command_direct] session_id: {session_id}")
    
    # Get session context
    session_context = sessions.get(session_id, {})
    
    # Route the command to the appropriate tool
    tool_name, parameters = command_router.route_command(command, session_context)
    
    print(f"[handle_command_direct] Routing to tool: {tool_name}")
    print(f"[handle_command_direct] Parameters: {parameters}")
    
    # Execute the tool directly
    try:
        if tool_name == "mutate_sequence":
            from tools.mutations import run_mutation_raw
            result = run_mutation_raw(parameters["sequence"], parameters["num_variants"])
            return {
                "text": result.get("text", "Sequence mutated successfully."),
                "input": parameters,
                "output": result.get("statistics", {}),
                "plot": result.get("plot", {})
            }
        
        elif tool_name == "sequence_alignment":
            from tools.alignment import run_alignment_tool
            result = run_alignment_tool(parameters["sequences"])
            return {
                "text": result.get("text", "Sequences aligned successfully."),
                "input": parameters["sequences"],
                "output": result.get("alignment", []),
                "statistics": result.get("statistics", {}),
                "plot": result.get("plot", {})
            }
        
        elif tool_name == "sequence_selection":
            from tools.sequence_selection import run_sequence_selection_raw
            result = run_sequence_selection_raw(
                parameters["aligned_sequences"], 
                parameters["selection_type"], 
                parameters["num_sequences"]
            )
            return {
                "text": result.get("text", "Sequence selection completed successfully."),
                "input": parameters,
                "output": result.get("selected_sequences", []),
                "selection_criteria": result.get("selection_criteria", {})
            }
        
        elif tool_name == "phylogenetic_tree":
            from tools.phylogenetic_tree import run_phylogenetic_tree_raw
            result = run_phylogenetic_tree_raw(parameters["aligned_sequences"])
            return {
                "text": result.get("text", "Phylogenetic tree created successfully."),
                "input": parameters["aligned_sequences"],
                "output": result.get("statistics", {}),
                "tree_newick": result.get("tree_newick", ""),
                "statistics": result.get("statistics", {})
            }
        
        elif tool_name == "dna_vendor_research":
            from tools.dna_vendor_research import run_dna_vendor_research_raw
            result = run_dna_vendor_research_raw(
                parameters["command"], 
                parameters["sequence_length"], 
                parameters["quantity"]
            )
            return {
                "text": f"DNA vendor research completed: {result.get('message', 'Research done')}",
                "input": parameters,
                "output": result,
                "plot": result.get("plot", {})
            }
        
        else:
            # Fallback to general command
            return {
                "text": f"Command '{command}' routed to {tool_name}",
                "input": parameters,
                "output": {"status": "processed", "tool": tool_name},
                "plot": {}
            }
    
    except Exception as e:
        print(f"[handle_command_direct] Error executing tool {tool_name}: {e}")
        return {
            "text": f"Error executing {tool_name}: {str(e)}",
            "input": parameters,
            "output": {"error": str(e)},
            "plot": {}
        }

@app.post("/create_session")
async def create_session():
    """Create a new session and return the session ID."""
    session_id = str(uuid.uuid4())
    sessions[session_id] = {
        "aligned_sequences": None,
        "selected_sequences": None,
        "mutated_sequences": None,
        "plasmid_data": None
    }
    return SessionResponse(session_id=session_id)

@app.post("/execute")
async def execute(req: CommandRequest):
    """Execute a command with session context."""
    session_id = req.session_id or "default"
    
    # Initialize command router
    command_router = CommandRouter()
    
    # Initialize session if it doesn't exist
    if session_id not in sessions:
        sessions[session_id] = {
            "aligned_sequences": None,
            "selected_sequences": None,
            "mutated_sequences": None,
            "plasmid_data": None
        }
    
    # Get session context
    session_context = sessions[session_id]
    
    # Enhance command with session context
    enhanced_command = req.command
    
    # If command mentions sequence selection (pick, select, choose) and we have aligned sequences in session
    selection_keywords = ['pick', 'select', 'choose']
    if any(keyword in req.command.lower() for keyword in selection_keywords) and session_context["aligned_sequences"]:
        enhanced_command = f"{req.command}\n\nAligned sequences from previous step:\n{session_context['aligned_sequences']}"
        print(f"🔧 Enhanced command with aligned sequences for selection")
    
    # If command mentions mutation and we have selected sequences in session
    mutation_keywords = ['mutate', 'mutation', 'variant']
    if any(keyword in req.command.lower() for keyword in mutation_keywords) and session_context["selected_sequences"]:
        sequences_text = "\n".join([f">{i+1}\n{seq}" for i, seq in enumerate(session_context["selected_sequences"])])
        enhanced_command = f"{req.command}\n\nSelected sequences from previous step:\n{sequences_text}"
        print(f"🔧 Enhanced command with selected sequences for mutation")
    
    # If command mentions analysis, phylogenetic tree, or visualization and we have mutated sequences in session
    analysis_keywords = ['analyze', 'statistics', 'stats', 'phylogenetic', 'tree', 'visualize', 'variant']
    if any(keyword in req.command.lower() for keyword in analysis_keywords) and session_context["mutated_sequences"]:
        sequences_text = "\n".join([f">mutant_{i+1}\n{seq}" for i, seq in enumerate(session_context["mutated_sequences"])])
        enhanced_command = f"{req.command}\n\nMutated sequences from previous step:\n{sequences_text}"
        print(f"🔧 Enhanced command with mutated sequences for analysis/visualization")
    
    print(f"Original command: {req.command}")
    print(f"Enhanced command: {enhanced_command}")
    print(f"Session context: {session_context}")
    print(f"Session ID: {session_id}")
    print(f"All sessions: {list(sessions.keys())}")
    
    # Execute the enhanced command using command router instead of ChatGPT
    result = await handle_command_direct(enhanced_command, session_id)
    
    # Clean up old session data to prevent context bloat
    if len(sessions[session_id]) > 10:  # If session has too many entries
        # Keep only the most recent data
        sessions[session_id] = {
            "aligned_sequences": sessions[session_id].get("aligned_sequences"),
            "selected_sequences": sessions[session_id].get("selected_sequences"),
            "mutated_sequences": sessions[session_id].get("mutated_sequences"),
            "plasmid_data": sessions[session_id].get("plasmid_data")
        }
    
    # Update session context based on the result
    print(f"🔧 Session update - result keys: {result.keys() if result else 'None'}")
    print(f"🔧 Session update - full result: {result}")
    
    # Update session based on the tool that was executed
    # Get the tool name from the command router
    tool_name, _ = command_router.route_command(req.command, session_context)
    
    print(f"🔧 Session update - detected tool: {tool_name}")
    
    if tool_name == "mutate_sequence" and result and "output" in result:
        # Store mutated sequences in session
        if "variants" in result["output"]:
            mutated_sequences = result["output"]["variants"]
            sessions[session_id]["mutated_sequences"] = mutated_sequences
            print(f"🔧 Updated session with mutated sequences: {len(mutated_sequences)} variants")
    
    elif tool_name == "sequence_alignment" and result and "output" in result:
        # Store aligned sequences in session
        if isinstance(result["output"], list):
            aligned_sequences = "\n".join([f">{seq['name']}\n{seq['sequence']}" for seq in result["output"]])
            sessions[session_id]["aligned_sequences"] = aligned_sequences
            print(f"🔧 Updated session with aligned sequences: {len(result['output'])} sequences")
    
    elif tool_name == "sequence_selection" and result and "output" in result:
        # Store selected sequences in session
        if isinstance(result["output"], list):
            selected_sequences = [seq["sequence"] for seq in result["output"]]
            sessions[session_id]["selected_sequences"] = selected_sequences
            print(f"🔧 Updated session with selected sequences: {len(selected_sequences)} sequences")
    
    return result

@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "healthy", "sessions": len(sessions)}

@app.get("/mcp/tools")
async def list_tools():
    """List available MCP tools."""
    tools = [
        {
            "name": "sequence_alignment",
            "description": "Align multiple DNA/RNA sequences using various algorithms",
            "parameters": {
                "sequences": "string",
                "algorithm": "string (optional)"
            }
        },
        {
            "name": "sequence_selection",
            "description": "Select sequences from aligned sequences based on criteria",
            "parameters": {
                "aligned_sequences": "string",
                "selection_type": "string (optional)",
                "num_sequences": "number (optional)"
            }
        },
        {
            "name": "mutate_sequence",
            "description": "Generate mutated variants of a DNA sequence",
            "parameters": {
                "sequence": "string",
                "num_variants": "number (optional)",
                "mutation_rate": "number (optional)"
            }
        },
        {
            "name": "analyze_sequence_data",
            "description": "Analyze sequence data for various properties",
            "parameters": {
                "data": "string",
                "analysis_type": "string (optional)"
            }
        },
        {
            "name": "visualize_alignment",
            "description": "Create visualizations of sequence alignments",
            "parameters": {
                "alignment_file": "string",
                "output_format": "string (optional)"
            }
        },
        {
            "name": "plasmid_visualization",
            "description": "Visualize plasmid constructs and cloning",
            "parameters": {
                "vector_name": "string",
                "cloning_sites": "string",
                "insert_sequence": "string"
            }
        },
        {
            "name": "synthesis_submission",
            "description": "Submit sequences for DNA synthesis",
            "parameters": {
                "sequences": "string",
                "vendor_preference": "string (optional)",
                "quantity": "string (optional)",
                "delivery_time": "string (optional)"
            }
        }
    ]
    return {"tools": tools}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main:app", host="0.0.0.0", port=8001, reload=True)
