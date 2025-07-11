from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from agent import handle_command
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
    
    # If command mentions aligned sequences and we have them in session
    if ("aligned sequences" in req.command.lower() or "select sequence" in req.command.lower()) and session_context["aligned_sequences"]:
        enhanced_command = f"{req.command}\n\nAligned sequences from previous step:\n{session_context['aligned_sequences']}"
    
    # If command mentions selected sequences and we have them in session
    if "selected sequences" in req.command.lower() and session_context["selected_sequences"]:
        sequences_text = "\n".join([f">{i+1}\n{seq}" for i, seq in enumerate(session_context["selected_sequences"])])
        enhanced_command = f"{req.command}\n\nSelected sequences from previous step:\n{sequences_text}"
    
    # If command mentions mutated sequences and we have them in session
    if "mutated sequences" in req.command.lower() and session_context["mutated_sequences"]:
        sequences_text = "\n".join([f">mutant_{i+1}\n{seq}" for i, seq in enumerate(session_context["mutated_sequences"])])
        enhanced_command = f"{req.command}\n\nMutated sequences from previous step:\n{sequences_text}"
    
    print(f"Original command: {req.command}")
    print(f"Enhanced command: {enhanced_command}")
    print(f"Session context: {session_context}")
    
    # Execute the enhanced command
    result = await handle_command(enhanced_command, session_id)
    
    # Update session context based on the result
    if result and "result" in result:
        messages = result["result"].get("messages", [])
        for message in messages:
            if message.get("type") == "tool":
                tool_name = message.get("name")
                if tool_name == "sequence_alignment":
                    try:
                        import json
                        tool_result = json.loads(message.get("content", "{}"))
                        if "output" in tool_result and isinstance(tool_result["output"], list):
                            aligned_sequences = "\n".join([f">{seq['name']}\n{seq['sequence']}" for seq in tool_result["output"]])
                            sessions[session_id]["aligned_sequences"] = aligned_sequences
                            print(f"Updated session with aligned sequences: {aligned_sequences}")
                    except Exception as e:
                        print(f"Error parsing sequence alignment result: {e}")
                
                elif tool_name == "sequence_selection":
                    try:
                        import json
                        tool_result = json.loads(message.get("content", "{}"))
                        if "output" in tool_result and isinstance(tool_result["output"], list):
                            selected_sequences = [seq["sequence"] for seq in tool_result["output"]]
                            sessions[session_id]["selected_sequences"] = selected_sequences
                            print(f"Updated session with selected sequences: {selected_sequences}")
                    except Exception as e:
                        print(f"Error parsing sequence selection result: {e}")
                
                elif tool_name == "mutate_sequence":
                    try:
                        import json
                        tool_result = json.loads(message.get("content", "{}"))
                        if "output" in tool_result and "variants" in tool_result["output"]:
                            mutated_sequences = tool_result["output"]["variants"]
                            sessions[session_id]["mutated_sequences"] = mutated_sequences
                            print(f"Updated session with mutated sequences: {mutated_sequences}")
                    except Exception as e:
                        print(f"Error parsing mutation result: {e}")
    
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
