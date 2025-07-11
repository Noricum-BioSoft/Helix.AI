from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Dict, Any, Optional, List
import asyncio
import json
import subprocess
import sys
from pathlib import Path

# Add the current directory to Python path for imports
sys.path.append(str(Path(__file__).parent))

from agent import handle_command
from history_manager import history_manager

app = FastAPI(title="DataBloom.AI Bioinformatics API", version="1.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class CommandRequest(BaseModel):
    command: str
    session_id: Optional[str] = None

class MCPToolRequest(BaseModel):
    tool_name: str
    arguments: Dict[str, Any]
    session_id: Optional[str] = None

class MCPResponse(BaseModel):
    success: bool
    result: Dict[str, Any]
    error: Optional[str] = None
    session_id: Optional[str] = None

class SequenceAlignmentRequest(BaseModel):
    sequences: str
    algorithm: str = "clustal"
    session_id: Optional[str] = None

class MutationRequest(BaseModel):
    sequence: str
    num_variants: int = 96
    mutation_rate: float = 0.1
    session_id: Optional[str] = None

class AnalysisRequest(BaseModel):
    data: str
    analysis_type: str = "alignment"
    session_id: Optional[str] = None

class VariantSelectionRequest(BaseModel):
    session_id: str
    selection_criteria: str = "diversity"
    num_variants: int = 10
    custom_filters: Optional[Dict[str, Any]] = None

class CommandParseRequest(BaseModel):
    command: str
    session_id: Optional[str] = None

class CommandExecuteRequest(BaseModel):
    parsed_command: Dict[str, Any]

class NaturalCommandRequest(BaseModel):
    command: str
    session_id: Optional[str] = None

class PlasmidVisualizationRequest(BaseModel):
    vector_name: str
    cloning_sites: str
    insert_sequence: str
    session_id: Optional[str] = None

class SessionRequest(BaseModel):
    user_id: Optional[str] = None

@app.post("/session/create")
async def create_session(req: SessionRequest):
    """Create a new session for tracking user interactions."""
    try:
        session_id = history_manager.create_session(req.user_id)
        return {
            "success": True,
            "session_id": session_id,
            "message": "Session created successfully"
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/session/{session_id}")
async def get_session_info(session_id: str):
    """Get session information and history."""
    try:
        session = history_manager.get_session(session_id)
        if not session:
            raise HTTPException(status_code=404, detail="Session not found")
        
        summary = history_manager.get_session_summary(session_id)
        return {
            "success": True,
            "session": session,
            "summary": summary
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/execute")
async def execute(req: CommandRequest):
    """Execute a general command using the existing agent with session tracking."""
    try:
        # Create session if not provided
        if not req.session_id:
            req.session_id = history_manager.create_session()
        
        # Execute command
        result = await handle_command(req.command)
        
        # Track in history
        history_manager.add_history_entry(
            req.session_id,
            req.command,
            "general_command",
            result
        )
        
        return {
            "success": True,
            "result": result,
            "session_id": req.session_id
        }
    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "session_id": req.session_id
        }

@app.post("/mcp/sequence-alignment")
async def sequence_alignment_mcp(req: SequenceAlignmentRequest):
    """Perform sequence alignment using MCP server with session tracking."""
    try:
        # Create session if not provided
        if not req.session_id:
            req.session_id = history_manager.create_session()
        
        result = await call_mcp_tool("sequence_alignment", {
            "sequences": req.sequences,
            "algorithm": req.algorithm
        })
        
        # Track in history
        history_manager.add_history_entry(
            req.session_id,
            f"Align sequences using {req.algorithm}",
            "sequence_alignment",
            result,
            {"algorithm": req.algorithm}
        )
        
        return MCPResponse(success=True, result=result, session_id=req.session_id)
    except Exception as e:
        return MCPResponse(success=False, result={}, error=str(e), session_id=req.session_id)

@app.post("/mcp/mutate-sequence")
async def mutate_sequence_mcp(req: MutationRequest):
    """Generate sequence mutations using MCP server with session tracking."""
    try:
        # Create session if not provided
        if not req.session_id:
            req.session_id = history_manager.create_session()
        
        result = await call_mcp_tool("mutate_sequence", {
            "sequence": req.sequence,
            "num_variants": req.num_variants,
            "mutation_rate": req.mutation_rate
        })
        
        # Track in history
        history_manager.add_history_entry(
            req.session_id,
            f"Generate {req.num_variants} variants with mutation rate {req.mutation_rate}",
            "mutate_sequence",
            result,
            {
                "num_variants": req.num_variants,
                "mutation_rate": req.mutation_rate
            }
        )
        
        return MCPResponse(success=True, result=result, session_id=req.session_id)
    except Exception as e:
        return MCPResponse(success=False, result={}, error=str(e), session_id=req.session_id)

@app.post("/mcp/analyze-sequence-data")
async def analyze_sequence_data_mcp(req: AnalysisRequest):
    """Analyze sequence data using MCP server with session tracking."""
    try:
        # Create session if not provided
        if not req.session_id:
            req.session_id = history_manager.create_session()
        
        result = await call_mcp_tool("analyze_sequence_data", {
            "data": req.data,
            "analysis_type": req.analysis_type
        })
        
        # Track in history
        history_manager.add_history_entry(
            req.session_id,
            f"Analyze sequence data with type {req.analysis_type}",
            "analyze_sequence_data",
            result,
            {"analysis_type": req.analysis_type}
        )
        
        return MCPResponse(success=True, result=result, session_id=req.session_id)
    except Exception as e:
        return MCPResponse(success=False, result={}, error=str(e), session_id=req.session_id)

@app.post("/mcp/select-variants")
async def select_variants_mcp(req: VariantSelectionRequest):
    """Select variants from previous mutation results."""
    try:
        # Add tools directory to path
        tools_path = str((Path(__file__).resolve().parent.parent / "tools").resolve())
        sys.path.insert(0, tools_path)
        
        import variant_selection
        
        result = variant_selection.run_variant_selection_raw(
            req.session_id,
            req.selection_criteria,
            req.num_variants,
            req.custom_filters
        )
        
        # Track in history
        history_manager.add_history_entry(
            req.session_id,
            f"Select {req.num_variants} variants using {req.selection_criteria} criteria",
            "select_variants",
            result,
            {
                "selection_criteria": req.selection_criteria,
                "num_variants": req.num_variants,
                "custom_filters": req.custom_filters
            }
        )
        
        return MCPResponse(success=True, result=result, session_id=req.session_id)
    except Exception as e:
        return MCPResponse(success=False, result={}, error=str(e), session_id=req.session_id)

@app.post("/mcp/parse-command")
async def parse_command_mcp(req: CommandParseRequest):
    """Parse a natural language command into structured parameters."""
    try:
        # Add tools directory to path
        tools_path = str((Path(__file__).resolve().parent.parent / "tools").resolve())
        sys.path.insert(0, tools_path)
        
        import command_parser
        
        result = command_parser.parse_command_raw(req.command, req.session_id)
        
        return MCPResponse(success=True, result=result, session_id=req.session_id)
    except Exception as e:
        return MCPResponse(success=False, result={}, error=str(e), session_id=req.session_id)

@app.post("/mcp/execute-command")
async def execute_command_mcp(req: CommandExecuteRequest):
    """Execute a parsed command using appropriate tools."""
    try:
        # Add tools directory to path
        tools_path = str((Path(__file__).resolve().parent.parent / "tools").resolve())
        sys.path.insert(0, tools_path)
        
        import command_executor
        
        result = command_executor.execute_command_raw(req.parsed_command)
        
        return MCPResponse(success=True, result=result, session_id=req.parsed_command.get("session_id"))
    except Exception as e:
        return MCPResponse(success=False, result={}, error=str(e))

@app.post("/mcp/handle-natural-command")
async def handle_natural_command_mcp(req: NaturalCommandRequest):
    """Handle a natural language command by parsing and executing it."""
    try:
        # Add tools directory to path
        tools_path = str((Path(__file__).resolve().parent.parent / "tools").resolve())
        sys.path.insert(0, tools_path)
        
        import command_handler
        
        result = command_handler.handle_command_raw(req.command, req.session_id)
        
        return MCPResponse(success=True, result=result, session_id=req.session_id)
    except Exception as e:
        return MCPResponse(success=False, result={}, error=str(e), session_id=req.session_id)

@app.post("/mcp/visualize-alignment")
async def visualize_alignment_mcp(alignment_file: str, output_format: str = "png"):
    """Visualize sequence alignment using MCP server."""
    try:
        result = await call_mcp_tool("visualize_alignment", {
            "alignment_file": alignment_file,
            "output_format": output_format
        })
        return MCPResponse(success=True, result=result)
    except Exception as e:
        return MCPResponse(success=False, result={}, error=str(e))

@app.post("/mcp/plasmid-visualization")
async def plasmid_visualization_mcp(req: PlasmidVisualizationRequest):
    """Generate plasmid visualization data from vector, cloning sites, and insert sequence."""
    try:
        # Create session if not provided
        if not req.session_id:
            req.session_id = history_manager.create_session()
        
        # Add tools directory to path
        tools_path = str((Path(__file__).resolve().parent.parent / "tools").resolve())
        sys.path.insert(0, tools_path)
        
        import plasmid_visualizer
        
        result = plasmid_visualizer.run_plasmid_visualization_raw(
            req.vector_name,
            req.cloning_sites,
            req.insert_sequence
        )
        
        # Track in history
        history_manager.add_history_entry(
            req.session_id,
            f"Visualize plasmid {req.vector_name} with {req.cloning_sites} and insert {req.insert_sequence}",
            "plasmid_visualization",
            result,
            {
                "vector_name": req.vector_name,
                "cloning_sites": req.cloning_sites,
                "insert_sequence": req.insert_sequence
            }
        )
        
        return MCPResponse(success=True, result=result, session_id=req.session_id)
    except Exception as e:
        return MCPResponse(success=False, result={}, error=str(e), session_id=req.session_id)

@app.get("/mcp/tools")
async def list_mcp_tools():
    """List available MCP tools."""
    try:
        # This would typically call the MCP server to list tools
        tools = [
            {
                "name": "sequence_alignment",
                "description": "Perform multiple sequence alignment on a set of DNA/RNA sequences",
                "parameters": {
                    "sequences": "string (FASTA format)",
                    "algorithm": "string (clustal, muscle, mafft)"
                }
            },
            {
                "name": "mutate_sequence",
                "description": "Generate mutations and variants of a given sequence",
                "parameters": {
                    "sequence": "string (DNA/RNA sequence)",
                    "num_variants": "integer (default: 96)",
                    "mutation_rate": "float (default: 0.1)"
                }
            },
            {
                "name": "analyze_sequence_data",
                "description": "Analyze sequence data and generate visualizations",
                "parameters": {
                    "data": "string (CSV or FASTA format)",
                    "analysis_type": "string (alignment, phylogeny, composition)"
                }
            },
            {
                "name": "select_variants",
                "description": "Select variants from previous mutation results based on criteria",
                "parameters": {
                    "session_id": "string (session ID)",
                    "selection_criteria": "string (diversity, random, length, custom)",
                    "num_variants": "integer (default: 10)",
                    "custom_filters": "object (optional custom filters)"
                }
            },
            {
                "name": "visualize_alignment",
                "description": "Create visual representation of sequence alignments",
                "parameters": {
                    "alignment_file": "string (path to FASTA file)",
                    "output_format": "string (png, svg, pdf)"
                }
            },
            {
                "name": "plasmid_visualization",
                "description": "Generate plasmid visualization data from vector, cloning sites, and insert sequence",
                "parameters": {
                    "vector_name": "string (e.g., pTet)",
                    "cloning_sites": "string (e.g., BsaI:123-456, EcoRI:789-1012)",
                    "insert_sequence": "string (DNA sequence to insert)"
                }
            }
        ]
        return {"tools": tools}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

async def call_mcp_tool(tool_name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Call an MCP tool and return the result."""
    # Add tools directory to path
    tools_path = str((Path(__file__).resolve().parent.parent / "tools").resolve())
    sys.path.insert(0, tools_path)
    
    if tool_name == "sequence_alignment":
        import alignment
        return alignment.run_alignment(arguments.get("sequences", ""))
    
    elif tool_name == "mutate_sequence":
        import mutations
        return mutations.run_mutation_raw(
            arguments.get("sequence", ""),
            arguments.get("num_variants", 96)
        )
    
    elif tool_name == "analyze_sequence_data":
        import bio
        import pandas as pd
        
        data = arguments.get("data", "")
        if isinstance(data, str):
            try:
                df = pd.read_csv(data)
            except:
                df = parse_fasta_to_dataframe(data)
        else:
            df = data
        
        return {"result": str(bio.align_and_visualize_fasta(df))}
    
    elif tool_name == "visualize_alignment":
        import bio
        return {"result": str(bio.align_and_visualize_fasta(None))}
    
    elif tool_name == "plasmid_visualization":
        import plasmid_visualizer
        return plasmid_visualizer.run_plasmid_visualization_raw(
            arguments.get("vector_name", ""),
            arguments.get("cloning_sites", ""),
            arguments.get("insert_sequence", "")
        )
    
    else:
        raise ValueError(f"Unknown tool: {tool_name}")

def parse_fasta_to_dataframe(fasta_content: str):
    """Parse FASTA content to DataFrame."""
    import pandas as pd
    
    sequences = []
    names = []
    current_seq = ""
    current_name = ""
    
    for line in fasta_content.split('\n'):
        if line.startswith('>'):
            if current_name and current_seq:
                names.append(current_name)
                sequences.append(current_seq)
            current_name = line[1:].strip()
            current_seq = ""
        else:
            current_seq += line.strip()
    
    if current_name and current_seq:
        names.append(current_name)
        sequences.append(current_seq)
    
    return pd.DataFrame({
        'name': names,
        'sequence': sequences
    })

@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "healthy", "service": "DataBloom.AI Bioinformatics API"}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main_with_mcp:app", host="0.0.0.0", port=8001, reload=True) 