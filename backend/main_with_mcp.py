from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Dict, Any, Optional
import asyncio
import json
import subprocess
import sys
from pathlib import Path

# Add the current directory to Python path for imports
sys.path.append(str(Path(__file__).parent))

from agent import handle_command

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

class MCPToolRequest(BaseModel):
    tool_name: str
    arguments: Dict[str, Any]

class MCPResponse(BaseModel):
    success: bool
    result: Dict[str, Any]
    error: Optional[str] = None

class SequenceAlignmentRequest(BaseModel):
    sequences: str
    algorithm: str = "clustal"

class MutationRequest(BaseModel):
    sequence: str
    num_variants: int = 96
    mutation_rate: float = 0.1

class AnalysisRequest(BaseModel):
    data: str
    analysis_type: str = "alignment"

@app.post("/execute")
async def execute(req: CommandRequest):
    """Execute a general command using the existing agent."""
    return await handle_command(req.command)

@app.post("/mcp/sequence-alignment")
async def sequence_alignment_mcp(req: SequenceAlignmentRequest):
    """Perform sequence alignment using MCP server."""
    try:
        result = await call_mcp_tool("sequence_alignment", {
            "sequences": req.sequences,
            "algorithm": req.algorithm
        })
        return MCPResponse(success=True, result=result)
    except Exception as e:
        return MCPResponse(success=False, error=str(e))

@app.post("/mcp/mutate-sequence")
async def mutate_sequence_mcp(req: MutationRequest):
    """Generate sequence mutations using MCP server."""
    try:
        result = await call_mcp_tool("mutate_sequence", {
            "sequence": req.sequence,
            "num_variants": req.num_variants,
            "mutation_rate": req.mutation_rate
        })
        return MCPResponse(success=True, result=result)
    except Exception as e:
        return MCPResponse(success=False, error=str(e))

@app.post("/mcp/analyze-sequence-data")
async def analyze_sequence_data_mcp(req: AnalysisRequest):
    """Analyze sequence data using MCP server."""
    try:
        result = await call_mcp_tool("analyze_sequence_data", {
            "data": req.data,
            "analysis_type": req.analysis_type
        })
        return MCPResponse(success=True, result=result)
    except Exception as e:
        return MCPResponse(success=False, error=str(e))

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
        return MCPResponse(success=False, error=str(e))

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
                "name": "visualize_alignment",
                "description": "Create visual representation of sequence alignments",
                "parameters": {
                    "alignment_file": "string (path to FASTA file)",
                    "output_format": "string (png, svg, pdf)"
                }
            }
        ]
        return {"tools": tools}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

async def call_mcp_tool(tool_name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Call an MCP tool and return the result."""
    # For now, we'll use the existing tools directly
    # In a full implementation, this would communicate with the MCP server
    
    if tool_name == "sequence_alignment":
        from tools.alignment import run_alignment
        return run_alignment(arguments.get("sequences", ""))
    
    elif tool_name == "mutate_sequence":
        from tools.mutations import mutate_sequence
        return mutate_sequence(
            arguments.get("sequence", ""),
            arguments.get("num_variants", 96)
        )
    
    elif tool_name == "analyze_sequence_data":
        from tools.bio import align_and_visualize_fasta
        import pandas as pd
        
        data = arguments.get("data", "")
        if isinstance(data, str):
            try:
                df = pd.read_csv(data)
            except:
                df = parse_fasta_to_dataframe(data)
        else:
            df = data
        
        return {"result": str(align_and_visualize_fasta(df))}
    
    elif tool_name == "visualize_alignment":
        from tools.bio import align_and_visualize_fasta
        return {"result": str(align_and_visualize_fasta(None))}
    
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