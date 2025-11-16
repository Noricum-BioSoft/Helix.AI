from fastapi import FastAPI, HTTPException
from fastapi import Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from typing import Dict, Any, Optional, List
import asyncio
import json
import subprocess
import sys
from pathlib import Path
import numpy as np
from datetime import datetime, timezone

# Add the current directory to Python path for imports
sys.path.append(str(Path(__file__).parent))

from agent import handle_command
from history_manager import history_manager

def serialize_for_json(obj):
    """Convert objects to JSON-serializable format, handling numpy types."""
    if isinstance(obj, dict):
        return {key: serialize_for_json(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [serialize_for_json(item) for item in obj]
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj

class CustomJSONResponse(JSONResponse):
    def render(self, content):
        serialized_content = serialize_for_json(content)
        return json.dumps(serialized_content, default=str).encode("utf-8")

app = FastAPI(title="Helix.AI Bioinformatics API", version="1.0.0")

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

class AgentCommandRequest(BaseModel):
    prompt: str
    session_id: Optional[str] = None
    files: Optional[List[Dict[str, Any]]] = None

# -------------------------------
# Limits and validation utilities
# -------------------------------

TEN_MB_BYTES = 10 * 1024 * 1024
MAX_PROMPT_TOKENS = 5000
MAX_PROMPTS_PER_DAY = 100

# In-memory per-identity counters. In production, move to Redis or DB.
_daily_prompt_counters: Dict[str, Dict[str, int]] = {}

def _get_request_identity(req: Request, session_id: Optional[str]) -> str:
    """
    Identify requester for rate limiting. Prefer session_id; fall back to client IP.
    """
    if session_id:
        return f"session:{session_id}"
    client_ip = req.client.host if req and req.client else "unknown"
    return f"ip:{client_ip}"

def _today_iso() -> str:
    return datetime.now(timezone.utc).date().isoformat()

def _check_and_increment_daily_counter(identity: str) -> None:
    day = _today_iso()
    by_identity = _daily_prompt_counters.setdefault(identity, {})
    count = by_identity.get(day, 0)
    if count >= MAX_PROMPTS_PER_DAY:
        raise HTTPException(status_code=429, detail=f"Daily prompt limit reached ({MAX_PROMPTS_PER_DAY}). Try again tomorrow.")
    by_identity[day] = count + 1

def _estimate_token_count(text: str) -> int:
    """
    Lightweight token estimator without external deps.
    Uses a conservative heuristic ~4 chars per token.
    """
    if not text:
        return 0
    char_tokens = max(1, round(len(text) / 4))
    word_tokens = len(text.split())
    return max(char_tokens, word_tokens)

def _validate_prompt_length(prompt: str) -> None:
    tokens = _estimate_token_count(prompt)
    if tokens > MAX_PROMPT_TOKENS:
        raise HTTPException(
            status_code=413,
            detail=f"Prompt too long: ~{tokens} tokens (max {MAX_PROMPT_TOKENS}). Please shorten your prompt."
        )

def _validate_files(files: Optional[List[Dict[str, Any]]]) -> None:
    if not files:
        return
    for idx, f in enumerate(files):
        # Accept 'size' in bytes if provided by the client.
        size = f.get("size")
        if isinstance(size, int) and size > TEN_MB_BYTES:
            raise HTTPException(
                status_code=413,
                detail=f"Uploaded file #{idx+1} exceeds 10 MB limit."
            )
        # If raw content is provided, guard enormous payloads
        content = f.get("content")
        if isinstance(content, str):
            # Rough size check: 1 char ~ 1 byte for ASCII payloads
            if len(content) > TEN_MB_BYTES:
                raise HTTPException(
                    status_code=413,
                    detail=f"Uploaded file #{idx+1} content exceeds 10 MB limit."
                )

class PlasmidVisualizationRequest(BaseModel):
    vector_name: str
    cloning_sites: str
    insert_sequence: str
    session_id: Optional[str] = None

class PlasmidForRepresentativesRequest(BaseModel):
    representatives: List[str]
    aligned_sequences: str
    vector_name: str = "pUC19"
    cloning_sites: str = "EcoRI, BamHI, HindIII"
    session_id: Optional[str] = None

class ReadTrimmingRequest(BaseModel):
    reads: str
    adapter: Optional[str] = None
    quality_threshold: int = 20
    session_id: Optional[str] = None

class ReadMergingRequest(BaseModel):
    forward_reads: str
    reverse_reads: str
    min_overlap: int = 12
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

@app.post("/create_session")
async def create_session_alias():
    """Alias for /session/create to support frontend compatibility."""
    try:
        session_id = history_manager.create_session()
        return {
            "session_id": session_id
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
async def execute(req: CommandRequest, request: Request):
    """Execute a general command using the existing agent with session tracking."""
    try:
        # Validate prompt length
        _validate_prompt_length(req.command)

        # Create session if not provided
        if not req.session_id:
            req.session_id = history_manager.create_session()
        
        # Rate limit per identity (session or IP)
        identity = _get_request_identity(request, req.session_id)
        _check_and_increment_daily_counter(identity)
        
        # Get session context from history manager
        session_context = {}
        if hasattr(history_manager, 'sessions') and req.session_id in history_manager.sessions:
            session_context = history_manager.sessions[req.session_id]
        
        # Use command router with session context
        from command_router import CommandRouter
        command_router = CommandRouter()
        
        # Route the command to the appropriate tool
        tool_name, parameters = command_router.route_command(req.command, session_context)
        print(f"🔧 Routed command '{req.command}' to tool '{tool_name}' with parameters: {parameters}")
        
        # DEBUG: Print session context before tool call
        print(f"🔧 [DEBUG] Session context before {tool_name} call:")
        print(f"  Session ID: {req.session_id}")
        print(f"  Session context keys: {list(session_context.keys()) if session_context else 'None'}")
        if session_context and "mutated_sequences" in session_context:
            print(f"  mutated_sequences count: {len(session_context['mutated_sequences'])}")
        if session_context and "mutation_results" in session_context:
            print(f"  mutation_results count: {len(session_context['mutation_results'])}")
        
        # Call the appropriate MCP tool
        try:
            result = await call_mcp_tool(tool_name, parameters)
            
            # Track in history with the correct tool name
            history_manager.add_history_entry(
                req.session_id,
                req.command,
                tool_name,  # Use the actual tool name instead of "general_command"
                result
            )
            
            return CustomJSONResponse({
                "success": True,
                "result": result,
                "session_id": req.session_id
            })
        except ValueError as e:
            # Handle unknown tool errors gracefully
            error_result = {
                "status": "error",
                "message": str(e),
                "tool": tool_name,
                "suggestion": "Try a different command or check available tools"
            }
            return CustomJSONResponse({
                "success": False,
                "result": error_result,
                "session_id": req.session_id
            })
    except Exception as e:
        return CustomJSONResponse({
            "success": False,
            "error": str(e),
            "session_id": req.session_id
        })

@app.post("/agent")
async def agent_command(req: AgentCommandRequest, request: Request):
    """Proxy endpoint for the bioinformatics agent orchestrator."""
    try:
        # Validate prompt length
        _validate_prompt_length(req.prompt)

        session_id = req.session_id or history_manager.create_session()

        # Rate limit per identity (session or IP)
        identity = _get_request_identity(request, session_id)
        _check_and_increment_daily_counter(identity)

        # Capture uploaded files in session context if provided
        if req.files:
            # Validate file sizes
            _validate_files(req.files)
            if hasattr(history_manager, "sessions"):
                session_store = history_manager.sessions.setdefault(session_id, {})
                uploaded = session_store.setdefault("uploaded_files", [])
                uploaded.extend(req.files)

        session_context = {}
        if hasattr(history_manager, "sessions"):
            session_context = history_manager.sessions.get(session_id, {})

        result = await handle_command(req.prompt, session_id=session_id, session_context=session_context)

        # Track agent interaction in history
        history_manager.add_history_entry(
            session_id,
            req.prompt,
            "agent",
            result,
        )

        return CustomJSONResponse({
            "success": True,
            "result": result,
            "session_id": session_id
        })
    except Exception as err:
        return CustomJSONResponse({
            "success": False,
            "error": str(err),
            "session_id": req.session_id
        })

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
        
        return CustomJSONResponse({
            "success": True,
            "result": result,
            "session_id": req.session_id
        })
    except Exception as e:
        return CustomJSONResponse({
            "success": False,
            "result": {},
            "error": str(e),
            "session_id": req.session_id
        })
    except Exception as e:
        return CustomJSONResponse({
            "success": False,
            "result": {},
            "error": str(e),
            "session_id": req.session_id
        })

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
        # Store mutated sequences in session context for downstream steps
        # Try to extract variants from result
        variants = None
        if isinstance(result, dict):
            if "variants" in result:
                variants = result["variants"]
            elif "output" in result and isinstance(result["output"], dict) and "variants" in result["output"]:
                variants = result["output"]["variants"]
        if variants:
            # Store in the in-memory session context
            if hasattr(history_manager, 'sessions') and req.session_id in history_manager.sessions:
                history_manager.sessions[req.session_id]["mutated_sequences"] = variants
                # Also store as mutation_results for select_variants
                history_manager.sessions[req.session_id]["mutation_results"] = variants
                
                # DEBUG: Print session contents after mutation
                print(f"🔧 [DEBUG] After mutation - Session {req.session_id} contents:")
                session_data = history_manager.sessions[req.session_id]
                for key, value in session_data.items():
                    if key in ["mutated_sequences", "mutation_results"]:
                        print(f"  {key}: {len(value) if isinstance(value, list) else type(value)} items")
                        if isinstance(value, list) and len(value) > 0:
                            print(f"    First item: {value[0]}")
                    else:
                        print(f"  {key}: {type(value)}")
        
        return CustomJSONResponse({
            "success": True,
            "result": result,
            "session_id": req.session_id
        })
    except Exception as e:
        return CustomJSONResponse({
            "success": False,
            "result": {},
            "error": str(e),
            "session_id": req.session_id
        })

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
        
        return CustomJSONResponse({
            "success": True,
            "result": result,
            "session_id": req.session_id
        })
    except Exception as e:
        return CustomJSONResponse({
            "success": False,
            "result": {},
            "error": str(e),
            "session_id": req.session_id
        })

@app.post("/mcp/select-variants")
async def select_variants_mcp(req: VariantSelectionRequest):
    """Select variants from previous mutation results."""
    try:
        # DEBUG: Print session contents before selection
        print(f"🔧 [DEBUG] Before selection - Session {req.session_id} contents:")
        if hasattr(history_manager, 'sessions') and req.session_id in history_manager.sessions:
            session_data = history_manager.sessions[req.session_id]
            for key, value in session_data.items():
                if key in ["mutated_sequences", "mutation_results"]:
                    print(f"  {key}: {len(value) if isinstance(value, list) else type(value)} items")
                    if isinstance(value, list) and len(value) > 0:
                        print(f"    First item: {value[0]}")
                else:
                    print(f"  {key}: {type(value)}")
        else:
            print(f"  Session {req.session_id} not found in history_manager.sessions")
        
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
        
        return CustomJSONResponse({
            "success": True,
            "result": result,
            "session_id": req.session_id
        })
    except Exception as e:
        return CustomJSONResponse({
            "success": False,
            "result": {},
            "error": str(e),
            "session_id": req.session_id
        })

@app.post("/mcp/parse-command")
async def parse_command_mcp(req: CommandParseRequest):
    """Parse a natural language command into structured parameters."""
    try:
        # Add tools directory to path
        tools_path = str((Path(__file__).resolve().parent.parent / "tools").resolve())
        sys.path.insert(0, tools_path)
        
        import command_parser
        
        result = command_parser.parse_command_raw(req.command, req.session_id)
        
        return CustomJSONResponse({
            "success": True,
            "result": result,
            "session_id": req.session_id
        })
    except Exception as e:
        return CustomJSONResponse({
            "success": False,
            "result": {},
            "error": str(e),
            "session_id": req.session_id
        })

@app.post("/mcp/execute-command")
async def execute_command_mcp(req: CommandExecuteRequest):
    """Execute a parsed command using appropriate tools."""
    try:
        # Add tools directory to path
        tools_path = str((Path(__file__).resolve().parent.parent / "tools").resolve())
        sys.path.insert(0, tools_path)
        
        import command_executor
        
        result = command_executor.execute_command_raw(req.parsed_command)
        
        return CustomJSONResponse({
            "success": True,
            "result": result,
            "session_id": req.parsed_command.get("session_id")
        })
    except Exception as e:
        return CustomJSONResponse({
            "success": False,
            "result": {},
            "error": str(e)
        })

@app.post("/mcp/handle-natural-command")
async def handle_natural_command_mcp(req: NaturalCommandRequest):
    """Handle a natural language command by parsing and executing it."""
    try:
        # Add tools directory to path
        tools_path = str((Path(__file__).resolve().parent.parent / "tools").resolve())
        sys.path.insert(0, tools_path)
        
        import command_handler
        
        result = command_handler.handle_command_raw(req.command, req.session_id)
        
        return CustomJSONResponse({
            "success": True,
            "result": result,
            "session_id": req.session_id
        })
    except Exception as e:
        return CustomJSONResponse({
            "success": False,
            "result": {},
            "error": str(e),
            "session_id": req.session_id
        })

@app.post("/mcp/visualize-alignment")
async def visualize_alignment_mcp(alignment_file: str, output_format: str = "png"):
    """Visualize sequence alignment using MCP server."""
    try:
        result = await call_mcp_tool("visualize_alignment", {
            "alignment_file": alignment_file,
            "output_format": output_format
        })
        return CustomJSONResponse({
            "success": True,
            "result": result
        })
    except Exception as e:
        return CustomJSONResponse({
            "success": False,
            "result": {},
            "error": str(e)
        })

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
        
        return CustomJSONResponse({
            "success": True,
            "result": result,
            "session_id": req.session_id
        })
    except Exception as e:
        return CustomJSONResponse({
            "success": False,
            "result": {},
            "error": str(e),
            "session_id": req.session_id
        })

@app.post("/mcp/plasmid-for-representatives")
async def plasmid_for_representatives_mcp(req: PlasmidForRepresentativesRequest):
    """Generate plasmid visualizations for representative sequences from clustering."""
    try:
        # Create session if not provided
        if not req.session_id:
            req.session_id = history_manager.create_session()
        
        # Add tools directory to path
        tools_path = str((Path(__file__).resolve().parent.parent / "tools").resolve())
        sys.path.insert(0, tools_path)
        
        import plasmid_visualizer
        
        result = plasmid_visualizer.create_plasmid_for_representatives(
            req.representatives,
            req.aligned_sequences,
            req.vector_name,
            req.cloning_sites
        )
        
        # Track in history
        history_manager.add_history_entry(
            req.session_id,
            f"Create plasmid visualizations for {len(req.representatives)} representatives in {req.vector_name}",
            "plasmid_for_representatives",
            result,
            {
                "representatives": req.representatives,
                "vector_name": req.vector_name,
                "cloning_sites": req.cloning_sites
            }
        )
        
        return CustomJSONResponse({
            "success": True,
            "result": result,
            "session_id": req.session_id
        })
    except Exception as e:
        return CustomJSONResponse({
            "success": False,
            "result": {},
            "error": str(e),
            "session_id": req.session_id
        })

@app.post("/mcp/read-trimming")
async def read_trimming_mcp(req: ReadTrimmingRequest):
    """Perform read trimming using MCP tooling."""
    try:
        session_id = req.session_id or history_manager.create_session()
        result = await call_mcp_tool("read_trimming", {
            "reads": req.reads,
            "adapter": req.adapter,
            "quality_threshold": req.quality_threshold,
        })
        history_manager.add_history_entry(
            session_id,
            "Trim sequencing reads",
            "read_trimming",
            result,
            {
                "adapter": req.adapter,
                "quality_threshold": req.quality_threshold,
            }
        )
        return CustomJSONResponse({
            "success": True,
            "result": result,
            "session_id": session_id
        })
    except Exception as e:
        return CustomJSONResponse({
            "success": False,
            "result": {},
            "error": str(e),
            "session_id": req.session_id
        })

@app.post("/mcp/read-merging")
async def read_merging_mcp(req: ReadMergingRequest):
    """Merge paired-end reads using MCP tooling."""
    try:
        session_id = req.session_id or history_manager.create_session()
        result = await call_mcp_tool("read_merging", {
            "forward_reads": req.forward_reads,
            "reverse_reads": req.reverse_reads,
            "min_overlap": req.min_overlap,
        })
        history_manager.add_history_entry(
            session_id,
            "Merge paired-end reads",
            "read_merging",
            result,
            {
                "min_overlap": req.min_overlap,
            }
        )
        return CustomJSONResponse({
            "success": True,
            "result": result,
            "session_id": session_id
        })
    except Exception as e:
        return CustomJSONResponse({
            "success": False,
            "result": {},
            "error": str(e),
            "session_id": req.session_id
        })

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
            },
            {
                "name": "plasmid_for_representatives",
                "description": "Create plasmid visualizations for representative sequences from clustering analysis",
                "parameters": {
                    "representatives": "list of strings (representative sequence names)",
                    "aligned_sequences": "string (FASTA format sequences)",
                    "vector_name": "string (e.g., pUC19)",
                    "cloning_sites": "string (e.g., EcoRI, BamHI, HindIII)"
                }
            },
            {
                "name": "read_trimming",
                "description": "Trim adapters and low-quality bases from sequencing reads",
                "parameters": {
                    "reads": "string (FASTQ format)",
                    "adapter": "string (adapter sequence to remove)",
                    "quality_threshold": "integer (minimum Phred score to retain)"
                }
            },
            {
                "name": "read_merging",
                "description": "Merge paired-end reads using overlap consensus",
                "parameters": {
                    "forward_reads": "string (FASTQ format with forward reads)",
                    "reverse_reads": "string (FASTQ format with reverse reads)",
                    "min_overlap": "integer (minimum overlap length for merging)"
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
    
    elif tool_name == "select_variants":
        import variant_selection
        return variant_selection.run_variant_selection_raw(
            arguments.get("session_id", ""),
            arguments.get("selection_criteria", "diversity"),
            arguments.get("num_variants", 10),
            arguments.get("custom_filters", None)
        )
    
    elif tool_name == "plasmid_visualization":
        import plasmid_visualizer
        return plasmid_visualizer.run_plasmid_visualization_raw(
            arguments.get("vector_name", ""),
            arguments.get("cloning_sites", ""),
            arguments.get("insert_sequence", "")
        )
    
    elif tool_name == "plasmid_for_representatives":
        import plasmid_visualizer
        return plasmid_visualizer.create_plasmid_for_representatives(
            arguments.get("representatives", []),
            arguments.get("aligned_sequences", ""),
            arguments.get("vector_name", "pUC19"),
            arguments.get("cloning_sites", "EcoRI, BamHI, HindIII")
        )

    elif tool_name == "read_trimming":
        import read_trimming
        return read_trimming.run_read_trimming_raw(
            arguments.get("reads", ""),
            arguments.get("adapter"),
            arguments.get("quality_threshold", 20),
        )

    elif tool_name == "read_merging":
        import read_merging
        return read_merging.run_read_merging_raw(
            arguments.get("forward_reads", ""),
            arguments.get("reverse_reads", ""),
            arguments.get("min_overlap", 12),
        )
    
    elif tool_name == "handle_natural_command":
        # Use the natural command handler
        import command_handler
        command = arguments.get("command", "")
        session_id = arguments.get("session_id", "")
        return command_handler.handle_command_raw(command, session_id)
    
    elif tool_name == "phylogenetic_tree":
        # Handle phylogenetic tree analysis
        import phylogenetic_tree
        aligned_sequences = arguments.get("aligned_sequences", "")
        return phylogenetic_tree.run_phylogenetic_tree_raw(aligned_sequences)
    
    elif tool_name == "clustering_analysis":
        # Handle clustering analysis
        import phylogenetic_tree
        aligned_sequences = arguments.get("aligned_sequences", "")
        num_clusters = arguments.get("num_clusters", 5)
        return phylogenetic_tree.run_clustering_from_tree(aligned_sequences, num_clusters)
    
    elif tool_name == "variant_selection":
        # Handle variant selection
        import phylogenetic_tree
        aligned_sequences = arguments.get("aligned_sequences", "")
        num_variants = arguments.get("num_variants", 10)
        return phylogenetic_tree.run_variant_selection_from_tree(aligned_sequences, num_variants)
    
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
    return {"status": "healthy", "service": "Helix.AI Bioinformatics API"}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main_with_mcp:app", host="0.0.0.0", port=8001, reload=True) 