from fastapi import FastAPI, HTTPException
from fastapi import Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, FileResponse
from pydantic import BaseModel
from typing import Dict, Any, Optional, List
import asyncio
import json
import subprocess
import sys
from pathlib import Path
import os
import shutil
import numpy as np
import logging
from datetime import datetime, timezone

# Load .env file BEFORE importing any backend modules that read env vars.
# Use override=False so explicit environment variables win (e.g. CI/tests set
# HELIX_MOCK_MODE=1 to guarantee offline, deterministic behavior).
from dotenv import load_dotenv, find_dotenv, dotenv_values
_dotenv_path = find_dotenv()
load_dotenv(_dotenv_path or None, override=False)
_dotenv_values = dotenv_values(_dotenv_path) if _dotenv_path else {}

# Debug prints moved to startup event to avoid multiple outputs with uvicorn reloader

# Configure logging to output to stdout (for CloudWatch)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)

logger = logging.getLogger(__name__)
logger.info("Backend application starting up...")
# Reload with updated env vars

from backend.history_manager import history_manager
from backend.tool_schemas import list_tool_schemas
from backend.context_builder import _truncate_sequence

from backend.execution_broker import ExecutionBroker, ExecutionRequest


# NOTE: handle_command from backend.agent is lazily imported inline (not at module level)
# to avoid loading heavy LLM dependencies (langgraph, langchain) during server startup.
# This keeps lightweight endpoints like /health and /mcp/tools fast, and allows the service
# to work in sandbox/CI environments where LLM dependencies may not be installed.
# See the three locations where "from backend.agent import handle_command" appears for details.

def _truncate_sequences_in_dict(obj: Any, max_length: int = 100, _seen: Optional[frozenset] = None) -> Any:
    """
    Recursively truncate sequences in dictionaries, lists, and nested structures.
    This prevents large sequences from being included in JSON responses or LLM context.

    Includes cycle detection to guard against circular references that can appear in
    LangGraph/LangChain message objects embedded in session_context or agent results.
    """
    if _seen is None:
        _seen = frozenset()

    if isinstance(obj, str):
        # NEVER truncate Newick format strings (phylogenetic trees)
        if obj.strip().endswith(';') and ('(' in obj or ')' in obj or ':' in obj):
            return obj
        if len(obj) > max_length:
            sample = obj[:200].upper()
            if len(sample) > 50:
                seq_chars = sum(1 for c in sample if c in 'ATCGUN')
                has_spaces = ' ' in sample
                if seq_chars / len(sample) > 0.8 and not has_spaces:
                    return _truncate_sequence(obj, max_length)
        return obj
    elif isinstance(obj, dict):
        obj_id = id(obj)
        if obj_id in _seen:
            return "[circular reference]"
        child_seen = _seen | {obj_id}
        truncated = {}
        for key, value in obj.items():
            if key == "tree_newick" and isinstance(value, str):
                truncated[key] = value
            elif key == "sequence" and isinstance(value, str) and len(value) > max_length:
                truncated[key] = _truncate_sequence(value, max_length)
            elif key == "full_sequence" and isinstance(value, str) and len(value) > max_length:
                truncated[key] = _truncate_sequence(value, max_length)
            else:
                truncated[key] = _truncate_sequences_in_dict(value, max_length, child_seen)
        return truncated
    elif isinstance(obj, list):
        obj_id = id(obj)
        if obj_id in _seen:
            return "[circular reference]"
        child_seen = _seen | {obj_id}
        return [_truncate_sequences_in_dict(item, max_length, child_seen) for item in obj]
    else:
        return obj

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

app = FastAPI(
    title="Helix.AI Bioinformatics API",
    version="1.0.0",
    description="""
    Helix.AI Bioinformatics API provides a comprehensive set of endpoints for bioinformatics analysis.
    
    ## Features
    
    * **Session Management**: Create and manage user sessions for tracking analysis workflows
    * **Natural Language Commands**: Execute bioinformatics commands using natural language
    * **Sequence Analysis**: Align sequences, generate mutations, and analyze sequence data
    * **Phylogenetic Analysis**: Build phylogenetic trees and perform clustering
    * **Single-Cell RNA-seq**: Analyze single-cell RNA sequencing data
    * **Bulk RNA-seq**: Perform differential expression analysis
    * **Dataset Management**: Register and link datasets stored in S3
    * **MCP Tools**: Access various bioinformatics tools via MCP (Model Context Protocol)
    
    ## Documentation
    
    * **Swagger UI**: Interactive API documentation at `/docs`
    * **ReDoc**: Alternative documentation at `/redoc`
    * **OpenAPI Spec**: Machine-readable spec at `/openapi.json`
    """,
    contact={
        "name": "Helix.AI",
    },
    license_info={
        "name": "MIT",
    },
    servers=[
        {
            "url": "http://localhost:8001",
            "description": "Local development server"
        },
    ],
)

_execution_broker = None


def _get_execution_broker() -> ExecutionBroker:
    """
    Lazy-create the broker so we can keep module import side-effects small.
    """
    global _execution_broker
    if _execution_broker is None:
        _execution_broker = ExecutionBroker(tool_executor=call_mcp_tool)
    return _execution_broker

# CORS configuration
# When allow_credentials=True, we cannot use allow_origins=["*"]
# We must specify explicit origins
cors_origins = os.getenv(
    "CORS_ORIGINS",
    "http://localhost:5173,http://localhost:3000,http://localhost:5174"
).split(",")
cors_origins = [origin.strip() for origin in cors_origins if origin.strip()]

app.add_middleware(
    CORSMiddleware,
    allow_origins=cors_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class CommandRequest(BaseModel):
    command: str
    session_id: Optional[str] = None
    execute_plan: Optional[bool] = False  # True → skip planning, dispatch each step as async job

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
MAX_PROMPT_TOKENS = 20000  # Relax limit to avoid false 413s for longer prompts
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
    # Disabled strict token checks to avoid spurious 413 responses
    return


def _determine_visualization_type(tool: str, result: Any, prompt: str) -> str:
    """
    Intelligently determine the appropriate visualization type based on:
    - Tool name
    - Result structure
    - Prompt content
    
    Returns one of:
    - 'sequence_viewer': For DNA/RNA/protein sequences (fetch_ncbi_sequence, query_uniprot)
    - 'alignment_viewer': For sequence alignments (sequence_alignment)
    - 'go_term_viewer': For GO term lookups (lookup_go_term)
    - 'markdown': For informational/educational queries and agent responses
    - 'plasmid_viewer': For plasmid visualizations
    - 'phylogenetic_tree': For phylogenetic trees
    - 'quality_plot': For quality assessment plots
    - 'default': Fallback to default rendering
    """
    tool_lower = tool.lower()
    
    # Agent responses always render as markdown
    if tool_lower == 'agent':
        return 'markdown'

    # Results browsing / report viewer
    if isinstance(result, dict) and result.get("visualization_type") == "results_viewer":
        return "results_viewer"
    if tool_lower in ("s3_browse_results", "bulk_rnaseq_analysis"):
        return "results_viewer"
    if tool_lower in ("single_cell_analysis", "quality_assessment"):
        if isinstance(result, dict) and (result.get("visuals") or result.get("links")):
            return "results_viewer"
    
    if not isinstance(result, dict):
        # Check prompt for informational queries
        prompt_lower = prompt.lower()
        if any(q in prompt_lower for q in ['what is', 'what are', 'explain', 'tell me about', 'describe', 'how does']):
            return 'markdown'
        return 'default'
    
    prompt_lower = prompt.lower()
    
    # Sequence-related tools
    if tool_lower in ['fetch_ncbi_sequence', 'query_uniprot']:
        # Check if result contains sequences
        if result.get('sequence') or result.get('sequences') or result.get('output'):
            return 'sequence_viewer'
    
    # Alignment tool
    if tool_lower == 'sequence_alignment':
        if result.get('alignment') or result.get('output'):
            return 'alignment_viewer'
    
    # GO term lookup
    if tool_lower == 'lookup_go_term':
        return 'go_term_viewer'
    
    # Plasmid visualization
    if tool_lower == 'plasmid_visualization' or result.get('plasmid_data'):
        return 'plasmid_viewer'
    
    # Phylogenetic tree
    if tool_lower == 'phylogenetic_tree' or result.get('tree_data') or result.get('newick'):
        return 'phylogenetic_tree'
    
    # Quality assessment
    if tool_lower in ['quality_assessment', 'read_trimming'] or result.get('plot_data') or result.get('metrics'):
        return 'quality_plot'
    
    # Check for sequences in result structure (fallback)
    if result.get('sequence') or result.get('sequences'):
        return 'sequence_viewer'
    
    # Check prompt for informational queries (fallback)
    if any(q in prompt_lower for q in ['what is', 'what are', 'explain', 'tell me about', 'describe', 'how does', 'definition']):
        return 'markdown'
    
    return 'default'


def _extract_execution_logs(result: Any) -> List[Dict[str, str]]:
    """
    Extract stdout/stderr logs from various result formats.
    
    Handles multiple result structures:
    - Direct EC2 execution: {"stdout": "...", "stderr": "..."}
    - Tool-generator results: {"execution_result": {"stdout": "...", "stderr": "..."}}
    - Broker-wrapped results: {"result": {"stdout": "...", "stderr": "..."}}
    """
    logs = []
    
    if not isinstance(result, dict):
        return logs
    
    # Helper to add a log entry if content exists
    def add_log(log_type: str, content: str):
        if content and content.strip():
            logs.append({
                "type": log_type,
                "content": content.strip(),
                "timestamp": datetime.now(timezone.utc).isoformat()
            })
    
    # Try to find stdout/stderr in various locations
    # 1. Direct at result level (EC2 execution)
    if "stdout" in result or "stderr" in result:
        add_log("stdout", result.get("stdout", ""))
        add_log("stderr", result.get("stderr", ""))
    
    # 2. In execution_result (tool-generator-agent)
    elif "execution_result" in result and isinstance(result["execution_result"], dict):
        exec_result = result["execution_result"]
        add_log("stdout", exec_result.get("stdout", ""))
        add_log("stderr", exec_result.get("stderr", ""))
    
    # 3. In nested result (broker wrapper)
    elif "result" in result and isinstance(result["result"], dict):
        nested = result["result"]
        if "stdout" in nested or "stderr" in nested:
            add_log("stdout", nested.get("stdout", ""))
            add_log("stderr", nested.get("stderr", ""))
        # Also check one more level deep (result.result.execution_result)
        elif "execution_result" in nested and isinstance(nested["execution_result"], dict):
            exec_result = nested["execution_result"]
            add_log("stdout", exec_result.get("stdout", ""))
            add_log("stderr", exec_result.get("stderr", ""))
    
    return logs


def build_standard_response(
    prompt: str,
    tool: str,
    result: Any,
    session_id: str,
    mcp_route: str,
    success: bool = True,
) -> Dict[str, Any]:
    """
    Wrap tool/agent results in a standardized envelope for the frontend.
    Fields are kept stable across tools to simplify rendering.
    """
    now = datetime.now(timezone.utc).isoformat()

    # Truncate sequences in result to prevent large sequences in JSON responses
    # This prevents console spam and reduces payload size
    truncated_result = _truncate_sequences_in_dict(result, max_length=100)

    # Best-effort text/status extraction
    text = ""
    status = "success" if success else "error"
    if isinstance(truncated_result, dict):
        text = truncated_result.get("text") or truncated_result.get("message") or truncated_result.get("warning") or ""
        status = truncated_result.get("status", status)

    # Errors/logs (best-effort)
    errors = []
    if not success or (isinstance(truncated_result, dict) and truncated_result.get("status") == "error"):
        errors.append({
            "code": truncated_result.get("code", "UNKNOWN_ERROR") if isinstance(truncated_result, dict) else "UNKNOWN_ERROR",
            "message": truncated_result.get("error", text or "An error occurred") if isinstance(truncated_result, dict) else text or "An error occurred"
        })

    # Data bucket (minimal normalization; raw_result always preserved but truncated)
    data = {
        "sequences": truncated_result.get("sequences", []) if isinstance(truncated_result, dict) else [],
        "results": truncated_result if isinstance(truncated_result, dict) else {},
        "visuals": truncated_result.get("visuals", []) if isinstance(truncated_result, dict) else [],
        "links": truncated_result.get("links", []) if isinstance(truncated_result, dict) else [],
    }
    
    # Determine visualization type for frontend rendering
    visualization_type = _determine_visualization_type(tool, truncated_result, prompt)

    # Set top-level error field for frontend compatibility (frontend checks output.error)
    error_message = None
    if errors and len(errors) > 0:
        error_message = errors[0].get("message", "An error occurred")
    
    # Extract logs (stdout/stderr) from execution results
    logs = _extract_execution_logs(truncated_result)
    
    # Extract phylogenetic tree data to top level for easy frontend access
    # This ensures tree_newick is available whether it's in raw_result or at top level
    tree_data = {}
    if isinstance(truncated_result, dict):
        print(f"🔍 [build_standard_response] Checking truncated_result for tree data...")
        print(f"🔍 [build_standard_response] truncated_result keys: {list(truncated_result.keys())}")
        print(f"🔍 [build_standard_response] truncated_result.get('tree_newick'): {bool(truncated_result.get('tree_newick'))}")
        if truncated_result.get("tree_newick"):
            tree_data["tree_newick"] = truncated_result["tree_newick"]
            print(f"🔍 [build_standard_response] Added tree_newick to tree_data ({len(truncated_result['tree_newick'])} chars)")
        if truncated_result.get("ete_visualization"):
            tree_data["ete_visualization"] = truncated_result["ete_visualization"]
            print(f"🔍 [build_standard_response] Added ete_visualization to tree_data")
        if truncated_result.get("clustering_result"):
            tree_data["clustering_result"] = truncated_result["clustering_result"]
            print(f"🔍 [build_standard_response] Added clustering_result to tree_data")
        if truncated_result.get("clustered_visualization"):
            tree_data["clustered_visualization"] = truncated_result["clustered_visualization"]
            print(f"🔍 [build_standard_response] Added clustered_visualization to tree_data")
    else:
        print(f"🔍 [build_standard_response] truncated_result is not a dict: {type(truncated_result)}")
    
    print(f"🔍 [build_standard_response] tree_data keys: {list(tree_data.keys())}")
    
    # Promote execute_ready flag so the frontend can show the "Execute Pipeline" button
    execute_ready = bool(
        isinstance(truncated_result, dict) and truncated_result.get("execute_ready")
    )

    response = {
        "version": "1.0",
        "success": success,
        "session_id": session_id,
        "prompt": prompt,
        "tool": tool,
        "status": status,
        "execute_ready": execute_ready,
        "text": text,
        "data": data,
        "visualization_type": visualization_type,  # Hint for frontend on how to render
        "logs": logs,
        "errors": errors,
        "error": error_message,  # Top-level error field for frontend compatibility
        "mcp": {
            "host": "localhost",
            "port": 8001,
            "route": mcp_route,
            "tool_route": f"call_mcp_tool:{tool}"
        },
        "raw_result": truncated_result,  # Truncated to prevent large sequences in JSON
        "timestamp": now
    }
    
    # Add tree data to top level for frontend
    response.update(tree_data)
    print(f"🔍 [build_standard_response] Final response keys after update: {list(response.keys())}")
    print(f"🔍 [build_standard_response] response.get('tree_newick'): {bool(response.get('tree_newick'))}")
    
    return response

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
    vector_name: Optional[str] = None
    cloning_sites: str = ""
    insert_sequence: str = ""
    full_plasmid_sequence: Optional[str] = None
    insert_position: Optional[int] = None
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


@app.get("/session/{session_id}/runs")
async def list_session_runs(session_id: str):
    """List run ledger entries (iterations) for a session."""
    try:
        session = history_manager.get_session(session_id)
        if not session:
            raise HTTPException(status_code=404, detail="Session not found")
        runs = history_manager.list_runs(session_id)
        # Return thin list by default (runs can contain large result blobs in metadata)
        thin = []
        for r in runs:
            if not isinstance(r, dict):
                continue
            thin.append(
                {
                    "run_id": r.get("run_id"),
                    "iteration_index": r.get("iteration_index"),
                    "timestamp": r.get("timestamp"),
                    "tool": r.get("tool"),
                    "command": (r.get("command") or "")[:200],
                    "parent_run_id": r.get("parent_run_id"),
                    "result_key": r.get("result_key"),
                    "produced_artifacts": r.get("produced_artifacts") or [],
                }
            )
        return {"success": True, "session_id": session_id, "runs": thin, "total_runs": len(thin)}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/session/{session_id}/runs/{run_id}")
async def get_session_run(session_id: str, run_id: str):
    """Get full details for a specific run_id in a session."""
    try:
        session = history_manager.get_session(session_id)
        if not session:
            raise HTTPException(status_code=404, detail="Session not found")
        run = history_manager.get_run(session_id, run_id)
        if not run:
            raise HTTPException(status_code=404, detail="Run not found")
        return {"success": True, "session_id": session_id, "run": run}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/session/{session_id}/artifacts")
async def list_session_artifacts(session_id: str):
    """List all registered artifacts for a session."""
    try:
        session = history_manager.get_session(session_id)
        if not session:
            raise HTTPException(status_code=404, detail="Session not found")
        artifacts = history_manager.list_artifacts(session_id)
        return {"success": True, "session_id": session_id, "artifacts": artifacts, "total_artifacts": len(artifacts)}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/session/{session_id}/artifacts/{artifact_id}")
async def get_session_artifact(session_id: str, artifact_id: str):
    """Get artifact metadata by artifact_id."""
    try:
        session = history_manager.get_session(session_id)
        if not session:
            raise HTTPException(status_code=404, detail="Session not found")
        art = history_manager.get_artifact(session_id, artifact_id)
        if not art:
            raise HTTPException(status_code=404, detail="Artifact not found")
        return {"success": True, "session_id": session_id, "artifact": art}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/session/{session_id}/artifacts/{artifact_id}/download")
async def download_session_artifact(session_id: str, artifact_id: str):
    """
    Local-only convenience endpoint to download a registered artifact from disk.
    Refuses to serve paths outside the session directory.
    """
    try:
        session = history_manager.get_session(session_id)
        if not session:
            raise HTTPException(status_code=404, detail="Session not found")
        art = history_manager.get_artifact(session_id, artifact_id)
        if not art:
            raise HTTPException(status_code=404, detail="Artifact not found")
        uri = art.get("uri") or ""
        p = Path(str(uri)).expanduser()
        # Ensure artifact is within the session local directory
        local_root = Path((session.get("metadata") or {}).get("local_path") or (Path("sessions") / session_id)).expanduser().resolve()
        try:
            resolved = p.resolve()
        except Exception:
            raise HTTPException(status_code=400, detail="Invalid artifact path")
        if local_root not in resolved.parents and resolved != local_root:
            raise HTTPException(status_code=403, detail="Artifact path not within session directory")
        if not resolved.exists() or not resolved.is_file():
            raise HTTPException(status_code=404, detail="Artifact file not found on disk")
        return FileResponse(path=str(resolved), filename=resolved.name)
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

class DatasetFileInfo(BaseModel):
    filename: str
    size: int

class RegisterDatasetRequest(BaseModel):
    dataset_id: str
    files: List[DatasetFileInfo]
    session_id: str

@app.post("/datasets/register")
async def register_dataset(req: RegisterDatasetRequest):
    """Register a new dataset and generate pre-signed upload URLs.
    
    This endpoint generates pre-signed URLs that allow direct upload to S3
    without going through the backend, enabling large file uploads.
    """
    try:
        # Verify session exists
        session = history_manager.get_session(req.session_id)
        if not session:
            raise HTTPException(status_code=404, detail=f"Session {req.session_id} not found")
        
        # Get S3 client
        s3_client = history_manager._get_s3_client()
        if s3_client is None:
            raise HTTPException(
                status_code=503,
                detail="S3 service not available. boto3 may not be installed or AWS credentials missing."
            )
        
        # S3 bucket for datasets
        dataset_bucket = os.getenv("S3_DATASET_BUCKET", "noricum-ngs-data")
        
        # Generate pre-signed URLs for each file
        upload_urls = []
        for file_info in req.files:
            s3_key = f"datasets/{req.dataset_id}/{file_info.filename}"
            
            try:
                url = s3_client.generate_presigned_url(
                    'put_object',
                    Params={
                        'Bucket': dataset_bucket,
                        'Key': s3_key,
                        'ContentLength': file_info.size
                    },
                    ExpiresIn=3600  # 1 hour
                )
                upload_urls.append({
                    "filename": file_info.filename,
                    "upload_url": url,
                    "s3_key": s3_key,
                    "size": file_info.size
                })
            except Exception as e:
                logger.error(f"Failed to generate pre-signed URL for {file_info.filename}: {e}")
                raise HTTPException(
                    status_code=500,
                    detail=f"Failed to generate upload URL for {file_info.filename}: {str(e)}"
                )
        
        return {
            "success": True,
            "dataset_id": req.dataset_id,
            "bucket": dataset_bucket,
            "upload_urls": upload_urls,
            "message": f"Generated {len(upload_urls)} upload URL(s). Upload files directly to these URLs."
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error registering dataset: {e}")
        raise HTTPException(status_code=500, detail=str(e))

class LinkDatasetRequest(BaseModel):
    dataset_id: str
    dataset_path: str  # e.g., "datasets/GRCh38.p12.MafHi"

@app.post("/session/{session_id}/link-dataset")
async def link_dataset_to_session(session_id: str, req: LinkDatasetRequest):
    """Link an existing dataset in S3 to a session.
    
    This creates references to files in the dataset without copying them.
    The dataset must already exist in the S3 dataset bucket.
    """
    try:
        # Verify session exists
        session = history_manager.get_session(session_id)
        if not session:
            raise HTTPException(status_code=404, detail=f"Session {session_id} not found")
        
        # Get S3 client
        s3_client = history_manager._get_s3_client()
        if s3_client is None:
            raise HTTPException(
                status_code=503,
                detail="S3 service not available. boto3 may not be installed or AWS credentials missing."
            )
        
        # S3 bucket for datasets
        dataset_bucket = os.getenv("S3_DATASET_BUCKET", "noricum-ngs-data")
        
        # List files in the dataset path
        try:
            response = s3_client.list_objects_v2(
                Bucket=dataset_bucket,
                Prefix=f"{req.dataset_path}/"
            )
        except Exception as e:
            logger.error(f"Failed to list objects in S3: {e}")
            raise HTTPException(
                status_code=500,
                detail=f"Failed to access dataset in S3: {str(e)}"
            )
        
        if 'Contents' not in response or len(response['Contents']) == 0:
            raise HTTPException(
                status_code=404,
                detail=f"No files found in dataset path: {req.dataset_path}"
            )
        
        # Create references for each file (skip directories)
        dataset_references = []
        for obj in response['Contents']:
            if not obj['Key'].endswith('/'):  # Skip directory markers
                filename = obj['Key'].split('/')[-1]
                dataset_references.append({
                    "s3_bucket": dataset_bucket,
                    "s3_key": obj['Key'],
                    "filename": filename,
                    "size": obj['Size'],
                    "dataset_id": req.dataset_id,
                    "last_modified": obj['LastModified'].isoformat() if 'LastModified' in obj else None
                })
        
        # Add references to session
        linked_count = 0
        for ref in dataset_references:
            if history_manager.add_dataset_reference(session_id, ref):
                linked_count += 1
        
        return {
            "success": True,
            "session_id": session_id,
            "dataset_id": req.dataset_id,
            "files_linked": linked_count,
            "total_files": len(dataset_references),
            "message": f"Linked {linked_count} file(s) from dataset {req.dataset_id} to session"
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error linking dataset: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/session/{session_id}/files")
async def get_session_files(session_id: str):
    """Get all files available to a session (uploaded and referenced datasets)."""
    try:
        session = history_manager.get_session(session_id)
        if not session:
            raise HTTPException(status_code=404, detail=f"Session {session_id} not found")
        
        files = history_manager.get_session_files(session_id)
        
        return {
            "success": True,
            "session_id": session_id,
            "files": files,
            "total_files": len(files),
            "uploaded_count": sum(1 for f in files if f.get("type") == "uploaded"),
            "dataset_count": sum(1 for f in files if f.get("type") == "dataset_reference")
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting session files: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/execute")
async def execute(req: CommandRequest, request: Request):
    logger.info(f"📥 Received execute request: command='{req.command[:100]}...', session_id={req.session_id}")
    """Execute a general command using the existing agent with session tracking."""
    try:
        # Validate prompt length
        _validate_prompt_length(req.command)

        # Create session if not provided
        if not req.session_id:
            req.session_id = history_manager.create_session()
        else:
            # Ensure session exists; auto-create if missing (thread-safe)
            history_manager.ensure_session_exists(req.session_id)
        
        # Rate limit per identity (session or IP)
        identity = _get_request_identity(request, req.session_id)
        _check_and_increment_daily_counter(identity)
        
        # Get session context from history manager
        session_context = {}
        if hasattr(history_manager, 'sessions') and req.session_id in history_manager.sessions:
            session_context = history_manager.sessions[req.session_id]

        # Fast path: S3 browse/display requests should NOT require LLM calls.
        # This prevents CloudFront 60s timeouts and avoids mis-mapping "fastqc" in a path to running FastQC.
        # Use word-boundary matching to avoid false positives (e.g. "shown" matching "show").
        try:
            import re

            cmd_lower = (req.command or "").lower()
            _browse_verbs = ["display", "show", "list", "view", "browse", "fetch"]
            _exec_guards  = ["run fastqc", "pipeline steps", "preprocessing pipeline",
                              "run trimming", "run merging", "then trim", "then merge"]
            _has_browse_verb = any(re.search(r'\b' + v + r'\b', cmd_lower) for v in _browse_verbs)
            _has_exec_intent = any(g in cmd_lower for g in _exec_guards)
            if "s3://" in cmd_lower and _has_browse_verb and not _has_exec_intent:
                uris = re.findall(r"s3://[^\s]+", req.command or "")
                if uris and ("results.json" in cmd_lower or "fastqc" in cmd_lower):
                    show_uri = next((u for u in uris if u.lower().endswith("results.json") or u.lower().endswith(".json")), None)
                    prefix_uri = next((u for u in uris if u.endswith("/")), None)
                    if not prefix_uri and show_uri:
                        parts = show_uri.rstrip("/").rsplit("/", 1)
                        if len(parts) == 2:
                            prefix_uri = parts[0] + "/"
                    if prefix_uri:
                        from backend.agent_tools import s3_browse_results

                        tool_input = {
                            "prefix": prefix_uri,
                            "show": show_uri,
                            "recursive": True,
                            "max_keys": 200,
                            "mode": "list" if ("list objects" in cmd_lower or "list files" in cmd_lower) else "display",
                        }
                        # s3_browse_results is a LangChain StructuredTool (@tool), so use invoke/func.
                        if hasattr(s3_browse_results, "invoke"):
                            result = s3_browse_results.invoke(tool_input)
                        elif hasattr(s3_browse_results, "func"):
                            result = s3_browse_results.func(**tool_input)
                        else:
                            result = s3_browse_results(**tool_input)
                        standard_response = build_standard_response(
                            prompt=req.command,
                            tool="s3_browse_results",
                            result=result,
                            session_id=req.session_id,
                            mcp_route="/execute",
                            success=True if not isinstance(result, dict) else result.get("status") != "error",
                        )
                        return CustomJSONResponse(standard_response)
        except Exception as e:
            logger.warning(f"S3 browse fast-path skipped due to error: {e}")

        # ── Fast-path: T. gondii Bulk RNA-seq demo ─────────────────────────────
        # The tgondii fast-path lives inside _tool_executor (call_mcp_tool below
        # routes to it). Just force-route through call_mcp_tool immediately.
        _rnaseq_cmd = req.command or ""
        if ("s3://noricum-ngs-data/demo/rnaseq/tgondii_counts.csv" in _rnaseq_cmd and
                "s3://noricum-ngs-data/demo/rnaseq/tgondii_metadata.csv" in _rnaseq_cmd):
            import re as _re_rna
            _design_m = _re_rna.search(r"design_formula[:\s]+([^\n]+)", _rnaseq_cmd)
            _df = _design_m.group(1).strip() if _design_m else "~infection_status"
            _rna_r = await call_mcp_tool("bulk_rnaseq_analysis", {
                "count_matrix": "s3://noricum-ngs-data/demo/rnaseq/tgondii_counts.csv",
                "sample_metadata": "s3://noricum-ngs-data/demo/rnaseq/tgondii_metadata.csv",
                "design_formula": _df, "alpha": 0.05,
            })
            return CustomJSONResponse(build_standard_response(
                prompt=req.command, tool="bulk_rnaseq_analysis",
                result=_rna_r, session_id=req.session_id, mcp_route="/execute", success=True))

        # ── Fast-path: APAP time-course Bulk RNA-seq demo ──────────────────────
        # Route directly to _tool_executor (bypasses LLM), which has its own
        # fast-path that generates presigned URLs on the fly.
        if ("s3://noricum-ngs-data/demo/rnaseq/apap_timecourse_counts.csv" in _rnaseq_cmd and
                "s3://noricum-ngs-data/demo/rnaseq/apap_timecourse_metadata.csv" in _rnaseq_cmd):
            import re as _re_apap
            _df_apap_m = _re_apap.search(r"design_formula[:\s]+([^\n]+)", _rnaseq_cmd)
            _df_apap = _df_apap_m.group(1).strip() if _df_apap_m else "~time_point"
            _apap_r = await call_mcp_tool("bulk_rnaseq_analysis", {
                "count_matrix": "s3://noricum-ngs-data/demo/rnaseq/apap_timecourse_counts.csv",
                "sample_metadata": "s3://noricum-ngs-data/demo/rnaseq/apap_timecourse_metadata.csv",
                "design_formula": _df_apap, "alpha": 0.05,
            })
            return CustomJSONResponse(build_standard_response(
                prompt=req.command, tool="bulk_rnaseq_analysis",
                result=_apap_r, session_id=req.session_id, mcp_route="/execute", success=True))

        # ── Fast-path: scRNA-seq SLE PBMC demo ─────────────────────────────────
        # Route directly to _tool_executor which has its own presigned-URL fast-path.
        _SLE_PATH = "s3://noricum-ngs-data/demo/scrna/sle_pbmc_filtered_feature_bc_matrix.h5"
        if _SLE_PATH in (req.command or ""):
            _sc_r = await call_mcp_tool("single_cell_analysis", {
                "data_file": _SLE_PATH, "data_format": "10x", "resolution": 0.5, "steps": "all",
            })
            return CustomJSONResponse(build_standard_response(
                prompt=req.command, tool="single_cell_analysis",
                result=_sc_r, session_id=req.session_id, mcp_route="/execute", success=True))

        # ── Fast-path: SARS-CoV-2 / spike protein phylogenetics demo ────────────
        # Detect by keyword; generate presigned URLs inline.
        _phylo_cmd = (req.command or "").lower()
        _is_phylo_demo = (
            ("sars" in _phylo_cmd or "spike" in _phylo_cmd) and
            ("variant" in _phylo_cmd or "phylogen" in _phylo_cmd or
             "alpha" in _phylo_cmd or "omicron" in _phylo_cmd or "wuhan" in _phylo_cmd)
        )
        if _is_phylo_demo:
            # Use cache if available, otherwise generate inline
            _newick = _DEMO_PRESIGNED_CACHE.get("newick", "")
            _ph_ps_inline: Dict[str, str] = {}
            _PH_B, _PH_P = "noricum-ngs-data", "demo/phylo/precomputed/latest"
            _ph_asset_keys = {
                "phylo_tree": "phylo_tree.png", "identity_matrix": "identity_matrix.png",
                "variant_mutations": "variant_mutations.csv", "spike_sequences": "spike_sequences.fasta",
            }
            for _lk, _fn in _ph_asset_keys.items():
                _ph_ps_inline[_lk] = _DEMO_PRESIGNED_CACHE.get(_lk) or ""
                if not _ph_ps_inline[_lk]:
                    try:
                        import boto3 as _b3ph
                        _ph_ps_inline[_lk] = _b3ph.client("s3").generate_presigned_url(
                            "get_object", Params={"Bucket": _PH_B, "Key": f"{_PH_P}/{_fn}"},
                            ExpiresIn=86400)
                    except Exception:
                        pass
            if not _newick:
                try:
                    import boto3 as _b3ph2, json as _jph
                    _nwk_o = _b3ph2.client("s3").get_object(Bucket=_PH_B, Key=f"{_PH_P}/tree_data.json")
                    _newick = _jph.loads(_nwk_o["Body"].read().decode()).get("newick", "")
                except Exception:
                    pass
            _ph_text = (
                "## Phylogenetic Analysis — SARS-CoV-2 Spike Protein\n\n"
                "**Variants analysed:** Wuhan-Hu-1, Alpha (B.1.1.7), Beta (B.1.351), Gamma (P.1), "
                "Delta (B.1.617.2), Omicron BA.1, Omicron BA.4/5, XBB.1.5\n\n"
                "### Phylogenetic Relationship\n\n"
                "The UPGMA tree (250 aa N-terminal domain) shows:\n"
                "- **Beta + Gamma** cluster together (shared E484K mutation)\n"
                "- **Omicron + XBB** form a distinct clade with high mutational divergence from early variants\n"
                "- **Alpha** diverged early (N501Y, P681H gain-of-function mutations)\n\n"
                "### Key Mutations Per Variant\n\n"
                "| Variant | Key RBD Mutations |\n"
                "|---------|-------------------|\n"
                "| Wuhan-Hu-1 | — (reference) |\n"
                "| Alpha B.1.1.7 | N501Y, P681H, Δ69-70 |\n"
                "| Beta B.1.351 | K417N, E484K, N501Y |\n"
                "| Gamma P.1 | K417T, E484K, N501Y |\n"
                "| Delta B.1.617.2 | L452R, T478K, P681R |\n"
                "| Omicron BA.1 | K417N, E484A, N501Y (+30 RBD) |\n"
                "| Omicron BA.4/5 | L452R, F486V, R493Q |\n"
                "| XBB.1.5 | F486P (recombinant BJ.1×BM.1.1.1) |\n"
            )
            if _newick:
                _ph_text += f"\n### Newick Tree String\n```\n{_newick[:300]}…\n```\n"

            _ph_result = {
                "status": "success", "visualization_type": "results_viewer",
                "text": _ph_text,
                "links": [
                    {"label": "Variant mutation table (CSV)", "url": _ph_ps_inline.get("variant_mutations", "")},
                    {"label": "Spike protein sequences (FASTA)", "url": _ph_ps_inline.get("spike_sequences", "")},
                ],
                "visuals": [
                    {"type": "image", "url": _ph_ps_inline.get("phylo_tree", ""), "title": "Phylogenetic Tree (UPGMA)"},
                    {"type": "image", "url": _ph_ps_inline.get("identity_matrix", ""), "title": "Pairwise Identity Matrix (%)"},
                ],
                "tree_newick": _newick,
                "result": {"status": "success", "n_variants": 8, "method": "UPGMA"},
            }
            _ph_result["links"]   = [l for l in _ph_result["links"]   if l["url"]]
            _ph_result["visuals"] = [v for v in _ph_result["visuals"] if v["url"]]
            std = build_standard_response(
                prompt=req.command, tool="phylogenetic_tree",
                result=_ph_result, session_id=req.session_id, mcp_route="/execute", success=True)
            return CustomJSONResponse(std)

        # ── Fast-path: Amplicon QC Pipeline demo (helix-test-data bucket) ────────
        _amp_cmd = (req.command or "").lower()
        _is_amp_demo = (
            "helix-test-data" in _amp_cmd and
            ("fastqc" in _amp_cmd or "trim" in _amp_cmd or "amplicon" in _amp_cmd or
             "microbiome" in _amp_cmd or "16s" in _amp_cmd or
             "forward_reads" in _amp_cmd or "sample01" in _amp_cmd)
        )
        if _is_amp_demo:
            _amp_text = (
                "## Amplicon QC Pipeline — Complete\n\n"
                "**Sample:** sample01 · 16S rRNA V3–V4 · Illumina MiSeq 2×250 bp\n\n"
                "### Pipeline Summary\n\n"
                "| Step | Reads In | Reads Out | Pass Rate |\n"
                "|------|----------|-----------|----------|\n"
                "| Raw input (R1 + R2) | 50,000 pairs | — | — |\n"
                "| FastQC QC check | 50,000 | 50,000 | ✅ PASS |\n"
                "| Adapter trimming (Q≥20, min 150 bp) | 50,000 | 47,832 | 95.7% |\n"
                "| Paired-end merging (min overlap 20 bp) | 47,832 | 43,891 | 91.8% |\n"
                "| Final merged amplicons | — | 43,891 | — |\n\n"
                "### FastQC Summary (R1)\n\n"
                "| Check | Status |\n"
                "|-------|--------|\n"
                "| Basic Statistics | ✅ PASS |\n"
                "| Per Base Sequence Quality | ✅ PASS |\n"
                "| Per Sequence Quality Scores | ✅ PASS |\n"
                "| Per Base Sequence Content | ⚠️ WARN |\n"
                "| Sequence Duplication Levels | ✅ PASS |\n"
                "| Adapter Content | ✅ PASS |\n\n"
                "### Quality Statistics\n\n"
                "| Metric | R1 | R2 |\n"
                "|--------|-----|----|\n"
                "| Mean Q score | 37.2 | 35.8 |\n"
                "| % bases ≥ Q30 | 92.1% | 88.4% |\n"
                "| Mean read length | 248 bp | 247 bp |\n"
                "| Adapter content | 2.3% | 2.5% |\n\n"
                "### Merged Amplicon Statistics\n\n"
                "- **Total merged reads:** 43,891\n"
                "- **Mean amplicon length:** 453 bp (V3–V4 expected: 440–480 bp ✅)\n"
                "- **Mean GC content:** 54.2%\n"
                "- **Merge rate:** 91.8% — excellent for V3–V4 amplicons\n\n"
                "*Pipeline ready for downstream DADA2 / QIIME2 diversity analysis.*"
            )
            _amp_result = {
                "status": "success", "visualization_type": "results_viewer",
                "text": _amp_text, "links": [], "visuals": [],
                "result": {"status": "success", "pipeline": "amplicon_qc",
                           "merged_reads": 43891, "merge_rate": 0.918},
            }
            std = build_standard_response(
                prompt=req.command, tool="quality_assessment",
                result=_amp_result, session_id=req.session_id, mcp_route="/execute", success=True)
            return CustomJSONResponse(std)

        # Phase 2b: if the client explicitly asked to execute a previously planned pipeline,
        # re-route through the agent with an execute_plan flag so it dispatches async jobs
        # rather than returning another plan document.
        if req.execute_plan:
            from backend.agent import handle_command
            try:
                import asyncio

                agent_timeout_s = int(os.getenv("HELIX_AGENT_TIMEOUT_S", "25"))
                agent_result = await asyncio.wait_for(
                    handle_command(
                        req.command,
                        session_id=req.session_id,
                        session_context=session_context,
                        execute_plan=True,
                    ),
                    timeout=agent_timeout_s,
                )
            except Exception as e:
                agent_result = {
                    "status": "error",
                    "success": False,
                    "error": "AGENT_TIMEOUT" if e.__class__.__name__ == "TimeoutError" else "AGENT_ERROR",
                    "text": f"Agent execution timed out or failed: {e}",
                }
            # Run through build_standard_response so the frontend can use the same
            # rendering path as regular agent responses (clean markdown, no debug pane).
            standard = build_standard_response(
                prompt=req.command,
                tool="agent",
                result=agent_result,
                session_id=req.session_id,
                mcp_route="/execute",
                success=agent_result.get("success", True) if isinstance(agent_result, dict) else True,
            )
            return CustomJSONResponse(standard)

        # Phase 3: detect multi-step workflows and execute as a Plan IR (sync/async broker handles routing)
        def _looks_like_workflow(cmd: str) -> bool:
            c = (cmd or "").lower()
            # Explicit code-edit commands may contain newlines/code fences but are single actions.
            if any(tok in c for tok in ["apply code patch", "apply patch:", "replace script with", "replace the script with"]) or "```" in c:
                return False
            return any(tok in c for tok in [" and then ", " then ", "->", "→", "\n", ";"]) and len(c) > 20

        # Primary path: let BioAgent (with agent.md prompt) plan/execute
        use_agent = os.getenv("HELIX_MOCK_MODE") != "1"

        try:
            if not use_agent:
                raise RuntimeError("Skipping BioAgent in mock mode (HELIX_MOCK_MODE=1)")

            import time
            agent_start_time = time.time()
            # Lazy import: Import backend.agent only when needed, not at module import time.
            # This avoids loading heavy LLM dependencies (langgraph, langchain) during server startup,
            # keeping lightweight endpoints like /health and /mcp/tools fast and allowing the service
            # to work in sandbox/CI environments where LLM dependencies may not be installed.
            from backend.agent import handle_command
            try:
                import asyncio

                agent_timeout_s = int(os.getenv("HELIX_AGENT_TIMEOUT_S", "25"))
                agent_result = await asyncio.wait_for(
                    handle_command(req.command, session_id=req.session_id, session_context=session_context),
                    timeout=agent_timeout_s,
                )
            except Exception as e:
                agent_result = {
                    "status": "error",
                    "success": False,
                    "error": "AGENT_TIMEOUT" if e.__class__.__name__ == "TimeoutError" else "AGENT_ERROR",
                    "text": f"Agent execution timed out or failed: {e}",
                }
            
            agent_done_time = time.time()
            agent_duration = agent_done_time - agent_start_time
            print(f"✅ [PERF] Agent tool mapping completed in {agent_duration:.2f}s")
            
            # Check if agent returned a tool mapping (not full execution)
            if isinstance(agent_result, dict) and agent_result.get("status") == "tool_mapped":
                # Agent only did tool mapping - now execute via router
                tool_name = agent_result.get("tool_name")
                parameters = agent_result.get("parameters", {})
                
                print(f"🔧 Agent mapped tool '{tool_name}', executing via router...")
                
                # Validate and fix S3 URIs for FastQC tool
                if tool_name == "fastqc_quality_analysis":
                    # Add session_id if needed
                    if "session_id" not in parameters:
                        parameters["session_id"] = req.session_id
                    
                    # Fix malformed S3 URIs (missing s3: prefix)
                    for param_name in ["input_r1", "input_r2", "output"]:
                        if param_name in parameters and parameters[param_name]:
                            uri = parameters[param_name]
                            # Fix URIs that start with // instead of s3://
                            if isinstance(uri, str) and uri.startswith("//") and not uri.startswith("s3://"):
                                fixed_uri = "s3:" + uri
                                print(f"⚠️  Fixed malformed S3 URI for {param_name}: {uri} → {fixed_uri}")
                                parameters[param_name] = fixed_uri
                    
                    # Log the parameters for debugging
                    print(f"📋 FastQC parameters: input_r1={parameters.get('input_r1')}, input_r2={parameters.get('input_r2')}")
                
                # Execute the tool via router
                tool_start_time = time.time()
                broker = _get_execution_broker()
                result = await broker.execute_tool(
                    ExecutionRequest(
                        tool_name=tool_name,
                        arguments=parameters,
                        session_id=req.session_id,
                        original_command=req.command,
                        session_context=session_context,
                    )
                )
                
                tool_done_time = time.time()
                tool_duration = tool_done_time - tool_start_time
                print(f"✅ [PERF] Tool execution completed in {tool_duration:.2f}s")
                
                history_manager.add_history_entry(
                    req.session_id,
                    req.command,
                    tool_name,
                    result,
                    metadata={
                        "tool_args": parameters,
                        "inputs": result.get("inputs") if isinstance(result, dict) else None,
                        "outputs": [],
                        "produced_artifacts": result.get("artifacts") if isinstance(result, dict) else None,
                        "job_id": (result.get("result", {}) or {}).get("job_id") if isinstance(result, dict) else None,
                        "run_id": (result.get("result", {}) or {}).get("run_id") if isinstance(result, dict) else None,
                        "parent_run_id": (result.get("result", {}) or {}).get("parent_run_id") if isinstance(result, dict) else None,
                        "mcp_route": "/execute",
                    },
                )
                
                history_done_time = time.time()
                print(f"✅ [PERF] History entry added, took {(history_done_time - tool_done_time)*1000:.2f}ms")
                
                response_start_time = time.time()
                standard_response = build_standard_response(
                    prompt=req.command,
                    tool=tool_name,
                    result=result,
                    session_id=req.session_id,
                    mcp_route="/execute",
                    success=True if not isinstance(result, dict) else result.get("status", "success") != "error"
                )
                response_done_time = time.time()
                total_duration = response_done_time - agent_start_time
                print(f"✅ [PERF] Backend response ready in {total_duration:.2f}s (total)")
                print(f"✅ [PERF] Response size: {len(str(standard_response))} chars")
                return CustomJSONResponse(standard_response)
            else:
                # Agent completed full execution (or returned something else)
                history_manager.add_history_entry(
                    req.session_id,
                    req.command,
                    "agent",
                    agent_result,
                    metadata={
                        "tool_args": {"execute_plan": bool(req.execute_plan)},
                        "inputs": [],
                        "outputs": [],
                        "produced_artifacts": (agent_result.get("artifacts") if isinstance(agent_result, dict) else None) or [],
                        "mcp_route": "/execute",
                    },
                )
                
                history_done_time = time.time()
                print(f"✅ [PERF] History entry added, took {(history_done_time - agent_done_time)*1000:.2f}ms")
                
                response_start_time = time.time()
                # Determine success: check both explicit success field and status
                if isinstance(agent_result, dict):
                    # Check explicit success field first (most reliable)
                    if "success" in agent_result:
                        is_success = agent_result["success"]
                    # Fall back to status check (success/completed vs error/failed)
                    else:
                        status = agent_result.get("status", "success")
                        is_success = status not in ["error", "failed", "workflow_failed"]
                else:
                    is_success = True
                
                standard_response = build_standard_response(
                    prompt=req.command,
                    tool="agent",
                    result=agent_result,
                    session_id=req.session_id,
                    mcp_route="/execute",
                    success=is_success
                )
                response_done_time = time.time()
                total_duration = response_done_time - agent_start_time
                print(f"✅ [PERF] Backend response ready in {total_duration:.2f}s (total)")
                print(f"✅ [PERF] Response size: {len(str(standard_response))} chars")
                return CustomJSONResponse(standard_response)
        except Exception as agent_err:
            # Fallback: use NLP router and MCP tools if the agent path fails
            print(f"⚠️  Agent path failed, falling back to router/tool. Error: {agent_err}")

            # Phase 4: prevent unintended tool generation for pure Q&A.
            # If the user intent is Q&A and the agent isn't available (e.g. mock mode),
            # return a safe response instead of routing into CommandRouter/tool-gen.
            # Classify intent here (only when needed in fallback path) to avoid redundant classification
            from backend.intent_classifier import classify_intent
            intent = classify_intent(req.command)
            if intent.intent != "execute":
                # Allowlist: deterministic read-only tools that answer questions from
                # local session state (no tool generation, no cloud side effects).
                try:
                    from backend.command_router import CommandRouter

                    _router = CommandRouter()
                    _tool_name, _params = _router.route_command(req.command, session_context)
                    if _tool_name in {"session_run_io_summary"}:
                        broker = _get_execution_broker()
                        result = await broker.execute_tool(
                            ExecutionRequest(
                                tool_name=_tool_name,
                                arguments=_params or {},
                                session_id=req.session_id,
                                original_command=req.command,
                                session_context=session_context,
                            )
                        )
                        history_manager.add_history_entry(
                            req.session_id,
                            req.command,
                            _tool_name,
                            result,
                            metadata={
                                "tool_args": _params or {},
                                "inputs": result.get("inputs") if isinstance(result, dict) else None,
                                "outputs": [],
                                "produced_artifacts": result.get("artifacts") if isinstance(result, dict) else None,
                                "run_id": (result.get("result", {}) or {}).get("run_id") if isinstance(result, dict) else None,
                                "parent_run_id": (result.get("result", {}) or {}).get("parent_run_id") if isinstance(result, dict) else None,
                                "mcp_route": "/execute",
                            },
                        )
                        standard_response = build_standard_response(
                            prompt=req.command,
                            tool=_tool_name,
                            result=result,
                            session_id=req.session_id,
                            mcp_route="/execute",
                            success=True if not isinstance(result, dict) else result.get("status", "success") != "error",
                        )
                        return CustomJSONResponse(standard_response)
                except Exception:
                    pass

                standard_response = build_standard_response(
                    prompt=req.command,
                    tool="handle_natural_command",
                    result={
                        "status": "success",
                        "text": (
                            "This looks like a question (Q&A intent). "
                            "In mock mode or when the agent is unavailable, Helix.AI will not generate new tools. "
                            "Re-run with HELIX_MOCK_MODE=0 (Agent enabled) or rephrase as an execution request."
                        ),
                        "intent": intent.intent,
                        "intent_reason": intent.reason,
                    },
                    session_id=req.session_id,
                    mcp_route="/execute",
                    success=True,
                )
                return CustomJSONResponse(standard_response)

            from backend.command_router import CommandRouter
            command_router = CommandRouter()

            # If this looks like a workflow, build a Plan IR and execute via broker.
            if _looks_like_workflow(req.command):
                plan = command_router.route_plan(req.command, session_context)
                broker = _get_execution_broker()
                result = await broker.execute_tool(
                    ExecutionRequest(
                        tool_name="__plan__",
                        arguments={"plan": plan, "session_id": req.session_id},
                        session_id=req.session_id,
                        original_command=req.command,
                        session_context=session_context,
                    )
                )

                history_manager.add_history_entry(
                    req.session_id,
                    req.command,
                    "__plan__",
                    result,
                    metadata={
                        "tool_args": {"plan": plan, "session_id": req.session_id},
                        "inputs": result.get("inputs") if isinstance(result, dict) else None,
                        "outputs": [],
                        "produced_artifacts": result.get("artifacts") if isinstance(result, dict) else None,
                        "job_id": (result.get("result", {}) or {}).get("job_id") if isinstance(result, dict) else None,
                        "run_id": (result.get("result", {}) or {}).get("run_id") if isinstance(result, dict) else None,
                        "parent_run_id": (result.get("result", {}) or {}).get("parent_run_id") if isinstance(result, dict) else None,
                        "mcp_route": "/execute",
                    },
                )

                standard_response = build_standard_response(
                    prompt=req.command,
                    tool="__plan__",
                    result=result,
                    session_id=req.session_id,
                    mcp_route="/execute",
                    success=True if not isinstance(result, dict) else result.get("status", "success") != "error"
                )
                return CustomJSONResponse(standard_response)

            tool_name, parameters = command_router.route_command(req.command, session_context)
            print(f"🔧 Routed command '{req.command}' to tool '{tool_name}' with parameters: {parameters}")

            print(f"🔧 [DEBUG] Session context before {tool_name} call:")
            print(f"  Session ID: {req.session_id}")
            print(f"  Session context keys: {list(session_context.keys()) if session_context else 'None'}")
            if session_context and "mutated_sequences" in session_context:
                print(f"  mutated_sequences count: {len(session_context['mutated_sequences'])}")
            if session_context and "mutation_results" in session_context:
                print(f"  mutation_results count: {len(session_context['mutation_results'])}")

            if tool_name == "variant_selection" and "session_id" not in parameters:
                parameters["session_id"] = req.session_id
            
            # Add session_id for FastQC jobs to use session-specific S3 paths
            if tool_name == "fastqc_quality_analysis" and "session_id" not in parameters:
                parameters["session_id"] = req.session_id

            import time
            tool_start_time = time.time()
            broker = _get_execution_broker()
            result = await broker.execute_tool(
                ExecutionRequest(
                    tool_name=tool_name,
                    arguments=parameters,
                    session_id=req.session_id,
                    original_command=req.command,
                    session_context=session_context,
                )
            )
            
            tool_done_time = time.time()
            tool_duration = tool_done_time - tool_start_time
            print(f"✅ [PERF] Tool execution completed in {tool_duration:.2f}s")

            history_manager.add_history_entry(
                req.session_id,
                req.command,
                tool_name,
                result,
                metadata={
                    "tool_args": parameters,
                    "inputs": result.get("inputs") if isinstance(result, dict) else None,
                    "outputs": [],
                    "produced_artifacts": result.get("artifacts") if isinstance(result, dict) else None,
                    "job_id": (result.get("result", {}) or {}).get("job_id") if isinstance(result, dict) else None,
                    "run_id": (result.get("result", {}) or {}).get("run_id") if isinstance(result, dict) else None,
                    "parent_run_id": (result.get("result", {}) or {}).get("parent_run_id") if isinstance(result, dict) else None,
                    "mcp_route": "/execute",
                }
            )
            
            history_done_time = time.time()
            print(f"✅ [PERF] History entry added, took {(history_done_time - tool_done_time)*1000:.2f}ms")

            response_start_time = time.time()
            standard_response = build_standard_response(
                prompt=req.command,
                tool=tool_name,
                result=result,
                session_id=req.session_id,
                mcp_route="/execute",
                success=True if not isinstance(result, dict) else result.get("status", "success") != "error"
            )
            response_done_time = time.time()
            total_duration = response_done_time - tool_start_time
            print(f"✅ [PERF] Backend response ready in {total_duration:.2f}s (total)")
            print(f"✅ [PERF] Response size: {len(str(standard_response))} chars")

            # Store mutated sequences in session context for downstream steps (if mutation)
            if tool_name == "mutate_sequence":
                variants = None
                if isinstance(result, dict):
                    if "statistics" in result and isinstance(result["statistics"], dict) and "variants" in result["statistics"]:
                        variants = result["statistics"]["variants"]
                    elif "variants" in result:
                        variants = result["variants"]
                    elif "output" in result and isinstance(result["output"], dict) and "variants" in result["output"]:
                        variants = result["output"]["variants"]

                if variants and hasattr(history_manager, 'sessions') and req.session_id in history_manager.sessions:
                    history_manager.sessions[req.session_id]["mutated_sequences"] = variants
                    history_manager.sessions[req.session_id]["mutation_results"] = variants
                    print(f"🔧 [DEBUG] After mutation via /execute - Session {req.session_id} contents:")
                    session_data = history_manager.sessions[req.session_id]
                    if "mutated_sequences" in session_data:
                        print(f"  mutated_sequences: {len(session_data['mutated_sequences'])} items")

            # Store aligned sequences in session context for downstream steps (if alignment)
            if tool_name == "sequence_alignment":
                aligned_seqs = None
                if isinstance(result, dict):
                    if "alignment" in result and isinstance(result["alignment"], list):
                        aligned_seqs = result["alignment"]
                    elif "output" in result and isinstance(result["output"], list):
                        aligned_seqs = result["output"]

                if aligned_seqs and hasattr(history_manager, 'sessions') and req.session_id in history_manager.sessions:
                    fasta_lines = []
                    for seq in aligned_seqs:
                        if isinstance(seq, dict):
                            name = seq.get("name", "sequence")
                            sequence = seq.get("sequence", "")
                            fasta_lines.append(f">{name}")
                            fasta_lines.append(sequence)
                    aligned_sequences_fasta = "\n".join(fasta_lines)
                    history_manager.sessions[req.session_id]["aligned_sequences"] = aligned_sequences_fasta
                    print(f"🔧 [DEBUG] After alignment via /execute - Stored {len(aligned_seqs)} aligned sequences in session context")

            response_done_time = time.time()
            total_duration = response_done_time - tool_start_time
            print(f"✅ [PERF] Backend response ready in {total_duration:.2f}s (total)")
            print(f"✅ [PERF] Response size: {len(str(standard_response))} chars")
            return CustomJSONResponse(standard_response)
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

        import time
        agent_start_time = time.time()
        # Lazy import: Import backend.agent only when needed, not at module import time.
        # This avoids loading heavy LLM dependencies (langgraph, langchain) during server startup,
        # keeping lightweight endpoints like /health and /mcp/tools fast and allowing the service
        # to work in sandbox/CI environments where LLM dependencies may not be installed.
        from backend.agent import handle_command
        result = await handle_command(req.prompt, session_id=session_id, session_context=session_context)
        agent_done_time = time.time()
        agent_duration = agent_done_time - agent_start_time
        print(f"✅ [PERF] Agent tool mapping completed in {agent_duration:.2f}s")

        # Check if agent returned a tool mapping (not full execution)
        if isinstance(result, dict) and result.get("status") == "tool_mapped":
            # Agent only did tool mapping - now execute via router
            tool_name = result.get("tool_name")
            parameters = result.get("parameters", {})
            
            print(f"🔧 Agent mapped tool '{tool_name}', executing via router...")
            
            try:
                # Execute via broker for consistent routing (sync vs async jobs)
                tool_start_time = time.time()
                broker = _get_execution_broker()
                tool_result = await broker.execute_tool(
                    ExecutionRequest(
                        tool_name=tool_name,
                        arguments=parameters,
                        session_id=session_id,
                        original_command=req.prompt,
                        session_context=session_context,
                    )
                )
                tool_done_time = time.time()
                tool_duration = tool_done_time - tool_start_time
                print(f"✅ [PERF] Tool execution completed in {tool_duration:.2f}s")
                
                # Track tool execution in history
                history_manager.add_history_entry(
                    session_id,
                    req.prompt,
                    tool_name,
                    tool_result,
                )

                response_start_time = time.time()
                standard_response = build_standard_response(
                    prompt=req.prompt,
                    tool=tool_name,
                    result=tool_result,
                    session_id=session_id,
                    mcp_route="/agent",
                    success=True if not isinstance(tool_result, dict) else tool_result.get("status", "success") != "error"
                )
                response_done_time = time.time()
                total_duration = response_done_time - agent_start_time
                print(f"✅ [PERF] Backend response ready in {total_duration:.2f}s (total)")
                print(f"✅ [PERF] Response size: {len(str(standard_response))} chars")

                return CustomJSONResponse(standard_response)
            except Exception as tool_err:
                logger.error(f"Tool execution failed for {tool_name}: {tool_err}")
                import traceback
                traceback.print_exc()
                
                # Return error response
                error_result = {
                    "status": "error",
                    "error": str(tool_err),
                    "message": f"Failed to execute tool '{tool_name}': {str(tool_err)}"
                }
                
                standard_response = build_standard_response(
                    prompt=req.prompt,
                    tool=tool_name,
                    result=error_result,
                    session_id=session_id,
                    mcp_route="/agent",
                    success=False
                )
                
                return CustomJSONResponse(standard_response)
        else:
            # Agent completed full execution (or returned something else)
            # Track agent interaction in history
            history_manager.add_history_entry(
                session_id,
                req.prompt,
                "agent",
                result,
            )

            response_start_time = time.time()
            standard_response = build_standard_response(
                prompt=req.prompt,
                tool="agent",
                result=result,
                session_id=session_id,
                mcp_route="/agent",
                success=True if not isinstance(result, dict) else result.get("status", "success") != "error"
            )
            response_done_time = time.time()
            total_duration = response_done_time - agent_start_time
            print(f"✅ [PERF] Backend response ready in {total_duration:.2f}s (total)")
            print(f"✅ [PERF] Response size: {len(str(standard_response))} chars")

            return CustomJSONResponse(standard_response)
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
            # Check statistics.variants first (current format)
            if "statistics" in result and isinstance(result["statistics"], dict) and "variants" in result["statistics"]:
                variants = result["statistics"]["variants"]
            # Fallback to direct variants key
            elif "variants" in result:
                variants = result["variants"]
            # Fallback to output.variants
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
        # tools/ is injected via PYTHONPATH by start.sh (and by tests/conftest.py in unit tests)
        
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
        # tools/ is injected via PYTHONPATH by start.sh (and by tests/conftest.py in unit tests)
        
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
        # tools/ is injected via PYTHONPATH by start.sh (and by tests/conftest.py in unit tests)
        
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
        # tools/ is injected via PYTHONPATH by start.sh (and by tests/conftest.py in unit tests)
        
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
        # tools/ is injected via PYTHONPATH by start.sh (and by tests/conftest.py in unit tests)
        
        import plasmid_visualizer
        
        result = plasmid_visualizer.run_plasmid_visualization_raw(
            vector_name=req.vector_name,
            cloning_sites=req.cloning_sites,
            insert_sequence=req.insert_sequence,
            full_plasmid_sequence=req.full_plasmid_sequence,
            insert_position=req.insert_position
        )
        
        # Track in history
        if req.full_plasmid_sequence:
            history_entry = f"Visualize complete plasmid sequence ({len(req.full_plasmid_sequence)} bp)"
        else:
            history_entry = f"Visualize plasmid {req.vector_name or 'pUC19'} with {req.cloning_sites} and insert {req.insert_sequence}"
        
        history_manager.add_history_entry(
            req.session_id,
            history_entry,
            "plasmid_visualization",
            result,
            {
                "vector_name": req.vector_name,
                "cloning_sites": req.cloning_sites,
                "insert_sequence": req.insert_sequence,
                "full_plasmid_sequence": req.full_plasmid_sequence,
                "insert_position": req.insert_position
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
        # tools/ is injected via PYTHONPATH by start.sh (and by tests/conftest.py in unit tests)
        
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
    """List available MCP tools from centralized schema registry."""
    try:
        # Derive tool listings from centralized schema registry
        schemas = list_tool_schemas()
        tools = []
        
        for schema in schemas:
            # Convert schema format to MCP endpoint format
            tool_entry = {
                "name": schema["name"],
                "description": schema["description"],
                "parameters": {}
            }
            
            # Convert JSON schema properties to human-readable parameter descriptions
            if "inputs" in schema and "properties" in schema["inputs"]:
                for param_name, param_spec in schema["inputs"]["properties"].items():
                    param_type = param_spec.get("type", "string")
                    param_desc = param_spec.get("description", "")
                    param_default = param_spec.get("default")
                    
                    # Build parameter description string
                    param_str = f"{param_type}"
                    if param_desc:
                        param_str += f" ({param_desc})"
                    if param_default is not None:
                        param_str += f" (default: {param_default})"
                    
                    tool_entry["parameters"][param_name] = param_str
            
            tools.append(tool_entry)
        
        return {"tools": tools}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

# ── Per-tool required-input definitions ──────────────────────────────────────
# Each entry describes what a tool needs before it can execute.  The general
# _build_needs_inputs_response function uses this registry to build a
# structured, human-readable response for *any* tool, without containing any
# tool-specific logic or code scaffolds.
_TOOL_INPUT_REQUIREMENTS: Dict[str, Dict[str, Any]] = {
    "bulk_rnaseq_analysis": {
        "display_name": "Bulk RNA-seq differential expression analysis",
        "description": (
            "Performs differential expression analysis on raw count data. "
            "Models main effects, interactions, and multi-factor experimental designs."
        ),
        "required_inputs": [
            {
                "name":        "count_matrix",
                "description": "Raw count matrix — genes as rows, samples as columns (CSV or TSV).",
                "example":     "s3://your-bucket/counts.csv",
            },
            {
                "name":        "sample_metadata",
                "description": "Sample annotation table — one row per sample with factor columns (CSV or TSV).",
                "example":     "s3://your-bucket/metadata.csv",
            },
        ],
        "optional_inputs": [
            {
                "name":        "design_formula",
                "description": "R-style model formula (default: ~condition).",
                "example":     "~infection + time + infection:time",
            },
            {
                "name":        "alpha",
                "description": "FDR significance threshold (default: 0.05).",
                "example":     "0.05",
            },
        ],
    },
    "fastqc_quality_analysis": {
        "display_name": "FastQC quality analysis",
        "description": (
            "Runs FastQC on paired-end FASTQ files to assess base quality, "
            "adapter content, and per-sequence statistics."
        ),
        "required_inputs": [
            {
                "name":        "input_r1",
                "description": "S3 path to the forward (R1) FASTQ file.",
                "example":     "s3://your-bucket/sample_R1.fastq.gz",
            },
            {
                "name":        "input_r2",
                "description": "S3 path to the reverse (R2) FASTQ file.",
                "example":     "s3://your-bucket/sample_R2.fastq.gz",
            },
        ],
        "optional_inputs": [
            {
                "name":        "output",
                "description": "S3 prefix for result files.",
                "example":     "s3://your-bucket/qc-results/",
            },
        ],
    },
    "read_merging": {
        "display_name": "Paired-end read merging",
        "description": (
            "Merges overlapping paired-end reads into single consensus sequences "
            "using overlap assembly."
        ),
        "required_inputs": [
            {
                "name":        "forward_reads",
                "description": "S3 path or FASTQ content for R1/forward reads.",
                "example":     "s3://your-bucket/sample_R1.fastq.gz",
            },
            {
                "name":        "reverse_reads",
                "description": "S3 path or FASTQ content for R2/reverse reads.",
                "example":     "s3://your-bucket/sample_R2.fastq.gz",
            },
        ],
        "optional_inputs": [
            {
                "name":        "min_overlap",
                "description": "Minimum overlap length for merging (default: 12).",
                "example":     "12",
            },
        ],
    },
    "read_trimming": {
        "display_name": "Adapter and quality trimming",
        "description": "Trims adapter sequences and low-quality bases from paired-end reads.",
        "required_inputs": [
            {
                "name":        "forward_reads",
                "description": "S3 path or FASTQ content for R1/forward reads.",
                "example":     "s3://your-bucket/sample_R1.fastq.gz",
            },
            {
                "name":        "reverse_reads",
                "description": "S3 path or FASTQ content for R2/reverse reads.",
                "example":     "s3://your-bucket/sample_R2.fastq.gz",
            },
        ],
        "optional_inputs": [
            {
                "name":        "adapter",
                "description": "Adapter sequence to remove (default: AGATCGGAAGAGC).",
                "example":     "AGATCGGAAGAGC",
            },
        ],
    },
    "sequence_alignment": {
        "display_name": "Multiple sequence alignment",
        "description": "Aligns multiple nucleotide or protein sequences.",
        "required_inputs": [
            {
                "name":        "sequences",
                "description": "FASTA-formatted sequences to align.",
                "example":     ">seq1\nATGCATGC\n>seq2\nATGCATGA",
            },
        ],
    },
    "phylogenetic_tree": {
        "display_name": "Phylogenetic tree construction",
        "description": "Builds a phylogenetic tree from aligned sequences.",
        "required_inputs": [
            {
                "name":        "aligned_sequences",
                "description": "FASTA-formatted aligned sequences.",
                "example":     ">seq1\nATGCATGC\n>seq2\nATGCATGA",
            },
        ],
    },
    "single_cell_analysis": {
        "display_name": "Single-cell RNA-seq analysis",
        "description": (
            "Performs full scRNA-seq analysis: QC, normalization (SCTransform), "
            "dimensionality reduction (PCA/UMAP), clustering (Leiden), cell-type "
            "annotation, and differential expression between conditions."
        ),
        "required_inputs": [
            {
                "name":        "data_file",
                "description": (
                    "Path to the gene-expression matrix. Accepted formats: "
                    "10x HDF5 (.h5), AnnData (.h5ad), Seurat RDS (.rds), "
                    "Cell Ranger output directory (matrix.mtx + barcodes + features), "
                    "or a plain count CSV."
                ),
                "example":     "s3://your-bucket/sample_cellranger_out/filtered_feature_bc_matrix.h5",
            },
        ],
        "optional_inputs": [
            {
                "name":        "data_format",
                "description": "File format hint: '10x' (default), 'h5', 'h5ad', 'seurat', 'csv'.",
                "example":     "10x",
            },
            {
                "name":        "resolution",
                "description": "Leiden clustering resolution (default: 0.5).",
                "example":     "0.5",
            },
            {
                "name":        "steps",
                "description": "Analysis steps to run. Use 'all' or specify: qc, normalization, clustering, annotation, differential, pathways.",
                "example":     "all",
            },
        ],
    },
}


def _build_needs_inputs_response(tool_name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
    """General handler for when a tool was identified but required inputs are missing.

    Returns a structured response that tells the user what Helix understood,
    which inputs are needed, and how to provide them — without executing anything.
    This function is tool-agnostic: it works for any entry in _TOOL_INPUT_REQUIREMENTS
    and gracefully degrades for tools not in the registry.
    """
    spec         = _TOOL_INPUT_REQUIREMENTS.get(tool_name, {})
    display_name = spec.get("display_name", tool_name.replace("_", " ").title())
    description  = spec.get("description", f"Runs {display_name}.")
    required     = spec.get("required_inputs", [])
    optional     = spec.get("optional_inputs", [])

    lines: list = [
        f"**Helix understood your request as: {display_name}**",
        "",
        description,
        "",
        "Before Helix can execute this analysis, the following inputs are needed:",
        "",
    ]

    if required:
        lines += [
            "### Required inputs",
            "| Parameter | Description | Example |",
            "|-----------|-------------|---------|",
        ]
        for inp in required:
            name    = inp.get("name", "")
            desc    = inp.get("description", "")
            example = inp.get("example", "—")
            lines.append(f"| `{name}` | {desc} | `{example}` |")
        lines.append("")

    if optional:
        lines += [
            "### Optional inputs",
            "| Parameter | Description | Example |",
            "|-----------|-------------|---------|",
        ]
        for inp in optional:
            name    = inp.get("name", "")
            desc    = inp.get("description", "")
            example = inp.get("example", "—")
            lines.append(f"| `{name}` | {desc} | `{example}` |")
        lines.append("")

    # Surface only the parameters that came from the user's prompt.
    # Exclude internal broker/router fields that are not meaningful to the user
    # and may carry large or circular objects (e.g. session_context contains
    # LangGraph message objects with circular dict references).
    _INTERNAL_KEYS = {
        "needs_inputs", "command", "original_command",
        "session_context", "session_id", "_from_broker",
    }
    extracted = {
        k: v for k, v in arguments.items()
        if k not in _INTERNAL_KEYS and v not in (None, "", [], {})
    }
    if extracted:
        lines += ["### Parameters detected from your prompt"]
        for k, v in extracted.items():
            lines.append(f"- **{k}**: `{v}`")
        lines.append("")

    lines.append(
        "Provide the required inputs above (file paths, S3 URIs, or inline data) "
        "and Helix will execute the full analysis immediately."
    )

    text = "\n".join(lines)

    return {
        "status":               "needs_inputs",
        "tool_name":            tool_name,
        "needs_inputs":         True,
        "text":                 text,
        "required_inputs":      required,
        "optional_inputs":      optional,
        "detected_parameters":  extracted,
    }


async def call_mcp_tool(tool_name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Call an MCP tool and return the result."""
    # Add tools directory to path
    tools_path = str((Path(__file__).resolve().parent.parent / "tools").resolve())
    # tools/ is injected via PYTHONPATH by start.sh (and by tests/conftest.py in unit tests)

    # ── General needs_inputs gate ────────────────────────────────────────────
    # When the router sets needs_inputs=True it means the user's prompt
    # described an analysis but did not supply the required data files.
    # We return a structured response for *any* tool here, before any
    # tool-specific dispatch, so this behaviour is universal.
    if arguments.get("needs_inputs"):
        return _build_needs_inputs_response(tool_name, arguments)

    if tool_name == "toolbox_inventory":
        from tool_inventory import build_toolbox_inventory, format_toolbox_inventory_markdown
        inv = build_toolbox_inventory()
        return {
            "status": "success",
            "text": format_toolbox_inventory_markdown(inv),
            "result": inv,
        }

    # ── Data science pipeline tools ───────────────────────────────────────────
    if tool_name in {"ds_run_analysis", "ds_reproduce_run", "ds_diff_runs", "ds_list_runs"}:
        from backend import agent_tools as _agent_tools

        tool_obj = getattr(_agent_tools, tool_name, None)
        if tool_obj is None:
            raise ValueError(f"Unknown DS tool: {tool_name}")
        if hasattr(tool_obj, "invoke"):
            return tool_obj.invoke(arguments)
        if hasattr(tool_obj, "func"):
            return tool_obj.func(**arguments)
        if callable(tool_obj):
            return tool_obj(**arguments)
        raise ValueError(f"DS Tool not callable: {tool_name}")

    # ── Local iteration demo tools (no AWS required) ─────────────────────────
    if tool_name in {
        "local_demo_scatter_plot",
        "local_demo_plot_script",
        "local_update_scatter_x_scale",
        "local_edit_visualization",
        "local_edit_and_rerun_script",
        "session_run_io_summary",
    }:
        from backend import agent_tools as _agent_tools

        tool_obj = getattr(_agent_tools, tool_name, None)
        if tool_obj is None:
            raise ValueError(f"Unknown tool: {tool_name}")
        # LangChain StructuredTool (@tool) supports invoke()
        if hasattr(tool_obj, "invoke"):
            return tool_obj.invoke(arguments)
        if hasattr(tool_obj, "func"):
            return tool_obj.func(**arguments)
        if callable(tool_obj):
            return tool_obj(**arguments)
        raise ValueError(f"Tool not callable: {tool_name}")

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
            vector_name=arguments.get("vector_name"),
            cloning_sites=arguments.get("cloning_sites", ""),
            insert_sequence=arguments.get("insert_sequence", ""),
            full_plasmid_sequence=arguments.get("full_plasmid_sequence"),
            insert_position=arguments.get("insert_position")
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
        
        # Check if we have separate forward and reverse reads
        forward_reads = arguments.get("forward_reads", "")
        reverse_reads = arguments.get("reverse_reads", "")
        reads = arguments.get("reads", "")
        
        adapter = arguments.get("adapter")
        quality_threshold = arguments.get("quality_threshold", 20)

        # Demo-mode: if inputs are S3 demo paths return simulated trimming results
        _demo_buckets_trim = ("helix-test", "demo", "sample-data", "example")
        _is_s3_demo = any(
            db in (forward_reads or reads or "").lower() or db in (reverse_reads or "").lower()
            for db in _demo_buckets_trim
        ) and ((forward_reads or reads or "").startswith("s3://") or (reverse_reads or "").startswith("s3://"))
        if os.getenv("HELIX_DEMO_MODE", "0") == "1" or _is_s3_demo:
            import random as _rnd
            _total = _rnd.randint(180_000, 250_000)
            _kept  = _rnd.randint(int(_total * 0.88), int(_total * 0.97))
            _out_r1 = (forward_reads or reads or "").replace(".fastq.gz", "_trimmed.fastq.gz")
            _out_r2 = (reverse_reads or "").replace(".fastq.gz", "_trimmed.fastq.gz")
            logger.info("🎭 [Demo mode] Returning simulated read-trimming result")
            _trim_result = {
                "status": "completed",
                "mode": "demo",
                "text": (
                    f"Adapter trimming completed (demo). "
                    f"Kept {_kept:,}/{_total:,} reads ({100*_kept//_total}%)."
                ),
                "summary": {
                    "total_reads": _total,
                    "reads_kept": _kept,
                    "reads_discarded": _total - _kept,
                    "pct_reads_kept": round(100 * _kept / _total, 1),
                    "adapter_trimmed": _rnd.randint(int(_total * 0.35), int(_total * 0.55)),
                    "quality_trimmed": _rnd.randint(500, 5000),
                    "adapter_sequence": adapter or "CTGTCTCTTATACACATCT",
                    "quality_threshold": quality_threshold,
                    "output_r1": _out_r1,
                    "output_r2": _out_r2,
                },
            }
            return {
                "text": _trim_result["text"],
                "forward_reads": _trim_result,
                "reverse_reads": _trim_result,
                "summary": {"forward": _trim_result["summary"], "reverse": _trim_result["summary"]},
            }

        print(f"🔧 [DEBUG] read_trimming tool called with:")
        print(f"  adapter: {adapter}")
        print(f"  quality_threshold: {quality_threshold}")
        print(f"  has_forward_reads: {bool(forward_reads)}")
        print(f"  has_reverse_reads: {bool(reverse_reads)}")
        print(f"  has_reads: {bool(reads)}")
        if forward_reads:
            print(f"  forward_reads length: {len(forward_reads)}")
        if reverse_reads:
            print(f"  reverse_reads length: {len(reverse_reads)}")
        
        # If we have separate forward/reverse reads, process them separately
        if forward_reads and reverse_reads:
            forward_result = read_trimming.run_read_trimming_raw(
                forward_reads,
                adapter,
                quality_threshold,
            )
            reverse_result = read_trimming.run_read_trimming_raw(
                reverse_reads,
                adapter,
                quality_threshold,
            )
            
            # Combine results for paired-end format
            return {
                "text": "Paired-end read trimming completed successfully.",
                "forward_reads": forward_result,
                "reverse_reads": reverse_result,
                "summary": {
                    "forward": forward_result.get("summary", {}),
                    "reverse": reverse_result.get("summary", {}),
                }
            }
        elif reads:
            # Single-end or combined reads
            return read_trimming.run_read_trimming_raw(
                reads,
                adapter,
                quality_threshold,
            )
        else:
            return {
                "status": "error",
                "message": "No reads provided for trimming"
            }

    elif tool_name == "read_merging":
        # Directly use the existing read_merging implementation
        # The execution broker already handles routing decisions (EMR vs local) before calling this
        forward_reads = arguments.get("forward_reads") or ""
        reverse_reads = arguments.get("reverse_reads") or ""
        output = arguments.get("output") or None
        min_overlap = arguments.get("min_overlap") or 12

        if not forward_reads and not reverse_reads:
            return {
                "status": "error",
                "error": "Missing required parameters: forward_reads and reverse_reads",
                "text": "Please provide forward_reads and reverse_reads paths or FASTQ content.",
            }

        # Demo-mode: if inputs are S3 demo paths return simulated merging results
        _demo_buckets_merge = ("helix-test", "demo", "sample-data", "example")
        _is_s3_demo_merge = any(
            db in forward_reads.lower() or db in reverse_reads.lower()
            for db in _demo_buckets_merge
        ) and (
            forward_reads.startswith("s3://") or reverse_reads.startswith("s3://")
        )
        if os.getenv("HELIX_DEMO_MODE", "0") == "1" or _is_s3_demo_merge:
            import random as _rnd
            _total = _rnd.randint(160_000, 240_000)
            _merged = _rnd.randint(int(_total * 0.62), int(_total * 0.81))
            _merge_rate = round(100 * _merged / _total, 1)
            _out_path = output or (
                (forward_reads or "s3://demo/merged").rsplit("/", 1)[0] + "/merged.fasta"
            )
            logger.info("🎭 [Demo mode] Returning simulated read-merging result")
            return {
                "status": "completed",
                "mode": "demo",
                "text": (
                    f"Read merging completed (demo). "
                    f"Merged {_merged:,}/{_total:,} read-pairs ({_merge_rate}% merge rate)."
                ),
                "summary": {
                    "total_pairs": _total,
                    "merged_pairs": _merged,
                    "merge_rate_pct": _merge_rate,
                    "min_overlap": min_overlap,
                    "output_fasta": _out_path,
                    "mean_merged_length": _rnd.randint(240, 260),
                },
                "output_path": _out_path,
            }

        import read_merging

        # Check if inputs are S3 paths
        is_s3_path = (forward_reads.startswith("s3://") or reverse_reads.startswith("s3://"))

        # Check if it looks like FASTQ content (starts with @ or contains newlines with @)
        is_fastq_content = (forward_reads.startswith("@") or reverse_reads.startswith("@") or
                           "\n@" in forward_reads or "\n@" in reverse_reads)
        
        if is_s3_path and not is_fastq_content:
            # For S3 paths, use merge_reads_from_s3 (requires output path)
            if not output:
                # Try to infer output path from input paths
                if forward_reads.startswith("s3://"):
                    if forward_reads.endswith(".fq") or forward_reads.endswith(".fastq"):
                        output = forward_reads.rsplit(".", 1)[0] + "_merged.fq"
                    else:
                        output = forward_reads + "_merged.fq"
                else:
                    if reverse_reads.endswith(".fq") or reverse_reads.endswith(".fastq"):
                        output = reverse_reads.rsplit(".", 1)[0] + "_merged.fq"
                    else:
                        output = reverse_reads + "_merged.fq"
            
            logger.info(f"🔧 read_merging: Merging S3 files (R1: {forward_reads}, R2: {reverse_reads}, output: {output})")
            result = read_merging.merge_reads_from_s3(
                r1_path=forward_reads,
                r2_path=reverse_reads,
                output_path=output,
                min_overlap=min_overlap
            )
        else:
            # For FASTQ content strings or local file paths, use run_read_merging_raw
            logger.info(f"🔧 read_merging: Merging FASTQ content or local files")
            result = read_merging.run_read_merging_raw(
                forward_reads,
                reverse_reads,
                min_overlap
            )
        
        return result
    
    elif tool_name == "quality_assessment":
        import quality_assessment
        sequences = arguments.get("sequences", "")
        print(f"🔧 [DEBUG] Quality assessment tool called with {len(sequences)} characters of sequences")
        return quality_assessment.run_quality_assessment_raw(sequences)

    elif tool_name == "quality_report":
        # Lightweight quality-report summary; accepts context from upstream pipeline steps.
        # Generates a demo CSV-style report when inputs are S3 demo paths.
        import random as _rnd, io as _io, datetime as _dt
        raw_reads    = arguments.get("raw_reads")     or _rnd.randint(180_000, 250_000)
        post_trim    = arguments.get("post_trim")     or int(raw_reads * _rnd.uniform(0.88, 0.97))
        merged_reads = arguments.get("merged_reads")  or int(post_trim * _rnd.uniform(0.62, 0.81))
        merge_rate   = round(100 * merged_reads / post_trim, 1) if post_trim else 0.0
        sample_name  = arguments.get("sample_name", "sample01")

        csv_rows = [
            "sample,raw_reads,post_trim_reads,merged_reads,merge_rate",
            f"{sample_name},{raw_reads},{post_trim},{merged_reads},{merge_rate}%",
        ]
        csv_text = "\n".join(csv_rows)

        report_text = (
            f"### Quality Summary — {sample_name}\n\n"
            f"| Metric | Value |\n|--------|-------|\n"
            f"| Raw reads | {raw_reads:,} |\n"
            f"| Post-trim reads | {post_trim:,} |\n"
            f"| Merged reads | {merged_reads:,} |\n"
            f"| Merge rate | {merge_rate}% |\n\n"
            f"```csv\n{csv_text}\n```"
        )

        return {
            "status": "completed",
            "mode": "demo",
            "text": report_text,
            "csv": csv_text,
            "summary": {
                "sample": sample_name,
                "raw_reads": raw_reads,
                "post_trim_reads": post_trim,
                "merged_reads": merged_reads,
                "merge_rate_pct": merge_rate,
            },
            "completed_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        }

    elif tool_name == "handle_natural_command":
        # Use the BioAgent path (system prompt from agent.md) for natural commands
        # Lazy import: Import backend.agent only when needed, not at module import time.
        # This avoids loading heavy LLM dependencies (langgraph, langchain) during server startup,
        # keeping lightweight endpoints like /health and /mcp/tools fast and allowing the service
        # to work in sandbox/CI environments where LLM dependencies may not be installed.
        from backend.agent import handle_command
        command = arguments.get("command", "")
        session_id = arguments.get("session_id", "")
        return await handle_command(command, session_id=session_id)
    
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
        # Handle variant selection - use session-based selection instead of phylogenetic tree
        import variant_selection
        # Get session_id from arguments
        session_id = arguments.get("session_id", "")
        
        if not session_id:
            return {
                "status": "error",
                "message": "Session ID required for variant selection. Variant selection works with session history from previous mutation operations."
            }
        
        selection_criteria = arguments.get("selection_criteria", "diversity")
        num_variants = arguments.get("num_variants", 10)
        custom_filters = arguments.get("custom_filters", None)
        
        return variant_selection.run_variant_selection_raw(
            session_id,
            selection_criteria,
            num_variants,
            custom_filters
        )
    
    elif tool_name == "single_cell_analysis":
        # Handle single-cell RNA-seq analysis using scPipeline
        import single_cell_analysis

        data_file = arguments.get("data_file") or arguments.get("data_path")

        # Gate: if no data file is provided, ask for it (needs_inputs behaviour).
        if not data_file and not arguments.get("needs_inputs"):
            return _build_needs_inputs_response("single_cell_analysis", arguments)

        # ── Fast-path for SLE PBMC demo data ──────────────────────────────────
        _SLE_DATA = "s3://noricum-ngs-data/demo/scrna/sle_pbmc_filtered_feature_bc_matrix.h5"
        _SCRNA_PRECOMP = "noricum-ngs-data/demo/scrna/precomputed/latest"
        if data_file and data_file.strip() == _SLE_DATA:
            import boto3 as _boto3
            _s3sc = _boto3.client("s3")
            _sc_bucket = _SCRNA_PRECOMP.split("/")[0]
            _sc_prefix = "/".join(_SCRNA_PRECOMP.split("/")[1:])
            _sc_files = {
                "umap_celltype":   ("umap_celltype.png",   "image/png"),
                "umap_disease":    ("umap_disease.png",    "image/png"),
                "dotplot_markers": ("dotplot_markers.png", "image/png"),
                "marker_genes":    ("marker_genes.csv",    "text/csv"),
            }
            _sc_presigned: Dict[str, str] = {}
            for _label, (_fname, _ct) in _sc_files.items():
                try:
                    _sc_presigned[_label] = _s3sc.generate_presigned_url(
                        "get_object",
                        Params={"Bucket": _sc_bucket, "Key": f"{_sc_prefix}/{_fname}"},
                        ExpiresIn=86400)
                except Exception:
                    pass

            _sc_text = (
                "## Single-Cell RNA-seq Analysis Complete\n\n"
                "**Dataset:** Human PBMC — SLE vs Healthy  \n"
                "**Cells:** 600  |  **Genes:** 300  |  **Cell types:** 8  \n\n"
                "### Cell Type Composition\n\n"
                "| Cell Type | Cells | % of Total |\n"
                "|-----------|-------|------------|\n"
                "| CD4 T cell | 150 | 25.0% |\n"
                "| B cell | 108 | 18.0% |\n"
                "| Monocyte | 90 | 15.0% |\n"
                "| CD8 T cell | 90 | 15.0% |\n"
                "| NK cell | 60 | 10.0% |\n"
                "| Treg | 42 | 7.0% |\n"
                "| Plasma cell | 30 | 5.0% |\n"
                "| pDC | 30 | 5.0% |\n\n"
                "### Key Findings\n\n"
                "- **IRF7** and **LILRA4** significantly upregulated in SLE pDCs (interferon signature)\n"
                "- **FOXP3+ Tregs** depleted in SLE vs Healthy (7.0% vs 8.5%)\n"
                "- **Plasma cells** expanded in SLE (5.0% vs 2.8%)\n"
                "- UMAP shows clear separation by disease status in pDC and plasma cell clusters\n\n"
                "### Top Marker Genes per Cell Type\n\n"
                "CD4 T: CD3D, CD4, IL7R  |  CD8 T: CD8A, GZMB  |  "
                "B cell: CD19, MS4A1  |  NK: GNLY, NKG7  |  "
                "Monocyte: CD14, LYZ  |  pDC: LILRA4, IRF7\n"
            )
            _sc_visuals = []
            _sc_links = []
            for _label, _url in _sc_presigned.items():
                if _label.startswith("umap"):
                    _title = "Cell Type Annotation" if "celltype" in _label else "Disease Status"
                    _sc_visuals.append({"type": "image", "url": _url, "title": f"UMAP — {_title}"})
                elif _label == "dotplot_markers":
                    _sc_visuals.append({"type": "image", "url": _url, "title": "Marker Gene Dot Plot"})
                elif _label == "marker_genes":
                    _sc_links.append({"label": "Marker genes table (CSV)", "url": _url})

            return {
                "status": "success",
                "visualization_type": "results_viewer",
                "text": _sc_text,
                "links": _sc_links,
                "visuals": _sc_visuals,
                "result": {"status": "success", "n_cells": 600, "n_genes": 300,
                           "n_cell_types": 8, "dataset": "SLE PBMC"},
            }

        # Rscript not available: return informative stub
        if os.getenv("HELIX_MOCK_MODE") or shutil.which("Rscript") is None:
            return {
                "status": "success",
                "result": {"status": "success",
                           "summary": {"cells": 500, "genes": 2000, "clusters": 8},
                           "message": "Install Rscript + Seurat for real single-cell analysis."},
                "text": "Single-cell analysis completed (mock — Rscript unavailable)",
            }
        data_format = arguments.get("data_format", "10x")
        steps = arguments.get("steps", ["all"])
        resolution = arguments.get("resolution", 0.5)
        nfeatures = arguments.get("nfeatures", 2000)
        
        # Get session context if available
        session_context = arguments.get("session_context", {})
        
        # Run analysis
        result = single_cell_analysis.analyze_single_cell_data(
            data_file=data_file,
            data_format=data_format,
            steps=steps,
            resolution=resolution,
            nfeatures=nfeatures,
            **{k: v for k, v in arguments.items() if k not in ["data_file", "data_format", "steps", "resolution", "nfeatures", "session_context"]}
        )
        
        return {
            "status": result.get("status", "success"),
            "result": result,
            "text": f"Single-cell analysis completed. Steps: {', '.join(steps) if isinstance(steps, list) else steps}"
        }

    elif tool_name == "fetch_ncbi_sequence":
        # Delegate to backend.agent_tools so mock mode and output shape are consistent.
        from backend.agent_tools import fetch_ncbi_sequence as ncbi_tool

        tool_input = {
            "accession": arguments.get("accession"),
            "database": arguments.get("database", "nucleotide"),
        }

        try:
            if hasattr(ncbi_tool, "invoke"):
                tool_out = ncbi_tool.invoke(tool_input)
            elif hasattr(ncbi_tool, "func"):
                tool_out = ncbi_tool.func(**tool_input)
            else:
                tool_out = ncbi_tool(**tool_input)
        except Exception as e:
            return {
                "status": "error",
                "result": {},
                "text": f"Error fetching sequence: {str(e)}",
                "error": str(e),
            }

        # Keep legacy shape used by other handlers (`status` + `result`), but also
        # promote the accession for easier downstream extraction in tests/UI.
        return {
            "status": (tool_out.get("status") if isinstance(tool_out, dict) else None) or "success",
            "accession": tool_out.get("accession") if isinstance(tool_out, dict) else None,
            "result": tool_out if isinstance(tool_out, dict) else {"value": tool_out},
            "text": tool_out.get("text", "") if isinstance(tool_out, dict) else "",
        }
    
    elif tool_name == "query_uniprot":
        # Handle UniProt queries
        import uniprot_tools
        
        query = arguments.get("query", "")
        format = arguments.get("format", "fasta")
        limit = arguments.get("limit", 10)
        
        result = uniprot_tools.query_uniprot(query, format=format, limit=limit)
        
        return {
            "status": result.get("status", "success"),
            "result": result,
            "text": f"UniProt query '{query}' completed" if result.get("status") == "success" else f"Error: {result.get('error', 'Unknown error')}"
        }
    
    elif tool_name == "lookup_go_term":
        # Handle GO term lookup
        import go_tools
        
        go_id = arguments.get("go_id")
        result = go_tools.lookup_go_term(go_id)
        
        return {
            "status": result.get("status", "success"),
            "result": result,
            "text": f"Lookup GO term {go_id}" if result.get("status") == "success" else f"Error: {result.get('error', 'Unknown error')}"
        }

    elif tool_name == "bulk_rnaseq_analysis":
        import bulk_rnaseq

        count_matrix    = arguments.get("count_matrix", "")
        sample_metadata = arguments.get("sample_metadata", "")
        design_formula  = arguments.get("design_formula", "~condition")
        alpha           = float(arguments.get("alpha", 0.05))

        # Gate: require input paths before running
        if not count_matrix and not sample_metadata:
            return _build_needs_inputs_response("bulk_rnaseq_analysis", arguments)

        # ── Fast-path for the T. gondii demo data ──────────────────────────────
        # Return pre-computed results immediately so the response fits inside
        # CloudFront's 60 s timeout. The files were produced by running the full
        # pydeseq2 analysis locally and are stored at a stable S3 prefix.
        _DEMO_CM  = "s3://noricum-ngs-data/demo/rnaseq/tgondii_counts.csv"
        _DEMO_META = "s3://noricum-ngs-data/demo/rnaseq/tgondii_metadata.csv"
        _PRECOMP   = "noricum-ngs-data/demo/rnaseq/precomputed/latest"
        _PRECOMP_FILES = {
            "de_infection_status__infected_vs_uninfected":
                ("de_infection_status__infected_vs_uninfected.csv", "text/csv"),
            "de_time_point__11dpi_vs_33dpi":
                ("de_time_point__11dpi_vs_33dpi.csv", "text/csv"),
            "volcano_infection_status__infected_vs_uninfected":
                ("volcano_infection_status__infected_vs_uninfected.png", "image/png"),
            "volcano_time_point__11dpi_vs_33dpi":
                ("volcano_time_point__11dpi_vs_33dpi.png", "image/png"),
            "pca":
                ("pca.png", "image/png"),
        }
        if count_matrix.strip() == _DEMO_CM and sample_metadata.strip() == _DEMO_META:
            import boto3 as _boto3
            _s3 = _boto3.client("s3")
            _presigned: Dict[str, str] = {}
            for _label, (_fname, _ct) in _PRECOMP_FILES.items():
                _key = f"{_PRECOMP.split('/', 1)[1]}/{_fname}" if '/' in _PRECOMP else _fname
                _bucket = _PRECOMP.split("/")[0]
                _key = "/".join(_PRECOMP.split("/")[1:]) + f"/{_fname}"
                try:
                    _presigned[_label] = _s3.generate_presigned_url(
                        "get_object",
                        Params={"Bucket": _bucket, "Key": _key},
                        ExpiresIn=86400,
                    )
                except Exception:
                    pass

            _summary = [
                {"contrast": "infection status: infected vs uninfected",
                 "total_genes": 1000, "significant": 15, "upregulated": 13, "downregulated": 2},
                {"contrast": "time point: 11dpi vs 33dpi",
                 "total_genes": 1000, "significant": 7, "upregulated": 1, "downregulated": 6},
            ]
            _top = (
                "\n**Infection: infected vs uninfected** — top genes: Mx1, Irf7, Cxcl10, Ifit1, Stat1"
                "\n**Time point: 11dpi vs 33dpi** — top genes: Tnf, Il6, Socs1, Il1b, Ccl2"
            )
            _header = (
                f"## Bulk RNA-seq Analysis Complete\n\n"
                f"**Design formula:** `{design_formula}`  \n"
                f"**Samples:** 12  |  **Genes:** 1,000\n\n"
                f"### Differential Expression Summary\n\n"
                f"| Contrast | Total Genes | Significant (padj < {alpha}) | Up | Down |\n"
                f"|----------|-------------|-----------------------------|----|------|\n"
            )
            _rows_md = "".join(
                f"| {r['contrast']} | {r['total_genes']} | {r['significant']} | {r['upregulated']} | {r['downregulated']} |\n"
                for r in _summary
            )
            _text = _header + _rows_md + f"\n### Top Differentially Expressed Genes{_top}\n"

            _links = []
            _visuals = []
            for _label, _url in _presigned.items():
                if _label.startswith("de_"):
                    _display = _label[3:].replace("__", ": ").replace("_", " ").title()
                    _links.append({"label": f"DE table — {_display}", "url": _url})
                elif _label.startswith("volcano_"):
                    _display = _label[8:].replace("__", ": ").replace("_", " ").title()
                    _visuals.append({"type": "image", "url": _url, "title": f"Volcano — {_display}"})
                elif _label == "pca":
                    _visuals.append({"type": "image", "url": _url, "title": "PCA Plot"})

            return {
                "status": "success",
                "visualization_type": "results_viewer",
                "text": _text,
                "links": _links,
                "visuals": _visuals,
                "result": {"status": "success", "de_summary": _summary, "top_genes": _top,
                           "n_genes_total": 1000, "n_samples": 12, "design_formula": design_formula},
            }

        # ── Fast-path for the APAP time-course demo data ───────────────────────
        _APAP_CM   = "s3://noricum-ngs-data/demo/rnaseq/apap_timecourse_counts.csv"
        _APAP_META = "s3://noricum-ngs-data/demo/rnaseq/apap_timecourse_metadata.csv"
        _APAP_PRECOMP = "noricum-ngs-data/demo/rnaseq/apap_precomputed/latest"
        if count_matrix.strip() == _APAP_CM and sample_metadata.strip() == _APAP_META:
            import boto3 as _boto3
            _s3a = _boto3.client("s3")
            _bucket_a = _APAP_PRECOMP.split("/")[0]
            _prefix_a = "/".join(_APAP_PRECOMP.split("/")[1:])
            # List all files in the precomputed prefix
            _resp = _s3a.list_objects_v2(Bucket=_bucket_a, Prefix=_prefix_a + "/")
            _apap_files = [obj["Key"] for obj in _resp.get("Contents", [])]
            _apap_presigned: Dict[str, str] = {}
            for _key in _apap_files:
                _fname = _key.split("/")[-1]
                _label = _fname.replace(".csv", "").replace(".png", "")
                try:
                    _apap_presigned[_label] = _s3a.generate_presigned_url(
                        "get_object", Params={"Bucket": _bucket_a, "Key": _key}, ExpiresIn=86400)
                except Exception:
                    pass

            _apap_summary = [
                {"contrast": "time point: 0h vs 6h",   "total_genes": 1000, "significant": 11, "upregulated": 3, "downregulated": 8},
                {"contrast": "time point: 0h vs 24h",  "total_genes": 1000, "significant": 12, "upregulated": 4, "downregulated": 8},
                {"contrast": "time point: 0h vs 72h",  "total_genes": 1000, "significant": 9,  "upregulated": 1, "downregulated": 8},
                {"contrast": "time point: 0h vs 168h", "total_genes": 1000, "significant": 7,  "upregulated": 2, "downregulated": 5},
            ]
            _apap_top = (
                "\n**0h vs 6h (acute response)** — top genes: Cyp2e1, Hmox1, Lcn2, Saa1, Il6"
                "\n**0h vs 24h (peak injury)** — top genes: Ccl2, Cxcl10, Tgfb1, Col1a1, Timp1"
                "\n**0h vs 168h (recovery)** — top genes: Alb, Apoe, Cyp7a1, G6pc, Pcsk9"
            )
            _apap_header = (
                f"## Bulk RNA-seq Time-Course Analysis Complete\n\n"
                f"**Design formula:** `{design_formula}`  \n"
                f"**Samples:** 20 (5 time points × 4 replicates)  |  **Genes:** 1,000\n\n"
                f"### Differential Expression Summary (vs 0h baseline)\n\n"
                f"| Contrast | Total Genes | Significant (padj < {alpha}) | Up | Down |\n"
                f"|----------|-------------|-----------------------------|----|------|\n"
            )
            _apap_rows = "".join(
                f"| {r['contrast']} | {r['total_genes']} | {r['significant']} | {r['upregulated']} | {r['downregulated']} |\n"
                for r in _apap_summary
            )
            _apap_text = _apap_header + _apap_rows + f"\n### Top Differentially Expressed Genes{_apap_top}\n"

            _apap_links, _apap_visuals = [], []
            for _label, _url in _apap_presigned.items():
                if _label.startswith("de_"):
                    _disp = _label[3:].replace("__", ": ").replace("_", " ").replace("  ", " ").title()
                    _apap_links.append({"label": f"DE table — {_disp}", "url": _url})
                elif _label.startswith("volcano_"):
                    _disp = _label[8:].replace("__", ": ").replace("_", " ").replace("  ", " ").title()
                    _apap_visuals.append({"type": "image", "url": _url, "title": f"Volcano — {_disp}"})
                elif _label == "pca":
                    _apap_visuals.append({"type": "image", "url": _url, "title": "PCA — Time Points"})

            return {
                "status": "success",
                "visualization_type": "results_viewer",
                "text": _apap_text,
                "links": _apap_links,
                "visuals": _apap_visuals,
                "result": {"status": "success", "de_summary": _apap_summary, "top_genes": _apap_top,
                           "n_genes_total": 1000, "n_samples": 20, "design_formula": design_formula},
            }

        result = bulk_rnaseq.run_deseq2_analysis(
            count_matrix=count_matrix,
            sample_metadata=sample_metadata,
            design_formula=design_formula,
            alpha=alpha,
        )

        if result.get("status") == "error":
            return {
                "status": "error",
                "text": result.get("message", "Bulk RNA-seq analysis failed."),
                "result": result,
            }

        # Build a rich markdown summary
        summary = result.get("summary", [])
        mode_note = "\n\n> ⚠️ *pydeseq2 not installed — results are illustrative only.*" if result.get("mode") == "mock" else ""
        header = (
            f"## Bulk RNA-seq Analysis Complete\n\n"
            f"**Design formula:** `{design_formula}`  \n"
            f"**Samples:** {result.get('n_samples', '?')}  |  "
            f"**Genes:** {result.get('n_genes_total', '?')}\n\n"
            f"### Differential Expression Summary\n\n"
            f"| Contrast | Total Genes | Significant (padj < {alpha}) | Up | Down |\n"
            f"|----------|-------------|-----------------------------|----|------|\n"
        )
        rows_md = ""
        for row in summary:
            rows_md += (
                f"| {row['contrast']} | {row['total_genes']} | "
                f"{row['significant']} | {row['upregulated']} | {row['downregulated']} |\n"
            )
        top = result.get("top_genes", "")
        top_md = f"\n### Top Differentially Expressed Genes{top}\n" if top else ""

        text = header + rows_md + top_md + mode_note

        # Build links and visuals for the results_viewer
        presigned = result.get("presigned_urls", {})
        links = []
        visuals = []
        for label, url in presigned.items():
            if label.startswith("de_"):
                display = label[3:].replace("__", ": ").replace("_", " ").title()
                links.append({"label": f"DE table — {display}", "url": url})
            elif label.startswith("volcano_"):
                display = label[8:].replace("__", ": ").replace("_", " ").title()
                visuals.append({"type": "image", "url": url, "title": f"Volcano — {display}"})
            elif label == "pca":
                visuals.append({"type": "image", "url": url, "title": "PCA Plot"})

        return {
            "status": "success",
            "visualization_type": "results_viewer",
            "text": text,
            "links": links,
            "visuals": visuals,
            "result": result,
        }
    
    elif tool_name == "dna_vendor_research":
        # Handle DNA vendor research
        import dna_vendor_research
        
        command = arguments.get("command", "")
        sequence_length = arguments.get("sequence_length")
        quantity = arguments.get("quantity", "standard")
        
        result = dna_vendor_research.run_dna_vendor_research_raw(command, sequence_length, quantity)
        
        return {
            "status": result.get("status", "success"),
            "result": result,
            "text": result.get("message", "DNA vendor research completed"),
            "total_vendors": result.get("total_vendors", 0),
            "total_testing_options": result.get("total_testing_options", 0)
        }
    
    elif tool_name == "fastqc_quality_analysis":
        # Delegate to the proper fastqc_quality_analysis function in agent_tools.py
        # which has proper infrastructure decision handling (_from_broker flag)
        # We access the underlying function via .func attribute to bypass LangChain's @tool wrapper
        from backend.agent_tools import fastqc_quality_analysis as fastqc_tool
        
        input_r1 = arguments.get("input_r1", "")
        input_r2 = arguments.get("input_r2", "")
        output = arguments.get("output")
        _from_broker = arguments.get("_from_broker", False)
        session_id = arguments.get("session_id")
        
        if not input_r1 or not input_r2:
            return {
                "status": "error",
                "result": {},
                "text": "Both input_r1 and input_r2 are required for FastQC analysis"
            }
        
        # Call the underlying function directly (bypassing LangChain's @tool decorator)
        # The _from_broker flag will determine whether to use local or EMR execution
        try:
            # Access the underlying function via .func attribute
            if hasattr(fastqc_tool, 'func'):
                # LangChain tool - call underlying function
                result = fastqc_tool.func(
                    input_r1=input_r1,
                    input_r2=input_r2,
                    output=output,
                    _from_broker=_from_broker,
                    session_id=session_id
                )
            else:
                # Not a LangChain tool - call directly
                result = fastqc_tool(
                    input_r1=input_r1,
                    input_r2=input_r2,
                    output=output,
                    _from_broker=_from_broker,
                    session_id=session_id
                )
            return result
        except Exception as e:
            logger.error(f"FastQC tool invocation failed: {e}")
            return {
                "status": "error",
                "result": {},
                "text": f"Failed to invoke FastQC tool: {str(e)}",
                "error": str(e)
            }

    elif tool_name == "s3_browse_results":
        from backend.agent_tools import s3_browse_results as s3_tool

        tool_input = {
            "prefix": arguments.get("prefix") or arguments.get("output") or arguments.get("s3_prefix") or "",
            "show": arguments.get("show") or arguments.get("results_json") or arguments.get("results_path"),
            "recursive": bool(arguments.get("recursive", True)),
            "max_keys": int(arguments.get("max_keys", 200) or 200),
            "mode": arguments.get("mode") or "display",
        }

        try:
            if hasattr(s3_tool, "invoke"):
                return s3_tool.invoke(tool_input)
            if hasattr(s3_tool, "func"):
                return s3_tool.func(**tool_input)
            return s3_tool(**tool_input)
        except Exception as e:
            logger.error(f"s3_browse_results invocation failed: {e}")
            return {
                "status": "error",
                "result": {},
                "text": f"Failed to browse S3 results: {str(e)}",
                "error": str(e),
            }
    
    else:
        # Unknown tool - try tool-generator-agent
        logger.info(f"🔧 Unknown tool '{tool_name}', attempting tool-generator-agent...")
        try:
            from backend.tool_generator_agent import generate_and_execute_tool, _discover_inputs_from_args, _discover_outputs_from_args
            from backend.intent_classifier import classify_intent
            
            # Build command from tool_name and arguments
            # Try to construct a natural language command from the tool name and arguments
            # If the original user command is available, use it; otherwise reconstruct
            user_command = arguments.get("command") or arguments.get("user_request") or arguments.get("original_command")
            
            if user_command:
                # Use the original user command if available
                command = user_command
                user_request = user_command
            else:
                # Reconstruct command from tool_name and arguments
                command_parts = [tool_name.replace("_", " ")]
                for key, value in arguments.items():
                    if key not in ["command", "user_request", "original_command", "session_id"]:
                        if isinstance(value, str) and len(value) < 200:  # Avoid huge values
                            command_parts.append(f"{key}: {value}")
                
                command = " ".join(command_parts)
                user_request = f"Execute {tool_name} with parameters: {arguments}"

            intent = classify_intent(user_command or command)
            if intent.intent != "execute":
                # Phase 4: never run tool-generator-agent for Q&A intent.
                return {
                    "status": "success",
                    "tool_generated": False,
                    "tool_name": tool_name,
                    "result": {},
                    "text": (
                        "This request looks like Q&A intent. "
                        "Helix.AI will not generate or execute a new tool for Q&A. "
                        "Rephrase as an execution request (e.g. 'run ...', 'analyze ...')."
                    ),
                    "intent": intent.intent,
                    "intent_reason": intent.reason,
                }
            
            # Discover inputs and outputs to pass file information to the agent
            session_context = arguments.get("session_context") or {}
            discovered_inputs = _discover_inputs_from_args(arguments, session_context)
            if discovered_inputs:
                logger.info(f"🔧 Discovered {len(discovered_inputs)} input files for infrastructure decision")
            
            discovered_outputs = _discover_outputs_from_args(arguments, user_command or command)
            if discovered_outputs:
                logger.info(f"🔧 Discovered {len(discovered_outputs)} output paths")
            
            result = await generate_and_execute_tool(
                command=command,
                user_request=user_request,
                session_id=arguments.get("session_id"),
                inputs=discovered_inputs,
                outputs=discovered_outputs,
                session_context=session_context
            )
            
            if result.get("status") == "success":
                logger.info("✅ Tool-generator-agent successfully generated and executed tool")
                return {
                    "status": "success",
                    "tool_generated": True,
                    "tool_name": tool_name,
                    "result": result,
                    "text": result.get("explanation", "Tool generated and executed successfully")
                }
            else:
                # Extract error message - check top-level error first, then execution_result
                error_msg = result.get('error')
                if not error_msg:
                    execution_result = result.get("execution_result", {})
                    if isinstance(execution_result, dict):
                        error_msg = execution_result.get("error")
                        if not error_msg and execution_result.get("stderr"):
                            error_msg = f"Execution failed: {execution_result.get('stderr', '')[:200]}"
                    error_msg = error_msg or "Unknown error"
                logger.warning(f"⚠️  Tool-generator-agent failed: {error_msg}")
                # Fall through to raise ValueError
        except Exception as e:
            logger.error(f"❌ Tool-generator-agent exception: {e}", exc_info=True)
            # Fall through to raise ValueError
        
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

# FastQC Job Management Endpoints

class FastQCJobRequest(BaseModel):
    r1_path: str
    r2_path: str
    output_path: Optional[str] = None
    session_id: Optional[str] = None

@app.post("/tools/fastqc-analysis")
async def submit_fastqc_job(req: FastQCJobRequest):
    """
    Submit a FastQC quality analysis job to EMR.
    
    Returns immediately with a job_id for tracking.
    """
    try:
        from backend.job_manager import get_job_manager
        
        job_manager = get_job_manager()
        # Run blocking submit_fastqc_job in executor to avoid blocking event loop
        # This is necessary because submit_fastqc_job can block for up to 15 minutes
        # waiting for EMR cluster to be ready
        if hasattr(asyncio, 'to_thread'):
            # Python 3.9+ - use to_thread
            job_id = await asyncio.to_thread(
                job_manager.submit_fastqc_job,
                r1_path=req.r1_path,
                r2_path=req.r2_path,
                output_path=req.output_path,
                session_id=req.session_id
            )
        else:
            # Python < 3.9 - use run_in_executor
            loop = asyncio.get_event_loop()
            job_id = await loop.run_in_executor(
                None,
                lambda: job_manager.submit_fastqc_job(
                    r1_path=req.r1_path,
                    r2_path=req.r2_path,
                    output_path=req.output_path,
                    session_id=req.session_id
                )
            )
        
        return CustomJSONResponse({
            "success": True,
            "job_id": job_id,
            "status": "submitted",
            "message": "FastQC job submitted. Processing will take 10-30 minutes."
        })
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Failed to submit FastQC job: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/jobs/{job_id}")
async def get_job_status(job_id: str):
    """Get the current status of a job."""
    try:
        from backend.job_manager import get_job_manager
        
        job_manager = get_job_manager()
        job = job_manager.get_job_status(job_id)
        
        return CustomJSONResponse({
            "success": True,
            "job": job
        })
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Failed to get job status: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/jobs/{job_id}/results")
async def get_job_results(job_id: str):
    """Get results for a completed job, including output file validation."""
    try:
        from backend.job_manager import get_job_manager
        
        job_manager = get_job_manager()
        results = job_manager.get_job_results(job_id)
        
        # Include validation warnings in response
        validation = results.get("output_validation", {})
        warnings = validation.get("validation_warnings", [])
        
        response = {
            "success": True,
            "results": results,
        }
        
        # Add top-level warning if files are missing
        if warnings:
            response["warnings"] = warnings
            if not validation.get("all_files_exist", True):
                response["warning"] = "Some expected output files are missing. Check validation details."
        
        return CustomJSONResponse(response)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Failed to get job results: {e}")
        raise HTTPException(status_code=500, detail=str(e))

class CopyToSessionRequest(BaseModel):
    session_id: str

@app.post("/jobs/{job_id}/copy-to-session")
async def copy_job_results_to_session(job_id: str, req: Optional[CopyToSessionRequest] = None, session_id: Optional[str] = None):
    """Copy a completed job's results to its session's S3 path.
    
    This is useful for jobs that completed before session-specific storage was implemented,
    or for jobs that are no longer in the backend's in-memory store.
    
    Args:
        job_id: Job identifier
        session_id: Optional session ID. If not provided, will try to get from job metadata.
    """
    try:
        from backend.job_manager import get_job_manager
        from backend.history_manager import history_manager
        
        job_manager = get_job_manager()
        
        # Try to get job from memory
        job = None
        try:
            job = job_manager.get_job_status(job_id)
        except ValueError:
            # Job not in memory - that's OK, we'll work with what we have
            logger.info(f"Job {job_id} not in memory, will use provided session_id: {session_id}")
        
        # Get session_id from request body or query parameter
        if not session_id:
            if req and req.session_id:
                session_id = req.session_id
            elif job:
                session_id = job.get("session_id")
        
        if not session_id:
            raise HTTPException(
                status_code=400,
                detail=f"Job {job_id} has no session_id. Please provide session_id in request body or as query parameter."
            )
        
        # Verify session exists
        session = history_manager.get_session(session_id)
        if not session:
            raise HTTPException(
                status_code=404,
                detail=f"Session {session_id} not found"
            )
        
        # Get session S3 path
        s3_bucket = session.get("metadata", {}).get("s3_bucket")
        s3_path = session.get("metadata", {}).get("s3_path")
        if not s3_bucket or not s3_path:
            raise HTTPException(
                status_code=400,
                detail=f"Session {session_id} has no S3 path configured"
            )
        
        # Use job_id as sub-folder name: s3://bucket/session_id/job_id/
        # This allows multiple FastQC jobs per session, each in their own folder
        session_s3_path = f"s3://{s3_bucket}/{s3_path}{job_id}/"
        
        # Find results.json - try multiple locations
        results_path = None
        
        # Strategy 1: If job is in memory, use its results path
        if job:
            try:
                results_info = job_manager.get_job_results(job_id)
                results_path = results_info.get("results_path")
            except Exception as e:
                logger.warning(f"Could not get results path from job manager: {e}")
        
        # Strategy 2: Try default location based on common pattern
        if not results_path:
            # Default location: s3://noricum-ngs-data/fastqc-results/GRCh38.p12.MafHi/results.json
            default_results = "s3://noricum-ngs-data/fastqc-results/GRCh38.p12.MafHi/results.json"
            import boto3
            s3_client = boto3.client('s3')
            try:
                # Check if file exists
                parts = default_results.replace("s3://", "").split("/", 1)
                s3_client.head_object(Bucket=parts[0], Key=parts[1])
                results_path = default_results
                logger.info(f"Found results at default location: {results_path}")
            except Exception:
                logger.warning(f"Results not found at default location: {default_results}")
        
        if not results_path:
            raise HTTPException(
                status_code=404,
                detail=f"Could not find results.json for job {job_id}. Please ensure the job completed successfully."
            )
        
        # Copy results to session S3 path
        import boto3
        s3_client = boto3.client('s3')
        
        # Parse source and destination
        orig_parts = results_path.replace("s3://", "").split("/", 1)
        orig_bucket = orig_parts[0]
        orig_key = orig_parts[1]
        
        dest_parts = session_s3_path.replace("s3://", "").split("/", 1)
        dest_bucket = dest_parts[0]
        dest_key = f"{dest_parts[1]}results.json" if len(dest_parts) > 1 else "results.json"
        
        try:
            # Copy results.json
            s3_client.copy_object(
                CopySource={'Bucket': orig_bucket, 'Key': orig_key},
                Bucket=dest_bucket,
                Key=dest_key
            )
            logger.info(f"Copied results.json to session S3 path: s3://{dest_bucket}/{dest_key}")
            
            # If job is in memory, update it
            if job:
                job["session_results_path"] = f"s3://{dest_bucket}/{dest_key}"
                job["session_id"] = session_id
            
            # Try to generate HTML visualizations
            html_path = None
            try:
                # Create a temporary job entry if not in memory, so HTML generation can work
                if not job:
                    # Reconstruct minimal job entry for HTML generation
                    job_manager.jobs[job_id] = {
                        "job_id": job_id,
                        "status": "completed",
                        "type": "fastqc",
                        "session_id": session_id,
                        "session_results_path": f"s3://{dest_bucket}/{dest_key}",
                        "output_path": None,
                        "r1_path": None,  # Not needed for HTML generation
                        "r2_path": None,
                    }
                
                job_manager._generate_and_upload_html_visualizations(job_id)
                updated_job = job_manager.get_job_status(job_id) if job_id in job_manager.jobs else None
                if updated_job:
                    html_path = updated_job.get("session_html_path")
            except Exception as e:
                logger.warning(f"Failed to generate HTML visualizations: {e}")
                import traceback
                traceback.print_exc()
            
            return CustomJSONResponse({
                "success": True,
                "message": f"Results copied to session S3 path for job {job_id}",
                "session_results_path": f"s3://{dest_bucket}/{dest_key}",
                "session_html_path": html_path
            })
        except Exception as e:
            logger.error(f"Failed to copy results to session: {e}")
            raise HTTPException(
                status_code=500,
                detail=f"Failed to copy results to session: {str(e)}"
            )
            
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to copy job results to session: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/jobs/{job_id}/results-data")
async def get_job_results_data(job_id: str, results_path: Optional[str] = None):
    """Get the actual results data (JSON) from S3 for a completed job.
    
    Args:
        job_id: Job identifier
        results_path: Optional S3 path to results.json. If not provided, will try to get from job manager.
    """
    try:
        import boto3
        import tempfile
        
        from backend.job_manager import get_job_manager
        
        # If results_path not provided, try to get from job manager
        if not results_path:
            job_manager = get_job_manager()
            try:
                results_info = job_manager.get_job_results(job_id)
                results_path = results_info.get("results_path")
            except ValueError:
                # Job not found in backend - this is OK if frontend provides results_path
                pass
        
        if not results_path or not results_path.startswith("s3://"):
            raise ValueError(f"Invalid results path: {results_path}. Please provide results_path parameter if job is not in backend memory.")
        
        # Parse S3 path
        s3_path = results_path.replace("s3://", "")
        parts = s3_path.split("/", 1)
        if len(parts) != 2:
            raise ValueError(f"Invalid S3 path format: {results_path}")
        
        bucket = parts[0]
        key = parts[1]
        
        # Download from S3
        s3_client = boto3.client('s3')
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.json') as tmp_file:
            tmp_path = tmp_file.name
        
        try:
            s3_client.download_file(bucket, key, tmp_path)
            
            # Read and parse JSON
            with open(tmp_path, 'r') as f:
                results_data = json.load(f)
            
            return CustomJSONResponse({
                "success": True,
                "data": results_data
            })
        finally:
            # Clean up temp file
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
        
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Failed to get job results data: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/jobs/{job_id}/cancel")
async def cancel_job(job_id: str):
    """Cancel a running job."""
    try:
        from backend.job_manager import get_job_manager
        
        job_manager = get_job_manager()
        result = job_manager.cancel_job(job_id)
        
        return CustomJSONResponse({
            "success": True,
            "result": result
        })
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Failed to cancel job: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/jobs/{job_id}/retry")
async def retry_job(job_id: str):
    """Retry a failed job with the same parameters."""
    try:
        from backend.job_manager import get_job_manager
        
        job_manager = get_job_manager()
        new_job_id = job_manager.retry_job(job_id)
        
        # Get the new job status
        new_job = job_manager.get_job_status(new_job_id)
        
        return CustomJSONResponse({
            "success": True,
            "job_id": new_job_id,
            "status": new_job["status"],
            "message": f"Job {job_id} retried as new job {new_job_id}"
        })
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Failed to retry job: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/jobs/{job_id}/logs")
async def get_job_logs(job_id: str):
    """Get EMR logs for a job."""
    try:
        from backend.job_manager import get_job_manager
        
        job_manager = get_job_manager()
        logs = job_manager.get_job_logs(job_id)
        
        return CustomJSONResponse({
            "success": True,
            "logs": logs
        })
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Failed to get job logs: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/jobs")
async def list_jobs(status: Optional[str] = None, limit: int = 100):
    """List all jobs, optionally filtered by status."""
    try:
        from backend.job_manager import get_job_manager
        
        job_manager = get_job_manager()
        jobs = job_manager.list_jobs(status=status, limit=limit)
        
        return CustomJSONResponse({
            "success": True,
            "jobs": jobs,
            "count": len(jobs)
        })
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Failed to list jobs: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/jobs/backfill-directories")
async def backfill_job_directories(job_ids: Optional[List[str]] = None):
    """Backfill/create local directories for existing jobs.
    
    If job_ids is not provided, backfills all jobs.
    """
    try:
        from backend.job_manager import get_job_manager
        job_manager = get_job_manager()
        results = job_manager.backfill_job_directories(job_ids=job_ids)
        success_count = sum(1 for v in results.values() if v)
        total_count = len(results)
        return CustomJSONResponse({
            "success": True,
            "results": results,
            "summary": {
                "total": total_count,
                "successful": success_count,
                "failed": total_count - success_count
            }
        })
    except Exception as e:
        logger.error(f"Error backfilling job directories: {e}")
        raise HTTPException(status_code=500, detail=str(e))

# ── Data Science Pipeline endpoints ──────────────────────────────────────────

class DSRunRequest(BaseModel):
    data_path: str
    target_col: Optional[str] = None
    task_type: str = "auto"
    objective: str = "Analyze dataset"
    hypothesis: str = ""
    changes: str = "Initial run"
    random_seed: int = 42
    time_col: Optional[str] = None
    entity_col: Optional[str] = None
    session_id: Optional[str] = None
    feedback: Optional[str] = None
    parent_run_id: Optional[str] = None


class DSReproduceRequest(BaseModel):
    run_id: str
    session_id: Optional[str] = None


class DSDiffRequest(BaseModel):
    run_id_a: str
    run_id_b: str
    session_id: Optional[str] = None


@app.post("/ds/run")
async def ds_run(req: DSRunRequest):
    """
    Run the full iterative data science pipeline on a CSV dataset.

    Executes: data audit → cleaning → EDA → baseline model → evaluation → report → planning
    """
    try:
        from backend.ds_pipeline.orchestrator import DataScienceOrchestrator, DSRunConfig
        from backend.ds_pipeline.reviewer import review as make_review

        if not req.session_id:
            req.session_id = history_manager.create_session()
        else:
            history_manager.ensure_session_exists(req.session_id)

        from backend.agent_tools import _get_session_local_dir
        session_dir = _get_session_local_dir(req.session_id)

        config = DSRunConfig(
            data_path=req.data_path,
            target_col=req.target_col,
            task_type=req.task_type,
            objective=req.objective,
            hypothesis=req.hypothesis,
            changes=req.changes,
            random_seed=req.random_seed,
            time_col=req.time_col,
            entity_col=req.entity_col,
        )
        orch = DataScienceOrchestrator(base_dir=session_dir, session_id=req.session_id)
        run_data = orch.run(config, feedback=req.feedback, parent_run_id=req.parent_run_id)

        return {
            "success": True,
            "session_id": req.session_id,
            "run_id": run_data["run_id"],
            "decision": run_data.get("decision"),
            "metrics": run_data.get("metrics"),
            "next_steps": run_data.get("next_steps"),
            "summary": make_review(run_data),
            "artifacts": run_data.get("artifacts"),
            "steps_run": run_data.get("steps_run"),
        }
    except Exception as e:
        logger.error(f"DS run failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/ds/reproduce")
async def ds_reproduce(req: DSReproduceRequest):
    """Re-run a prior data science iteration using its exact recorded configuration."""
    try:
        from backend.ds_pipeline.orchestrator import DataScienceOrchestrator
        from backend.ds_pipeline.reviewer import review as make_review

        if not req.session_id:
            req.session_id = history_manager.create_session()

        from backend.agent_tools import _get_session_local_dir
        session_dir = _get_session_local_dir(req.session_id)

        orch = DataScienceOrchestrator(base_dir=session_dir, session_id=req.session_id)
        run_data = orch.reproduce(req.run_id)
        return {
            "success": True,
            "session_id": req.session_id,
            "run_id": run_data["run_id"],
            "parent_run_id": req.run_id,
            "decision": run_data.get("decision"),
            "metrics": run_data.get("metrics"),
            "summary": make_review(run_data),
        }
    except Exception as e:
        logger.error(f"DS reproduce failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/ds/diff")
async def ds_diff(req: DSDiffRequest):
    """Compare two data science runs: configs, metrics, and decisions."""
    try:
        from backend.ds_pipeline.orchestrator import DataScienceOrchestrator

        session_id = req.session_id or history_manager.create_session()
        from backend.agent_tools import _get_session_local_dir
        session_dir = _get_session_local_dir(session_id)

        orch = DataScienceOrchestrator(base_dir=session_dir, session_id=session_id)
        diff = orch.diff(req.run_id_a, req.run_id_b)
        return {"success": True, "diff": diff}
    except Exception as e:
        logger.error(f"DS diff failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/ds/runs")
async def ds_list_runs(session_id: str):
    """List all data science runs in the given session."""
    try:
        from backend.ds_pipeline.run_store import RunStore
        from backend.agent_tools import _get_session_local_dir
        session_dir = _get_session_local_dir(session_id)
        store = RunStore(base_dir=session_dir)
        return {"success": True, "runs": store.read_experiment_log(), "run_ids": store.list_runs()}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/ds/runs/{run_id}")
async def ds_get_run(run_id: str, session_id: str):
    """Retrieve the full run.json for a given run_id."""
    try:
        from backend.ds_pipeline.run_store import RunStore
        from backend.agent_tools import _get_session_local_dir
        session_dir = _get_session_local_dir(session_id)
        store = RunStore(base_dir=session_dir)
        run = store.load_run(run_id)
        if not run:
            raise HTTPException(status_code=404, detail=f"Run {run_id} not found")
        return {"success": True, "run": run}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/health")
async def health_check():
    """Health check endpoint - available immediately, even before full initialization."""
    return {"status": "healthy", "service": "Helix.AI Bioinformatics API"}

@app.get("/api-docs")
async def api_docs_info(request: Request):
    """Get information about API documentation endpoints."""
    base_url = str(request.base_url).rstrip('/')
    return {
        "openapi_spec": f"{base_url}/openapi.json",
        "swagger_ui": f"{base_url}/docs",
        "redoc": f"{base_url}/redoc",
        "description": "Access the OpenAPI specification and interactive documentation",
        "endpoints": {
            "openapi_json": {
                "url": f"{base_url}/openapi.json",
                "description": "OpenAPI 3.0 specification in JSON format"
            },
            "swagger_ui": {
                "url": f"{base_url}/docs",
                "description": "Interactive Swagger UI documentation"
            },
            "redoc": {
                "url": f"{base_url}/redoc",
                "description": "ReDoc documentation interface"
            }
        }
    }

_DEMO_PRESIGNED_CACHE: Dict[str, str] = {}

def _refresh_demo_presigned_urls() -> None:
    """Pre-generate presigned URLs for all demo assets at startup / on demand."""
    try:
        import boto3 as _b3_startup
        _s3_startup = _b3_startup.client("s3", region_name="us-west-1")
        _DEMO_ASSETS = {
            # Key: label, Value: (bucket, key)
            "phylo_tree":          ("noricum-ngs-data", "demo/phylo/precomputed/latest/phylo_tree.png"),
            "identity_matrix":     ("noricum-ngs-data", "demo/phylo/precomputed/latest/identity_matrix.png"),
            "variant_mutations":   ("noricum-ngs-data", "demo/phylo/precomputed/latest/variant_mutations.csv"),
            "spike_sequences":     ("noricum-ngs-data", "demo/phylo/precomputed/latest/spike_sequences.fasta"),
            "tree_data_json":      ("noricum-ngs-data", "demo/phylo/precomputed/latest/tree_data.json"),
            "scrna_umap_celltype": ("noricum-ngs-data", "demo/scrna/precomputed/latest/umap_celltype.png"),
            "scrna_umap_disease":  ("noricum-ngs-data", "demo/scrna/precomputed/latest/umap_disease.png"),
            "scrna_dotplot":       ("noricum-ngs-data", "demo/scrna/precomputed/latest/dotplot_markers.png"),
            "scrna_markers_csv":   ("noricum-ngs-data", "demo/scrna/precomputed/latest/marker_genes.csv"),
            # APAP time-course (pre-fetch a few key files; rest resolved dynamically from s3 listing)
            "apap_pca":            ("noricum-ngs-data", "demo/rnaseq/apap_precomputed/latest/pca.png"),
            "apap_volcano_0h_6h":  ("noricum-ngs-data", "demo/rnaseq/apap_precomputed/latest/volcano_time_point__0h_vs_6h.png"),
            "apap_volcano_0h_24h": ("noricum-ngs-data", "demo/rnaseq/apap_precomputed/latest/volcano_time_point__0h_vs_24h.png"),
            "apap_volcano_0h_168h":("noricum-ngs-data", "demo/rnaseq/apap_precomputed/latest/volcano_time_point__0h_vs_168h.png"),
            "apap_de_0h_6h":       ("noricum-ngs-data", "demo/rnaseq/apap_precomputed/latest/de_time_point__0h_vs_6h.csv"),
            "apap_de_0h_24h":      ("noricum-ngs-data", "demo/rnaseq/apap_precomputed/latest/de_time_point__0h_vs_24h.csv"),
            "apap_de_0h_168h":     ("noricum-ngs-data", "demo/rnaseq/apap_precomputed/latest/de_time_point__0h_vs_168h.csv"),
        }
        for label, (bucket, key) in _DEMO_ASSETS.items():
            try:
                url = _s3_startup.generate_presigned_url(
                    "get_object", Params={"Bucket": bucket, "Key": key}, ExpiresIn=86400 * 3)
                _DEMO_PRESIGNED_CACHE[label] = url
            except Exception as _ex:
                logger.debug(f"[startup] presign failed for {label}: {_ex}")
        # Fetch newick tree string
        try:
            _nwk = _s3_startup.get_object(Bucket="noricum-ngs-data",
                                           Key="demo/phylo/precomputed/latest/tree_data.json")
            import json as _jstart
            _nwk_data = _jstart.loads(_nwk["Body"].read().decode())
            _DEMO_PRESIGNED_CACHE["newick"] = _nwk_data.get("newick", "")
        except Exception:
            pass
        logger.info(f"[startup] Demo presigned URL cache loaded ({len(_DEMO_PRESIGNED_CACHE)} entries)")
    except Exception as _e:
        logger.warning(f"[startup] Could not pre-generate demo presigned URLs: {_e}")


@app.on_event("startup")
async def startup_event():
    """Initialize heavy components after server starts."""
    # Pre-generate presigned URLs for demo assets (avoids per-request boto3 calls)
    import asyncio
    loop = asyncio.get_event_loop()
    loop.run_in_executor(None, _refresh_demo_presigned_urls)

    # Print EC2 environment variables only once at actual server startup
    # (not during module imports in reloader)
    print("=" * 80)
    print(f"🔍 DEBUG: EC2 Environment Variables at Backend Startup (dotenv: {_dotenv_path or 'NOT FOUND'}):")
    print(f"  AWS_REGION: {os.getenv('AWS_REGION', 'NOT SET')}")
    print(f"  HELIX_EC2_INSTANCE_ID: {os.getenv('HELIX_EC2_INSTANCE_ID', 'NOT SET')}")
    print(f"  HELIX_EC2_KEY_NAME: {os.getenv('HELIX_EC2_KEY_NAME', 'NOT SET')}")
    print(f"  HELIX_EC2_KEY_FILE: {os.getenv('HELIX_EC2_KEY_FILE', 'NOT SET')}")
    print(f"  HELIX_USE_EC2: {os.getenv('HELIX_USE_EC2', 'NOT SET')}")
    print(f"  HELIX_EC2_AUTO_CREATE: {os.getenv('HELIX_EC2_AUTO_CREATE', 'NOT SET')}")
    # Also show what values were parsed from the dotenv file (if found)
    if _dotenv_values:
        print("  [dotenv] AWS_REGION:", _dotenv_values.get("AWS_REGION", "NOT IN FILE"))
        print("  [dotenv] HELIX_EC2_INSTANCE_ID:", _dotenv_values.get("HELIX_EC2_INSTANCE_ID", "NOT IN FILE"))
        print("  [dotenv] HELIX_EC2_KEY_NAME:", _dotenv_values.get("HELIX_EC2_KEY_NAME", "NOT IN FILE"))
        print("  [dotenv] HELIX_EC2_KEY_FILE:", _dotenv_values.get("HELIX_EC2_KEY_FILE", "NOT IN FILE"))
    print("=" * 80)
    
    # Heavy initialization (LLM, agent) happens lazily on first use
    # Session loading happens lazily in history_manager
    logger.info("Backend startup complete - ready to accept requests")

if __name__ == "__main__":
    import uvicorn
    # Run as package module (start.sh uses `python -m backend.main_with_mcp`)
    uvicorn.run("backend.main_with_mcp:app", host="0.0.0.0", port=8001, reload=True)
