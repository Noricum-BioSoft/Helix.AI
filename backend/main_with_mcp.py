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
import os
import shutil
import numpy as np
import logging
from datetime import datetime, timezone

logger = logging.getLogger(__name__)

# Add the current directory to Python path for imports
sys.path.append(str(Path(__file__).parent))

from history_manager import history_manager
from tool_schemas import list_tool_schemas
from context_builder import _truncate_sequence

from execution_broker import ExecutionBroker, ExecutionRequest


def _get_bioagent_handle_command():
    """
    Lazy import of BioAgent to avoid importing heavy LLM dependencies at module import time.
    This keeps lightweight endpoints (like /health, /mcp/tools) working in sandbox/CI.
    """
    from agent import handle_command  # local import by design
    return handle_command

def _truncate_sequences_in_dict(obj: Any, max_length: int = 100) -> Any:
    """
    Recursively truncate sequences in dictionaries, lists, and nested structures.
    This prevents large sequences from being included in JSON responses or LLM context.
    """
    if isinstance(obj, str):
        # Truncate if it's longer than max_length
        if len(obj) > max_length:
            # Check if it looks like a sequence (mostly ATCGUN characters, not text with spaces/punctuation)
            # Sample first 200 chars to check
            sample = obj[:200].upper()
            if len(sample) > 50:
                # If >80% of characters are ATCGUN and no spaces, it's likely a sequence
                seq_chars = sum(1 for c in sample if c in 'ATCGUN')
                has_spaces = ' ' in sample
                if seq_chars / len(sample) > 0.8 and not has_spaces:
                    return _truncate_sequence(obj, max_length)
        return obj
    elif isinstance(obj, dict):
        truncated = {}
        for key, value in obj.items():
            # Always truncate "sequence" fields if they're strings
            if key == "sequence" and isinstance(value, str) and len(value) > max_length:
                truncated[key] = _truncate_sequence(value, max_length)
            # Truncate full_sequence fields too - they shouldn't be in JSON responses
            # Full sequences should be stored separately or accessed via download, not in API responses
            elif key == "full_sequence" and isinstance(value, str) and len(value) > max_length:
                truncated[key] = _truncate_sequence(value, max_length)
            else:
                truncated[key] = _truncate_sequences_in_dict(value, max_length)
        return truncated
    elif isinstance(obj, list):
        return [_truncate_sequences_in_dict(item, max_length) for item in obj]
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
        "visuals": [],
        "links": [],
    }
    
    # Determine visualization type for frontend rendering
    visualization_type = _determine_visualization_type(tool, truncated_result, prompt)

    # Set top-level error field for frontend compatibility (frontend checks output.error)
    error_message = None
    if errors and len(errors) > 0:
        error_message = errors[0].get("message", "An error occurred")
    
    return {
        "version": "1.0",
        "success": success,
        "session_id": session_id,
        "prompt": prompt,
        "tool": tool,
        "status": status,
        "text": text,
        "data": data,
        "visualization_type": visualization_type,  # Hint for frontend on how to render
        "logs": [],
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
    """Execute a general command using the existing agent with session tracking."""
    try:
        # Validate prompt length
        _validate_prompt_length(req.command)

        # Create session if not provided
        if not req.session_id:
            req.session_id = history_manager.create_session()
        else:
            # Ensure session exists; auto-create if missing
            if not history_manager.get_session(req.session_id):
                # Create local directory for the session
                session_dir = history_manager.storage_dir / req.session_id
                session_dir.mkdir(exist_ok=True)
                logger.info(f"Auto-created session directory: {session_dir}")
                
                # Create S3 path for the session
                s3_path = history_manager._create_s3_session_path(req.session_id)
                history_manager.sessions[req.session_id] = {
                    "session_id": req.session_id,
                    "user_id": None,
                    "created_at": datetime.now().isoformat(),
                    "updated_at": datetime.now().isoformat(),
                    "history": [],
                    "results": {},
                    "metadata": {
                        "s3_path": s3_path,
                        "s3_bucket": history_manager.s3_bucket_name if s3_path else None,
                        "local_path": str(session_dir)
                    }
                }
                history_manager._save_session(req.session_id)
        
        # Rate limit per identity (session or IP)
        identity = _get_request_identity(request, req.session_id)
        _check_and_increment_daily_counter(identity)
        
        # Get session context from history manager
        session_context = {}
        if hasattr(history_manager, 'sessions') and req.session_id in history_manager.sessions:
            session_context = history_manager.sessions[req.session_id]

        # Phase 3: detect multi-step workflows and execute as a Plan IR (sync/async broker handles routing)
        def _looks_like_workflow(cmd: str) -> bool:
            c = (cmd or "").lower()
            return any(tok in c for tok in [" and then ", " then ", "->", "→", "\n", ";"]) and len(c) > 20

        # Primary path: let BioAgent (with agent.md prompt) plan/execute
        use_agent = os.getenv("HELIX_MOCK_MODE") != "1"

        try:
            if not use_agent:
                raise RuntimeError("Skipping BioAgent in mock mode (HELIX_MOCK_MODE=1)")

            import time
            agent_start_time = time.time()
            handle_command = _get_bioagent_handle_command()
            agent_result = await handle_command(req.command, session_id=req.session_id, session_context=session_context)
            
            agent_done_time = time.time()
            agent_duration = agent_done_time - agent_start_time
            print(f"✅ [PERF] Agent tool mapping completed in {agent_duration:.2f}s")
            
            # Check if agent returned a tool mapping (not full execution)
            if isinstance(agent_result, dict) and agent_result.get("status") == "tool_mapped":
                # Agent only did tool mapping - now execute via router
                tool_name = agent_result.get("tool_name")
                parameters = agent_result.get("parameters", {})
                
                print(f"🔧 Agent mapped tool '{tool_name}', executing via router...")
                
                # Add session_id for tools that need it
                if tool_name == "fastqc_quality_analysis" and "session_id" not in parameters:
                    parameters["session_id"] = req.session_id
                
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
                    result
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
                )
                
                history_done_time = time.time()
                print(f"✅ [PERF] History entry added, took {(history_done_time - agent_done_time)*1000:.2f}ms")
                
                response_start_time = time.time()
                standard_response = build_standard_response(
                    prompt=req.command,
                    tool="agent",
                    result=agent_result,
                    session_id=req.session_id,
                    mcp_route="/execute",
                    success=True if not isinstance(agent_result, dict) else agent_result.get("status", "success") != "error"
                )
                response_done_time = time.time()
                total_duration = response_done_time - agent_start_time
                print(f"✅ [PERF] Backend response ready in {total_duration:.2f}s (total)")
                print(f"✅ [PERF] Response size: {len(str(standard_response))} chars")
                return CustomJSONResponse(standard_response)
        except Exception as agent_err:
            # Fallback: use NLP router and MCP tools if the agent path fails
            print(f"⚠️  Agent path failed, falling back to router/tool. Error: {agent_err}")

            from command_router import CommandRouter
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
                result
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
        handle_command = _get_bioagent_handle_command()
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
                # Execute the tool via router
                tool_start_time = time.time()
                tool_result = await call_mcp_tool(tool_name, parameters)
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

async def call_mcp_tool(tool_name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Call an MCP tool and return the result."""
    # Add tools directory to path
    tools_path = str((Path(__file__).resolve().parent.parent / "tools").resolve())
    sys.path.insert(0, tools_path)
    
    if tool_name == "toolbox_inventory":
        from tool_inventory import build_toolbox_inventory, format_toolbox_inventory_markdown
        inv = build_toolbox_inventory()
        return {
            "status": "success",
            "text": format_toolbox_inventory_markdown(inv),
            "result": inv,
        }

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
        # Check if inputs are S3 paths - if so, use tool-generator-agent
        forward_reads = arguments.get("forward_reads", "")
        reverse_reads = arguments.get("reverse_reads", "")
        
        # Check if either input is an S3 path or file path (not FASTQ content)
        # FASTQ content typically starts with "@" (read header) or contains newlines
        is_s3_path = (forward_reads.startswith("s3://") or reverse_reads.startswith("s3://"))
        
        # Also check for file paths (local paths starting with /)
        is_file_path = ((forward_reads.startswith("/") and not forward_reads.startswith("@")) or
                       (reverse_reads.startswith("/") and not reverse_reads.startswith("@")))
        
        # Check if it looks like FASTQ content (starts with @ or contains newlines with @)
        is_fastq_content = (forward_reads.startswith("@") or reverse_reads.startswith("@") or
                           "\n@" in forward_reads or "\n@" in reverse_reads)
        
        # Also check if the original command contains S3 paths
        original_command = arguments.get("command") or arguments.get("user_request") or arguments.get("original_command", "")
        if original_command and "s3://" in original_command:
            is_s3_path = True
        
        # Use tool-generator-agent if S3 paths or file paths (but not FASTQ content)
        use_tool_generator = is_s3_path or (is_file_path and not is_fastq_content)
        
        if use_tool_generator:
            # Use tool-generator-agent for S3 paths
            logger.info(f"🔧 read_merging tool detected S3 paths, routing to tool-generator-agent...")
            from tool_generator_agent import generate_and_execute_tool
            
            # Use original command if available, otherwise reconstruct
            if original_command:
                command = original_command
                user_request = original_command
            else:
                command = f"merge forward R1 and reverse R2 reads: R1: {forward_reads} R2: {reverse_reads}"
                user_request = command
            
            result = await generate_and_execute_tool(
                command=command,
                user_request=user_request,
                session_id=arguments.get("session_id")
            )
            
            if result.get("status") == "success":
                logger.info("✅ Tool-generator-agent successfully generated and executed read merging tool")
                return {
                    "status": "success",
                    "tool_generated": True,
                    "tool_name": tool_name,
                    "result": result,
                    "text": result.get("explanation", "Read merging completed successfully")
                }
            else:
                logger.warning(f"⚠️  Tool-generator-agent failed: {result.get('error', 'Unknown error')}")
                # Fall through to try existing tool as fallback
                pass
        
        # Use existing tool for local FASTQ content
        import read_merging
        return read_merging.run_read_merging_raw(
            forward_reads,
            reverse_reads,
            arguments.get("min_overlap", 12),
        )
    
    elif tool_name == "quality_assessment":
        import quality_assessment
        sequences = arguments.get("sequences", "")
        print(f"🔧 [DEBUG] Quality assessment tool called with {len(sequences)} characters of sequences")
        return quality_assessment.run_quality_assessment_raw(sequences)

    elif tool_name == "handle_natural_command":
        # Use the BioAgent path (system prompt from agent.md) for natural commands
        bioagent_handle_command = _get_bioagent_handle_command()
        command = arguments.get("command", "")
        session_id = arguments.get("session_id", "")
        return await bioagent_handle_command(command, session_id=session_id)
    
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
        
        # Mock mode or missing Rscript: return a stub success so UI renders
        if os.getenv("HELIX_MOCK_MODE") or shutil.which("Rscript") is None:
            return {
                "status": "success",
                "result": {
                    "status": "success",
                    "summary": {
                        "cells": 500,
                        "genes": 2000,
                        "clusters": 8,
                        "top_markers": ["GeneA", "GeneB", "GeneC"]
                    },
                    "plots": {
                        "umap": "mock_umap.png",
                        "qc": "mock_qc.png"
                    },
                    "message": "Mock single-cell analysis (HELIX_MOCK_MODE or Rscript unavailable)."
                },
                "text": "Single-cell analysis completed (mock)"
            }
        
        data_file = arguments.get("data_file")
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
        # Handle NCBI sequence fetching
        import ncbi_tools
        
        accession = arguments.get("accession")
        database = arguments.get("database", "nucleotide")
        
        result = ncbi_tools.fetch_sequence_from_ncbi(accession, database)
        
        # Truncate sequence for LLM consumption to avoid sending very long sequences
        if result.get("status") == "success" and "sequence" in result:
            full_sequence = result.get("sequence", "")
            sequence_length = len(full_sequence)
            # For very large sequences, truncate aggressively to prevent:
            # 1. LLM timeout (agent path)
            # 2. Huge JSON responses (500MB+)
            # Full sequence should be stored in session/history, not in API response
            if sequence_length > 50:
                truncated_sequence = full_sequence[:50] + f"... (truncated, full length: {sequence_length:,} bp)"
                # Create a modified result with truncated sequence
                # Don't include full_sequence in response - it's too large for JSON
                truncated_result = result.copy()
                truncated_result["sequence"] = truncated_sequence
                # Remove full_sequence from response to prevent 500MB+ payloads
                truncated_result.pop("full_sequence", None)
                result = truncated_result
        
        length = result.get("length")
        description = result.get("description", "")
        length_text = f" ({length} bp)" if length else ""
        desc_text = f": {description}" if description else ""
        return {
            "status": result.get("status", "success"),
            "result": result,
            "text": f"Fetched sequence {accession}{length_text} from {database}{desc_text}" if result.get("status") == "success" else f"Error: {result.get('error', 'Unknown error')}"
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
        # Handle bulk RNA-seq analysis
        import bulk_rnaseq
        
        # In mock mode or when Rscript missing, supply defaults so the stub returns success
        if os.getenv("HELIX_MOCK_MODE") or shutil.which("Rscript") is None:
            if not arguments.get("count_matrix"):
                arguments["count_matrix"] = str(Path(__file__).resolve().parent.parent / "tests" / "data" / "counts.csv")
            if not arguments.get("sample_metadata"):
                arguments["sample_metadata"] = str(Path(__file__).resolve().parent.parent / "tests" / "data" / "metadata.csv")
        
        count_matrix = arguments.get("count_matrix")
        sample_metadata = arguments.get("sample_metadata")
        design_formula = arguments.get("design_formula", "~condition")
        alpha = arguments.get("alpha", 0.05)
        if not count_matrix or not sample_metadata:
            return {
                "status": "error",
                "result": {},
                "text": "Count matrix and sample metadata paths are required for bulk RNA-seq analysis"
            }
        
        result = bulk_rnaseq.run_deseq2_analysis(
            count_matrix=count_matrix,
            sample_metadata=sample_metadata,
            design_formula=design_formula,
            alpha=alpha
        )
        
        return {
            "status": result.get("status", "success"),
            "result": result,
            "text": "Bulk RNA-seq analysis completed" if result.get("status") == "success" else f"Error: {result.get('message', 'Unknown error')}"
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
        # Handle FastQC quality analysis
        from job_manager import get_job_manager
        
        input_r1 = arguments.get("input_r1", "")
        input_r2 = arguments.get("input_r2", "")
        output = arguments.get("output")
        session_id = arguments.get("session_id")  # Get session_id from arguments if provided
        
        if not input_r1 or not input_r2:
            return {
                "status": "error",
                "result": {},
                "text": "Both input_r1 and input_r2 are required for FastQC analysis"
            }
        
        job_manager = get_job_manager()
        try:
            # Run blocking submit_fastqc_job in executor to avoid blocking event loop
            # This is necessary because submit_fastqc_job can block for up to 15 minutes
            # waiting for EMR cluster to be ready
            if hasattr(asyncio, 'to_thread'):
                # Python 3.9+ - use to_thread
                job_id = await asyncio.to_thread(
                    job_manager.submit_fastqc_job,
                    r1_path=input_r1,
                    r2_path=input_r2,
                    output_path=output,
                    session_id=session_id
                )
            else:
                # Python < 3.9 - use run_in_executor
                loop = asyncio.get_event_loop()
                job_id = await loop.run_in_executor(
                    None,
                    lambda: job_manager.submit_fastqc_job(
                        r1_path=input_r1,
                        r2_path=input_r2,
                        output_path=output,
                        session_id=session_id
                    )
                )
            
            return {
                "type": "job",
                "status": "submitted",
                "job_id": job_id,
                "message": "FastQC job submitted. Processing will take 10-30 minutes.",
                "input_r1": input_r1,
                "input_r2": input_r2,
                "output": output
            }
        except Exception as e:
            return {
                "status": "error",
                "result": {},
                "text": f"Failed to submit FastQC job: {str(e)}",
                "error": str(e)
            }
    
    else:
        # Unknown tool - try tool-generator-agent
        logger.info(f"🔧 Unknown tool '{tool_name}', attempting tool-generator-agent...")
        try:
            from tool_generator_agent import generate_and_execute_tool
            
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
            
            result = await generate_and_execute_tool(
                command=command,
                user_request=user_request,
                session_id=arguments.get("session_id")
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
                logger.warning(f"⚠️  Tool-generator-agent failed: {result.get('error', 'Unknown error')}")
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
        from job_manager import get_job_manager
        
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
        from job_manager import get_job_manager
        
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
    """Get results for a completed job."""
    try:
        from job_manager import get_job_manager
        
        job_manager = get_job_manager()
        results = job_manager.get_job_results(job_id)
        
        return CustomJSONResponse({
            "success": True,
            "results": results
        })
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
        from job_manager import get_job_manager
        from history_manager import history_manager
        
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
        
        from job_manager import get_job_manager
        
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
        from job_manager import get_job_manager
        
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
        from job_manager import get_job_manager
        
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
        from job_manager import get_job_manager
        
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
        from job_manager import get_job_manager
        
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
        from job_manager import get_job_manager
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

@app.on_event("startup")
async def startup_event():
    """Initialize heavy components after server starts."""
    # This runs after the server is ready to accept connections
    # Heavy initialization (LLM, agent) happens at import time, but at least
    # the health endpoint is available quickly
    pass

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main_with_mcp:app", host="0.0.0.0", port=8001, reload=True) 
