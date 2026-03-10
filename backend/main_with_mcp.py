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
        _execution_broker = ExecutionBroker(tool_executor=dispatch_tool)
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


# ──────────────────────────────────────────────────────────────────────────────
# Dispatch helpers — shared by every branch of execute()
# ──────────────────────────────────────────────────────────────────────────────

def _is_success(result: Any) -> bool:
    """Return True unless the result carries an explicit error/failed status."""
    if not isinstance(result, dict):
        return True
    if "success" in result:
        return bool(result["success"])
    return result.get("status", "success") not in ("error", "failed", "workflow_failed")


def _extract_metadata(result: Any, tool_args: Optional[dict] = None) -> dict:
    """
    Extract run-linkage fields from a tool result dict.

    Checks both the top-level result dict and the optional nested ``result``
    sub-dict so we handle all tool shapes (direct returns and broker-wrapped).
    """
    if not isinstance(result, dict):
        return {
            "tool_args": tool_args,
            "inputs": None,
            "outputs": [],
            "produced_artifacts": None,
            "mcp_route": "/execute",
        }
    inner: dict = (result.get("result") or {}) if isinstance(result.get("result"), dict) else {}
    return {
        "tool_args": tool_args,
        "inputs": result.get("inputs"),
        "outputs": [],
        "produced_artifacts": (
            result.get("produced_artifacts")
            or inner.get("produced_artifacts")
            or result.get("artifacts")
            or inner.get("artifacts")
        ),
        "job_id":       result.get("job_id")       or inner.get("job_id"),
        "run_id":       result.get("run_id")       or inner.get("run_id"),
        "parent_run_id": result.get("parent_run_id") or inner.get("parent_run_id"),
        "mcp_route": "/execute",
    }


async def _dispatch_result(
    req: "CommandRequest",
    tool: str,
    result: Any,
    *,
    tool_args: Optional[dict] = None,
    record_history: bool = True,
    execution_path: Optional[str] = None,
) -> "CustomJSONResponse":
    """
    Single exit point shared by every dispatch branch.

    Optionally records a history entry, then builds and returns the standard
    JSON response.  Pass ``record_history=False`` for paths (e.g. S3 browse)
    that intentionally skip the run ledger.
    """
    if record_history:
        history_manager.add_history_entry(
            req.session_id,
            req.command,
            tool,
            result,
            metadata=_extract_metadata(result, tool_args=tool_args),
        )
    if execution_path:
        logger.info("Routing path=%s tool=%s session=%s", execution_path, tool, req.session_id)
    std = build_standard_response(
        prompt=req.command,
        tool=tool,
        result=result,
        session_id=req.session_id,
        mcp_route="/execute",
        success=_is_success(result),
        execution_path=execution_path,
    )
    return CustomJSONResponse(std)


def _build_pipeline_execution_storage_result(result: Any, session_id: str) -> Any:
    """
    Build a trimmed plan result for session storage: pipeline execution info only,
    no system prompts, full_response, or huge code blocks. Adds pipeline_entry_point,
    storage paths, and per-step inputs/outputs summary.
    """
    if not isinstance(result, dict):
        return result
    inner = result.get("result") or result
    if inner.get("type") != "plan_result" or not isinstance(inner.get("steps"), list):
        return result

    paths = history_manager.get_session_storage_paths(session_id)
    steps_summary: List[Dict[str, Any]] = []
    trimmed_steps: List[Dict[str, Any]] = []

    for idx, step in enumerate(inner["steps"], 1):
        step_id = step.get("id") or f"step{idx}"
        tool_name = step.get("tool_name") or ""
        args = step.get("arguments") or {}
        step_result = step.get("result") or {}

        # Inputs summary (no previous_plan_steps blob)
        input_keys = [k for k in args.keys() if k != "previous_plan_steps"]
        inputs_summary = {k: _summarize_for_storage(args[k]) for k in input_keys}
        if "previous_plan_steps" in args:
            inputs_summary["previous_plan_steps"] = "(from previous step(s))"

        # Outputs summary
        outputs_summary: Dict[str, Any] = {}
        if isinstance(step_result, dict):
            status = step_result.get("status") or ("success" if step_result.get("alignment") or step_result.get("stdout") else "unknown")
            if step_result.get("alignment"):
                al = step_result["alignment"]
                outputs_summary["alignment"] = f"{len(al)} sequences" if isinstance(al, list) else str(al)[:80]
            if step_result.get("execution_result"):
                er = step_result["execution_result"]
                if isinstance(er, dict) and er.get("stdout"):
                    outputs_summary["stdout"] = (er["stdout"] or "").strip()[:200]
            if step_result.get("text"):
                outputs_summary["text"] = (step_result["text"] or "").strip()[:200]
        else:
            status = "unknown"

        steps_summary.append({
            "step_index": idx,
            "id": step_id,
            "tool_name": tool_name,
            "inputs": inputs_summary,
            "outputs": outputs_summary,
            "status": status,
            "step_dir": None,  # reserved: per-step working dir when we add it
        })

        # Trim step for storage: drop full_response, long code_preview, previous_plan_steps
        trimmed_step = {
            "id": step_id,
            "tool_name": tool_name,
            "arguments": {k: v for k, v in args.items() if k != "previous_plan_steps"},
        }
        if "previous_plan_steps" in args:
            trimmed_step["arguments"]["_previous_plan_steps"] = "(included at execution time)"
        res = dict(step_result) if isinstance(step_result, dict) else {}
        for key in ("full_response", "code_preview"):
            res.pop(key, None)
        if res.get("result") and isinstance(res["result"], dict):
            res["result"] = {k: v for k, v in res["result"].items() if k not in ("full_response", "code_preview")}
        trimmed_step["result"] = res
        trimmed_steps.append(trimmed_step)

    pipeline_execution = {
        "entry_point": "backend.execution_broker.ExecutionBroker._execute_plan_sync (invoked from POST /execute when _looks_like_workflow)",
        "storage": {
            "storage_root": paths["storage_root"],
            "session_file": paths["session_file"],
            "session_dir": paths["session_dir"],
        },
        "steps": steps_summary,
    }

    trimmed_inner = {
        "status": inner.get("status"),
        "type": "plan_result",
        "plan_version": inner.get("plan_version"),
        "steps": trimmed_steps,
        "result": trimmed_steps[-1]["result"] if trimmed_steps else {},
        "pipeline_execution": pipeline_execution,
    }
    return {
        **result,
        "result": trimmed_inner,
    }


def _summarize_for_storage(val: Any, max_len: int = 120) -> Any:
    """Shorten values for pipeline execution storage."""
    if val is None:
        return None
    if isinstance(val, str):
        return val[:max_len] + ("..." if len(val) > max_len else "")
    if isinstance(val, (int, float, bool)):
        return val
    if isinstance(val, dict):
        return "(object)"
    if isinstance(val, list):
        return f"({len(val)} items)"
    return str(val)[:max_len]


def _apply_session_context_side_effects(
    session_id: str, tool_name: str, result: Any
) -> None:
    """
    Persist mutation / alignment outputs into the in-memory session dict so
    downstream tool calls (e.g. variant_selection, sequence_alignment) can
    read them back via session_context. For __plan__, apply side effects from
    each step so session handling matches single-tool flow.
    """
    if not isinstance(result, dict):
        return
    if not (hasattr(history_manager, "sessions") and session_id in history_manager.sessions):
        return

    # Plan result: apply side effects from each step so session has aligned_sequences etc.
    if tool_name == "__plan__":
        inner = result.get("result") or result
        if inner.get("type") == "plan_result" and isinstance(inner.get("steps"), list):
            for step in inner["steps"]:
                step_tool = step.get("tool_name")
                step_result = step.get("result") or {}
                if step_tool and isinstance(step_result, dict):
                    _apply_session_context_side_effects(session_id, step_tool, step_result)
        return

    if tool_name == "mutate_sequence":
        variants = (
            (result.get("statistics") or {}).get("variants")
            or result.get("variants")
            or (result.get("output") or {}).get("variants")
        )
        if variants:
            history_manager.sessions[session_id]["mutated_sequences"] = variants
            history_manager.sessions[session_id]["mutation_results"] = variants
            print(
                f"🔧 [DEBUG] After mutation via /execute - Session {session_id} "
                f"mutated_sequences: {len(variants)} items"
            )

    if tool_name == "sequence_alignment":
        aligned_seqs = result.get("alignment") if isinstance(result.get("alignment"), list) else None
        if aligned_seqs is None and isinstance(result.get("output"), list):
            aligned_seqs = result["output"]
        if aligned_seqs:
            fasta_lines: list = []
            for seq in aligned_seqs:
                if isinstance(seq, dict):
                    fasta_lines.append(f">{seq.get('name', 'sequence')}")
                    fasta_lines.append(seq.get("sequence", ""))
            history_manager.sessions[session_id]["aligned_sequences"] = "\n".join(fasta_lines)
            print(
                f"🔧 [DEBUG] After alignment via /execute - Stored "
                f"{len(aligned_seqs)} aligned sequences in session context"
            )


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
    execution_path: Optional[str] = None,
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
        if truncated_result.get("tree_newick"):
            tree_data["tree_newick"] = truncated_result["tree_newick"]
        if truncated_result.get("ete_visualization"):
            tree_data["ete_visualization"] = truncated_result["ete_visualization"]
        if truncated_result.get("clustering_result"):
            tree_data["clustering_result"] = truncated_result["clustering_result"]
        if truncated_result.get("clustered_visualization"):
            tree_data["clustered_visualization"] = truncated_result["clustered_visualization"]
    
    # Promote execute_ready flag so the frontend can show the "Execute Pipeline" button
    execute_ready = bool(
        isinstance(truncated_result, dict) and truncated_result.get("execute_ready")
    )

    # Unwrap ExecutionBroker envelope (type == "execution_result") so that
    # patch_and_rerun / iteration tools' inner fields (run_id, script_path, links …)
    # are reachable at the top level of the result we inspect below.
    _inner_result: Dict = truncated_result if isinstance(truncated_result, dict) else {}
    if _inner_result.get("type") == "execution_result" and isinstance(_inner_result.get("result"), dict):
        _inner_result = _inner_result["result"]

    # Inject analysis.py download link when the result carries a script_path.
    # Works for both patch_and_rerun responses (explicit script_path) and
    # BioOrchestrator runs (script saved under sessions/<sid>/runs/<run_id>/).
    _script_path_str: Optional[str] = _inner_result.get("script_path") or None
    _run_id_hint: Optional[str] = (
        _inner_result.get("run_id")
        or (isinstance(truncated_result, dict) and truncated_result.get("run_id"))
        or None
    )
    if not _script_path_str and _run_id_hint and session_id:
        # BioOrchestrator saves scripts under sessions/<sid>/runs/<run_id>/
        _candidate = (
            Path(__file__).parent.parent
            / "sessions" / session_id / "runs" / _run_id_hint / "analysis.py"
        )
        if _candidate.exists():
            _script_path_str = str(_candidate)

    if _script_path_str:
        _script_url = f"/download/script?path={_script_path_str}"
        _bundle_url = (
            f"/download/bundle?session_id={session_id}"
            + (f"&run_id={_run_id_hint}" if _run_id_hint else "")
        )
        # Start from data.links, then merge any links from the inner result
        _existing_links: List[Dict] = list(data.get("links") or [])
        for _lnk in (_inner_result.get("links") or []):
            if isinstance(_lnk, dict) and _lnk not in _existing_links:
                _existing_links.append(_lnk)

        # Ensure analysis.py and bundle.zip are always present
        if not any(
            isinstance(lnk, dict) and "analysis.py" in (lnk.get("label") or "")
            for lnk in _existing_links
        ):
            _existing_links.append({"label": "analysis.py", "url": _script_url})
        if not any(
            isinstance(lnk, dict) and "bundle" in (lnk.get("label") or "").lower()
            for lnk in _existing_links
        ):
            _existing_links.append({"label": "bundle.zip", "url": _bundle_url})
        data["links"] = _existing_links

    response = {
        "version": "1.0",
        "success": success,
        "session_id": session_id,
        # Promote run_id to top-level so the frontend can read it without unwrapping envelopes
        "run_id": _run_id_hint or None,
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
            "tool_route": f"dispatch_tool:{tool}"
        },
        "raw_result": truncated_result,  # Truncated to prevent large sequences in JSON
        "timestamp": now
    }

    if os.getenv("HELIX_DEBUG_ROUTING", "0").lower() in ("1", "true", "yes"):
        response["execution_path"] = execution_path or "unknown"
        if isinstance(truncated_result, dict):
            if "intent" in truncated_result:
                response["intent"] = truncated_result.get("intent")
            if "intent_reason" in truncated_result:
                response["intent_reason"] = truncated_result.get("intent_reason")

    # Add tree data to top level for frontend
    response.update(tree_data)
    
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


@app.get("/download/bundle")
async def download_bundle(session_id: str, run_id: Optional[str] = None):
    """Build and stream a reproducibility ZIP for a session run.

    The bundle contains:
    - ``README.md``           human-readable walkthrough (prompts, I/O, steps)
    - ``run_manifest.json``   machine-readable metadata for the full run chain
    - ``analysis.py``         re-runnable script for the target (latest) run
    - ``plots/``              PNG plots (≤ 5 MB each)
    - ``tables/``             JSON / CSV tables (≤ 5 MB each)
    - ``iteration_history/``  scripts + params for every ancestor run
    - ``large_files.txt``     references for files that exceeded the size limit

    Parameters
    ----------
    session_id : str
        Active session ID.
    run_id : str, optional
        Specific run to package.  Defaults to the most recent scriptable run.
    """
    try:
        from backend.bundle_generator import build_bundle
        buf, filename = build_bundle(session_id=session_id, run_id=run_id)
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc))
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Bundle generation failed: {exc}")

    from fastapi.responses import StreamingResponse
    return StreamingResponse(
        buf,
        media_type="application/zip",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


@app.get("/download/script")
async def download_script(path: str):
    """Serve any analysis.py (or other text artifact) that lives under the sessions/ tree.

    The ``path`` query parameter must be an absolute filesystem path that resolves
    inside the project's ``sessions/`` directory — any attempt to escape is rejected
    with 403.

    Example::

        GET /download/script?path=/abs/path/to/sessions/xxx/runs/yyy/analysis.py
    """
    sessions_root = (Path(__file__).parent.parent / "sessions").resolve()
    try:
        resolved = Path(path).resolve()
    except Exception:
        raise HTTPException(status_code=400, detail="Invalid path")

    # Security: must be inside sessions/
    try:
        resolved.relative_to(sessions_root)
    except ValueError:
        raise HTTPException(status_code=403, detail="Path is outside the sessions directory")

    if not resolved.exists() or not resolved.is_file():
        raise HTTPException(status_code=404, detail="File not found")

    return FileResponse(
        path=str(resolved),
        filename=resolved.name,
        media_type="text/plain",
    )


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
                            execution_path="fast_path_s3_browse",
                        )
                        return CustomJSONResponse(standard_response)
        except Exception as e:
            logger.warning(f"S3 browse fast-path skipped due to error: {e}")

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
                execution_path="agent_execute_plan",
            )
            return CustomJSONResponse(standard)

        # Phase 2c: deterministic router allowlist fast path (agent-first strategy).
        # Keep this path narrow and deterministic. Everything else should flow to
        # the agent path below for semantic tool selection.
        _phase2c_allowlist = {"toolbox_inventory", "session_run_io_summary"}
        try:
            import time as _t
            from backend.command_router import CommandRouter as _PreRouter
            _pre_router = _PreRouter()
            _pre_tool, _pre_params = _pre_router.route_command(req.command, session_context)
            logger.info(f"[Phase2c] route_command → tool='{_pre_tool}'")
            if _pre_tool in _phase2c_allowlist:
                _pre_start = _t.time()
                if _pre_params is None:
                    _pre_params = {}
                _pre_params.setdefault("session_id", req.session_id)
                # Direct tool call — bypasses broker & infra-decision-agent
                _pre_result = await dispatch_tool(_pre_tool, _pre_params)
                logger.info(f"[Phase2c] '{_pre_tool}' done in {_t.time()-_pre_start:.2f}s")
                return await _dispatch_result(
                    req, _pre_tool, _pre_result, tool_args=_pre_params, execution_path="phase2c_router"
                )
        except Exception as _pre_err:
            logger.warning(f"[Phase2c] fast path raised exception ({_pre_err}), falling through to agent")

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
            logger.info("Agent tool mapping completed in %.2fs", agent_duration)

            # Check if agent returned a tool mapping (not full execution)
            if isinstance(agent_result, dict) and agent_result.get("status") == "tool_mapped":
                # Agent only did tool mapping - now execute via router
                tool_name = agent_result.get("tool_name")
                parameters = agent_result.get("parameters", {})
                
                logger.debug("Agent mapped tool '%s', executing via router...", tool_name)
                
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
                                logger.warning("Fixed malformed S3 URI for %s: %s → %s", param_name, uri, fixed_uri)
                                parameters[param_name] = fixed_uri
                    
                    logger.debug("FastQC parameters: input_r1=%s, input_r2=%s", parameters.get('input_r1'), parameters.get('input_r2'))
                
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
                logger.info("Tool '%s' completed in %.2fs", tool_name, tool_done_time - tool_start_time)
                return await _dispatch_result(
                    req, tool_name, result, tool_args=parameters, execution_path="agent_tool_mapped"
                )
            else:
                # Agent completed full execution (or returned something else)
                return await _dispatch_result(
                    req,
                    "agent",
                    agent_result,
                    tool_args={"execute_plan": bool(req.execute_plan)},
                    execution_path="agent",
                )
        except Exception as agent_err:
            # Fallback: use NLP router and MCP tools if the agent path fails
            logger.warning("Agent path failed, falling back to router/tool: %s", agent_err)

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
                        return await _dispatch_result(
                            req, _tool_name, result, tool_args=_params or {}, execution_path="fallback_router"
                        )
                except Exception:
                    pass

                standard_response = build_standard_response(
                    prompt=req.command,
                    tool="handle_natural_command",
                    result={
                        "status": "success",
                        "text": (
                            "This looks like a question. Enable the agent (HELIX_MOCK_MODE=0) "
                            "or use the /chat endpoint for Q&A. If you intended execution, "
                            "rephrase as: 'Run <analysis> on <inputs>'."
                        ),
                        "intent": intent.intent,
                        "intent_reason": intent.reason,
                    },
                    session_id=req.session_id,
                    mcp_route="/execute",
                    success=True,
                    execution_path="qa_safe_fallback",
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
                _apply_session_context_side_effects(req.session_id, "__plan__", result)
                storage_result = _build_pipeline_execution_storage_result(result, req.session_id)
                return await _dispatch_result(
                    req, "__plan__", storage_result,
                    tool_args={"plan": plan, "session_id": req.session_id},
                    execution_path="fallback_router_plan",
                )

            tool_name, parameters = command_router.route_command(req.command, session_context)
            logger.debug("Routed '%s' → tool='%s'", req.command[:60], tool_name)

            if tool_name == "variant_selection" and "session_id" not in parameters:
                parameters["session_id"] = req.session_id
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
            logger.info("Tool '%s' completed in %.2fs", tool_name, time.time() - tool_start_time)

            _apply_session_context_side_effects(req.session_id, tool_name, result)
            return await _dispatch_result(
                req, tool_name, result, tool_args=parameters, execution_path="fallback_router"
            )
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
        logger.info("Agent tool mapping completed in %.2fs", agent_duration)

        # Check if agent returned a tool mapping (not full execution)
        if isinstance(result, dict) and result.get("status") == "tool_mapped":
            # Agent only did tool mapping - now execute via router
            tool_name = result.get("tool_name")
            parameters = result.get("parameters", {})

            logger.debug("Agent mapped tool '%s', executing via router...", tool_name)
            
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
                logger.info("Tool '%s' completed in %.2fs", tool_name, tool_duration)
                
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
                    success=True if not isinstance(tool_result, dict) else tool_result.get("status", "success") != "error",
                    execution_path="agent_endpoint_tool",
                )
                response_done_time = time.time()
                total_duration = response_done_time - agent_start_time
                logger.info("Response ready in %.2fs", total_duration)
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
                    success=False,
                    execution_path="agent_endpoint_tool_error",
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
                success=True if not isinstance(result, dict) else result.get("status", "success") != "error",
                execution_path="agent_endpoint",
            )
            response_done_time = time.time()
            total_duration = response_done_time - agent_start_time
            logger.info("Response ready in %.2fs", total_duration)

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
        
        result = await dispatch_tool("sequence_alignment", {
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
        
        result = await dispatch_tool("mutate_sequence", {
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
                
                logger.debug("Session %s: stored %d mutated sequences", req.session_id, len(variants) if isinstance(variants, list) else 0)
        
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
        
        result = await dispatch_tool("analyze_sequence_data", {
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
        logger.debug("select_variants called for session %s", req.session_id)
        
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
        result = await dispatch_tool("visualize_alignment", {
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
        result = await dispatch_tool("read_trimming", {
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
        result = await dispatch_tool("read_merging", {
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


def _save_analysis_script(
    tool_name: str,
    arguments: Dict[str, Any],
    result: Dict[str, Any],
    session_id: str = "",
) -> Optional[str]:
    """Generate and persist an analysis.py for *tool_name* + *arguments*.

    Called after the tool result is available so we can embed the run_id that
    the BioOrchestrator assigned.  Returns the script path string on success,
    or None if generation is not supported for this tool.
    """
    try:
        from backend import code_generator as _cg
        from backend.history_manager import history_manager

        # Derive the run_id from the result (BioOrchestrator embeds it)
        run_id = (
            result.get("run_id")
            or result.get("result", {}).get("run_id", "")
            or ""
        )
        if not run_id:
            # Mint a simple timestamp-based id when orchestrator run_id is absent
            import time
            run_id = f"run-{int(time.time() * 1000)}"

        if not session_id:
            session_id = "default"

        # Store scripts under sessions/<session_id>/runs/<run_id>/
        sessions_root = Path(__file__).parent.parent / "sessions"
        run_dir = sessions_root / session_id / "runs" / run_id
        run_dir.mkdir(parents=True, exist_ok=True)

        script_text = _cg.generate(
            tool_name,
            arguments,
            run_id=run_id,
            session_id=session_id,
            run_dir=str(run_dir),
        )
        if script_text is None:
            return None

        script_path = run_dir / "analysis.py"
        script_path.write_text(script_text)

        # Register the script as an artifact on the matching run record
        try:
            runs = history_manager.list_runs(session_id)
            for r in reversed(runs):
                if isinstance(r, dict) and r.get("run_id") == run_id:
                    arts = r.setdefault("produced_artifacts", [])
                    # Avoid duplicates
                    if not any(
                        a.get("type") == "script" and a.get("uri") == str(script_path)
                        for a in arts if isinstance(a, dict)
                    ):
                        arts.append({
                            "type":  "script",
                            "uri":   str(script_path),
                            "title": "analysis.py",
                        })
                    break
        except Exception:
            pass  # non-fatal

        return str(script_path)
    except Exception as exc:
        import logging
        logging.getLogger(__name__).debug("Script generation failed: %s", exc)
        return None


# Tools for which we generate and save an analysis.py after each execution
_SCRIPTABLE_TOOLS = frozenset({
    "bulk_rnaseq_analysis",
    "single_cell_analysis",
    "phylogenetic_tree",
    "fastqc_quality_analysis",
})


async def dispatch_tool(tool_name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Dispatch to internal Helix tools and return the result."""
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
        "patch_and_rerun",
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
        import boto3
        import gzip
        
        # Check if we have separate forward and reverse reads
        forward_reads = arguments.get("forward_reads", "")
        reverse_reads = arguments.get("reverse_reads", "")
        reads = arguments.get("reads", "")
        
        adapter = arguments.get("adapter")
        quality_threshold = arguments.get("quality_threshold", 20)

        logger.debug("read_trimming: adapter=%s q=%s fw=%s rv=%s reads=%s",
                     adapter, quality_threshold, bool(forward_reads), bool(reverse_reads), bool(reads))

        def _load_fastq_content(value: str) -> str:
            if not value:
                return ""
            if value.startswith("s3://"):
                s3 = boto3.client("s3")
                bucket, key = value.replace("s3://", "").split("/", 1)
                obj = s3.get_object(Bucket=bucket, Key=key)
                body = obj["Body"].read()
                if key.endswith(".gz"):
                    return gzip.decompress(body).decode("utf-8", errors="replace")
                return body.decode("utf-8", errors="replace")
            # Treat filesystem path-like values as file inputs
            if (
                value.startswith("/")
                or value.endswith(".fastq")
                or value.endswith(".fastq.gz")
                or value.endswith(".fq")
                or value.endswith(".fq.gz")
            ):
                p = Path(value)
                if not p.exists():
                    raise FileNotFoundError(f"FASTQ file not found: {value}")
                if value.endswith(".gz"):
                    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
                        return fh.read()
                return p.read_text(encoding="utf-8", errors="replace")
            # Otherwise assume the string already contains FASTQ content
            return value
        
        # If we have separate forward/reverse reads, process them separately
        if forward_reads and reverse_reads:
            try:
                forward_content = _load_fastq_content(forward_reads)
                reverse_content = _load_fastq_content(reverse_reads)
                forward_result = read_trimming.run_read_trimming_raw(
                    forward_content,
                    adapter,
                    quality_threshold,
                )
                reverse_result = read_trimming.run_read_trimming_raw(
                    reverse_content,
                    adapter,
                    quality_threshold,
                )
            except Exception as e:
                return {
                    "status": "error",
                    "error": str(e),
                    "text": f"Failed to trim reads: {e}",
                }
            
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
            try:
                content = _load_fastq_content(reads)
                return read_trimming.run_read_trimming_raw(
                    content,
                    adapter,
                    quality_threshold,
                )
            except Exception as e:
                return {
                    "status": "error",
                    "error": str(e),
                    "text": f"Failed to trim reads: {e}",
                }
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
        logger.debug("quality_assessment called with %d chars", len(sequences))
        return quality_assessment.run_quality_assessment_raw(sequences)

    elif tool_name == "quality_report":
        # Deterministic quality-report summary from upstream pipeline metrics.
        import datetime as _dt
        raw_reads    = arguments.get("raw_reads")
        post_trim    = arguments.get("post_trim")
        merged_reads = arguments.get("merged_reads")
        if raw_reads is None or post_trim is None or merged_reads is None:
            return {
                "status": "error",
                "text": (
                    "quality_report requires `raw_reads`, `post_trim`, and `merged_reads` "
                    "from upstream pipeline steps."
                ),
            }
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

    elif tool_name == "bio_rerun":
        # Re-run the last bioinformatics tool with optional parameter overrides.
        from backend.bio_pipeline import BioOrchestrator as _BioOrch
        _orch_rerun = _BioOrch(tool_executor=dispatch_tool, history_manager=history_manager)
        _session_id_rr = arguments.get("session_id", "")
        _changes_rr = arguments.get("changes", {})
        _target_run = arguments.get("target_run", "latest")

        if _target_run == "latest":
            # Look up the most recent run for this session from the artifacts dir
            from pathlib import Path as _Path
            _arts_root = _Path(__file__).parent.parent / "artifacts"
            _run_dirs = sorted(
                _arts_root.glob("*/run.json"),
                key=lambda p: p.stat().st_mtime if p.exists() else 0,
                reverse=True,
            )
            if not _run_dirs:
                return {
                    "status": "error",
                    "text": "No prior bioinformatics run found. Run an analysis first.",
                }
            _target_run = _run_dirs[0].parent.name

        result = await _orch_rerun.rerun(_target_run, _changes_rr, session_id=_session_id_rr)
        return {
            "status": result.get("status", "success"),
            "visualization_type": "results_viewer",
            "text": result.get("text", result.get("summary_text", "Re-run complete.")),
            "visuals": result.get("visuals", []),
            "links": result.get("links", []),
            "run_id": result.get("run_id"),
            "parent_run_id": result.get("parent_run_id"),
            "delta": result.get("delta", {}),
            "result": result,
        }

    elif tool_name == "bio_diff_runs":
        # Return a structured comparison of two bioinformatics runs.
        from backend.bio_pipeline import BioOrchestrator as _BioOrch
        from pathlib import Path as _Path
        _orch_diff = _BioOrch(tool_executor=dispatch_tool)
        _run_a = arguments.get("run_id_a", "latest")
        _run_b = arguments.get("run_id_b", "prior")

        # Resolve "latest" / "prior" pseudo-IDs
        _arts_root_diff = _Path(__file__).parent.parent / "artifacts"
        _run_dirs_diff = sorted(
            _arts_root_diff.glob("*/run.json"),
            key=lambda p: p.stat().st_mtime if p.exists() else 0,
            reverse=True,
        )
        _resolved_ids = [p.parent.name for p in _run_dirs_diff]
        if _run_a == "latest":
            _run_a = _resolved_ids[0] if _resolved_ids else ""
        if _run_b == "prior":
            _run_b = _resolved_ids[1] if len(_resolved_ids) > 1 else ""

        if not _run_a or not _run_b:
            return {
                "status": "error",
                "text": "Could not resolve run IDs. Please run an analysis first.",
            }

        diff = await _orch_diff.diff_runs(_run_a, _run_b)
        if diff.get("status") == "error":
            return {"status": "error", "text": diff.get("message", "Diff failed.")}

        _param_changes = diff.get("param_changes", {})
        _change_lines = "\n".join(
            f"- **{k}**: `{v['run_a']}` → `{v['run_b']}`"
            for k, v in _param_changes.items()
        ) or "- No parameter changes detected."

        text = (
            f"## Run Comparison\n\n"
            f"**Run A:** `{_run_a[:8]}…`  \n"
            f"**Run B:** `{_run_b[:8]}…`\n\n"
            f"### Parameter Changes\n{_change_lines}\n\n"
            f"### Metrics (Run B)\n{diff.get('delta_a_to_b', {}).get('narrative', 'N/A')}"
        )
        return {
            "status": "success",
            "visualization_type": "text",
            "text": text,
            "result": diff,
        }

    elif tool_name == "handle_natural_command":
        # Use the BioAgent path (system prompt from agent.md) for natural commands
        # Lazy import: Import backend.agent only when needed, not at module import time.
        # This avoids loading heavy LLM dependencies (langgraph, langchain) during server startup,
        # keeping lightweight endpoints like /health and /mcp/tools fast and allowing the service
        # to work in sandbox/CI environments where LLM dependencies may not be installed.
        from backend.agent import handle_command
        command = arguments.get("command", "")
        session_id = arguments.get("session_id", "") or ""
        session_context = dict(arguments.get("session_context") or {})
        if arguments.get("previous_plan_steps"):
            session_context["previous_plan_steps"] = arguments["previous_plan_steps"]
        return await handle_command(command, session_id=session_id, session_context=session_context)
    
    elif tool_name == "phylogenetic_tree":
        # Handle phylogenetic tree analysis — use backend Python-native module
        from backend import phylogenetic_tree as _phylo_mod
        aligned_sequences = arguments.get("aligned_sequences", "")
        raw = _phylo_mod.run_phylogenetic_tree_raw(aligned_sequences)

        # Build a visuals list from the various plot fields the module returns
        _phylo_visuals = []
        if raw.get("tree_plot_b64"):
            _phylo_visuals.append({
                "type": "image_b64",
                "data": raw["tree_plot_b64"],
                "title": "Phylogenetic Tree",
            })
        if raw.get("dist_matrix_b64"):
            _phylo_visuals.append({
                "type": "image_b64",
                "data": raw["dist_matrix_b64"],
                "title": "Pairwise Distance Matrix",
            })
        if raw.get("ete_visualization"):
            _phylo_visuals.append({
                "type": "image_b64",
                "data": raw["ete_visualization"],
                "title": "Annotated Tree (ETE3)",
            })
        if not _phylo_visuals and raw.get("plot"):
            _phylo_visuals.append({
                "type": "image_b64",
                "data": raw["plot"],
                "title": "Distance Heatmap",
            })

        phylo_response = {
            **raw,
            "visuals": _phylo_visuals,
            "visualization_type": "results_viewer",
        }
        _save_analysis_script(
            "phylogenetic_tree",
            {"aligned_sequences": aligned_sequences},
            phylo_response,
            session_id=arguments.get("session_id", ""),
        )
        return phylo_response

    elif tool_name == "clustering_analysis":
        # Handle clustering analysis
        from backend import phylogenetic_tree as _phylo_mod
        aligned_sequences = arguments.get("aligned_sequences", "")
        num_clusters = arguments.get("num_clusters", 5)
        return _phylo_mod.run_clustering_from_tree(aligned_sequences, num_clusters)
    
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
        # Handle single-cell RNA-seq analysis — use backend Python-native module
        from backend import single_cell_analysis

        data_file = arguments.get("data_file") or arguments.get("data_path")

        # Gate: if no data file is provided, ask for it (needs_inputs behaviour).
        if not data_file and not arguments.get("needs_inputs"):
            return _build_needs_inputs_response("single_cell_analysis", arguments)

        data_format = arguments.get("data_format", "10x")
        steps = arguments.get("steps", "all")
        if isinstance(steps, list):
            steps = ",".join(steps)
        resolution = float(arguments.get("resolution", 0.5))

        result = single_cell_analysis.analyze_single_cell_data(
            data_file=data_file,
            data_format=data_format,
            steps=steps,
            resolution=resolution,
            question=arguments.get("question"),
        )

        # Build visuals from base64 plots
        plots = result.get("plots", {})
        _sc_visuals = []
        for label, b64 in plots.items():
            titles = {
                "umap_clusters":   "t-SNE — Cluster Overview",
                "umap_celltype":   "t-SNE — Cell Type Annotation",
                "umap_condition":  "t-SNE — Disease Status",
                "dotplot_markers": "Marker Gene Dot Plot",
            }
            _sc_visuals.append({
                "type": "image_b64",
                "data": b64,
                "title": titles.get(label, label.replace("_", " ").title()),
            })

        sc_response = {
            "status": result.get("status", "success"),
            "visualization_type": "results_viewer",
            "text": result.get("text", "Single-cell analysis complete."),
            "links": [],
            "visuals": _sc_visuals,
            "result": result,
        }
        _save_analysis_script(
            "single_cell_analysis",
            {
                "data_file": data_file,
                "data_format": data_format,
                "steps": steps,
                "resolution": resolution,
            },
            sc_response,
            session_id=arguments.get("session_id", ""),
        )
        return sc_response

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
        from backend import bulk_rnaseq

        count_matrix    = arguments.get("count_matrix", "")
        sample_metadata = arguments.get("sample_metadata", "")
        design_formula  = arguments.get("design_formula", "~condition")
        alpha           = float(arguments.get("alpha", 0.05))

        # Gate: require input paths before running
        if not count_matrix and not sample_metadata:
            return _build_needs_inputs_response("bulk_rnaseq_analysis", arguments)

        x_scale = arguments.get("x_scale", "log2")
        result = bulk_rnaseq.run_deseq2_analysis(
            count_matrix=count_matrix,
            sample_metadata=sample_metadata,
            design_formula=design_formula,
            alpha=alpha,
            x_scale=x_scale,
        )

        if result.get("status") == "error":
            return {
                "status": "error",
                "text": result.get("message", "Bulk RNA-seq analysis failed."),
                "result": result,
            }

        # Build a rich markdown summary
        summary = result.get("summary", [])
        mode = result.get("mode", "real")
        mode_note = (
            "\n\n> ⚠️ *Synthetic demo data used — real S3 data unavailable.*"
            if mode == "synthetic" else ""
        )
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

        # Build visuals from base64 plots returned by bulk_rnaseq module
        plots = result.get("plots", {})
        visuals = []
        for label, b64 in plots.items():
            if label.startswith("volcano_"):
                display = label[8:].replace("__", ": ").replace("_vs_", " vs ").replace("_", " ").title()
                visuals.append({"type": "image_b64", "data": b64, "title": f"Volcano — {display}"})
            elif label == "pca":
                visuals.append({"type": "image_b64", "data": b64, "title": "PCA — Sample Overview"})

        response = {
            "status": "success",
            "visualization_type": "results_viewer",
            "text": text,
            "links": [],
            "visuals": visuals,
            "result": result,
        }
        # Persist a re-runnable analysis.py alongside this result
        _save_analysis_script(
            "bulk_rnaseq_analysis",
            {
                "count_matrix": count_matrix,
                "sample_metadata": sample_metadata,
                "design_formula": design_formula,
                "alpha": alpha,
                "x_scale": x_scale,
            },
            response,
            session_id=arguments.get("session_id", ""),
        )
        return response
    
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
                        "Helix couldn't match your request to a known tool and this looks like Q&A intent. "
                        "Helix will not generate or execute a new tool for Q&A. "
                        "Try examples in /docs/user/CAPABILITIES_GUIDE.md or rephrase as "
                        "'Run <analysis> on <inputs>'."
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


# Backward compatibility alias (to be removed in a future major release).
call_mcp_tool = dispatch_tool

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

@app.on_event("startup")
async def startup_event():
    """Initialize heavy components after server starts."""

    # Log EC2 / AWS config at startup so operators can confirm env is wired correctly
    logger.info("Backend startup — dotenv: %s", _dotenv_path or "NOT FOUND")
    logger.info("  AWS_REGION=%s  HELIX_USE_EC2=%s  HELIX_EC2_AUTO_CREATE=%s",
                os.getenv("AWS_REGION", "NOT SET"),
                os.getenv("HELIX_USE_EC2", "NOT SET"),
                os.getenv("HELIX_EC2_AUTO_CREATE", "NOT SET"))
    
    # Heavy initialization (LLM, agent) happens lazily on first use
    # Session loading happens lazily in history_manager
    logger.info("Backend startup complete - ready to accept requests")

if __name__ == "__main__":
    import uvicorn
    # Run as package module (start.sh uses `python -m backend.main_with_mcp`)
    uvicorn.run("backend.main_with_mcp:app", host="0.0.0.0", port=8001, reload=True)
