from fastapi import FastAPI, HTTPException
from fastapi import File, Request, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, FileResponse
from pydantic import BaseModel
from typing import Dict, Any, Optional, List
import asyncio
import base64
import json
import subprocess
import sys
from pathlib import Path
import os
import shutil
import numpy as np
import logging
import re
from datetime import datetime, timezone
import uuid

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

from backend.history_manager import history_manager, sanitize_command_for_storage
from backend.tool_schemas import list_tool_schemas
from backend.context_builder import _truncate_sequence
from backend.orchestration.action_planner import (
    build_single_step_plan as build_single_step_action_plan,
    infer_action,
)
from backend.orchestration.approval_policy import (
    READ_ONLY_ROUTER_TOOLS,
    has_explicit_execute_intent as approval_has_explicit_execute_intent,
    is_approval_command as approval_is_approval_command,
    requires_approval_semantics as approval_requires_approval_semantics,
    should_stage_for_approval as approval_should_stage_for_approval,
)
from backend.orchestration.execution_router import (
    normalize_tool_selection,
    normalize_fastqc_parameters,
    preflight_tool_bindings as router_preflight_tool_bindings,
)
from backend.orchestration.upload_intake_policy import (
    evaluate_upload_intake,
    get_policy_profile,
    policy_profile_defaults,
)
from backend.orchestration.artifact_resolver import resolve_semantic_reference
from backend.orchestration.tool_registry import dispatch_via_registry
from backend.orchestration.visualization_resolver import determine_visualization_type

from backend.execution_broker import ExecutionBroker, ExecutionRequest


# NOTE: handle_command from backend.agent is lazily imported inline (not at module level)
# to avoid loading heavy LLM dependencies (langgraph, langchain) during server startup.
# This keeps lightweight endpoints like /health and /tools/list fast, and allows the service
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


def _agent_disabled() -> bool:
    """
    True when the LLM agent should be skipped (deterministic router only).
    Use for CI, offline, or testing. Production typically runs with agent enabled.
    Reads HELIX_AGENT_DISABLED=1 or legacy HELIX_MOCK_MODE=1.
    """
    return os.getenv("HELIX_AGENT_DISABLED") == "1" or os.getenv("HELIX_MOCK_MODE") == "1"


async def _run_agent_with_retry(coro_factory, timeout_s: int, retries: int = 1, backoff_s: float = 0.5):
    """
    Run an async agent call with timeout + simple retry/backoff.
    Returns (result, diagnostics).
    """
    diagnostics = {
        "attempts": 0,
        "timeouts": 0,
        "errors": [],
        "fallback_used": False,
    }
    last_exc: Optional[Exception] = None
    for attempt in range(retries + 1):
        diagnostics["attempts"] = attempt + 1
        try:
            result = await asyncio.wait_for(coro_factory(), timeout=timeout_s)
            return result, diagnostics
        except Exception as e:
            last_exc = e
            name = e.__class__.__name__
            if name == "TimeoutError":
                diagnostics["timeouts"] += 1
            diagnostics["errors"].append({"type": name, "message": str(e)})
            if attempt < retries:
                await asyncio.sleep(backoff_s * (2 ** attempt))
    diagnostics["fallback_used"] = True
    err = "AGENT_TIMEOUT" if (last_exc and last_exc.__class__.__name__ == "TimeoutError") else "AGENT_ERROR"
    return {
        "status": "error",
        "success": False,
        "error": err,
        "text": f"Agent execution timed out or failed: {last_exc}",
        "diagnostics": diagnostics,
    }, diagnostics


_READ_ONLY_ROUTER_TOOLS = READ_ONLY_ROUTER_TOOLS


def _is_approval_command(command: str) -> bool:
    return approval_is_approval_command(command)


def _has_explicit_execute_intent(command: str) -> bool:
    return approval_has_explicit_execute_intent(command)


def _requires_approval_semantics(command: str) -> bool:
    return approval_requires_approval_semantics(command)


def _should_stage_for_approval(tool_name: str, command: str, params: Optional[Dict[str, Any]] = None) -> bool:
    return approval_should_stage_for_approval(
        tool_name,
        command,
        params,
        action_type=infer_action(command, tool_name),
    )


def _historical_recreation_ready_for_execution(
    command: str,
    tool_name: str,
    session_context: Optional[Dict[str, Any]],
) -> bool:
    c = (command or "").lower()
    if tool_name not in {"bio_diff_runs", "bio_rerun"}:
        return False
    if not any(k in c for k in ("recreate", "reconstruct", "reproduce")):
        return False
    if not any(k in c for k in ("before", "first", "original", "prior", "version", "state")):
        return False
    if not isinstance(session_context, dict):
        return False
    resolved = resolve_semantic_reference(session_context, command)
    return str(resolved.get("status") or "") == "resolved"


def _preflight_tool_bindings(tool_name: str, parameters: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    return router_preflight_tool_bindings(tool_name, parameters, _build_needs_inputs_response)


def _build_single_step_plan(command: str, tool_name: str, params: Dict[str, Any]) -> Dict[str, Any]:
    return build_single_step_action_plan(command, tool_name, params)


def _autobind_plan_inputs(plan: Dict[str, Any], session_context: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Best-effort input binding for approved plans.
    Uses session context first, then per-tool requirement examples as demo fallbacks.
    """
    if not isinstance(plan, dict):
        return plan
    steps = plan.get("steps")
    if not isinstance(steps, list):
        return plan

    req_registry = globals().get("_TOOL_INPUT_REQUIREMENTS", {}) or {}
    sess = session_context or {}

    for step in steps:
        if not isinstance(step, dict):
            continue
        tool_name = step.get("tool_name")
        args = step.get("arguments")
        if not isinstance(args, dict):
            args = {}
            step["arguments"] = args

        spec = req_registry.get(tool_name) if isinstance(req_registry, dict) else None
        required = spec.get("required_inputs", []) if isinstance(spec, dict) else []
        for req_inp in required:
            if not isinstance(req_inp, dict):
                continue
            key = req_inp.get("name")
            if not key:
                continue
            cur = args.get(key)
            if cur not in (None, "", [], {}):
                continue
            # 1) session-derived value (if present)
            sess_val = sess.get(key)
            if sess_val not in (None, "", [], {}):
                args[key] = sess_val
                continue
            # 2) requirement example fallback (lets synthetic/demo-capable tools run)
            example = req_inp.get("example")
            if example not in (None, "", [], {}):
                args[key] = example

        # Only clear routing flag when required fields are now concretely bound.
        if args.get("needs_inputs") is True:
            missing_after_bind = []
            for req_inp in required:
                if not isinstance(req_inp, dict):
                    continue
                key = req_inp.get("name")
                if not key:
                    continue
                if args.get(key) in (None, "", [], {}):
                    missing_after_bind.append(key)
            if not missing_after_bind:
                args["needs_inputs"] = False

    return plan


def _validate_plan_bindings(plan: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """
    Validate required per-step input bindings before execution.
    Returns a structured diagnostic payload when validation fails, else None.
    """
    if not isinstance(plan, dict):
        return {
            "status": "error",
            "error_type": "invalid_plan",
            "message": "Plan payload is not a valid object.",
            "issues": [{"type": "invalid_plan", "detail": "Expected dict plan payload"}],
        }
    steps = plan.get("steps")
    if not isinstance(steps, list):
        return {
            "status": "error",
            "error_type": "invalid_plan",
            "message": "Plan has no executable steps.",
            "issues": [{"type": "invalid_plan", "detail": "Missing steps list"}],
        }

    req_registry = globals().get("_TOOL_INPUT_REQUIREMENTS", {}) or {}
    issues: List[Dict[str, Any]] = []
    for idx, step in enumerate(steps, start=1):
        if not isinstance(step, dict):
            continue
        tool_name = step.get("tool_name")
        args = step.get("arguments") if isinstance(step.get("arguments"), dict) else {}
        spec = req_registry.get(tool_name) if isinstance(req_registry, dict) else None
        required = spec.get("required_inputs", []) if isinstance(spec, dict) else []
        missing = []
        for req_inp in required:
            if not isinstance(req_inp, dict):
                continue
            key = req_inp.get("name")
            if not key:
                continue
            if args.get(key) in (None, "", [], {}):
                missing.append(key)
        if missing:
            issues.append(
                {
                    "step_index": idx,
                    "step_id": step.get("id") or f"step{idx}",
                    "tool_name": tool_name,
                    "type": "missing_artifact",
                    "missing_inputs": missing,
                }
            )

    if not issues:
        return None
    return {
        "status": "error",
        "error_type": "artifact_binding_failed",
        "message": "Plan execution blocked due to unresolved required inputs.",
        "issues": issues,
    }


def _analyze_plan_execution(steps: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Classify executed-plan step outcomes for top-level response status/text."""
    step_results = []
    for step in steps or []:
        res = step.get("result")
        if isinstance(res, dict) and res:
            step_results.append(res)
    if not step_results:
        return {"executed": False, "status": "workflow_planned", "has_needs_inputs": False}
    statuses = [str(r.get("status", "")).lower() for r in step_results]
    has_error = any(s in {"error", "failed"} for s in statuses)
    has_needs_inputs = any((r.get("needs_inputs") is True) or (str(r.get("status", "")).lower() == "needs_inputs") for r in step_results)
    if has_error:
        return {"executed": True, "status": "error", "has_needs_inputs": has_needs_inputs}
    if has_needs_inputs:
        return {"executed": True, "status": "needs_inputs", "has_needs_inputs": True}
    return {"executed": True, "status": "success", "has_needs_inputs": False}


def _store_pending_plan(session_id: str, command: str, plan: Dict[str, Any]) -> None:
    if not (hasattr(history_manager, "sessions") and session_id in history_manager.sessions):
        return
    command_storage = sanitize_command_for_storage(command)
    history_manager.sessions[session_id]["pending_plan"] = {
        "command": command_storage["command"],
        "plan": plan,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "command_storage": command_storage["meta"],
    }


def _get_pending_plan(session_id: str) -> Optional[Dict[str, Any]]:
    if not (hasattr(history_manager, "sessions") and session_id in history_manager.sessions):
        return None
    pending = history_manager.sessions[session_id].get("pending_plan")
    return pending if isinstance(pending, dict) else None


def _clear_pending_plan(session_id: str) -> None:
    if not (hasattr(history_manager, "sessions") and session_id in history_manager.sessions):
        return
    history_manager.sessions[session_id].pop("pending_plan", None)


def _pop_pending_plan(session_id: str) -> Optional[Dict[str, Any]]:
    pending = _get_pending_plan(session_id)
    _clear_pending_plan(session_id)
    return pending


def _should_clear_pending_plan_after_execution(result: Any) -> bool:
    """
    Keep pending plans around when execution still needs inputs, but clear them
    once execution was actually submitted/completed or hard-failed.
    """
    if not isinstance(result, dict):
        return False
    inner = result.get("result") or result
    if isinstance(inner, dict) and inner.get("type") == "execution_result" and isinstance(inner.get("result"), dict):
        inner = inner["result"]
    if isinstance(inner, dict) and inner.get("type") == "job":
        return True
    if isinstance(inner, dict) and inner.get("type") == "plan_result" and isinstance(inner.get("steps"), list):
        plan_exec = _analyze_plan_execution(inner.get("steps") or [])
        return bool(plan_exec.get("executed")) and str(plan_exec.get("status", "")).lower() in {"success", "error"}
    status = str(
        (inner.get("status") if isinstance(inner, dict) else "")
        or result.get("status")
        or ""
    ).lower()
    return status in {"success", "completed", "pipeline_submitted", "submitted"}


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

class HelixToolRequest(BaseModel):
    tool_name: str
    arguments: Dict[str, Any]
    session_id: Optional[str] = None

class HelixToolResponse(BaseModel):
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
# Max upload size for hosted/small-dataset mode. Enforced on prompt attachments and dataset register.
# Set HELIX_MAX_UPLOAD_BYTES (bytes) or HELIX_MAX_UPLOAD_MB (megabytes). Default 20 MB; set to 0 to disable limit.
TWENTY_MB_BYTES = 20 * 1024 * 1024
UPLOAD_CHUNK_BYTES = 1024 * 1024  # 1MB
ALLOWED_UPLOAD_EXTENSIONS = {
    # Sequence / FASTA / FASTQ
    ".fasta",
    ".fa",
    ".fas",
    ".fastq",
    ".fq",
    ".gz",
    # Tabular / structured data
    ".csv",
    ".tsv",
    ".xlsx",
    ".xls",
    # Variants
    ".vcf",
    ".bcf",
    # Genomic intervals / annotation
    ".bed",
    ".gff",
    ".gff3",
    ".gtf",
    # Single-cell
    ".h5ad",
    ".loom",
    ".h5",
    # Alignment
    ".sam",
    ".bam",
    ".cram",
    # Generic text
    ".txt",
}

def _get_max_upload_bytes() -> int:
    """Return max allowed upload size in bytes. 0 means no limit (backend behavior unchanged)."""
    val = os.environ.get("HELIX_MAX_UPLOAD_BYTES")
    if val is not None and val != "":
        try:
            return int(val)
        except ValueError:
            pass
    val = os.environ.get("HELIX_MAX_UPLOAD_MB")
    if val is not None and val != "":
        try:
            return int(float(val) * 1024 * 1024)
        except ValueError:
            pass
    return TWENTY_MB_BYTES


def _sanitize_uploaded_filename(filename: str) -> str:
    """Normalize a user-provided filename to a safe local file name."""
    base = Path(filename or "").name.strip()
    if not base:
        base = f"upload-{uuid.uuid4().hex[:8]}.dat"
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", base)
    sanitized = sanitized.lstrip(".")
    if not sanitized:
        sanitized = f"upload-{uuid.uuid4().hex[:8]}.dat"
    return sanitized


def _is_allowed_uploaded_filename(filename: str) -> bool:
    lower_name = filename.lower()
    if lower_name.endswith(".fastq.gz") or lower_name.endswith(".fq.gz"):
        return True
    return any(lower_name.endswith(ext) for ext in ALLOWED_UPLOAD_EXTENSIONS)


def _emit_policy_audit_event(event_type: str, session_id: str, payload: Dict[str, Any]) -> None:
    event = {
        "event_type": event_type,
        "session_id": session_id,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "payload": payload,
    }
    logger.info("POLICY_AUDIT %s", json.dumps(event, default=str))


def _get_policy_pending_uploads(session_id: str) -> List[Dict[str, Any]]:
    session = history_manager.get_session(session_id) or {}
    metadata = session.get("metadata", {})
    uploaded_files = metadata.get("uploaded_files", []) or []
    pending: List[Dict[str, Any]] = []
    for entry in uploaded_files:
        if bool(entry.get("requires_policy_approval")) and entry.get("policy_state") == "approval_required":
            pending.append(entry)
    return pending


def _approve_policy_pending_uploads(session_id: str, approver_note: str) -> int:
    session = history_manager.get_session(session_id)
    if not session:
        return 0
    metadata = session.setdefault("metadata", {})
    uploaded_files = metadata.setdefault("uploaded_files", [])
    approved_count = 0
    for entry in uploaded_files:
        if bool(entry.get("requires_policy_approval")) and entry.get("policy_state") == "approval_required":
            entry["policy_state"] = "approved_for_execution"
            entry["policy_approved_at"] = datetime.now(timezone.utc).isoformat()
            entry["policy_approval_note"] = approver_note
            approved_count += 1
    if approved_count > 0:
        session["updated_at"] = datetime.now().isoformat()
        history_manager._save_session(session_id)  # noqa: SLF001 - central session persistence hook
    return approved_count


def _make_unique_path(path: Path) -> Path:
    if not path.exists():
        return path
    stem = path.stem
    suffix = path.suffix
    counter = 1
    while True:
        candidate = path.with_name(f"{stem}_{counter}{suffix}")
        if not candidate.exists():
            return candidate
        counter += 1

MAX_PROMPT_TOKENS = 20000  # Relax limit to avoid false 413s for longer prompts
MAX_PROMPTS_PER_DAY = 100
MAX_SESSIONS_PER_HOUR = int(os.getenv("HELIX_MAX_SESSIONS_PER_HOUR", "60"))

# In-memory per-identity counters. In production, move to Redis or DB.
_daily_prompt_counters: Dict[str, Dict[str, int]] = {}
_hourly_session_create_counters: Dict[str, Dict[str, int]] = {}

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


def _current_hour_bucket() -> str:
    now = datetime.now(timezone.utc)
    return now.strftime("%Y-%m-%dT%H:00Z")


def _check_and_increment_session_create_counter(request: Request) -> None:
    """Limit session creation bursts by source IP to reduce bot-driven session churn."""
    client_ip = request.client.host if request and request.client else "unknown"
    # Keep local development flows unrestricted.
    if client_ip in {"127.0.0.1", "::1"}:
        return
    by_ip = _hourly_session_create_counters.setdefault(client_ip, {})
    hour = _current_hour_bucket()
    count = by_ip.get(hour, 0)
    if count >= MAX_SESSIONS_PER_HOUR:
        raise HTTPException(
            status_code=429,
            detail=(
                "Session creation rate limit exceeded for this source. "
                f"Try again next hour (limit: {MAX_SESSIONS_PER_HOUR}/hour)."
            ),
        )
    by_ip[hour] = count + 1

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


def _materialize_run_artifacts(
    *,
    session_id: str,
    tool: str,
    result: Any,
    tool_args: Optional[Dict[str, Any]] = None,
    command: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Persist executable/downloadable artifacts for a run and return metadata
    additions for history_manager.add_history_entry().

    This runs at execution time (not on download click) so bundles are complete
    and reproducibility is deterministic.
    """
    if not session_id or not isinstance(result, dict):
        return {}

    inner = result.get("result") if isinstance(result.get("result"), dict) else {}
    run_id = (
        result.get("run_id")
        or inner.get("run_id")
        or f"run_{uuid.uuid4().hex[:12]}"
    )

    sessions_root = Path(__file__).parent.parent / "sessions"
    run_dir = sessions_root / session_id / "runs" / run_id
    plots_dir = run_dir / "plots"
    tables_dir = run_dir / "tables"
    plots_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    artifacts: List[Dict[str, Any]] = []

    def _safe_name(text: str, fallback: str) -> str:
        cleaned = "".join(ch if (ch.isalnum() or ch in ("-", "_")) else "_" for ch in (text or ""))
        cleaned = cleaned.strip("_")
        return cleaned[:80] or fallback

    # Persist script for scriptable tools.
    if tool in _SCRIPTABLE_TOOLS:
        try:
            from backend import code_generator as _cg
            params = dict(tool_args or {})
            params.pop("_from_broker", None)
            script_text = _cg.generate(
                tool,
                params,
                run_id=run_id,
                session_id=session_id,
                run_dir=str(run_dir),
            )
            if script_text:
                script_path = run_dir / "analysis.py"
                script_path.write_text(script_text)
                artifacts.append(
                    {
                        "type": "script",
                        "title": "analysis.py",
                        "uri": str(script_path),
                        "format": "py",
                    }
                )
                result["script_path"] = str(script_path)
        except Exception as exc:
            logger.debug("Could not persist analysis.py for run %s: %s", run_id, exc)

    visuals = result.get("visuals")
    if not isinstance(visuals, list) and isinstance(inner, dict):
        visuals = inner.get("visuals")
    if not isinstance(visuals, list):
        visuals = []

    # Persist image_b64 visuals as PNG files.
    for idx, visual in enumerate(visuals, start=1):
        if not isinstance(visual, dict):
            continue
        if visual.get("type") != "image_b64":
            continue
        data_b64 = visual.get("data")
        if not isinstance(data_b64, str) or not data_b64.strip():
            continue
        if "," in data_b64 and data_b64.strip().startswith("data:"):
            data_b64 = data_b64.split(",", 1)[1]
        title = visual.get("title") or f"plot_{idx}"
        fname = f"{idx:02d}_{_safe_name(title, f'plot_{idx}')}.png"
        path = plots_dir / fname
        try:
            path.write_bytes(base64.b64decode(data_b64))
            artifacts.append(
                {
                    "type": "plot",
                    "title": title,
                    "uri": str(path),
                    "format": "png",
                }
            )
        except Exception as exc:
            logger.debug("Could not persist b64 visual '%s': %s", title, exc)

    # Persist Newick tree if available.
    tree_newick = (
        result.get("tree_newick")
        or inner.get("tree_newick")
        or result.get("newick_tree")
        or inner.get("newick_tree")
    )
    if isinstance(tree_newick, str) and tree_newick.strip():
        newick_path = tables_dir / "tree.newick"
        try:
            newick_path.write_text(tree_newick)
            artifacts.append(
                {
                    "type": "table",
                    "title": "tree.newick",
                    "uri": str(newick_path),
                    "format": "newick",
                }
            )
        except Exception as exc:
            logger.debug("Could not persist Newick tree: %s", exc)

    # Persist raw result snapshot for reproducibility.
    try:
        snapshot_path = tables_dir / "result.json"
        snapshot_path.write_text(json.dumps(result, indent=2, default=str))
        artifacts.append(
            {
                "type": "table",
                "title": "result.json",
                "uri": str(snapshot_path),
                "format": "json",
            }
        )
    except Exception as exc:
        logger.debug("Could not persist result snapshot for run %s: %s", run_id, exc)

    # Persist a normalized run export layout for easier navigation.
    try:
        export_root = run_dir / "export_v1"
        plan_dir = export_root / "01_plan_input"
        execute_dir = export_root / "02_execute"
        jobs_dir = export_root / "03_jobs"
        artifacts_dir = export_root / "04_artifacts" / "downloads"
        bundles_dir = export_root / "05_bundles"
        verification_dir = export_root / "06_verification"
        for d in (plan_dir, execute_dir, jobs_dir, artifacts_dir, bundles_dir, verification_dir):
            d.mkdir(parents=True, exist_ok=True)

        request_payload = {
            "session_id": session_id,
            "run_id": run_id,
            "tool": tool,
            "command": command,
            "tool_args": dict(tool_args or {}),
        }
        (plan_dir / "request.json").write_text(json.dumps(request_payload, indent=2, default=str))
        (execute_dir / "stage2.json").write_text(json.dumps(result, indent=2, default=str))

        discovered_jobs: List[Dict[str, Any]] = []
        data = result.get("data") if isinstance(result.get("data"), dict) else {}
        data_results = data.get("results") if isinstance(data.get("results"), dict) else {}
        raw_result = result.get("raw_result") if isinstance(result.get("raw_result"), dict) else {}
        inner_result = result.get("result") if isinstance(result.get("result"), dict) else {}

        for candidate in (
            result.get("jobs"),
            raw_result.get("jobs"),
            inner_result.get("jobs"),
            data_results.get("jobs"),
        ):
            if isinstance(candidate, list):
                for job in candidate:
                    if isinstance(job, dict) and job.get("job_id"):
                        discovered_jobs.append(job)

        # Also support single-job async response shapes.
        for candidate in (
            result.get("job"),
            raw_result.get("job"),
            inner_result.get("job"),
            result.get("result"),
            raw_result.get("result"),
            data_results.get("result"),
        ):
            if isinstance(candidate, dict) and candidate.get("job_id"):
                discovered_jobs.append(candidate)

        deduped_jobs: List[Dict[str, Any]] = []
        seen_job_ids = set()
        for job in discovered_jobs:
            job_id = str(job.get("job_id") or "")
            if not job_id or job_id in seen_job_ids:
                continue
            seen_job_ids.add(job_id)
            deduped_jobs.append(job)

        # Ensure a consistent run-level job exists for every execution.
        run_job_id = f"{run_id}__run"
        run_job_status = "completed" if _is_success(result) else "failed"
        run_job = {
            "job_id": run_job_id,
            "job_type": "run",
            "tool_name": tool,
            "status": run_job_status,
            "session_id": session_id,
            "run_id": run_id,
            "message": "Run-level execution summary record.",
            "generated_at": datetime.now(timezone.utc).isoformat(),
        }
        deduped_jobs.insert(0, run_job)

        (jobs_dir / "jobs_index.json").write_text(
            json.dumps({"jobs": deduped_jobs}, indent=2, default=str)
        )
        (jobs_dir / f"job_{run_job_id}.json").write_text(
            json.dumps(run_job, indent=2, default=str)
        )

        copied_files: List[str] = []
        for idx, art in enumerate(artifacts, start=1):
            if not isinstance(art, dict):
                continue
            uri = art.get("uri")
            if not isinstance(uri, str):
                continue
            src = Path(uri)
            if not src.exists() or not src.is_file():
                continue
            title = str(art.get("title") or src.name)
            safe = "".join(ch if (ch.isalnum() or ch in ("-", "_", ".")) else "_" for ch in title).strip("_")
            dst_name = f"{idx:02d}_{safe or src.name}"
            dst = artifacts_dir / dst_name
            shutil.copy2(src, dst)
            copied_files.append(dst_name)

        (bundles_dir / "README.txt").write_text(
            "Bundle ZIPs are generated on demand via /download/bundle for this run.\n"
        )

        verification_payload = {
            "session_id": session_id,
            "run_id": run_id,
            "tool": tool,
            "success": _is_success(result),
            "download_count": len(copied_files),
            "has_jobs": len(deduped_jobs) > 0,
            "generated_at": datetime.now(timezone.utc).isoformat(),
        }
        (verification_dir / "verification.json").write_text(
            json.dumps(verification_payload, indent=2, default=str)
        )
    except Exception as exc:
        logger.debug("Could not persist export_v1 layout for run %s: %s", run_id, exc)

    return {
        "run_id": run_id,
        "produced_artifacts": artifacts,
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
    effective_tool_args = dict(tool_args) if isinstance(tool_args, dict) else tool_args
    if isinstance(effective_tool_args, dict) and "action_type" not in effective_tool_args:
        try:
            from backend.action_plan import infer_action_type

            effective_tool_args["action_type"] = infer_action_type(req.command, tool)
        except Exception:
            pass

    if record_history:
        metadata = _extract_metadata(result, tool_args=effective_tool_args)
        persisted = _materialize_run_artifacts(
            session_id=req.session_id,
            tool=tool,
            result=result,
            tool_args=effective_tool_args,
            command=req.command,
        )
        if isinstance(persisted, dict):
            persisted_run_id = persisted.get("run_id")
            if persisted_run_id and isinstance(result, dict) and not result.get("run_id"):
                result["run_id"] = persisted_run_id
            if persisted_run_id and not metadata.get("run_id"):
                metadata["run_id"] = persisted_run_id
            persisted_arts = persisted.get("produced_artifacts") or []
            existing_arts = metadata.get("produced_artifacts") or []
            if isinstance(existing_arts, list) and isinstance(persisted_arts, list):
                metadata["produced_artifacts"] = existing_arts + persisted_arts
        history_manager.add_history_entry(
            req.session_id,
            req.command,
            tool,
            result,
            metadata=metadata,
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


async def _execute_routed_tool(
    req: CommandRequest,
    *,
    tool_name: str,
    parameters: Dict[str, Any],
    session_context: Dict[str, Any],
    execution_path: str,
    validation_execution_path: str,
) -> "CustomJSONResponse":
    """
    Shared execution path for both agent-mapped and fallback-routed tool calls.
    """
    selection = normalize_tool_selection(
        command=req.command,
        tool_name=tool_name,
        parameters=parameters,
        session_context=session_context or {},
    )
    selected_tool = selection.tool_name
    selected_args = selection.arguments

    norm = normalize_fastqc_parameters(
        tool_name=selected_tool,
        arguments=selected_args,
        command=req.command,
        session_context=session_context or {},
    )
    params = norm.arguments
    if not norm.ok and isinstance(norm.result, dict):
        return await _dispatch_result(
            req,
            selected_tool,
            norm.result,
            tool_args=params,
            execution_path=validation_execution_path,
        )

    preflight = _preflight_tool_bindings(selected_tool, params)
    if isinstance(preflight, dict):
        return await _dispatch_result(
            req,
            selected_tool,
            preflight,
            tool_args=params,
            execution_path=f"{validation_execution_path}_binding",
        )

    if selected_tool == "fastqc_quality_analysis":
        if "session_id" not in params:
            params["session_id"] = req.session_id
        for param_name in ["input_r1", "input_r2", "output"]:
            if param_name in params and params[param_name]:
                uri = params[param_name]
                if isinstance(uri, str) and uri.startswith("//") and not uri.startswith("s3://"):
                    params[param_name] = "s3:" + uri

    broker = _get_execution_broker()
    result = await broker.execute_tool(
        ExecutionRequest(
            tool_name=selected_tool,
            arguments=params,
            session_id=req.session_id,
            original_command=req.command,
            session_context=session_context,
        )
    )
    _apply_session_context_side_effects(req.session_id, selected_tool, result)
    return await _dispatch_result(req, selected_tool, result, tool_args=params, execution_path=execution_path)


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
    return determine_visualization_type(tool, result, prompt)


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


def _build_actionable_fallback_text(
    tool: str,
    status: str,
    text: str,
    truncated_result: Any,
) -> str:
    """
    Guard against empty/non-actionable success text in benchmark-critical paths.
    """
    if isinstance(text, str) and text.strip():
        return text
    status_l = str(status or "").lower()
    if status_l in {"error", "failed"}:
        return text or "Execution failed before an actionable response could be produced."

    diagnostics = {}
    if isinstance(truncated_result, dict):
        diagnostics = truncated_result.get("diagnostics") if isinstance(truncated_result.get("diagnostics"), dict) else {}
        if not diagnostics:
            maybe_result = truncated_result.get("result")
            if isinstance(maybe_result, dict) and isinstance(maybe_result.get("selector_diagnostics"), dict):
                diagnostics = maybe_result.get("selector_diagnostics") or {}

    if tool in {"bio_rerun", "bio_diff_runs", "patch_and_rerun"}:
        unresolved = diagnostics.get("unresolved") if isinstance(diagnostics.get("unresolved"), list) else []
        unresolved_msg = ""
        if unresolved:
            unresolved_msg = f" Unresolved selectors: {', '.join(str(x) for x in unresolved[:3])}."
        return (
            "The request did not produce an executable historical state yet. "
            "Please specify an explicit reference such as `first DEG results` vs `current DEG results`, "
            "or provide the concrete `run_id` to continue."
            + unresolved_msg
        )

    if tool == "go_enrichment_analysis":
        return (
            "GO enrichment was requested but no actionable gene list was available. "
            "Provide a `gene_list` directly or reference a DEG state (for example `current DEG results`)."
        )

    return (
        "The request returned without actionable output. "
        "Please provide the target artifact/state or required inputs so Helix can execute deterministically."
    )


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

    # __plan__ / workflow_plan: ensure frontend always gets markdown (never raw key-value dump)
    if tool == "__plan__" and isinstance(truncated_result, dict):
        # Unwrap broker envelope so we see plan_result or job
        inner = truncated_result.get("result") or truncated_result
        _was_execution_envelope = False
        if inner.get("type") == "execution_result" and isinstance(inner.get("result"), dict):
            _was_execution_envelope = True
            inner = inner["result"]
        if isinstance(inner, dict) and inner.get("type") == "plan_result" and isinstance(inner.get("steps"), list):
            steps_list = inner["steps"]
            plan_exec = _analyze_plan_execution(steps_list)
            lines = []
            for i, step in enumerate(steps_list, 1):
                step_id = step.get("id") or f"step{i}"
                tool_name_step = step.get("tool_name") or ""
                lines.append(f"{i}. **{step_id}** (`{tool_name_step}`)")
            steps_md = "\n".join(lines) if lines else "(no steps)"
            inner_status = str(inner.get("status") or "").lower()
            execute_ready_flag = bool(inner.get("execute_ready")) if "execute_ready" in inner else True
            approval_required_flag = bool(inner.get("approval_required")) if "approval_required" in inner else False
            binding_diag = inner.get("binding_diagnostics")
            if _was_execution_envelope or plan_exec.get("executed"):
                if plan_exec.get("status") == "needs_inputs":
                    text = (
                        "## Pipeline Execution Needs Inputs\n\n"
                        "The approved plan was executed, but one or more steps still need additional inputs.\n\n"
                        "**Steps:**\n\n" + steps_md
                    )
                elif plan_exec.get("status") == "error":
                    text = (
                        "## Pipeline Execution Failed\n\n"
                        "One or more steps failed during plan execution.\n\n"
                        "**Steps:**\n\n" + steps_md
                    )
                else:
                    text = (
                        "## Pipeline Execution Complete\n\n"
                        "The approved plan has been executed.\n\n"
                        "**Steps:**\n\n" + steps_md
                    )
                status = plan_exec.get("status") or "success"
            else:
                if inner_status == "needs_inputs" or not execute_ready_flag:
                    missing_hint = ""
                    if isinstance(binding_diag, dict):
                        issues = binding_diag.get("issues")
                        if isinstance(issues, list) and issues:
                            first = issues[0] if isinstance(issues[0], dict) else {}
                            missing = first.get("missing_inputs")
                            if isinstance(missing, list) and missing:
                                missing_hint = (
                                    "\n\n**Missing inputs:** "
                                    + ", ".join(f"`{str(k)}`" for k in missing if str(k).strip())
                                )
                    text = (
                        "## Pipeline Plan (Needs Inputs)\n\n"
                        "**Steps:**\n\n" + steps_md + "\n\n"
                        "Required inputs are still missing before this pipeline can execute."
                        + missing_hint
                    )
                    status = "needs_inputs"
                else:
                    if approval_required_flag:
                        # When binding_diagnostics is present, there are missing file inputs.
                        # The plan is staged for approval but execution needs those files.
                        if isinstance(binding_diag, dict) and binding_diag.get("issues"):
                            _missing_for_text: list = []
                            for _issue in (binding_diag.get("issues") or []):
                                if isinstance(_issue, dict):
                                    _missing_for_text.extend(_issue.get("missing_inputs") or [])
                            _missing_note = (
                                "\n\n*Provide `" + "`, `".join(str(m) for m in _missing_for_text if str(m).strip()) + "` then click **I approve** to run.*"
                                if _missing_for_text
                                else "\n\n*Approve this plan to proceed.*"
                            )
                            text = (
                                "## Pipeline Plan\n\n"
                                "**Steps:**\n\n" + steps_md + _missing_note
                            )
                        else:
                            text = (
                                "## Pipeline Plan\n\n"
                                "**Steps:**\n\n" + steps_md + "\n\n"
                                "*All inputs are available. Click **I approve** to execute these reviewed steps.*"
                            )
                    else:
                        text = (
                            "## Pipeline Plan\n\n"
                            "**Steps:**\n\n" + steps_md + "\n\n"
                            "*All inputs are available. Click **Execute Pipeline** to run these steps.*"
                        )
                    status = "workflow_planned"
        elif isinstance(inner, dict) and inner.get("type") == "job":
            text = (
                "## Pipeline Execution Submitted\n\n"
                "The approved workflow is now running asynchronously. Track progress in the **Jobs** panel."
            )
            status = "pipeline_submitted"
        elif not text:
            text = (
                "## Pipeline Plan\n\n"
                "Your workflow has been planned. Use **Execute Pipeline** to run it, or **Load & run** to use example data."
            )
            status = "workflow_planned"

    text = _build_actionable_fallback_text(tool, status, text, truncated_result)

    # Download links should only be surfaced for completed/executed outputs.
    # Planning, needs-inputs, and in-progress states should not expose downloads.
    _non_executed_statuses = {
        "workflow_planned",
        "needs_inputs",
        "tool_mapped",
        "pipeline_submitted",
        "workflow_needs_clarification",
        "submitted",
        "running",
        "queued",
        "pending",
    }
    _allow_download_links = bool(
        success and status not in {"error", "failed", "workflow_failed"} and status not in _non_executed_statuses
    )

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
    approval_required = bool(
        isinstance(truncated_result, dict) and truncated_result.get("approval_required")
    )
    # __plan__ with plan_result: frontend expects execute_ready so it shows the Execute Pipeline button
    if tool == "__plan__" and isinstance(truncated_result, dict):
        inner = truncated_result.get("result") or truncated_result
        _was_execution_envelope = False
        if isinstance(inner, dict) and inner.get("type") == "execution_result" and isinstance(inner.get("result"), dict):
            _was_execution_envelope = True
            inner = inner["result"]
        if isinstance(inner, dict) and inner.get("type") == "plan_result" and isinstance(inner.get("steps"), list):
            plan_exec = _analyze_plan_execution(inner.get("steps") or [])
            if "execute_ready" in inner:
                execute_ready = bool(inner.get("execute_ready"))
            else:
                execute_ready = False if (_was_execution_envelope or plan_exec.get("executed")) else True
            if "approval_required" in inner:
                approval_required = bool(inner.get("approval_required"))

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

    # Start from data.links and always merge links present on the unwrapped result.
    _existing_links: List[Dict] = list(data.get("links") or [])
    for _lnk in (_inner_result.get("links") or []):
        if isinstance(_lnk, dict) and _lnk not in _existing_links:
            _existing_links.append(_lnk)

    # Enrich links from the run-ledger artifact registry so every completed
    # demo/pipeline response can expose concrete downloadable outputs.
    _effective_run_id: Optional[str] = _run_id_hint
    if not _effective_run_id and session_id:
        try:
            _runs = history_manager.list_runs(session_id)
            if isinstance(_runs, list) and _runs:
                _latest = _runs[-1] if isinstance(_runs[-1], dict) else None
                if _latest:
                    _effective_run_id = _latest.get("run_id")
        except Exception:
            pass

    if _allow_download_links and session_id and _effective_run_id:
        try:
            _run = history_manager.get_run(session_id, _effective_run_id)
            _arts = _run.get("produced_artifacts") if isinstance(_run, dict) else None
            if isinstance(_arts, list):
                _seen_urls = {
                    lnk.get("url")
                    for lnk in _existing_links
                    if isinstance(lnk, dict) and isinstance(lnk.get("url"), str)
                }
                for _art in _arts:
                    if not isinstance(_art, dict):
                        continue
                    _artifact_id = _art.get("artifact_id")
                    if not _artifact_id:
                        continue
                    _url = f"/session/{session_id}/artifacts/{_artifact_id}/download"
                    if _url in _seen_urls:
                        continue
                    _existing_links.append(
                        {
                            "label": _art.get("title") or _art.get("name") or "artifact",
                            "url": _url,
                            "artifact_id": _artifact_id,
                            "type": _art.get("type"),
                            "format": _art.get("format"),
                        }
                    )
                    _seen_urls.add(_url)
        except Exception:
            pass

    if _allow_download_links and _script_path_str:
        _script_url = f"/download/script?path={_script_path_str}"
        _bundle_url = (
            f"/download/bundle?session_id={session_id}"
            + (f"&run_id={_run_id_hint}" if _run_id_hint else "")
        )
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

    # For successful responses with a session, always provide a reproducibility
    # bundle link so demos and pipelines consistently expose downloadable output.
    if _allow_download_links and session_id and not any(
        isinstance(lnk, dict) and "bundle" in (lnk.get("label") or "").lower()
        for lnk in _existing_links
    ):
        _bundle_url = (
            f"/download/bundle?session_id={session_id}"
            + (f"&run_id={_effective_run_id}" if _effective_run_id else "")
        )
        _existing_links.append({"label": "bundle.zip", "url": _bundle_url})

    # If a Newick tree was generated, expose it as a direct downloadable artifact.
    _tree_newick = (
        _inner_result.get("tree_newick")
        or _inner_result.get("newick_tree")
        or tree_data.get("tree_newick")
    )
    if _allow_download_links and isinstance(_tree_newick, str) and _tree_newick.strip():
        from urllib.parse import quote as _url_quote
        _newick_url = f"data:text/plain;charset=utf-8,{_url_quote(_tree_newick)}"
        if not any(
            isinstance(lnk, dict) and "newick" in (lnk.get("label") or "").lower()
            for lnk in _existing_links
        ):
            _existing_links.append({"label": "tree.newick", "url": _newick_url})

    # Keep links/downloads strictly for executed outputs only.
    # Planning / needs-inputs / in-progress states must not expose result artifacts.
    data["links"] = _existing_links if _allow_download_links else []
    if _allow_download_links:
        data["downloadable_artifacts"] = [
            lnk for lnk in _existing_links
            if isinstance(lnk, dict) and isinstance(lnk.get("url"), str) and lnk.get("url")
        ]
    else:
        data["downloadable_artifacts"] = []

    # Derive canonical workflow_state — status takes precedence over raw_result field
    # so that needs_inputs always resolves to WAITING_FOR_INPUTS regardless of what
    # the inner plan_preview carried.
    _raw_workflow_state: str = (
        (isinstance(truncated_result, dict) and truncated_result.get("workflow_state")) or ""
    )
    _workflow_state: str = ""
    # Always re-derive from status first; only fall through to raw_result field if no override needed
    if not _workflow_state:
        from backend.workflow_checkpoint import WorkflowState as _WS
        # Order matters: needs_inputs takes precedence over approval_required
        if status == "needs_inputs":
            _workflow_state = _WS.WAITING_FOR_INPUTS.value
        elif approval_required and status == "workflow_planned":
            _workflow_state = _WS.WAITING_FOR_APPROVAL.value
        elif status == "workflow_planned":
            _workflow_state = _WS.PLANNING.value
        elif status in {"success", "completed"}:
            _workflow_state = _WS.COMPLETED.value
        elif status in {"error", "failed"}:
            _workflow_state = _WS.FAILED.value
        else:
            _workflow_state = _WS.IDLE.value

    response = {
        "version": "1.0",
        "success": success,
        "session_id": session_id,
        # Promote run_id to top-level so the frontend can read it without unwrapping envelopes
        "run_id": _run_id_hint or None,
        "prompt": prompt,
        "tool": tool,
        "status": status,
        "workflow_state": _workflow_state,
        "execute_ready": execute_ready,
        "approval_required": approval_required,
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
    """Reject uploads that exceed HELIX_MAX_UPLOAD_BYTES / HELIX_MAX_UPLOAD_MB (default 20 MB). No limit if set to 0."""
    if not files:
        return
    max_bytes = _get_max_upload_bytes()
    if max_bytes <= 0:
        return
    for idx, f in enumerate(files):
        size = f.get("size")
        if isinstance(size, int) and size > max_bytes:
            raise HTTPException(
                status_code=413,
                detail=f"Uploaded file #{idx+1} exceeds the allowed size limit ({max_bytes // (1024*1024)} MB)."
            )
        content = f.get("content")
        if isinstance(content, str) and len(content) > max_bytes:
            raise HTTPException(
                status_code=413,
                detail=f"Uploaded file #{idx+1} content exceeds the allowed size limit ({max_bytes // (1024*1024)} MB)."
            )


async def _persist_upload_to_session(
    session_id: str,
    upload: UploadFile,
    *,
    max_bytes: int,
) -> Dict[str, Any]:
    """Persist one uploaded file under session uploads/raw and return metadata."""
    original_name = upload.filename or "upload.dat"
    safe_name = _sanitize_uploaded_filename(original_name)
    if not _is_allowed_uploaded_filename(safe_name):
        raise HTTPException(
            status_code=415,
            detail=(
                f"Unsupported file type for '{original_name}'. "
                "Allowed: FASTA/FASTQ (.fa .fas .fasta .fastq .fq .gz), "
                "tabular (.csv .tsv .xlsx .xls), "
                "variants (.vcf .bcf), intervals (.bed .gff .gff3 .gtf), "
                "single-cell (.h5ad .loom .h5), alignment (.sam .bam .cram), "
                "plain text (.txt)."
            ),
        )

    upload_paths = history_manager.ensure_session_upload_directories(session_id)
    raw_dir = Path(upload_paths["raw"])
    dest_path = _make_unique_path(raw_dir / safe_name)

    first_chunk = await upload.read(UPLOAD_CHUNK_BYTES)
    intake_decision = evaluate_upload_intake(
        filename=safe_name,
        content_type=upload.content_type or "application/octet-stream",
        first_chunk=first_chunk,
    )
    if intake_decision.decision == "block":
        await upload.close()
        _emit_policy_audit_event(
            "upload_blocked",
            session_id,
            {
                "filename": original_name,
                "sanitized_filename": safe_name,
                "content_type": upload.content_type or "application/octet-stream",
                "intake_policy": intake_decision.to_dict(),
            },
        )
        raise HTTPException(
            status_code=400,
            detail=(
                f"Upload blocked by intake policy for '{original_name}'. "
                f"Reason: {intake_decision.reason}"
            ),
        )

    size = len(first_chunk or b"")
    try:
        with open(dest_path, "wb") as out_f:
            if first_chunk:
                if max_bytes > 0 and size > max_bytes:
                    out_f.close()
                    try:
                        dest_path.unlink()
                    except Exception:
                        pass
                    raise HTTPException(
                        status_code=413,
                        detail=(
                            f"File '{original_name}' exceeds the allowed upload size "
                            f"({max_bytes // (1024*1024)} MB)."
                        ),
                    )
                out_f.write(first_chunk)
            while True:
                chunk = await upload.read(UPLOAD_CHUNK_BYTES)
                if not chunk:
                    break
                size += len(chunk)
                if max_bytes > 0 and size > max_bytes:
                    out_f.close()
                    try:
                        dest_path.unlink()
                    except Exception:
                        pass
                    raise HTTPException(
                        status_code=413,
                        detail=(
                            f"File '{original_name}' exceeds the allowed upload size "
                            f"({max_bytes // (1024*1024)} MB)."
                        ),
                    )
                out_f.write(chunk)
    finally:
        await upload.close()

    # ── Upload-time file profiling ────────────────────────────────────────────
    # Run the format-appropriate profiler synchronously in a thread pool so we
    # don't block the event loop.  Profiling errors are non-fatal: we log them
    # and continue with an empty schema_preview.
    schema_preview: Dict[str, Any] = {}
    try:
        from backend.file_intelligence.profiler import profile_file, SUPPORTED_EXTENSIONS
        suffix = Path(safe_name).suffix.lower()
        if suffix in SUPPORTED_EXTENSIONS or suffix == ".gz":
            loop = asyncio.get_running_loop()
            schema_preview = await loop.run_in_executor(
                None, profile_file, str(dest_path.resolve())
            )
    except Exception as _prof_exc:
        logger.warning("upload profiling failed for %s: %s", safe_name, _prof_exc)

    file_id = f"upload_{uuid.uuid4().hex[:12]}"
    metadata = {
        "file_id": file_id,
        "name": original_name,
        "filename": dest_path.name,
        "original_filename": original_name,
        "size": size,
        "content_type": upload.content_type or "application/octet-stream",
        "local_path": str(dest_path.resolve()),
        "relative_path": str(dest_path.resolve().relative_to(history_manager.storage_dir.resolve())),
        "uploaded_at": datetime.now(timezone.utc).isoformat(),
        "s3_bucket": None,
        "s3_key": "",
        "intake_policy": intake_decision.to_dict(),
        "requires_policy_approval": bool(intake_decision.approval_required),
        "policy_state": "approval_required" if intake_decision.approval_required else "cleared",
        "schema_preview": schema_preview,
    }
    if not history_manager.add_uploaded_file(session_id, metadata):
        raise HTTPException(status_code=500, detail=f"Failed to persist uploaded file '{original_name}'")
    _emit_policy_audit_event(
        "upload_accepted",
        session_id,
        {
            "filename": original_name,
            "sanitized_filename": safe_name,
            "size": size,
            "intake_policy": intake_decision.to_dict(),
        },
    )
    return metadata

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
async def create_session(req: SessionRequest, request: Request):
    """Create a new session for tracking user interactions."""
    try:
        _check_and_increment_session_create_counter(request)
        session_id = history_manager.create_session(req.user_id)
        return {
            "success": True,
            "session_id": session_id,
            "message": "Session created successfully"
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/create_session")
async def create_session_alias(request: Request):
    """Alias for /session/create to support frontend compatibility."""
    try:
        _check_and_increment_session_create_counter(request)
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
        # Sanitize NaN/Inf that may be present in legacy session files written
        # before the serialize_for_json fix.  Without this, json.dumps raises
        # ValueError which propagates past the CORS middleware (no CORS headers).
        from backend.file_intelligence.tabular import _sanitize_json
        return _sanitize_json({
            "success": True,
            "session": session,
            "summary": summary,
        })
    except HTTPException:
        raise
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


@app.get("/session/{session_id}/artifact-aliases")
async def get_session_artifact_aliases(session_id: str):
    """Return semantic artifact aliases for evaluator/debugger use."""
    try:
        session = history_manager.get_session(session_id)
        if not session:
            raise HTTPException(status_code=404, detail="Session not found")
        aliases = history_manager.get_artifact_aliases(session_id)
        return {"success": True, "session_id": session_id, "aliases": aliases}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/session/{session_id}/lineage")
async def get_session_lineage(session_id: str):
    """Return artifact lineage edges for a session."""
    try:
        session = history_manager.get_session(session_id)
        if not session:
            raise HTTPException(status_code=404, detail="Session not found")
        edges = history_manager.get_lineage_edges(session_id)
        return {"success": True, "session_id": session_id, "edges": edges}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/session/{session_id}/historical-states")
async def get_session_historical_states(session_id: str):
    """Return semantically indexed historical artifact states."""
    try:
        session = history_manager.get_session(session_id)
        if not session:
            raise HTTPException(status_code=404, detail="Session not found")
        states = history_manager.get_historical_states(session_id)
        return {"success": True, "session_id": session_id, "states": states}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/benchmark/score/{session_id}")
async def score_benchmark_session(
    session_id: str,
    run_file: Optional[str] = None,
    benchmark_yaml: str = "tests/bioinformatics_benchmark.yaml",
):
    """
    Score a benchmark run using the standalone scorer.
    By default, uses the latest run file in benchmarks/runs.
    """
    try:
        from backend.orchestration.benchmark_evaluator import evaluate_benchmark_run_file

        if run_file:
            run_path = Path(run_file)
        else:
            runs_dir = Path("benchmarks/runs")
            candidates = sorted(runs_dir.glob("iteration_*_raw.json"), key=lambda p: p.stat().st_mtime, reverse=True)
            if not candidates:
                # Backward-compatible fallback in tests/
                candidates = sorted(Path("tests").glob("_latest_benchmark_turns_run_iter*.json"), key=lambda p: p.stat().st_mtime, reverse=True)
            if not candidates:
                raise HTTPException(status_code=404, detail="No benchmark run artifacts found")
            run_path = candidates[0]

        if not run_path.exists():
            raise HTTPException(status_code=404, detail=f"Run file not found: {run_path}")
        scored = evaluate_benchmark_run_file(run_path, Path(benchmark_yaml))
        if scored.get("session_id") and scored.get("session_id") != session_id:
            # Session-specific endpoint should not accidentally score another session silently.
            raise HTTPException(
                status_code=400,
                detail=f"Run file session_id={scored.get('session_id')} does not match requested {session_id}",
            )
        return {"success": True, "session_id": session_id, "run_file": str(run_path), "score": scored}
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
    When HELIX_MAX_UPLOAD_BYTES or HELIX_MAX_UPLOAD_MB is set (e.g. 10 MB for hosted),
    any file over that size is rejected so only small datasets are supported.
    """
    try:
        max_bytes = _get_max_upload_bytes()
        if max_bytes > 0:
            for file_info in req.files:
                if file_info.size > max_bytes:
                    limit_mb = max_bytes // (1024 * 1024)
                    raise HTTPException(
                        status_code=413,
                        detail=f"File '{file_info.filename}' exceeds the allowed upload size ({limit_mb} MB). This deployment supports small datasets only."
                    )
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


@app.post("/session/{session_id}/uploads")
async def upload_session_files(session_id: str, files: List[UploadFile] = File(...)):
    """Upload one or more files into session-local upload directories."""
    try:
        session = history_manager.get_session(session_id)
        if not session:
            raise HTTPException(status_code=404, detail=f"Session {session_id} not found")
        if not files:
            raise HTTPException(status_code=400, detail="No files provided")

        max_bytes = _get_max_upload_bytes()
        uploaded: List[Dict[str, Any]] = []
        for upload in files:
            file_info = await _persist_upload_to_session(session_id, upload, max_bytes=max_bytes)
            uploaded.append(file_info)

        profile = get_policy_profile()
        profile_defaults = policy_profile_defaults(profile)
        approval_required_count = sum(
            1 for item in uploaded if bool(item.get("requires_policy_approval"))
        )

        return {
            "success": True,
            "session_id": session_id,
            "files": uploaded,
            "uploaded_count": len(uploaded),
            "max_upload_mb": (max_bytes // (1024 * 1024)) if max_bytes > 0 else None,
            "policy_profile": profile,
            "policy_defaults": profile_defaults,
            "approval_required_uploads": approval_required_count,
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error uploading files to session {session_id}: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/session/{session_id}/uploads/approve")
async def approve_policy_uploads(session_id: str, request: Request):
    """Approve policy-gated uploaded files for downstream execution."""
    session = history_manager.get_session(session_id)
    if not session:
        raise HTTPException(status_code=404, detail=f"Session {session_id} not found")

    pending = _get_policy_pending_uploads(session_id)
    if not pending:
        return {
            "success": True,
            "session_id": session_id,
            "approved_count": 0,
            "pending_before": 0,
        }

    approved_count = _approve_policy_pending_uploads(session_id, "approved_via_upload_approval_endpoint")
    _emit_policy_audit_event(
        "upload_policy_approved",
        session_id,
        {
            "approved_count": approved_count,
            "pending_before": len(pending),
            "source_ip": request.client.host if request and request.client else "unknown",
        },
    )
    return {
        "success": True,
        "session_id": session_id,
        "approved_count": approved_count,
        "pending_before": len(pending),
    }

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

        # ── Load workflow checkpoint ──────────────────────────────────────────
        from backend.workflow_checkpoint import WorkflowCheckpoint, WorkflowState
        _checkpoint = history_manager.load_checkpoint(req.session_id)
        logger.info(
            "[Checkpoint] session=%s state=%s resume_node=%s",
            req.session_id, _checkpoint.state.value, _checkpoint.resume_node,
        )

        pending_policy_uploads = _get_policy_pending_uploads(req.session_id)
        if _is_approval_command(req.command) and pending_policy_uploads:
            approved_count = _approve_policy_pending_uploads(
                req.session_id, "approved_via_execute_approval_command"
            )
            _emit_policy_audit_event(
                "upload_policy_approved",
                req.session_id,
                {
                    "approved_count": approved_count,
                    "pending_before": len(pending_policy_uploads),
                    "via": "execute_approval_command",
                },
            )
            pending_policy_uploads = _get_policy_pending_uploads(req.session_id)

        if req.execute_plan and pending_policy_uploads:
            _emit_policy_audit_event(
                "execute_blocked_policy_pending_uploads",
                req.session_id,
                {
                    "pending_upload_count": len(pending_policy_uploads),
                    "command": req.command[:200],
                },
            )
            raise HTTPException(
                status_code=409,
                detail=(
                    "Execution blocked: one or more uploaded files require policy approval "
                    "before high-risk execution. Approve them via /session/{session_id}/uploads/approve "
                    "or issue an approval command first."
                ),
            )

        # ── Route by checkpoint state ─────────────────────────────────────────
        # Case D: user is approving a pending plan
        if _is_approval_command(req.command):
            # Prefer checkpoint-stored plan; fall back to legacy pending_plan for backward compat
            _cp_plan = (
                _checkpoint.pending_plan
                if _checkpoint.state == WorkflowState.WAITING_FOR_APPROVAL
                   and isinstance(_checkpoint.pending_plan, dict)
                else None
            )
            pending_plan = _cp_plan or _get_pending_plan(req.session_id)
            if isinstance(pending_plan, dict) and isinstance(pending_plan.get("plan"), dict):
                plan = _autobind_plan_inputs(pending_plan["plan"], session_context)
                binding_check = _validate_plan_bindings(plan)
                if binding_check:
                    # Still missing inputs after approval — remain in WAITING_FOR_INPUTS
                    _new_cp = WorkflowCheckpoint.waiting_for_inputs(
                        pending_plan=pending_plan,
                        missing_inputs=list(binding_check.get("missing_inputs", [])),
                    )
                    history_manager.save_checkpoint(req.session_id, _new_cp)
                    std = build_standard_response(
                        prompt=req.command,
                        tool="__plan__",
                        result={
                            "status": "needs_inputs",
                            "text": binding_check.get("message"),
                            "binding_diagnostics": binding_check,
                            "workflow_state": WorkflowState.WAITING_FOR_INPUTS.value,
                        },
                        session_id=req.session_id,
                        mcp_route="/execute",
                        success=True,
                        execution_path="approval_binding_validation_failed",
                    )
                    return CustomJSONResponse(std)

                # All inputs bound — transition to EXECUTING
                _exec_cp = WorkflowCheckpoint(state=WorkflowState.EXECUTING)
                history_manager.save_checkpoint(req.session_id, _exec_cp)

                broker = _get_execution_broker()
                result = await broker.execute_tool(
                    ExecutionRequest(
                        tool_name="__plan__",
                        arguments={"plan": plan, "session_id": req.session_id},
                        session_id=req.session_id,
                        original_command=pending_plan.get("command") or req.command,
                        session_context=session_context,
                    )
                )
                _apply_session_context_side_effects(req.session_id, "__plan__", result)
                storage_result = _build_pipeline_execution_storage_result(result, req.session_id)
                if _should_clear_pending_plan_after_execution(storage_result):
                    _clear_pending_plan(req.session_id)
                    history_manager.save_checkpoint(req.session_id, WorkflowCheckpoint.completed(
                        run_id=storage_result.get("run_id") if isinstance(storage_result, dict) else None
                    ))
                return await _dispatch_result(
                    req, "__plan__", storage_result,
                    tool_args={"plan": plan, "session_id": req.session_id},
                    execution_path="approval_execute_pending_plan",
                )

            no_pending = build_standard_response(
                prompt=req.command,
                tool="__plan__",
                result={
                    "status": "workflow_needs_clarification",
                    "text": (
                        "There is no pending workflow to approve. "
                        "Please ask for a workflow plan first, then approve it."
                    ),
                    "workflow_state": WorkflowState.IDLE.value,
                },
                session_id=req.session_id,
                mcp_route="/execute",
                success=True,
                execution_path="approval_no_pending_plan",
            )
            return CustomJSONResponse(no_pending)

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

        # Universal approval gate (action-based):
        # return a pending plan preview unless user explicitly asked to execute.
        if not req.execute_plan and not _is_approval_command(req.command):
            try:
                from backend.command_router import CommandRouter as _ApprovalRouter

                _approval_router = _ApprovalRouter()
                # Route with a timeout so slow LLM routing calls can't hang the endpoint.
                _gate_timeout_s = int(os.getenv("HELIX_GATE_ROUTE_TIMEOUT_S", "20"))
                try:
                    _approval_tool, _approval_params = await asyncio.wait_for(
                        asyncio.get_running_loop().run_in_executor(
                            None, lambda: _approval_router.route_command_with_shadow(req.command, session_context)
                        ),
                        timeout=_gate_timeout_s,
                    )
                except asyncio.TimeoutError:
                    logger.warning("Approval gate routing timed out after %ss — skipping gate", _gate_timeout_s)
                    _approval_tool, _approval_params = "handle_natural_command", {}
                _ready_hist_recreate = _historical_recreation_ready_for_execution(
                    req.command,
                    _approval_tool,
                    session_context,
                )
                if _should_stage_for_approval(_approval_tool, req.command, _approval_params) and not _ready_hist_recreate:
                    _approval_params = _approval_params or {}
                    _approval_params.setdefault("session_id", req.session_id)
                    pending_plan = _build_single_step_plan(req.command, _approval_tool, _approval_params)
                    _store_pending_plan(req.session_id, req.command, pending_plan)
                    binding_check = _validate_plan_bindings(pending_plan)
                    _missing = list(binding_check.get("missing_inputs", [])) if binding_check else []
                    # Always stage for WAITING_FOR_APPROVAL when intent is clear.
                    # Missing file bindings don't prevent planning — they just mean execution
                    # needs those files when the plan is eventually approved and run.
                    history_manager.save_checkpoint(
                        req.session_id,
                        WorkflowCheckpoint.waiting_for_approval(
                            pending_plan={"plan": pending_plan, "command": req.command},
                        ),
                    )
                    _workflow_state_for_preview = WorkflowState.WAITING_FOR_APPROVAL.value
                    # Status is always workflow_planned so the scorer (and UI) treat this as
                    # a plan awaiting approval, regardless of whether file inputs are bound.
                    preview_status = "workflow_planned"
                    plan_preview = {
                        "status": preview_status,
                        "type": "plan_result",
                        "plan_version": pending_plan.get("version", "v1"),
                        "steps": pending_plan.get("steps", []),
                        "execute_ready": True,
                        "approval_required": True,
                        "workflow_state": _workflow_state_for_preview,
                    }
                    if binding_check:
                        plan_preview["binding_diagnostics"] = binding_check
                    return await _dispatch_result(
                        req,
                        "__plan__",
                        plan_preview,
                        tool_args={"plan": pending_plan, "session_id": req.session_id},
                        execution_path="approval_gated_plan",
                    )
            except Exception as _approval_err:
                logger.warning("Approval pre-gate skipped due to error: %s", _approval_err)

        # Phase 2b: if the client explicitly asked to execute a previously planned pipeline,
        # re-route through the agent with an execute_plan flag so it dispatches async jobs
        # rather than returning another plan document.
        if req.execute_plan:
            from backend.agent import handle_command
            agent_timeout_s = int(os.getenv("HELIX_AGENT_TIMEOUT_S", "25"))
            agent_result, agent_diag = await _run_agent_with_retry(
                lambda: handle_command(
                    req.command,
                    session_id=req.session_id,
                    session_context=session_context,
                    execute_plan=True,
                ),
                timeout_s=agent_timeout_s,
                retries=int(os.getenv("HELIX_AGENT_RETRIES", "1")),
                backoff_s=float(os.getenv("HELIX_AGENT_BACKOFF_S", "0.5")),
            )
            if isinstance(agent_result, dict):
                agent_result.setdefault("diagnostics", agent_diag)
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
        # Also handle FastQC validation and MultiQC unsupported BEFORE agent to ensure
        # clear, immediate responses for these common follow-up requests.
        _phase2c_allowlist = {"toolbox_inventory", "session_run_io_summary", "visualize_job_results"}
        _phase2c_semantic_allowlist = {"bio_diff_runs", "bio_rerun", "patch_and_rerun"}
        _phase2c_validation_tools = {"fastqc_quality_analysis", "unsupported_tool"}
        try:
            import time as _t
            from backend.command_router import CommandRouter as _PreRouter
            _pre_router = _PreRouter()
            _pre_tool, _pre_params = _pre_router.route_command_with_shadow(req.command, session_context)
            logger.info(f"[Phase2c] route_command → tool='{_pre_tool}'")
            if _pre_tool in _phase2c_validation_tools:
                # FastQC with session_resolution_error or validation failure, or MultiQC unsupported
                _pre_params = _pre_params or {}
                _pre_params.setdefault("session_id", req.session_id)
                _pre_params.setdefault("original_command", req.command)
                _pre_params.setdefault("session_context", session_context)
                _pre_result = await dispatch_tool(_pre_tool, _pre_params)
                return await _dispatch_result(
                    req, _pre_tool, _pre_result, tool_args=_pre_params, execution_path="phase2c_validation"
                )
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
            if _pre_tool in _phase2c_semantic_allowlist:
                _action = infer_action(req.command, _pre_tool)
                if _action in {"compare_versions", "generate_plot", "rerun_downstream_steps", "subset_data"}:
                    _pre_start = _t.time()
                    if _pre_params is None:
                        _pre_params = {}
                    _pre_params.setdefault("session_id", req.session_id)
                    _pre_result = await dispatch_tool(_pre_tool, _pre_params)
                    logger.info(f"[Phase2c] '{_pre_tool}' semantic fast-path in {_t.time()-_pre_start:.2f}s")
                    return await _dispatch_result(
                        req, _pre_tool, _pre_result, tool_args=_pre_params, execution_path="phase2c_semantic_router"
                    )
        except Exception as _pre_err:
            logger.warning(f"[Phase2c] fast path raised exception ({_pre_err}), falling through to agent")

        # Phase 3: detect multi-step workflows and execute as a Plan IR (sync/async broker handles routing)
        def _looks_like_workflow(cmd: str) -> bool:
            c = (cmd or "").lower()
            # Explicit code-edit commands may contain newlines/code fences but are single actions.
            if any(tok in c for tok in ["apply code patch", "apply patch:", "replace script with", "replace the script with"]) or "```" in c:
                return False
            # Treat only explicit workflow delimiters as multi-step signals.
            # Newlines alone are common in data-rich single commands (e.g., key-value follow-ups).
            return any(tok in c for tok in [" and then ", " then ", "->", "→", ";"]) and len(c) > 20

        # Primary path: let BioAgent (with agent.md prompt) plan/execute
        use_agent = not _agent_disabled()

        try:
            if not use_agent:
                raise RuntimeError("Agent disabled (HELIX_AGENT_DISABLED=1 or HELIX_MOCK_MODE=1)")

            import time
            agent_start_time = time.time()
            # Lazy import: Import backend.agent only when needed, not at module import time.
            # This avoids loading heavy LLM dependencies (langgraph, langchain) during server startup,
            # keeping lightweight endpoints like /health and /tools/list fast and allowing the service
            # to work in sandbox/CI environments where LLM dependencies may not be installed.
            from backend.agent import handle_command
            agent_timeout_s = int(os.getenv("HELIX_AGENT_TIMEOUT_S", "25"))
            agent_result, agent_diag = await _run_agent_with_retry(
                lambda: handle_command(req.command, session_id=req.session_id, session_context=session_context),
                timeout_s=agent_timeout_s,
                retries=int(os.getenv("HELIX_AGENT_RETRIES", "1")),
                backoff_s=float(os.getenv("HELIX_AGENT_BACKOFF_S", "0.5")),
            )
            if isinstance(agent_result, dict):
                agent_result.setdefault("diagnostics", agent_diag)
            
            agent_done_time = time.time()
            agent_duration = agent_done_time - agent_start_time
            logger.info("Agent tool mapping completed in %.2fs", agent_duration)

            # Reliability fallback: if agent call timed out/failed, fall through to
            # deterministic router instead of returning a hard error immediately.
            if isinstance(agent_result, dict) and agent_result.get("error") in {"AGENT_TIMEOUT", "AGENT_ERROR"}:
                raise RuntimeError(f"Agent unavailable: {agent_result.get('error')}")

            # Check if agent returned a tool mapping (not full execution)
            if isinstance(agent_result, dict) and agent_result.get("status") == "tool_mapped":
                # Agent only did tool mapping - now execute via broker with same validation as router
                tool_name = agent_result.get("tool_name")
                parameters = agent_result.get("parameters", {})

                logger.debug("Agent mapped tool '%s', executing via router...", tool_name)

                # Apply same validation as deterministic router: unsupported tools and session-aware FastQC
                from backend.command_router import CommandRouter as _AgentCheckRouter
                _check_router = _AgentCheckRouter()
                _check_tool, _check_params = _check_router.route_command_with_shadow(req.command, session_context)
                if _check_tool == "unsupported_tool":
                    _check_params = _check_params or {}
                    _check_params.setdefault("session_id", req.session_id)
                    _pre_result = await dispatch_tool("unsupported_tool", _check_params)
                    return await _dispatch_result(
                        req, "unsupported_tool", _pre_result,
                        tool_args=_check_params, execution_path="agent_validation_unsupported",
                    )
                if _check_tool == "fastqc_quality_analysis" and (_check_params or {}).get("session_resolution_error"):
                    parameters["session_resolution_error"] = _check_params["session_resolution_error"]

                return await _execute_routed_tool(
                    req,
                    tool_name=tool_name,
                    parameters=parameters,
                    session_context=session_context,
                    execution_path="agent_tool_mapped",
                    validation_execution_path="agent_validation_message",
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
            try:
                intent = classify_intent(
                    req.command,
                    session_context=session_context,
                    workflow_state=_checkpoint.state.value,
                )
            except TypeError:
                # Backward-compatible path for tests/mocks or legacy classifier call shapes.
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
                            "This looks like a question. Enable the agent (HELIX_AGENT_DISABLED=0) "
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

            # If this looks like a workflow, execute Plan IR directly via broker.
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
                    req,
                    "__plan__",
                    storage_result,
                    tool_args={"plan": plan, "session_id": req.session_id},
                    execution_path="fallback_router_plan_execute",
                )

            tool_name, parameters = command_router.route_command(req.command, session_context)
            logger.debug("Routed '%s' → tool='%s'", req.command[:60], tool_name)

            if tool_name == "variant_selection" and "session_id" not in parameters:
                parameters["session_id"] = req.session_id
            if tool_name == "fastqc_quality_analysis" and "session_id" not in parameters:
                parameters["session_id"] = req.session_id

            return await _execute_routed_tool(
                req,
                tool_name=tool_name,
                parameters=parameters,
                session_context=session_context,
                execution_path="fallback_router",
                validation_execution_path="validation_message",
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
        # keeping lightweight endpoints like /health and /tools/list fast and allowing the service
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

@app.post("/tools/sequence-alignment")
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

@app.post("/tools/mutate-sequence")
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

@app.post("/tools/analyze-sequence-data")
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

@app.post("/tools/select-variants")
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

@app.post("/tools/parse-command")
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

@app.post("/tools/execute-command")
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

@app.post("/tools/handle-natural-command")
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

@app.post("/tools/visualize-alignment")
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

@app.post("/tools/plasmid-visualization")
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

@app.post("/tools/plasmid-for-representatives")
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

@app.post("/tools/read-trimming")
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

@app.post("/tools/read-merging")
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

@app.get("/tools/list")
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
    "go_enrichment_analysis": {
        "display_name": "GO/pathway enrichment analysis",
        "description": (
            "Runs Gene Ontology/pathway enrichment on a significant gene list "
            "from a prior DEG result state."
        ),
        "required_inputs": [
            {
                "name": "gene_list",
                "description": "List of gene symbols to test for enrichment.",
                "example": "STAT1, IRF7, CXCL10, ISG15",
            },
        ],
        "optional_inputs": [
            {
                "name": "source_selector",
                "description": "Historical selector for DEG source (for example 'current DEG results').",
                "example": "current DEG results",
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

    handled, registry_result = await dispatch_via_registry(
        tool_name,
        arguments,
        needs_inputs_builder=_build_needs_inputs_response,
    )
    if handled:
        return registry_result

    # ── Tabular QA (Code Interpreter) ─────────────────────────────────────────
    if tool_name == "tabular_qa":
        from backend.tabular_qa.agent import run_tabular_qa
        session_id = arguments.get("session_id", "")
        question = arguments.get("question") or arguments.get("command") or arguments.get("objective", "")
        file_path = arguments.get("file_path", "")
        sheet = arguments.get("sheet")
        profile = arguments.get("profile") or {}

        # If no explicit file_path, find first tabular file in session uploads
        if not file_path and session_id:
            session_data = history_manager.get_session(session_id) or {}
            for up in (session_data.get("uploaded_files") or []):
                name = (up.get("filename") or up.get("name") or "").lower()
                if any(name.endswith(ext) for ext in {".csv", ".tsv", ".xlsx", ".xls"}):
                    file_path = up.get("local_path") or ""
                    profile = up.get("schema_preview") or profile
                    break

        if not file_path:
            return {
                "status": "error",
                "text": "No tabular file found in session. Please upload a CSV, TSV, or Excel file first.",
            }

        qa_result = run_tabular_qa(
            question=question,
            session_id=session_id,
            file_path=file_path,
            profile=profile,
            sheet=sheet,
        )
        return {
            "status": "success" if qa_result["success"] else "error",
            "text": qa_result["answer"],
            "result": qa_result.get("result"),
            "code": qa_result.get("code"),
            "attempts": qa_result.get("attempts"),
        }

    # ── Tabular analysis (deterministic operations) ───────────────────────────
    if tool_name == "tabular_analysis":
        from backend.ds_pipeline.pipelines.ingest import ingest_tabular
        session_id = arguments.get("session_id", "")
        data_path = arguments.get("data_path", "")
        sheet = arguments.get("sheet")

        if not data_path and session_id:
            session_data = history_manager.get_session(session_id) or {}
            for up in (session_data.get("uploaded_files") or []):
                name = (up.get("filename") or up.get("name") or "").lower()
                if any(name.endswith(ext) for ext in {".csv", ".tsv", ".xlsx", ".xls"}):
                    data_path = up.get("local_path") or ""
                    break

        if not data_path:
            return {
                "status": "error",
                "text": "No tabular file found. Please upload a CSV, TSV, or Excel file first.",
            }

        try:
            _conn, df, meta = ingest_tabular(data_path, sheet=sheet)
            return {
                "status": "success",
                "text": (
                    f"Loaded {meta['n_rows']} rows × {meta['n_cols']} columns"
                    + (f" from sheet '{meta['source_sheet']}'" if meta.get('source_sheet') else "")
                    + ". Ready for analysis."
                ),
                "schema_preview": {
                    "n_rows": meta["n_rows"],
                    "n_cols": meta["n_cols"],
                    "columns": [c["name"] for c in meta["columns"]],
                    "source_sheet": meta.get("source_sheet"),
                    "available_sheets": meta.get("available_sheets"),
                    "sample": meta["sample"],
                },
                "objective": arguments.get("objective", ""),
                "operations": arguments.get("operations", []),
            }
        except Exception as exc:
            return {"status": "error", "text": str(exc)}

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
        from backend.orchestration.state_selection import select_rerun_anchor
        _orch_rerun = _BioOrch(tool_executor=dispatch_tool, history_manager=history_manager)
        _session_id_rr = arguments.get("session_id", "")
        _changes_rr = arguments.get("changes", {})
        _target_run = arguments.get("target_run", "latest")
        _orig_cmd_rr = str(arguments.get("original_command") or arguments.get("command") or "")
        _rerunnable_tools = {
            "bulk_rnaseq_analysis",
            "single_cell_analysis",
            "quality_assessment",
            "fastqc_quality_analysis",
        }

        def _session_run_by_id(_sid: str, _rid: str) -> Optional[Dict[str, Any]]:
            if not _sid or not _rid:
                return None
            _runs = history_manager.list_runs(_sid)
            for _r in _runs:
                if isinstance(_r, dict) and _r.get("run_id") == _rid:
                    return _r
            return None

        def _extract_rerunnable_step_from_plan_run(_run: Dict[str, Any]) -> Optional[Dict[str, Any]]:
            if not isinstance(_run, dict):
                return None
            if _run.get("tool") != "__plan__":
                return None
            _result = _run.get("result")
            if not isinstance(_result, dict):
                return None
            _steps = _result.get("steps")
            if not isinstance(_steps, list):
                return None
            for _step in reversed(_steps):
                if not isinstance(_step, dict):
                    continue
                _step_tool = _step.get("tool_name")
                _step_args = _step.get("arguments") if isinstance(_step.get("arguments"), dict) else {}
                if _step_tool in _rerunnable_tools and _step_args:
                    return {
                        "run_id": _run.get("run_id"),
                        "tool": _step_tool,
                        "tool_args": _step_args,
                        "source": "plan_step",
                    }
            return None

        def _resolve_session_rerun_candidate(_sid: str, _target: str) -> Optional[Dict[str, Any]]:
            if not _sid:
                return None
            _runs = history_manager.list_runs(_sid)
            _candidates: List[Dict[str, Any]] = []
            for _r in _runs:
                if not isinstance(_r, dict):
                    continue
                if (
                    _r.get("tool") in _rerunnable_tools
                    and isinstance(_r.get("tool_args"), dict)
                    and len(_r.get("tool_args") or {}) > 0
                ):
                    _candidates.append(_r)
                    continue
                _plan_candidate = _extract_rerunnable_step_from_plan_run(_r)
                if _plan_candidate:
                    _candidates.append(_plan_candidate)
            if not _candidates:
                return None
            if _target == "latest":
                return _candidates[-1]
            if _target == "prior":
                return _candidates[-2] if len(_candidates) > 1 else None
            return _session_run_by_id(_sid, _target)

        _session_obj_rr = history_manager.get_session(_session_id_rr) if _session_id_rr else {}
        if _target_run in {"latest", "prior"}:
            _target_run = select_rerun_anchor(_session_obj_rr or {}, _target_run, _orig_cmd_rr)
            _session_candidate = _resolve_session_rerun_candidate(_session_id_rr, _target_run)
            if isinstance(_session_candidate, dict) and _session_candidate.get("run_id"):
                _target_run = _session_candidate["run_id"]
            else:
            # Prefer session-local run ledger to avoid cross-session run leakage.
                _session_runs = history_manager.list_runs(_session_id_rr) if _session_id_rr else []
                _session_run_ids_raw = [
                    r.get("run_id") for r in _session_runs
                    if isinstance(r, dict) and isinstance(r.get("run_id"), str)
                ]
                from pathlib import Path as _Path
                _arts_root = _Path(__file__).parent.parent / "artifacts"
                _session_run_ids = [
                    rid for rid in _session_run_ids_raw
                    if (_arts_root / rid / "run.json").exists()
                ]
                if _target_run == "latest" and _session_run_ids:
                    _target_run = _session_run_ids[-1]
                elif _target_run == "prior" and len(_session_run_ids) >= 2:
                    _target_run = _session_run_ids[-2]
                else:
                    if _session_id_rr:
                        return {
                            "status": "needs_inputs",
                            "text": (
                                "Rerun target was not executable from current session history. "
                                "I found no prior run with reusable parameters or artifact-backed config. "
                                "Provide an explicit `run_id`, or rerun with concrete inputs "
                                "(for example `count_matrix` and `sample_metadata`)."
                            ),
                            "diagnostics": {
                                "issue": "missing_rerunnable_anchor",
                                "requested_target": arguments.get("target_run", "latest"),
                                "resolved_target": _target_run,
                                "available_session_runs": _session_run_ids_raw[-5:],
                            },
                        }
                    # Fallback for older runs that may only exist in artifacts/.
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
                    _idx = 1 if _target_run == "prior" else 0
                    _target_run = _run_dirs[_idx].parent.name if len(_run_dirs) > _idx else _run_dirs[0].parent.name

        # If target run is a session-ledger run (run_...) without a BioOrchestrator artifact,
        # rerun directly from its stored tool + tool_args.
        from pathlib import Path as _Path
        _art_path = _Path(__file__).parent.parent / "artifacts" / _target_run / "run.json"
        if isinstance(_target_run, str) and _target_run.startswith("run_") and not _art_path.exists():
            _candidate = _resolve_session_rerun_candidate(_session_id_rr, _target_run)
            if not _candidate and _target_run in {"latest", "prior"}:
                _candidate = _resolve_session_rerun_candidate(_session_id_rr, _target_run)
            if _candidate and _candidate.get("tool") in _rerunnable_tools:
                _base_args = dict(_candidate.get("tool_args") or {})
                _new_args = {**_base_args, **(_changes_rr or {})}
                result = await _orch_rerun.run(
                    tool_name=_candidate.get("tool"),
                    params=_new_args,
                    parent_run_id=None,
                    session_id=_session_id_rr or None,
                    objective=f"rerun:{_candidate.get('tool')}",
                )
            else:
                result = {
                    "status": "needs_inputs",
                    "message": (
                        "No deterministic rerun anchor with reusable parameters was found in session history. "
                        "Provide an explicit run reference (for example `first DEG results`) "
                        "or rerun with concrete parameters."
                    ),
                    "diagnostics": {
                        "issue": "missing_rerunnable_parameters",
                        "target_run": _target_run,
                    },
                }
        else:
            result = await _orch_rerun.rerun(_target_run, _changes_rr, session_id=_session_id_rr)
        _result_status = result.get("status", "success")
        _text_default = result.get("text", result.get("summary_text", "Re-run complete."))
        _message_default = result.get("message") or _text_default
        return {
            "status": _result_status,
            "visualization_type": "results_viewer",
            "text": _message_default if _result_status == "error" else _text_default,
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
        from backend.orchestration.state_selection import select_diff_anchors
        from pathlib import Path as _Path
        _orch_diff = _BioOrch(tool_executor=dispatch_tool)
        _run_a = arguments.get("run_id_a", "latest")
        _run_b = arguments.get("run_id_b", "prior")
        _orig_cmd = str(arguments.get("original_command") or arguments.get("command") or "")

        # Resolve "latest" / "prior" pseudo-IDs with session-local preference.
        _sid_diff = arguments.get("session_id", "") or ""
        _session_runs_diff = history_manager.list_runs(_sid_diff) if _sid_diff else []
        _resolved_ids_raw = [
            r.get("run_id") for r in _session_runs_diff
            if isinstance(r, dict) and isinstance(r.get("run_id"), str)
        ]
        _arts_root_diff = _Path(__file__).parent.parent / "artifacts"
        _resolved_ids = [
            rid for rid in _resolved_ids_raw
            if (_arts_root_diff / rid / "run.json").exists()
        ]

        _session_obj_diff = history_manager.get_session(_sid_diff) if _sid_diff else {}
        _run_a, _run_b, _selector_diag = select_diff_anchors(
            _session_obj_diff or {},
            str(_run_a or ""),
            str(_run_b or ""),
            _orig_cmd,
        )

        if (_run_a == "latest" or _run_b == "prior") and len(_resolved_ids) < 2:
            if _sid_diff:
                return {
                    "status": "success",
                    "text": (
                        "I could not resolve two comparable historical runs in this session yet. "
                        "Run at least two comparable analyses first (for example, initial and corrected runs)."
                    ),
                }
            _run_dirs_diff = sorted(
                _arts_root_diff.glob("*/run.json"),
                key=lambda p: p.stat().st_mtime if p.exists() else 0,
                reverse=True,
            )
            # Convert newest-first artifact sort to oldest->newest ordering.
            _resolved_ids = list(reversed([p.parent.name for p in _run_dirs_diff]))
        if _run_a == "latest":
            _run_a = _resolved_ids[-1] if _resolved_ids else ""
        if _run_b == "prior":
            _run_b = _resolved_ids[-2] if len(_resolved_ids) > 1 else ""

        if not _run_a or not _run_b:
            return {
                "status": "needs_inputs",
                "text": (
                    "Historical comparison requires two concrete analysis states. "
                    "Specify explicit references (for example `first DEG results` vs `current DEG results`) "
                    "or provide `run_id_a` and `run_id_b`."
                ),
                "result": {"selector_diagnostics": _selector_diag, "run_id_a": _run_a, "run_id_b": _run_b},
            }

        # Ensure run anchors exist in artifact-backed run store expected by diff engine.
        def _exists_artifact_run(_rid: str) -> bool:
            if not _rid:
                return False
            return (_arts_root_diff / _rid / "run.json").exists()

        def _session_run_by_id(_rid: str) -> Optional[Dict[str, Any]]:
            for _r in _session_runs_diff:
                if isinstance(_r, dict) and _r.get("run_id") == _rid:
                    return _r
            return None

        def _effective_run_payload(_run: Optional[Dict[str, Any]]) -> Dict[str, Any]:
            if not isinstance(_run, dict):
                return {"tool_name": None, "params": {}, "summary": ""}
            _tool_name = _run.get("tool")
            _params = _run.get("tool_args") if isinstance(_run.get("tool_args"), dict) else {}
            _summary = ""
            _result = _run.get("result")
            if isinstance(_result, dict):
                _summary = str(_result.get("text") or _result.get("summary_text") or _result.get("message") or "")
                if _run.get("tool") == "__plan__" and isinstance(_result.get("steps"), list):
                    for _step in reversed(_result.get("steps") or []):
                        if not isinstance(_step, dict):
                            continue
                        _step_tool = _step.get("tool_name")
                        _step_args = _step.get("arguments") if isinstance(_step.get("arguments"), dict) else {}
                        if _step_tool and _step_args:
                            _tool_name = _step_tool
                            _params = _step_args
                            _step_result = _step.get("result")
                            if isinstance(_step_result, dict):
                                _summary = str(
                                    _step_result.get("text")
                                    or _step_result.get("summary_text")
                                    or _step_result.get("message")
                                    or _summary
                                )
                            break
            return {"tool_name": _tool_name, "params": _params, "summary": _summary}

        if not _exists_artifact_run(_run_a) and _resolved_ids:
            _run_a = _resolved_ids[-1]
        if not _exists_artifact_run(_run_b):
            if len(_resolved_ids) > 1:
                _run_b = _resolved_ids[-2]
            elif _resolved_ids:
                _run_b = _resolved_ids[-1]

        if not _exists_artifact_run(_run_a) or not _exists_artifact_run(_run_b):
            _run_a_session = _session_run_by_id(_run_a)
            _run_b_session = _session_run_by_id(_run_b)
            _payload_a = _effective_run_payload(_run_a_session)
            _payload_b = _effective_run_payload(_run_b_session)
            if _payload_a.get("tool_name") and _payload_b.get("tool_name"):
                _params_a = _payload_a.get("params") if isinstance(_payload_a.get("params"), dict) else {}
                _params_b = _payload_b.get("params") if isinstance(_payload_b.get("params"), dict) else {}
                _all_keys = sorted(set(_params_a.keys()) | set(_params_b.keys()))
                _param_changes = {
                    _k: {"run_a": _params_a.get(_k), "run_b": _params_b.get(_k)}
                    for _k in _all_keys
                    if _params_a.get(_k) != _params_b.get(_k)
                }
                _change_lines = "\n".join(
                    f"- **{k}**: `{v['run_a']}` -> `{v['run_b']}`"
                    for k, v in _param_changes.items()
                ) or "- No parameter changes detected from session-ledger params."
                _text = (
                    "## Run Comparison (Session Ledger)\n\n"
                    f"**Run A:** `{_run_a[:8]}...` (`{_payload_a.get('tool_name')}`)  \n"
                    f"**Run B:** `{_run_b[:8]}...` (`{_payload_b.get('tool_name')}`)\n\n"
                    f"### Parameter Changes\n{_change_lines}\n\n"
                    "### Notes\n"
                    "Artifact-backed run manifests were unavailable, so this comparison uses recorded session parameters/results."
                )
                return {
                    "status": "success",
                    "visualization_type": "text",
                    "text": _text,
                    "result": {
                        "status": "success",
                        "run_id_a": _run_a,
                        "run_id_b": _run_b,
                        "tool_name": _payload_a.get("tool_name"),
                        "param_changes": _param_changes,
                        "narrative_a": _payload_a.get("summary", "")[:300],
                        "narrative_b": _payload_b.get("summary", "")[:300],
                        "selector_diagnostics": _selector_diag,
                        "source": "session_ledger_fallback",
                    },
                }
            return {
                "status": "needs_inputs",
                "text": (
                    "Comparison anchors were identified, but neither state has enough persisted metadata to compute a reliable diff yet. "
                    "Re-run both target analyses or provide explicit run IDs that include stored parameters."
                ),
                "result": {"selector_diagnostics": _selector_diag, "run_id_a": _run_a, "run_id_b": _run_b},
            }

        diff = await _orch_diff.diff_runs(_run_a, _run_b)
        if diff.get("status") == "error":
            return {
                "status": "needs_inputs",
                "text": (
                    "Comparison anchors were selected, but persisted run manifests are incomplete for strict artifact diffing. "
                    "Please rerun the compared analyses so both states include saved run manifests."
                ),
                "result": {"selector_diagnostics": _selector_diag, "run_id_a": _run_a, "run_id_b": _run_b, "diff_error": diff.get("message")},
            }

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
            "result": {**diff, "selector_diagnostics": _selector_diag},
        }

    elif tool_name == "handle_natural_command":
        # Use the BioAgent path (system prompt from agent.md) for natural commands
        # Lazy import: Import backend.agent only when needed, not at module import time.
        # This avoids loading heavy LLM dependencies (langgraph, langchain) during server startup,
        # keeping lightweight endpoints like /health and /tools/list fast and allowing the service
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
        # Support both single-accession and multi-accession fetches.
        from tools.ncbi_tools import fetch_sequence_from_ncbi

        database = arguments.get("database", "nucleotide")
        single_accession = arguments.get("accession")
        multi_accessions = arguments.get("accessions")
        targets = []
        if isinstance(multi_accessions, list):
            targets = [str(a).strip() for a in multi_accessions if str(a).strip()]
        elif single_accession:
            targets = [str(single_accession).strip()]

        if not targets:
            return {
                "status": "error",
                "result": {},
                "text": "Error fetching sequence: no accession(s) provided",
                "error": "no accession(s) provided",
            }

        fetched = []
        failed = []
        for accession in targets:
            row = None
            tried_dbs = [database]
            # Fallback: if accession isn't available in requested DB, try the other one.
            if database == "protein":
                tried_dbs.append("nucleotide")
            elif database == "nucleotide":
                tried_dbs.append("protein")
            for db_choice in tried_dbs:
                try:
                    row = fetch_sequence_from_ncbi(accession=accession, database=db_choice)
                except Exception as e:
                    row = {"status": "error", "accession": accession, "error": str(e), "database": db_choice}
                if isinstance(row, dict) and row.get("status") == "success" and row.get("sequence"):
                    row["database"] = db_choice
                    break
            if isinstance(row, dict) and row.get("status") == "success" and row.get("sequence"):
                fetched.append(row)
            else:
                failed.append(
                    {
                        "accession": accession,
                        "error": (row.get("error") if isinstance(row, dict) else str(row)) or "Unknown error",
                    }
                )

        if len(targets) == 1:
            # Preserve legacy response shape for single fetch.
            if fetched:
                tool_out = fetched[0]
                return {
                    "status": "success",
                    "accession": tool_out.get("accession"),
                    "result": tool_out,
                    "text": f"Fetched sequence {tool_out.get('accession')} ({tool_out.get('length', 0)} aa/bp)",
                }
            err = failed[0].get("error") if failed else "Unknown error"
            return {
                "status": "error",
                "result": {},
                "text": f"Error fetching sequence: {err}",
                "error": err,
            }

        # Keep FASTA headers minimal (accession only) so downstream parsers that
        # expect strict FASTA tokenization do not treat description text as bases.
        combined_fasta = "\n".join(
            f">{row.get('accession')}\n{row.get('sequence', '')}".strip()
            for row in fetched
        ).strip()
        summary_text = (
            f"Fetched {len(fetched)}/{len(targets)} NCBI sequences."
            + (f" Failed: {len(failed)}." if failed else "")
        )
        return {
            "status": "success" if fetched else "error",
            "accessions": [row.get("accession") for row in fetched],
            "result": {
                "database": database,
                "fetched": fetched,
                "failed": failed,
                "combined_fasta": combined_fasta,
                "count_requested": len(targets),
                "count_fetched": len(fetched),
            },
            "combined_fasta": combined_fasta,
            "text": summary_text,
            "error": None if fetched else "; ".join(f.get("error", "Unknown error") for f in failed),
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

    elif tool_name == "go_enrichment_analysis":
        from backend.orchestration.artifact_resolver import resolve_semantic_reference

        _sid = arguments.get("session_id", "") or ""
        _gene_list = arguments.get("gene_list")
        _genes = _gene_list if isinstance(_gene_list, list) else []
        _source_selector = str(arguments.get("source_selector") or "current DEG results")
        _source_target = None
        if _sid:
            _session = history_manager.get_session(_sid) or {}
            _resolved = resolve_semantic_reference(_session, _source_selector)
            if _resolved.get("status") == "resolved":
                _source_target = _resolved.get("target")

        if not _genes:
            _base = _build_needs_inputs_response("go_enrichment_analysis", arguments)
            _base["status"] = "needs_inputs"
            _base["text"] = (
                "GO/pathway enrichment was recognized, but Helix does not yet have an extracted actionable gene list for this request. "
                "Provide `gene_list` explicitly, or ask Helix to export significant upregulated genes first and then run enrichment."
            )
            _base["diagnostics"] = {
                "intent": "run_enrichment",
                "source_selector": _source_selector,
                "resolved_source_target": _source_target,
                "missing": ["gene_list", "enrichment_executor"],
            }
            return _base

        return {
            "status": "needs_inputs",
            "text": (
                "GO/pathway enrichment intent is confirmed, but a dedicated enrichment executor is not configured in this runtime. "
                "Gene list was captured successfully; enable/install an enrichment backend (for example g:Profiler/clusterProfiler) to execute."
            ),
            "diagnostics": {
                "intent": "run_enrichment",
                "source_selector": _source_selector,
                "resolved_source_target": _source_target,
                "captured_gene_count": len(_genes),
                "missing_executor": "go_enrichment_backend",
            },
            "result": {
                "gene_list": _genes,
                "source_target": _source_target,
            },
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
        text = header + rows_md + top_md

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
        # Session-aware: check for resolution error (e.g. user asked for merged FASTA)
        session_resolution_error = arguments.get("session_resolution_error")
        if session_resolution_error:
            return {
                "status": "success",
                "result": {},
                "text": session_resolution_error,
                "validation_message": True,
            }

        input_r1 = arguments.get("input_r1", "")
        input_r2 = arguments.get("input_r2", "")
        # Enrich from session when agent path provided incomplete params
        if (not input_r1 or not input_r2) and arguments.get("session_context"):
            from backend.session_param_extractor import get_fastqc_inputs_from_session
            sess_r1, sess_r2, sess_out, sess_err = get_fastqc_inputs_from_session(
                arguments["session_context"],
                arguments.get("original_command") or arguments.get("command") or "",
            )
            if sess_err:
                return {
                    "status": "success",
                    "result": {},
                    "text": sess_err,
                    "validation_message": True,
                }
            if sess_r1 and sess_r2:
                input_r1 = input_r1 or sess_r1
                input_r2 = input_r2 or sess_r2
                if sess_out and not arguments.get("output"):
                    arguments["output"] = sess_out

        # Format and structure validation with clear error messages
        from backend.input_validation import validate_fastqc_inputs
        output = arguments.get("output")
        _from_broker = arguments.get("_from_broker", False)
        session_id = arguments.get("session_id")

        validation = validate_fastqc_inputs(input_r1, input_r2, arguments.get("original_command"))
        if not validation.valid:
            msg = validation.message or "Invalid FastQC inputs."
            if validation.suggestion:
                msg += "\n\n**What to do:**\n" + validation.suggestion
            return {
                "status": "success",
                "result": {},
                "text": msg,
                "validation_message": True,
            }

        # Delegate to the proper fastqc_quality_analysis function in agent_tools.py
        from backend.agent_tools import fastqc_quality_analysis as fastqc_tool
        
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

    elif tool_name == "unsupported_tool":
        from backend.unsupported_tools import get_unsupported_response
        requested = (arguments.get("requested_tool") or "").strip().lower().replace(" ", "_")
        unsupported = get_unsupported_response(requested) if requested else None
        if unsupported:
            text = f"**{requested.replace('_', ' ').title()} is not supported**\n\n"
            text += unsupported.get("reason", "") + "\n\n**Alternatives:**\n" + unsupported.get("alternatives", "")
            return {
                "status": "success",
                "result": {},
                "text": text,
                "tool_name": "unsupported_tool",
                "requested_tool": requested,
                "validation_message": True,
            }
        return {
            "status": "success",
            "result": {},
            "text": "This tool is not supported. Please try a different request.",
            "tool_name": "unsupported_tool",
        }

    else:
        # Check for known unsupported tools before tool-generator-agent
        from backend.unsupported_tools import get_unsupported_response
        unsupported = get_unsupported_response(tool_name)
        if unsupported:
            text = f"**{tool_name.replace('_', ' ').title()} is not supported**\n\n"
            text += unsupported.get("reason", "") + "\n\n**Alternatives:**\n" + unsupported.get("alternatives", "")
            return {
                "status": "success",
                "tool_name": tool_name,
                "text": text,
                "requested_tool": tool_name,
                "alternatives_provided": True,
            }

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
    return {
        "status": "healthy",
        "service": "Helix.AI Bioinformatics API",
        "agent_disabled": _agent_disabled(),
        "mock_mode": _agent_disabled(),  # backward compatibility
    }

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
    # Run as package module (start.sh uses `python -m backend.main`).
    # Use reload_dirs to restrict watching to Python source directories only.
    # Without this, the default behaviour watches the entire CWD, so every
    # session JSON write triggers a worker reload, dropping in-flight browser
    # requests (manifests as 500 + missing CORS headers in Chrome).
    uvicorn.run(
        "backend.main:app",
        host="0.0.0.0",
        port=8001,
        reload=True,
        reload_dirs=["backend", "tools"],
    )
