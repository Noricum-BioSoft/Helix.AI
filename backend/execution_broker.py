from __future__ import annotations

import asyncio
import os
import re
import logging
from dataclasses import dataclass
from typing import Any, Awaitable, Callable, Dict, List, Optional, Tuple
from pathlib import Path

from backend.plan_ir import Plan

logger = logging.getLogger(__name__)


ToolExecutor = Callable[[str, Dict[str, Any]], Awaitable[Dict[str, Any]]]

S3_URI_RE = re.compile(r"^s3://([^/]+)/(.+)$")


@dataclass
class ExecutionRequest:
    """
    Broker input describing what the user wants to run.

    This is intentionally minimal for Phase 1. Later phases will add:
    - plan steps (Plan IR)
    - discovered input assets + estimated bytes
    - routing decision (sync vs async) + infra target (EC2 vs EMR)
    """

    tool_name: str
    arguments: Dict[str, Any]
    session_id: str
    original_command: str
    session_context: Optional[Dict[str, Any]] = None


@dataclass
class InputAsset:
    """
    A discovered input asset (S3 object or local path) that may contribute to routing.
    """

    uri: str
    size_bytes: Optional[int]
    source: str  # e.g. "args", "session:uploaded", "session:dataset_reference"


@dataclass
class RoutingDecision:
    mode: str  # "sync" | "async"
    reason: str
    threshold_bytes: int
    estimated_bytes: int
    unknown_inputs: int
    override: Optional[str] = None


class ExecutionBroker:
    """
    Central entrypoint for executing a tool/workflow.

    Phase 1 goal:
    - Provide a single call site for tool execution so we can later attach:
      - routing policy (sync vs async)
      - input discovery + size estimation
      - generalized job submission and tracking

    For now (Phase 1), this is a thin wrapper over the existing tool executor.
    """

    # Small curated override list (Phase 1.2.B)
    FORCE_ASYNC_TOOLS = {"fastqc_quality_analysis"}
    FORCE_SYNC_TOOLS: set[str] = set()

    def __init__(self, tool_executor: ToolExecutor):
        self._tool_executor = tool_executor

    async def execute_tool(self, req: ExecutionRequest) -> Dict[str, Any]:
        # Phase 3: Plan execution path (multi-step workflow)
        if req.tool_name == "__plan__" or (isinstance(req.arguments, dict) and "plan" in req.arguments):
            plan_dict = req.arguments.get("plan") if isinstance(req.arguments, dict) else None
            plan = Plan.parse_obj(plan_dict)
            inputs = self._discover_inputs(plan.dict(), req.session_context or {})
            estimated_bytes, unknown = self._estimate_total_bytes(inputs)
            decision = self._evaluate_routing_policy(
                tool_name="__plan__",
                estimated_bytes=estimated_bytes,
                unknown_inputs=unknown,
            )
            if decision.mode == "async":
                output = await self._submit_plan_emr_job(req, plan)
                return self._wrap_result(
                    tool_name="__plan__",
                    mode="async",
                    decision=decision,
                    inputs=inputs,
                    output=output,
                )

            output = await self._execute_plan_sync(plan)
            return self._wrap_result(
                tool_name="__plan__",
                mode="sync",
                decision=decision,
                inputs=inputs,
                output=output,
            )

        inputs = self._discover_inputs(req.arguments or {}, req.session_context or {})
        estimated_bytes, unknown = self._estimate_total_bytes(inputs)
        decision = self._evaluate_routing_policy(
            tool_name=req.tool_name,
            estimated_bytes=estimated_bytes,
            unknown_inputs=unknown,
        )

        if decision.mode == "async" and req.tool_name == "fastqc_quality_analysis":
            result = await self._submit_fastqc_job(req)
            return self._wrap_result(
                tool_name=req.tool_name,
                mode="async",
                decision=decision,
                inputs=inputs,
                output=result,
            )

        if decision.mode == "async":
            result = await self._submit_universal_emr_job(req)
            return self._wrap_result(
                tool_name=req.tool_name,
                mode="async",
                decision=decision,
                inputs=inputs,
                output=result,
            )

        # Default: sync execution via existing tool executor (local/EC2/tool-generator).
        # Mark that we're calling from the broker to prevent infinite loops
        tool_args = dict(req.arguments) if req.arguments else {}
        tool_args["_from_broker"] = True
        output = await self._tool_executor(req.tool_name, tool_args)
        return self._wrap_result(
            tool_name=req.tool_name,
            mode="sync",
            decision=decision,
            inputs=inputs,
            output=output,
        )

    def _wrap_result(
        self,
        *,
        tool_name: str,
        mode: str,
        decision: RoutingDecision,
        inputs: List[InputAsset],
        output: Any,
    ) -> Dict[str, Any]:
        """
        Produce a consistent broker envelope while preserving original tool output under `result`.
        """
        status = "success"
        text = ""
        if isinstance(output, dict):
            status = output.get("status", "success") or "success"
            text = output.get("text") or output.get("message") or ""

        return {
            "status": status,
            "type": "execution_result",
            "mode": mode,
            "tool_name": tool_name,
            "text": text,
            # Preserve original tool output for frontend + compatibility.
            "result": output if isinstance(output, dict) else {"value": output},
            "artifacts": (output.get("artifacts") if isinstance(output, dict) else None) or [],
            "routing": {
                "mode": decision.mode,
                "reason": decision.reason,
                "threshold_bytes": decision.threshold_bytes,
                "estimated_bytes": decision.estimated_bytes,
                "unknown_inputs": decision.unknown_inputs,
                "override": decision.override,
            },
            "inputs": [
                {"uri": i.uri, "size_bytes": i.size_bytes, "source": i.source}
                for i in inputs
            ],
        }

    def _discover_inputs(self, arguments: Dict[str, Any], session_context: Dict[str, Any]) -> List[InputAsset]:
        """
        Phase 1.3: Discover inputs from tool arguments + session uploads/dataset references.

        Limitations:
        - Inline sequences / text blobs have unknown "bytes" cost (treated as unknown).
        - If boto3/AWS credentials are unavailable, S3 object sizes may be unknown.
        """
        assets: List[InputAsset] = []

        # 1) From arguments (recursive scan)
        for s in self._iter_strings(arguments):
            if isinstance(s, str) and s.startswith("s3://"):
                assets.append(InputAsset(uri=s, size_bytes=self._try_get_s3_size(s), source="args"))
            elif self._looks_like_local_path(s):
                assets.append(InputAsset(uri=s, size_bytes=self._try_get_local_size(s), source="args"))

        # 2) From session uploads + dataset references
        metadata = (session_context or {}).get("metadata", {}) if isinstance(session_context, dict) else {}

        for f in metadata.get("uploaded_files", []) or []:
            uri = self._session_file_to_s3_uri(f, session_context)
            if uri:
                assets.append(
                    InputAsset(
                        uri=uri,
                        size_bytes=self._safe_int(f.get("size")),
                        source="session:uploaded",
                    )
                )

        for ref in metadata.get("dataset_references", []) or []:
            uri = self._session_file_to_s3_uri(ref, session_context)
            if uri:
                assets.append(
                    InputAsset(
                        uri=uri,
                        size_bytes=self._safe_int(ref.get("size")),
                        source="session:dataset_reference",
                    )
                )

        # De-dupe by uri (keep first non-null size if available)
        dedup: Dict[str, InputAsset] = {}
        for a in assets:
            existing = dedup.get(a.uri)
            if existing is None:
                dedup[a.uri] = a
            else:
                if existing.size_bytes is None and a.size_bytes is not None:
                    dedup[a.uri] = a
        return list(dedup.values())

    def _estimate_total_bytes(self, inputs: List[InputAsset]) -> Tuple[int, int]:
        total = 0
        unknown = 0
        for i in inputs:
            if isinstance(i.size_bytes, int) and i.size_bytes >= 0:
                total += i.size_bytes
            else:
                unknown += 1
        return total, unknown

    def _evaluate_routing_policy(self, *, tool_name: str, estimated_bytes: int, unknown_inputs: int) -> RoutingDecision:
        """
        Phase 1.2: Routing policy v1 (sync vs async).

        Evaluation order (per checklist):
        (A) bytes threshold (default 100MB, env-configurable)
        (B) tool overrides (curated list)
        (C) timeout promotion hook (stub)
        """
        threshold = self._get_async_bytes_threshold()
        mode = "sync"
        reason = "below_threshold"
        override = None

        # (A) bytes threshold
        if estimated_bytes > threshold:
            mode = "async"
            reason = "bytes_threshold_exceeded"

        # (B) tool overrides
        if tool_name in self.FORCE_ASYNC_TOOLS:
            mode = "async"
            override = "force_async"
            reason = "tool_override_force_async"
        elif tool_name in self.FORCE_SYNC_TOOLS:
            mode = "sync"
            override = "force_sync"
            reason = "tool_override_force_sync"

        # (C) timeout promotion hook (stub)
        mode, reason = self._timeout_promotion_hook(mode=mode, reason=reason)

        return RoutingDecision(
            mode=mode,
            reason=reason,
            threshold_bytes=threshold,
            estimated_bytes=estimated_bytes,
            unknown_inputs=unknown_inputs,
            override=override,
        )

    def _timeout_promotion_hook(self, *, mode: str, reason: str) -> Tuple[str, str]:
        """
        Phase 1.2.C: stub for promoting long-running sync tasks to async later.
        """
        return mode, reason

    def _get_async_bytes_threshold(self) -> int:
        default = 100 * 1024 * 1024  # 100MB
        raw = os.getenv("HELIX_ASYNC_BYTES_THRESHOLD")
        if not raw:
            return default
        try:
            v = int(raw)
            return v if v > 0 else default
        except Exception:
            return default

    async def _submit_fastqc_job(self, req: ExecutionRequest) -> Dict[str, Any]:
        """
        Async (EMR) FastQC submission. Returns immediately with job_id.
        """
        from job_manager import get_job_manager

        args = req.arguments or {}
        input_r1 = args.get("input_r1") or args.get("r1_path") or ""
        input_r2 = args.get("input_r2") or args.get("r2_path") or ""
        output = args.get("output") or args.get("output_path")
        session_id = args.get("session_id") or req.session_id

        if not input_r1 or not input_r2:
            return {
                "status": "error",
                "message": "Both input_r1 and input_r2 are required for FastQC analysis",
            }

        jm = get_job_manager()
        if hasattr(asyncio, "to_thread"):
            job_id = await asyncio.to_thread(
                jm.submit_fastqc_job,
                r1_path=input_r1,
                r2_path=input_r2,
                output_path=output,
                session_id=session_id,
            )
        else:
            loop = asyncio.get_event_loop()
            job_id = await loop.run_in_executor(
                None,
                lambda: jm.submit_fastqc_job(
                    r1_path=input_r1,
                    r2_path=input_r2,
                    output_path=output,
                    session_id=session_id,
                ),
            )

        return {
            "type": "job",
            "status": "submitted",
            "job_id": job_id,
            "message": "FastQC job submitted. Processing will take 10-30 minutes.",
            "input_r1": input_r1,
            "input_r2": input_r2,
            "output": output,
        }

    async def _submit_universal_emr_job(self, req: ExecutionRequest) -> Dict[str, Any]:
        """
        Phase 2: generic EMR runner submission for non-FastQC tools.
        """
        from job_manager import get_job_manager

        jm = get_job_manager()
        tool_args = req.arguments or {}
        session_id = tool_args.get("session_id") or req.session_id

        if hasattr(asyncio, "to_thread"):
            job_id = await asyncio.to_thread(
                jm.submit_universal_emr_job,
                req.tool_name,
                tool_args,
                session_id,
            )
        else:
            loop = asyncio.get_event_loop()
            job_id = await loop.run_in_executor(
                None,
                lambda: jm.submit_universal_emr_job(req.tool_name, tool_args, session_id),
            )

        return {
            "type": "job",
            "status": "submitted",
            "job_id": job_id,
            "message": f"EMR job submitted for tool '{req.tool_name}'.",
        }

    async def _submit_plan_emr_job(self, req: ExecutionRequest, plan: Plan) -> Dict[str, Any]:
        from job_manager import get_job_manager

        jm = get_job_manager()
        session_id = (req.arguments or {}).get("session_id") or req.session_id

        if hasattr(asyncio, "to_thread"):
            job_id = await asyncio.to_thread(jm.submit_plan_emr_job, plan.dict(), session_id)
        else:
            loop = asyncio.get_event_loop()
            job_id = await loop.run_in_executor(None, lambda: jm.submit_plan_emr_job(plan.dict(), session_id))

        return {
            "type": "job",
            "status": "submitted",
            "job_id": job_id,
            "message": "EMR plan job submitted.",
        }

    async def _execute_plan_sync(self, plan: Plan) -> Dict[str, Any]:
        """
        Execute plan steps sequentially (sync mode).
        Supports minimal references using {"$ref": "steps.<id>.result.<path>"} in arguments.
        """
        context: Dict[str, Any] = {"steps": {}}
        step_outputs: List[Dict[str, Any]] = []

        for step in plan.steps:
            resolved_args = self._resolve_refs(step.arguments, context)
            out = await self._tool_executor(step.tool_name, resolved_args)
            step_record = {
                "id": step.id,
                "tool_name": step.tool_name,
                "arguments": resolved_args,
                "result": out,
            }
            context["steps"][step.id] = step_record
            step_outputs.append(step_record)

        return {
            "status": "success",
            "type": "plan_result",
            "plan_version": plan.version,
            "steps": step_outputs,
            "result": step_outputs[-1]["result"] if step_outputs else {},
        }

    def _resolve_refs(self, obj: Any, context: Dict[str, Any]) -> Any:
        if isinstance(obj, dict) and "$ref" in obj and isinstance(obj["$ref"], str):
            return self._get_ref_value(obj["$ref"], context)
        if isinstance(obj, dict):
            return {k: self._resolve_refs(v, context) for k, v in obj.items()}
        if isinstance(obj, list):
            return [self._resolve_refs(v, context) for v in obj]
        return obj

    def _get_ref_value(self, ref: str, context: Dict[str, Any]) -> Any:
        # ref format: steps.<step_id>.result.<path...>
        parts = ref.split(".")
        cur: Any = context
        for p in parts:
            if isinstance(cur, dict):
                cur = cur.get(p)
            elif isinstance(cur, list) and p.isdigit():
                cur = cur[int(p)]
            else:
                return None
        return cur

    def _iter_strings(self, obj: Any) -> List[str]:
        out: List[str] = []
        if isinstance(obj, str):
            out.append(obj)
        elif isinstance(obj, dict):
            for v in obj.values():
                out.extend(self._iter_strings(v))
        elif isinstance(obj, list) or isinstance(obj, tuple):
            for v in obj:
                out.extend(self._iter_strings(v))
        return out

    def _looks_like_local_path(self, s: str) -> bool:
        return isinstance(s, str) and (s.startswith("/") or s.startswith("./") or s.startswith("../"))

    def _try_get_local_size(self, p: str) -> Optional[int]:
        try:
            return int(Path(p).expanduser().resolve().stat().st_size)
        except Exception:
            return None

    def _try_get_s3_size(self, uri: str) -> Optional[int]:
        """
        Attempt to get S3 object size by calling head_object.
        
        Returns None if:
        - URI doesn't match S3 pattern
        - boto3 is not available
        - AWS credentials are not configured
        - IAM role lacks permissions
        - Object doesn't exist
        - Network/other errors occur
        
        Note: This method uses the default region. For cross-region buckets,
        the operation may fail silently. Consider enhancing with region detection.
        """
        m = S3_URI_RE.match(uri or "")
        if not m:
            return None
        bucket, key = m.group(1), m.group(2)
        try:
            import boto3  # type: ignore
            from botocore.exceptions import ClientError
            
            # Try default region first (works for most cases)
            default_region = os.getenv('AWS_REGION', os.getenv('AWS_DEFAULT_REGION', 'us-east-1'))
            client = boto3.client("s3", region_name=default_region)
            
            try:
                resp = client.head_object(Bucket=bucket, Key=key)
                size = int(resp.get("ContentLength")) if resp and "ContentLength" in resp else None
                if size is not None:
                    logger.debug(f"Retrieved S3 object size for {uri}: {size} bytes")
                return size
            except ClientError as e:
                error_code = e.response.get('Error', {}).get('Code', '')
                # If region error, try to detect bucket region
                if error_code in ['PermanentRedirect', '301']:
                    try:
                        # Get bucket region
                        bucket_location = client.get_bucket_location(Bucket=bucket)
                        region = bucket_location.get('LocationConstraint') or 'us-east-1'
                        if region != default_region:
                            logger.debug(f"Bucket {bucket} is in {region}, retrying with correct region")
                            client = boto3.client("s3", region_name=region)
                            resp = client.head_object(Bucket=bucket, Key=key)
                            size = int(resp.get("ContentLength")) if resp and "ContentLength" in resp else None
                            if size is not None:
                                logger.debug(f"Retrieved S3 object size for {uri}: {size} bytes")
                            return size
                    except Exception:
                        pass  # Fall through to error handling
                
                # Log specific error types for debugging
                if error_code == '403' or 'AccessDenied' in error_code or 'Forbidden' in error_code:
                    logger.debug(f"Access denied to {uri}, cannot get size (check IAM permissions)")
                elif error_code == '404' or 'NoSuchKey' in error_code:
                    logger.debug(f"Object not found: {uri}")
                else:
                    logger.debug(f"Failed to get size for {uri}: {error_code}: {str(e)[:100]}")
                return None
                
        except ImportError:
            logger.debug(f"boto3 not available, cannot get size for {uri}")
            return None
        except Exception as e:
            # Log other exceptions
            error_type = type(e).__name__
            if "NoCredentialsError" in error_type or "Credentials" in str(e):
                logger.debug(f"AWS credentials not configured, cannot get size for {uri}")
            else:
                logger.debug(f"Failed to get size for {uri}: {error_type}: {str(e)[:100]}")
            return None

    def _session_file_to_s3_uri(self, file_info: Dict[str, Any], session_context: Dict[str, Any]) -> Optional[str]:
        if not isinstance(file_info, dict):
            return None
        bucket = file_info.get("s3_bucket") or (session_context.get("metadata", {}) if isinstance(session_context, dict) else {}).get("s3_bucket")
        key = file_info.get("s3_key")
        if bucket and key:
            return f"s3://{bucket}/{key}"
        return None

    def _safe_int(self, v: Any) -> Optional[int]:
        try:
            if v is None:
                return None
            i = int(v)
            return i if i >= 0 else None
        except Exception:
            return None


