from __future__ import annotations

import asyncio
import os
import re
import logging
import tempfile
import threading
import uuid
from datetime import datetime, timezone
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
    # FastQC: Use sync (local) for small files, async (EMR) for large files
    # The infrastructure agent will decide based on file size
    FORCE_ASYNC_TOOLS: set[str] = set()  # Removed fastqc_quality_analysis
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
            decision = await self._evaluate_routing_policy(
                tool_name="__plan__",
                estimated_bytes=estimated_bytes,
                unknown_inputs=unknown,
                inputs=inputs,
                command=req.original_command,
                session_context=req.session_context,
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
        decision = await self._evaluate_routing_policy(
            tool_name=req.tool_name,
            estimated_bytes=estimated_bytes,
            unknown_inputs=unknown,
            inputs=inputs,
            command=req.original_command,
            session_context=req.session_context,
        )

        # CloudFront has an origin response timeout (~60s). Some "local sync" tools can exceed that.
        # For FastQC local execution, prefer returning a job_id immediately and running the tool in
        # the background on the backend host.
        if req.tool_name == "fastqc_quality_analysis" and decision.mode == "sync":
            promoted = RoutingDecision(
                mode="async",
                reason="promoted_to_async_local_job_for_cloudfront_timeout",
                threshold_bytes=decision.threshold_bytes,
                estimated_bytes=decision.estimated_bytes,
                unknown_inputs=decision.unknown_inputs,
                override=decision.override or "local_job",
            )
            result = await self._submit_local_tool_job(req)
            return self._wrap_result(
                tool_name=req.tool_name,
                mode="async",
                decision=promoted,
                inputs=inputs,
                output=result,
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
        # Pass through original_command, session_context, and session_id so tools can access them
        if req.original_command:
            tool_args["original_command"] = req.original_command
        if req.session_context:
            tool_args["session_context"] = req.session_context
        if req.session_id:
            tool_args["session_id"] = req.session_id
        output = await self._tool_executor(req.tool_name, tool_args)
        return self._wrap_result(
            tool_name=req.tool_name,
            mode="sync",
            decision=decision,
            inputs=inputs,
            output=output,
        )

    async def _submit_local_tool_job(self, req: ExecutionRequest) -> Dict[str, Any]:
        """
        Run a tool locally in a background thread, returning immediately with a job_id.
        """
        from backend.job_manager import get_job_manager

        jm = get_job_manager()
        tool_args = dict(req.arguments) if req.arguments else {}
        tool_args["_from_broker"] = True
        if req.original_command:
            tool_args["original_command"] = req.original_command
        if req.session_context:
            tool_args["session_context"] = req.session_context
        if req.session_id:
            tool_args["session_id"] = req.session_id

        job_id = jm.create_local_tool_job(
            tool_name=req.tool_name,
            tool_args=tool_args,
            session_id=req.session_id,
            original_command=req.original_command,
        )

        def _runner():
            jm.set_local_job_running(job_id)
            try:
                loop = asyncio.new_event_loop()
                asyncio.set_event_loop(loop)
                output = loop.run_until_complete(self._tool_executor(req.tool_name, tool_args))
                jm.set_local_job_completed(job_id, output)
            except Exception as e:
                jm.set_local_job_failed(job_id, str(e))
            finally:
                try:
                    loop.close()
                except Exception:
                    pass

        t = threading.Thread(target=_runner, name=f"helix-local-job-{job_id}", daemon=True)
        t.start()

        return {
            "type": "job",
            "status": "submitted",
            "job_id": job_id,
            "tool_name": req.tool_name,
            "message": "Job submitted for local execution. Track progress via /jobs/{job_id}.",
        }

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

        # Promote display-critical fields to the broker envelope so that
        # build_standard_response can find them without drilling into `result`.
        promoted: Dict[str, Any] = {}
        if isinstance(output, dict):
            for _k in ("links", "visuals", "visualization_type"):
                if output.get(_k):
                    promoted[_k] = output[_k]

        return {
            "status": status,
            "type": "execution_result",
            "mode": mode,
            "tool_name": tool_name,
            "text": text,
            **promoted,
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

        # 3) From previous plan steps (e.g. step 1 alignment → input for step 2 consensus)
        prev_steps = (session_context or {}).get("previous_plan_steps") or {}
        if isinstance(prev_steps, dict):
            for key in sorted(prev_steps.keys()):
                step = prev_steps.get(key)
                if not isinstance(step, dict):
                    continue
                result = step.get("result") or {}
                alignment = result.get("alignment")
                if isinstance(alignment, list) and len(alignment) > 0:
                    lines = []
                    for rec in alignment:
                        name = (rec.get("name") or rec.get("id") or "seq")
                        seq = (rec.get("sequence") or "").strip()
                        if seq:
                            lines.append(f">{name}\n{seq}")
                    if lines:
                        fasta_content = "\n".join(lines)
                        fd, path = tempfile.mkstemp(suffix=".fasta", prefix="alignment_from_plan_", text=True)
                        try:
                            os.write(fd, fasta_content.encode("utf-8"))
                        finally:
                            os.close(fd)
                        size = len(fasta_content.encode("utf-8"))
                        assets.append(InputAsset(uri=path, size_bytes=size, source="previous_plan_step"))
                        logger.info("Discovered alignment from previous_plan_steps as input for generated tool")
                    break
                aligned_str = result.get("aligned_sequences")
                if isinstance(aligned_str, str) and aligned_str.strip():
                    fd, path = tempfile.mkstemp(suffix=".fasta", prefix="alignment_from_plan_", text=True)
                    try:
                        os.write(fd, aligned_str.strip().encode("utf-8"))
                    finally:
                        os.close(fd)
                    size = len(aligned_str.strip().encode("utf-8"))
                    assets.append(InputAsset(uri=path, size_bytes=size, source="previous_plan_step"))
                    logger.info("Discovered aligned_sequences from previous_plan_steps as input for generated tool")
                    break

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

    async def _evaluate_routing_policy(
        self, 
        *, 
        tool_name: str, 
        estimated_bytes: int, 
        unknown_inputs: int,
        inputs: Optional[List[InputAsset]] = None,
        command: Optional[str] = None,
        session_context: Optional[Dict[str, Any]] = None
    ) -> RoutingDecision:
        """
        Phase 1.2: Routing policy v1 (sync vs async).

        Evaluation order (per checklist):
        (A) Optional: Infrastructure Decision Agent (if enabled and available)
        (B) bytes threshold (default 100MB, env-configurable)
        (C) tool overrides (curated list)
        (D) timeout promotion hook (stub)
        """
        threshold = self._get_async_bytes_threshold()
        mode = "sync"
        reason = "below_threshold"
        override = None

        # (A) Use Infrastructure Decision Agent for all commands (primary decision mechanism)
        # Infrastructure agent decides optimal execution environment (Local/EC2/EMR/Batch/Lambda)
        # and we map that to sync/async routing:
        # - EMR and Batch → async (long-running jobs, distributed processing)
        # - EC2, Local, Lambda → sync (faster execution, suitable for sync)
        infra_decision = None
        if inputs and command:
            try:
                from backend.infrastructure_decision_agent import (
                    decide_infrastructure,
                    InputAsset as InfraInputAsset,
                )
                
                # Convert InputAsset to InfraInputAsset
                infra_inputs = [
                    InfraInputAsset(uri=inp.uri, size_bytes=inp.size_bytes, source=inp.source)
                    for inp in inputs
                ]
                
                # Call infrastructure decision agent (async)
                infra_decision = await decide_infrastructure(
                    command=command or tool_name,
                    inputs=infra_inputs,
                    outputs=[],
                    session_context=session_context or {}
                )
                
                # Map infrastructure decision to sync/async mode
                infra = infra_decision.infrastructure.upper()
                if infra == "EMR":
                    mode = "async"
                    reason = f"infrastructure_agent_recommended_emr: {infra_decision.reasoning[:100]}"
                elif infra == "BATCH":
                    # Batch is typically async for medium/large jobs
                    mode = "async"
                    reason = f"infrastructure_agent_recommended_batch: {infra_decision.reasoning[:100]}"
                elif infra == "EC2":
                    mode = "sync"
                    reason = f"infrastructure_agent_recommended_ec2: {infra_decision.reasoning[:100]}"
                elif infra == "LOCAL":
                    mode = "sync"
                    reason = f"infrastructure_agent_recommended_local: {infra_decision.reasoning[:100]}"
                elif infra == "LAMBDA":
                    # Lambda is typically sync (serverless, fast)
                    mode = "sync"
                    reason = f"infrastructure_agent_recommended_lambda: {infra_decision.reasoning[:100]}"
                
                logger.info(f"🔧 ExecutionBroker: Infrastructure Decision Agent recommended {infra} -> {mode} mode")
                
            except Exception as e:
                logger.warning(f"⚠️  Infrastructure Decision Agent failed in ExecutionBroker: {e}, falling back to threshold-based routing")
                # Fall through to threshold-based logic

        # Special-case: FastQC has a local implementation for small inputs.
        #
        # The Infrastructure Decision Agent prompt suggests EMR if S3 sizes are unknown
        # (conservative assumption). In practice, for FastQC we prefer attempting local
        # execution when inputs are known-small OR size is unknown, and only require EMR
        # when size clearly exceeds the async threshold.
        if tool_name == "fastqc_quality_analysis":
            if estimated_bytes > 0 and estimated_bytes <= threshold:
                if mode == "async":
                    mode = "sync"
                    override = "fastqc_prefer_local_small_inputs"
                    reason = "fastqc_override_sync_small_inputs"
            elif estimated_bytes == 0 and unknown_inputs > 0:
                # Unknown-size FastQ inputs: prefer local attempt over failing EMR submission.
                if mode == "async":
                    mode = "sync"
                    override = "fastqc_prefer_local_unknown_size"
                    reason = "fastqc_override_sync_unknown_input_sizes"
        
        # (B) Fallback: bytes threshold (if infrastructure agent unavailable or didn't decide)
        # Per infrastructure-decision-agent.md: Medium (100MB-10GB) and Large (>10GB) files use async
        # Default threshold is 100MB, so medium and large files should route to async
        if mode == "sync" and reason == "below_threshold":
            # Route to async if files exceed medium threshold (100MB)
            # This covers both medium (100MB-10GB) and large (>10GB) files
            if estimated_bytes > threshold:
                mode = "async"
                reason = "bytes_threshold_exceeded_medium_or_large_files"

        # (C) tool overrides
        # FORCE_ASYNC_TOOLS always take precedence (e.g., FastQC has no local implementation yet)
        # FORCE_SYNC_TOOLS can be overridden by infrastructure decisions
        if tool_name in self.FORCE_ASYNC_TOOLS:
            # Tool MUST be async (e.g., FastQC - no local execution implemented)
            mode = "async"
            override = "force_async"
            reason = "tool_override_force_async"
            if infra_decision and infra_decision.infrastructure.upper() == "LOCAL":
                logger.warning(f"⚠️  Infrastructure recommended LOCAL but {tool_name} requires async (no local implementation)")
        elif infra_decision and infra_decision.infrastructure.upper() == "LOCAL":
            # Infrastructure decision recommends LOCAL - respect it for tools that support local execution
            logger.info(f"🔧 ExecutionBroker: Respecting LOCAL recommendation for {tool_name}")
            # mode already set to "sync" by infrastructure decision
        elif tool_name in self.FORCE_SYNC_TOOLS:
            mode = "sync"
            override = "force_sync"
            reason = "tool_override_force_sync"

        # (D) timeout promotion hook (stub)
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
        from backend.job_manager import get_job_manager

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
        from backend.job_manager import get_job_manager

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
        from backend.job_manager import get_job_manager

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
            # So step 2+ can use step 1 output (e.g. consensus from alignment); see CONSENSUS_FROM_ALIGNMENT_GAP.md
            if step.tool_name == "handle_natural_command" and context.get("steps"):
                resolved_args = dict(resolved_args)
                resolved_args["previous_plan_steps"] = dict(context["steps"])
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


