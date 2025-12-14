from __future__ import annotations

import asyncio
import os
import re
from dataclasses import dataclass
from typing import Any, Awaitable, Callable, Dict, List, Optional, Tuple
from pathlib import Path


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
        inputs = self._discover_inputs(req.arguments or {}, req.session_context or {})
        estimated_bytes, unknown = self._estimate_total_bytes(inputs)
        decision = self._evaluate_routing_policy(
            tool_name=req.tool_name,
            estimated_bytes=estimated_bytes,
            unknown_inputs=unknown,
        )

        if decision.mode == "async" and req.tool_name == "fastqc_quality_analysis":
            # Phase 1: only FastQC is wired for async (EMR job submission).
            result = await self._submit_fastqc_job(req)
            return self._wrap_result(
                tool_name=req.tool_name,
                mode="async",
                decision=decision,
                inputs=inputs,
                output=result,
            )

        # Default: sync execution via existing tool executor (local/EC2/tool-generator).
        output = await self._tool_executor(req.tool_name, req.arguments)
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
        m = S3_URI_RE.match(uri or "")
        if not m:
            return None
        bucket, key = m.group(1), m.group(2)
        try:
            import boto3  # type: ignore
            client = boto3.client("s3")
            resp = client.head_object(Bucket=bucket, Key=key)
            return int(resp.get("ContentLength")) if resp and "ContentLength" in resp else None
        except Exception:
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


