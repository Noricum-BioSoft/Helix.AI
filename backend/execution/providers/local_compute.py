"""
LocalComputeProvider (Phase 3A) — runs Helix's in-process / sandbox tools via
``dispatch_tool``. Synchronous execution: ``execute`` runs the tool and caches
the result against the handle, which ``get_status``/``retrieve_results`` read.

Idempotency: ``client_dedup``. A local tool run is not transactional, so the
*fabric* dedups by looking up ``idempotency_key`` in the ledger before calling
``execute``; this provider does not re-run a key it has already seen in-process.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Callable, Dict, List, Optional

from backend.contracts.execution_request import ExecutionRequest, RunStatus
from backend.contracts.provenance import Actor, ProvenanceRecord
from backend.contracts.scientific_objective import ArtifactRef
from backend.execution.providers.base import (
    Estimate,
    ExecutionHandle,
    ExecutionResult,
    ExecutionStatus,
    PreparedExecution,
    ValidationResult,
)
from backend.execution.providers.common import ProviderBase, run_coroutine_sync
from shared.capability_registry import CapabilityDescriptor

ToolRunner = Callable[[str, Dict[str, Any]], Dict[str, Any]]


def _default_tool_runner(tool_name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
    from backend.main import dispatch_tool

    return run_coroutine_sync(dispatch_tool(tool_name, arguments))


_SUCCESS_STATUSES = {"success", "succeeded", "ok", "complete", "completed"}


class LocalComputeProvider(ProviderBase):
    provider_id = "local_compute"

    def __init__(
        self,
        descriptors: List[CapabilityDescriptor],
        config: Optional[Dict[str, Any]] = None,
        *,
        tool_runner: Optional[ToolRunner] = None,
    ):
        super().__init__(descriptors, config)
        self._tool_runner = tool_runner or _default_tool_runner
        self._results: Dict[str, ExecutionResult] = {}

    def validate_request(self, request: ExecutionRequest) -> ValidationResult:
        errors: List[str] = []
        if self._descriptor(request.capability_id) is None:
            errors.append(f"unknown capability for {self.provider_id}: {request.capability_id}")
        if not request.capability_id.startswith("local_compute:"):
            errors.append(f"capability_id {request.capability_id} is not a local_compute capability")
        return ValidationResult(ok=not errors, errors=errors)

    def estimate(self, request: ExecutionRequest) -> Estimate:
        return Estimate(cost_range_usd=[0.0, 0.0], turnaround_minutes=5.0, confidence=0.5, assumptions=["local run, no cloud cost"])

    def prepare_execution(self, request: ExecutionRequest) -> PreparedExecution:
        self._remember_trace(request.idempotency_key, request.trace_id)
        tool_name = request.capability_id.split(":", 1)[1]
        return PreparedExecution(
            execution_request_id=request.execution_request_id,
            idempotency_key=request.idempotency_key,
            provider_id=self.provider_id,
            capability_id=request.capability_id,
            payload={"tool_name": tool_name, "arguments": dict(request.parameters)},
            staged_inputs=[i.uri for i in request.inputs],
        )

    def execute(self, prepared: PreparedExecution) -> ExecutionHandle:
        tool_name = prepared.payload["tool_name"]
        arguments = prepared.payload.get("arguments") or {}
        raw = self._tool_runner(tool_name, arguments)
        status: RunStatus = "succeeded" if str((raw or {}).get("status", "success")).lower() in _SUCCESS_STATUSES else "failed"
        outputs: List[ArtifactRef] = []
        for art in (raw or {}).get("downloadable_artifacts", []) or []:
            if isinstance(art, dict) and art.get("url"):
                outputs.append(ArtifactRef(uri=str(art["url"]), kind=str(art.get("type") or "file")))
        handle_id = f"local:{prepared.idempotency_key[:16]}"
        self._results[handle_id] = ExecutionResult(
            status=status,
            outputs=outputs,
            metrics={"tool_name": tool_name},
            cost_actual_usd=0.0,
            error=None if status == "succeeded" else str((raw or {}).get("error") or (raw or {}).get("text") or "tool failed"),
        )
        return ExecutionHandle(
            provider_id=self.provider_id,
            provider_handle=handle_id,
            idempotency_key=prepared.idempotency_key,
            execution_request_id=prepared.execution_request_id,
            submitted_at=datetime.now(timezone.utc),
        )

    def get_status(self, handle: ExecutionHandle) -> ExecutionStatus:
        result = self._results.get(handle.provider_handle)
        if result is None:
            return ExecutionStatus(status="unknown")
        return ExecutionStatus(status=result.status, progress=1.0 if result.status == "succeeded" else None)

    def retrieve_results(self, handle: ExecutionHandle) -> ExecutionResult:
        result = self._results.get(handle.provider_handle)
        if result is None:
            return ExecutionResult(status="unknown", error="no result for handle")
        return result

    def cancel(self, handle: ExecutionHandle) -> None:  # local runs are synchronous; nothing to cancel
        return None

    def provenance(self, handle: ExecutionHandle) -> ProvenanceRecord:
        result = self._results.get(handle.provider_handle)
        return ProvenanceRecord(
            event_id=f"prov:{handle.provider_handle}",
            trace_id=self._trace_for(handle.idempotency_key),
            event_type="completed",
            actor=Actor(kind="system", id=self.provider_id),
            provider_id=self.provider_id,
            execution_backend="local",
            outputs=[o.uri for o in (result.outputs if result else [])],
        )
