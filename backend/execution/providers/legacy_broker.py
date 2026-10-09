"""
LegacyBrokerProvider (Phase 3A) — the compatibility adapter for capabilities
that are not (yet) first-class providers: data-retrieval tools and any tool
still dispatched through Helix's existing path. Wraps ``dispatch_tool`` so the
fabric can treat these uniformly.

Idempotency: ``none`` — these are read-only/external fetches; the fabric still
records the request, but there is no cross-run dedup guarantee here.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Callable, Dict, List, Optional

from backend.contracts.execution_request import ExecutionRequest, RunStatus
from backend.contracts.provenance import Actor, ProvenanceRecord
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

_SUCCESS_STATUSES = {"success", "succeeded", "ok", "complete", "completed"}


def _default_tool_runner(tool_name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
    from backend.main import dispatch_tool

    return run_coroutine_sync(dispatch_tool(tool_name, arguments))


class LegacyBrokerProvider(ProviderBase):
    provider_id = "legacy_broker"

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
        return ValidationResult(ok=True)

    def estimate(self, request: ExecutionRequest) -> Estimate:
        return Estimate(cost_range_usd=[0.0, 0.0], turnaround_minutes=1.0, confidence=0.3)

    def prepare_execution(self, request: ExecutionRequest) -> PreparedExecution:
        self._remember_trace(request.idempotency_key, request.trace_id)
        tool_name = request.capability_id.split(":", 1)[1] if ":" in request.capability_id else request.capability_id
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
        raw = self._tool_runner(tool_name, prepared.payload.get("arguments") or {})
        status: RunStatus = "succeeded" if str((raw or {}).get("status", "success")).lower() in _SUCCESS_STATUSES else "failed"
        handle_id = f"legacy:{prepared.idempotency_key[:16]}"
        self._results[handle_id] = ExecutionResult(
            status=status,
            metrics={"tool_name": tool_name},
            error=None if status == "succeeded" else str((raw or {}).get("error") or "tool failed"),
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
        return ExecutionStatus(status=result.status if result else "unknown")

    def retrieve_results(self, handle: ExecutionHandle) -> ExecutionResult:
        return self._results.get(handle.provider_handle) or ExecutionResult(status="unknown", error="no result for handle")

    def cancel(self, handle: ExecutionHandle) -> None:
        return None

    def provenance(self, handle: ExecutionHandle) -> ProvenanceRecord:
        return ProvenanceRecord(
            event_id=f"prov:{handle.provider_handle}",
            trace_id=self._trace_for(handle.idempotency_key),
            event_type="completed",
            actor=Actor(kind="system", id=self.provider_id),
            provider_id=self.provider_id,
            execution_backend="legacy_broker",
        )
