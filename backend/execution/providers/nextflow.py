"""
NextflowProvider (Phase 3A) — wraps Nextflow/nf-core pipeline submission.

Async execution: ``execute`` submits a run and returns a handle;
``get_status``/``retrieve_results`` poll the job. Executor, queue, workdir and
config come from the ``ProviderConfig`` (execution profile), never from here.

Idempotency: ``client_dedup``. Nextflow's ``-name`` is derived from the
idempotency key (``-name helix_{key[:12]}``); the fabric also looks the key up
in the ledger before re-submitting.

The submit/status/results backends are injectable so the conformance suite can
drive the provider with a fake (no Nextflow binary required in CI).
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
    ProviderNotConfigured,
    ValidationResult,
)
from backend.execution.providers.common import ProviderBase
from shared.capability_registry import CapabilityDescriptor

# submit(run_name, pipeline, revision, params, config) -> provider_handle (str)
SubmitFn = Callable[[str, str, Optional[str], Dict[str, Any], Dict[str, Any]], str]
# status(provider_handle) -> RunStatus
StatusFn = Callable[[str], str]
# results(provider_handle) -> {"outputs": [uri...], "cost_actual_usd": float|None, "error": str|None}
ResultsFn = Callable[[str], Dict[str, Any]]


def _run_name(idempotency_key: str) -> str:
    return f"helix_{idempotency_key[:12]}"


class NextflowProvider(ProviderBase):
    provider_id = "nextflow"

    def __init__(
        self,
        descriptors: List[CapabilityDescriptor],
        config: Optional[Dict[str, Any]] = None,
        *,
        submit_fn: Optional[SubmitFn] = None,
        status_fn: Optional[StatusFn] = None,
        results_fn: Optional[ResultsFn] = None,
    ):
        super().__init__(descriptors, config)
        self._submit_fn = submit_fn
        self._status_fn = status_fn
        self._results_fn = results_fn

    def validate_request(self, request: ExecutionRequest) -> ValidationResult:
        errors: List[str] = []
        if self._descriptor(request.capability_id) is None:
            errors.append(f"unknown capability for {self.provider_id}: {request.capability_id}")
        if not request.capability_id.startswith("nextflow:"):
            errors.append(f"capability_id {request.capability_id} is not a nextflow capability")
        return ValidationResult(ok=not errors, errors=errors)

    def estimate(self, request: ExecutionRequest) -> Estimate:
        descriptor = self._descriptor(request.capability_id)
        rng = descriptor.cost.unit_usd_range if descriptor and descriptor.cost.unit_usd_range else [0.0, 1.0]
        minutes = descriptor.turnaround.estimated_minutes if descriptor and descriptor.turnaround.estimated_minutes else 120.0
        return Estimate(cost_range_usd=[float(rng[0]), float(rng[1])], turnaround_minutes=float(minutes), confidence=0.3)

    def prepare_execution(self, request: ExecutionRequest) -> PreparedExecution:
        self._remember_trace(request.idempotency_key, request.trace_id)
        descriptor = self._descriptor(request.capability_id)
        pipeline = request.capability_id.split(":", 1)[1]
        revision = None
        if descriptor and descriptor.provenance.pinned_release:
            revision = descriptor.provenance.pinned_release.get("revision")
        return PreparedExecution(
            execution_request_id=request.execution_request_id,
            idempotency_key=request.idempotency_key,
            provider_id=self.provider_id,
            capability_id=request.capability_id,
            payload={
                "run_name": _run_name(request.idempotency_key),
                "pipeline": pipeline,
                "revision": revision,
                "params": dict(request.parameters),
                "config": dict(self._config),
            },
            staged_inputs=[i.uri for i in request.inputs],
        )

    def execute(self, prepared: PreparedExecution) -> ExecutionHandle:
        if self._submit_fn is None:
            raise ProviderNotConfigured("NextflowProvider has no submit backend (inject submit_fn or configure the executor)")
        p = prepared.payload
        provider_handle = self._submit_fn(p["run_name"], p["pipeline"], p.get("revision"), p.get("params") or {}, p.get("config") or {})
        return ExecutionHandle(
            provider_id=self.provider_id,
            provider_handle=str(provider_handle),
            idempotency_key=prepared.idempotency_key,
            execution_request_id=prepared.execution_request_id,
            submitted_at=datetime.now(timezone.utc),
        )

    def get_status(self, handle: ExecutionHandle) -> ExecutionStatus:
        if self._status_fn is None:
            return ExecutionStatus(status="unknown")
        return ExecutionStatus(status=self._status_fn(handle.provider_handle))  # type: ignore[arg-type]

    def retrieve_results(self, handle: ExecutionHandle) -> ExecutionResult:
        if self._results_fn is None:
            return ExecutionResult(status="unknown", error="no results backend")
        raw = self._results_fn(handle.provider_handle) or {}
        status: RunStatus = raw.get("status") or ("succeeded" if not raw.get("error") else "failed")
        return ExecutionResult(
            status=status,
            outputs=[ArtifactRef(uri=str(u), kind="directory") for u in (raw.get("outputs") or [])],
            cost_actual_usd=raw.get("cost_actual_usd"),
            error=raw.get("error"),
        )

    def cancel(self, handle: ExecutionHandle) -> None:
        return None

    def provenance(self, handle: ExecutionHandle) -> ProvenanceRecord:
        return ProvenanceRecord(
            event_id=f"prov:{handle.provider_handle}",
            trace_id=self._trace_for(handle.idempotency_key),
            event_type="completed",
            actor=Actor(kind="system", id=self.provider_id),
            provider_id=self.provider_id,
            execution_backend="nextflow",
        )
