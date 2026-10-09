"""
MockExperimentalProvider (Phase 3A) — an in-memory "wet lab" that proves the
experimental execution path without any instrument. Every run is synthetic and
flagged ``provider="mock"``.

- Refuses to run without an approval bound to the request
  (``ExecutionRequest.approval_id``); the fabric also enforces this, so it is a
  defence in depth.
- Idempotency: ``native``. The same ``idempotency_key`` always maps to the same
  fake run (same handle, same synthetic measurements) — re-submitting after a
  timeout never creates a second lab run.
"""

from __future__ import annotations

import hashlib
import random
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

from backend.contracts.execution_request import ExecutionRequest
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
from backend.execution.providers.common import ProviderBase
from shared.capability_registry import CapabilityDescriptor


class MockExperimentalProvider(ProviderBase):
    provider_id = "mock_experimental_provider"

    def __init__(self, descriptors: List[CapabilityDescriptor], config: Optional[Dict[str, Any]] = None):
        super().__init__(descriptors, config)
        # Keyed by idempotency_key → the one true run for that key.
        self._runs: Dict[str, Dict[str, Any]] = {}

    def validate_request(self, request: ExecutionRequest) -> ValidationResult:
        errors: List[str] = []
        if self._descriptor(request.capability_id) is None:
            errors.append(f"unknown capability for {self.provider_id}: {request.capability_id}")
        if not request.approval_id:
            errors.append("mock_experimental refuses to run without an approval bound to the request")
        return ValidationResult(ok=not errors, errors=errors)

    def estimate(self, request: ExecutionRequest) -> Estimate:
        descriptor = self._descriptor(request.capability_id)
        rng = (descriptor.cost.unit_usd_range if descriptor and descriptor.cost.unit_usd_range else [50.0, 200.0])
        days = descriptor.turnaround.estimated_days if descriptor and descriptor.turnaround.estimated_days else 10
        return Estimate(
            cost_range_usd=[float(rng[0]), float(rng[1])],
            turnaround_minutes=float(days) * 24 * 60,
            confidence=0.4,
            assumptions=["synthetic mock lab; not a real measurement"],
        )

    def prepare_execution(self, request: ExecutionRequest) -> PreparedExecution:
        self._remember_trace(request.idempotency_key, request.trace_id)
        return PreparedExecution(
            execution_request_id=request.execution_request_id,
            idempotency_key=request.idempotency_key,
            provider_id=self.provider_id,
            capability_id=request.capability_id,
            payload={"capability_id": request.capability_id, "parameters": dict(request.parameters)},
            staged_inputs=[i.uri for i in request.inputs],
        )

    def _synthesize(self, prepared: PreparedExecution) -> Dict[str, Any]:
        # Deterministic synthetic measurements seeded by the idempotency key.
        seed = int(hashlib.sha256(prepared.idempotency_key.encode("utf-8")).hexdigest()[:8], 16)
        rng = random.Random(seed)
        measurements = [
            {"replicate": i + 1, "value": round(rng.uniform(0.1, 10.0), 3), "unit": "a.u.", "provider": "mock"}
            for i in range(3)
        ]
        handle_id = f"mock:{prepared.idempotency_key[:16]}"
        return {
            "handle_id": handle_id,
            "measurements": measurements,
            "submitted_at": datetime.now(timezone.utc),
        }

    def execute(self, prepared: PreparedExecution) -> ExecutionHandle:
        run = self._runs.get(prepared.idempotency_key)
        if run is None:
            run = self._synthesize(prepared)
            self._runs[prepared.idempotency_key] = run  # native idempotency: first write wins
        return ExecutionHandle(
            provider_id=self.provider_id,
            provider_handle=run["handle_id"],
            idempotency_key=prepared.idempotency_key,
            execution_request_id=prepared.execution_request_id,
            submitted_at=run["submitted_at"],
        )

    def _run_for_handle(self, handle: ExecutionHandle) -> Optional[Dict[str, Any]]:
        return self._runs.get(handle.idempotency_key)

    def get_status(self, handle: ExecutionHandle) -> ExecutionStatus:
        run = self._run_for_handle(handle)
        return ExecutionStatus(status="succeeded" if run else "unknown", progress=1.0 if run else None)

    def retrieve_results(self, handle: ExecutionHandle) -> ExecutionResult:
        run = self._run_for_handle(handle)
        if run is None:
            return ExecutionResult(status="unknown", error="no mock run for handle")
        return ExecutionResult(
            status="succeeded",
            outputs=[ArtifactRef(uri=f"mock://{run['handle_id']}/measurements.json", kind="dataset")],
            metrics={"measurements": run["measurements"], "provider": "mock"},
            cost_actual_usd=None,
        )

    def cancel(self, handle: ExecutionHandle) -> None:
        self._runs.pop(handle.idempotency_key, None)

    def provenance(self, handle: ExecutionHandle) -> ProvenanceRecord:
        run = self._run_for_handle(handle)
        return ProvenanceRecord(
            event_id=f"prov:{handle.provider_handle}",
            trace_id=self._trace_for(handle.idempotency_key),
            event_type="completed",
            actor=Actor(kind="system", id=self.provider_id),
            provider_id=self.provider_id,
            execution_backend="mock_experimental",
            outputs=[f"mock://{run['handle_id']}/measurements.json"] if run else [],
        )
