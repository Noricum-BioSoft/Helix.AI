"""
ExecutionFabric (Phase 3A) — the single place a staged, approved intent becomes
a real provider submission.

``run(session_id, intent, approval, request)``:

1. Invariants (fail closed): request binds to the intent; non read-only work
   needs an approved HumanApproval bound to the intent; experimental work always
   needs approval.
2. Idempotency: if a submission for ``request.idempotency_key`` already exists in
   the ledger, re-attach to its run instead of resubmitting (money for compute,
   physical duplicates for wet-lab). Otherwise persist the request *first*, then
   ``validate → prepare → execute``; if ``execute`` times out, re-query by key
   before any retry.
3. Record the ``ExecutionRun`` and a ``ProvenanceRecord`` carrying the trace id.

Providers are resolved from the capability descriptor's ``provider`` (the
registry owns that mapping); the intent's ``provider_id`` is still validated
against the request.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Optional

from backend.contracts.execution_intent import ExecutionIntent
from backend.contracts.execution_request import ExecutionRequest, ExecutionRun
from backend.contracts.human_approval import HumanApproval
from backend.contracts.ids import new_id
from backend.contracts.provenance import Actor, ProvenanceRecord
from backend.execution.providers.base import ExecutionHandle, PreparedExecution
from backend.execution.registry import CapabilityRegistry
from backend.orchestration import invariants
from backend.orchestration.ledger import LocalLedger, get_ledger

logger = logging.getLogger(__name__)

READ_ONLY_CATEGORIES = frozenset({"data", "advisory"})


class FabricError(RuntimeError):
    pass


class ProviderTimeout(RuntimeError):
    """Raised by a provider's execute() to signal an ambiguous timeout (submission may or may not have landed)."""


class ExecutionFabric:
    def __init__(
        self,
        registry: CapabilityRegistry,
        providers: Dict[str, Any],
        ledger: Optional[LocalLedger] = None,
    ):
        self.registry = registry
        self.providers = providers
        self._ledger = ledger

    def _ledger_for(self) -> LocalLedger:
        return self._ledger or get_ledger()

    def _provider_for(self, capability_id: str):
        descriptor = self.registry.get(capability_id)
        if descriptor is None:
            raise FabricError(f"unknown capability: {capability_id}")
        provider = self.providers.get(descriptor.provider)
        if provider is None:
            raise FabricError(f"no enabled provider for {descriptor.provider} (capability {capability_id})")
        return descriptor, provider

    def run(
        self,
        session_id: str,
        intent: ExecutionIntent,
        approval: Optional[HumanApproval],
        request: ExecutionRequest,
    ) -> ExecutionRun:
        ledger = self._ledger_for()
        descriptor, provider = self._provider_for(request.capability_id)

        # ── invariants (fail closed) ─────────────────────────────────────────
        invariants.check_request_matches_intent(request, intent)
        read_only = descriptor.category in READ_ONLY_CATEGORIES
        if not read_only:
            invariants.check_approval_is_approved(approval)
            invariants.check_approval_matches_intent(approval, intent)  # type: ignore[arg-type]
        invariants.check_experimental_requires_approval(descriptor.category, approval)

        # ── idempotency: re-attach to an existing submission ─────────────────
        existing = ledger.request_by_idempotency_key(session_id, request.idempotency_key)
        if existing is not None:
            run = ledger.run_for_request(session_id, existing.execution_request_id)
            if run is not None:
                logger.info(
                    "[fabric] idempotent re-attach trace=%s key=%s run=%s (no resubmit)",
                    intent.trace_id, request.idempotency_key[:12], run.execution_run_id,
                )
                return run
            # Request persisted but no run yet (crash between persist and execute):
            # re-query the provider by key is not possible generically, so continue
            # to execute using the persisted request (safe for client_dedup/native).
            request = existing

        # First submission: persist the request BEFORE any provider call.
        ledger.record_request(session_id, request)

        validation = provider.validate_request(request)
        if not validation.ok:
            run = self._failed_run(intent, request, "; ".join(validation.errors) or "validation failed")
            ledger.record_run(session_id, run)
            return run

        prepared: PreparedExecution = provider.prepare_execution(request)
        try:
            handle: ExecutionHandle = provider.execute(prepared)
        except ProviderTimeout:
            # Ambiguous: the submission may have landed. Re-query by key before any retry.
            run = ledger.run_for_request(session_id, request.execution_request_id)
            if run is not None:
                return run
            logger.warning("[fabric] execute timed out with no recorded run; leaving request pending (no blind retry)")
            run = self._pending_run(intent, request)
            ledger.record_run(session_id, run)
            return run

        result = provider.retrieve_results(handle)
        run = ExecutionRun(
            execution_run_id=new_id(),
            trace_id=intent.trace_id,
            execution_request_id=request.execution_request_id,
            status=result.status,
            provider_id=request.provider_id,
            provider_handle=handle.provider_handle,
            finished_at=None,
            cost_actual_usd=result.cost_actual_usd,
            outputs=list(result.outputs),
            error=result.error,
        )
        ledger.record_run(session_id, run)
        self._record_provenance(session_id, intent, descriptor, handle, result)
        logger.info(
            "[fabric] executed trace=%s provider=%s capability=%s status=%s run=%s",
            intent.trace_id, descriptor.provider, request.capability_id, run.status, run.execution_run_id,
        )
        return run

    # ── helpers ──────────────────────────────────────────────────────────────

    def _failed_run(self, intent: ExecutionIntent, request: ExecutionRequest, error: str) -> ExecutionRun:
        return ExecutionRun(
            execution_run_id=new_id(),
            trace_id=intent.trace_id,
            execution_request_id=request.execution_request_id,
            status="failed",
            provider_id=request.provider_id,
            error=error,
        )

    def _pending_run(self, intent: ExecutionIntent, request: ExecutionRequest) -> ExecutionRun:
        return ExecutionRun(
            execution_run_id=new_id(),
            trace_id=intent.trace_id,
            execution_request_id=request.execution_request_id,
            status="unknown",
            provider_id=request.provider_id,
            error="submission ambiguous after timeout; re-query before retry",
        )

    def _record_provenance(self, session_id: str, intent, descriptor, handle, result) -> None:
        record = ProvenanceRecord(
            event_id=new_id(),
            trace_id=intent.trace_id,
            event_type="completed" if result.status == "succeeded" else "submitted",
            actor=Actor(kind="system", id=descriptor.provider),
            provider_id=descriptor.provider,
            execution_backend=descriptor.provider,
            inputs=[],
            outputs=[o.uri for o in result.outputs if o.uri],
        )
        try:
            self._ledger_for().record_audit_event(
                session_id,
                {
                    "event_type": "execution_run",
                    "trace_id": intent.trace_id,
                    "provider_id": descriptor.provider,
                    "capability_id": descriptor.capability_id,
                    "status": result.status,
                    "provider_handle": handle.provider_handle,
                },
            )
        except Exception:  # provenance is best-effort; never fail a run on audit write
            logger.debug("[fabric] audit write failed", exc_info=True)
        return record
