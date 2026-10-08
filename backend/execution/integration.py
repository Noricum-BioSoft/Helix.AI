"""
Execution Fabric integration seam (Phase 3A).

``build_default_fabric`` wires the active profile → registry → providers →
fabric. ``try_fabric_execution`` is the fallback-safe hook the ExecutionBroker
calls: it delegates a single-tool run to the fabric only when

* ``HELIX_EXECUTION_FABRIC_V1`` is on,
* the session checkpoint is ``READY_TO_EXECUTE`` with a staged, approved intent,
* and the intent's capability matches the tool being run.

In every other case it returns ``None`` and the caller keeps its legacy path,
so behaviour is unchanged with the flag off (or when no platform intent exists).
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Optional

from backend.config.feature_flags import execution_fabric_enabled
from backend.contracts.execution_request import ExecutionRequest, ExecutionRun, make_idempotency_key
from backend.contracts.ids import new_id
from backend.execution.fabric import ExecutionFabric
from backend.execution.providers.factory import build_providers
from backend.execution.registry import CapabilityRegistry, get_registry

logger = logging.getLogger(__name__)


def build_default_fabric(**adapter_overrides: Any) -> ExecutionFabric:
    from backend.config.execution_profile import load_execution_profile

    profile = load_execution_profile(check_adapters=False)
    registry = CapabilityRegistry.from_config(profile)
    providers = build_providers(profile, registry, **adapter_overrides)
    return ExecutionFabric(registry, providers)


def _run_to_result(run: ExecutionRun) -> Dict[str, Any]:
    return {
        "status": "success" if run.status == "succeeded" else run.status,
        "run_id": run.execution_run_id,
        "provider_id": run.provider_id,
        "provider_handle": run.provider_handle,
        "outputs": [o.model_dump(mode="json") for o in run.outputs],
        "error": run.error,
        "via": "execution_fabric",
        "text": f"Executed via {run.provider_id} ({run.status}).",
    }


def try_fabric_execution(
    tool_name: str,
    arguments: Dict[str, Any],
    session_context: Optional[Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
    """Delegate ``tool_name`` to the fabric when a matching approved intent is staged; else ``None``."""
    if not execution_fabric_enabled():
        return None
    session_id = (session_context or {}).get("session_id")
    if not session_id:
        return None
    try:
        from backend.history_manager import history_manager
        from backend.orchestration.approval_service import verify_approved_for_execution
        from backend.orchestration.ledger import get_ledger
        from backend.workflow_checkpoint import WorkflowState

        checkpoint = history_manager.load_checkpoint(session_id)
        if checkpoint.state not in (WorkflowState.READY_TO_EXECUTE, WorkflowState.EXECUTING):
            return None
        if not checkpoint.pending_execution_intent_id:
            return None

        ledger = get_ledger()
        intent = ledger.load_intent(session_id, checkpoint.pending_execution_intent_id)
        if intent is None:
            return None

        registry = get_registry()
        # Only delegate when the capability resolves to the same tool being run.
        if registry.resolve_for_tool(tool_name) != intent.capability_id and intent.capability_id != f"local_compute:{tool_name}":
            return None
        if registry.get(intent.capability_id) is None:
            return None

        approval = verify_approved_for_execution(session_id, checkpoint, ledger=ledger)

        request = ExecutionRequest(
            execution_request_id=new_id(),
            trace_id=intent.trace_id,
            execution_intent_id=intent.execution_intent_id,
            execution_intent_hash=intent.execution_intent_hash,
            approval_id=approval.approval_id,
            idempotency_key=make_idempotency_key(intent.execution_intent_hash),
            provider_id=intent.provider_id,
            capability_id=intent.capability_id,
            parameters=dict(arguments or {}),
        )
        fabric = build_default_fabric()
        run = fabric.run(session_id, intent, approval, request)
        logger.info("[fabric] broker delegated tool=%s → run=%s status=%s", tool_name, run.execution_run_id, run.status)
        return _run_to_result(run)
    except Exception:
        # The fabric seam must never break the legacy path; fall back on any error.
        logger.warning("[fabric] delegation failed; falling back to legacy path", exc_info=True)
        return None
