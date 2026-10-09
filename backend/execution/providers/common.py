"""
Shared helpers for ExecutionProvider implementations (Phase 3A).

Providers are synchronous per the ``ExecutionProvider`` protocol. Some wrap an
async callable (e.g. Helix's ``dispatch_tool``); ``run_coroutine_sync`` runs
such a coroutine safely whether or not an event loop is already running.
"""

from __future__ import annotations

import asyncio
import threading
from typing import Any, Awaitable, Dict, List, TypeVar

from shared.capability_registry import CapabilityDescriptor, IdempotencySupport

T = TypeVar("T")


def run_coroutine_sync(coro: Awaitable[T]) -> T:
    """Run ``coro`` to completion from sync code, even inside a running loop."""
    try:
        asyncio.get_running_loop()
    except RuntimeError:
        return asyncio.run(coro)  # no loop running

    # A loop is already running in this thread: run the coroutine in a worker
    # thread with its own loop so we don't deadlock the caller's loop.
    result: Dict[str, Any] = {}

    def _runner() -> None:
        new_loop = asyncio.new_event_loop()
        try:
            result["value"] = new_loop.run_until_complete(coro)
        except BaseException as exc:  # noqa: BLE001 — re-raised below
            result["error"] = exc
        finally:
            new_loop.close()

    worker = threading.Thread(target=_runner, daemon=True)
    worker.start()
    worker.join()
    if "error" in result:
        raise result["error"]
    return result["value"]


class ProviderBase:
    """Common descriptor plumbing shared by the concrete providers.

    Subclasses set ``provider_id`` and implement the protocol methods. The
    descriptors a provider serves are injected (the registry owns the YAML), so
    providers never read config files themselves.
    """

    provider_id: str = "base"

    def __init__(self, descriptors: List[CapabilityDescriptor], config: Dict[str, Any] | None = None):
        self._descriptors = list(descriptors)
        self._config = dict(config or {})
        self._trace_by_key: Dict[str, str] = {}

    def _remember_trace(self, idempotency_key: str, trace_id: str) -> None:
        if trace_id:
            self._trace_by_key[idempotency_key] = trace_id

    def _trace_for(self, idempotency_key: str) -> str:
        from backend.contracts.ids import new_trace_id

        return self._trace_by_key.get(idempotency_key) or new_trace_id()

    def describe_capabilities(self) -> List[CapabilityDescriptor]:
        return list(self._descriptors)

    def _descriptor(self, capability_id: str) -> CapabilityDescriptor | None:
        return next((d for d in self._descriptors if d.capability_id == capability_id), None)

    @property
    def idempotency_support(self) -> IdempotencySupport:
        if self._descriptors:
            return self._descriptors[0].integration.idempotency_support
        return "none"
