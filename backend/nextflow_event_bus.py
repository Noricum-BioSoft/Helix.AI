"""Asyncio pub/sub event bus for Nextflow job lifecycle events.

Nextflow posts weblog events to POST /internal/nextflow/events.
The SSE endpoint GET /jobs/{job_id}/stream subscribes to receive them.

Each job gets a list of asyncio.Queue instances — one per connected
SSE client.  Events are broadcast to all subscribers simultaneously.

Usage::

    bus = get_event_bus()

    # Producer (weblog receiver endpoint)
    await bus.publish(job_id, {"type": "process_completed", "trace": {...}})

    # Consumer (SSE endpoint)
    q = bus.subscribe(job_id)
    try:
        event = await asyncio.wait_for(q.get(), timeout=30.0)
    finally:
        bus.unsubscribe(job_id, q)
"""

from __future__ import annotations

import asyncio
import logging
from typing import Dict, List

logger = logging.getLogger(__name__)


class NextflowEventBus:
    """Thread-safe asyncio pub/sub for Nextflow job events."""

    def __init__(self) -> None:
        self._queues: Dict[str, List[asyncio.Queue]] = {}
        # Keep a short replay buffer per job so a reconnecting client
        # catches up on events it missed.
        self._history: Dict[str, List[dict]] = {}
        self._max_history = 200

    # ------------------------------------------------------------------
    # Producer side
    # ------------------------------------------------------------------

    async def publish(self, job_id: str, event: dict) -> None:
        """Broadcast *event* to all active subscribers of *job_id*."""
        if not job_id:
            logger.warning("publish called with empty job_id — event dropped: %s", event)
            return

        # Append to replay buffer
        buf = self._history.setdefault(job_id, [])
        buf.append(event)
        if len(buf) > self._max_history:
            buf.pop(0)

        queues = list(self._queues.get(job_id, []))
        logger.debug("EventBus publish job=%s event=%s subscribers=%d", job_id, event.get("type"), len(queues))
        for q in queues:
            await q.put(event)

    # ------------------------------------------------------------------
    # Consumer side
    # ------------------------------------------------------------------

    def subscribe(self, job_id: str, replay: bool = True) -> asyncio.Queue:
        """Return a new Queue pre-loaded with any buffered events.

        Set *replay=True* (default) to replay buffered events to a
        reconnecting client so it catches up immediately.
        """
        q: asyncio.Queue = asyncio.Queue()
        self._queues.setdefault(job_id, []).append(q)

        if replay:
            for evt in self._history.get(job_id, []):
                q.put_nowait(evt)

        logger.debug("EventBus subscribe job=%s queue_id=%d", job_id, id(q))
        return q

    def unsubscribe(self, job_id: str, q: asyncio.Queue) -> None:
        """Remove *q* from the subscriber list for *job_id*."""
        subscribers = self._queues.get(job_id, [])
        try:
            subscribers.remove(q)
        except ValueError:
            pass
        if not subscribers:
            self._queues.pop(job_id, None)
        logger.debug("EventBus unsubscribe job=%s", job_id)

    def clear_history(self, job_id: str) -> None:
        """Remove replay buffer for a finished job to free memory."""
        self._history.pop(job_id, None)


# ---------------------------------------------------------------------------
# Module-level singleton
# ---------------------------------------------------------------------------

_event_bus: NextflowEventBus = NextflowEventBus()


def get_event_bus() -> NextflowEventBus:
    """Return the process-wide NextflowEventBus singleton."""
    return _event_bus
