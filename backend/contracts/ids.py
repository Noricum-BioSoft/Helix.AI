"""
Domain identifiers, trace ids and stable hashing for platform contracts.

Rules (plan rev 2, principles 5 and 6):
- IDs for shared domain objects are minted here and nowhere else. Stores
  (local ledger, DataWeaver) never generate ids for shared objects.
- One ``trace_id`` per scientific loop; every contract carries it.
- Hashes are computed over canonical JSON so they are stable across
  processes, Python versions and re-serialisation.
"""

from __future__ import annotations

import hashlib
import json
import secrets
import time
import uuid
from typing import Any

TRACE_PREFIX = "orch_"


def uuid7() -> uuid.UUID:
    """RFC 9562 UUIDv7 (time-ordered). Python 3.9 has no ``uuid.uuid7``."""
    ms = int(time.time() * 1000) & ((1 << 48) - 1)
    rand = secrets.randbits(74)
    rand_a = rand >> 62  # 12 bits
    rand_b = rand & ((1 << 62) - 1)  # 62 bits
    value = (ms << 80) | (0x7 << 76) | (rand_a << 64) | (0b10 << 62) | rand_b
    return uuid.UUID(int=value)


def new_id() -> str:
    """Mint a new domain id (uuid7, canonical string form)."""
    return str(uuid7())


def new_trace_id() -> str:
    """Mint a new trace id for one scientific loop: ``orch_<uuid7 hex>``."""
    return TRACE_PREFIX + uuid7().hex


def is_trace_id(value: str) -> bool:
    return isinstance(value, str) and value.startswith(TRACE_PREFIX) and len(value) == len(TRACE_PREFIX) + 32


def canonical_json(obj: Any) -> str:
    """Deterministic JSON: sorted keys, compact separators, non-JSON types via str()."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str, ensure_ascii=False)


def stable_hash(obj: Any) -> str:
    """sha256 hex of ``canonical_json(obj)``. Pydantic models are dumped in JSON mode first."""
    if hasattr(obj, "model_dump"):
        obj = obj.model_dump(mode="json")
    return hashlib.sha256(canonical_json(obj).encode("utf-8")).hexdigest()
