"""
Typed accessors for platform feature flags. All default OFF.

Flags are read at call time (not import time) so tests can monkeypatch the
environment. Truthy values: 1, true, yes, on (case-insensitive).
"""

from __future__ import annotations

import os
from typing import Literal

KnowledgeStoreMode = Literal["off", "local", "dataweaver", "dual"]

SCIENCE_GATE_V1 = "HELIX_SCIENCE_GATE_V1"
EXECUTION_FABRIC_V1 = "HELIX_EXECUTION_FABRIC_V1"
EXECUTION_RECOMMENDER_V1 = "HELIX_EXECUTION_RECOMMENDER_V1"
KNOWLEDGE_STORE = "HELIX_KNOWLEDGE_STORE"

_TRUTHY = {"1", "true", "yes", "on"}


def _flag(name: str) -> bool:
    return os.getenv(name, "").strip().lower() in _TRUTHY


def science_gate_enabled() -> bool:
    return _flag(SCIENCE_GATE_V1)


def execution_fabric_enabled() -> bool:
    return _flag(EXECUTION_FABRIC_V1)


def execution_recommender_enabled() -> bool:
    return _flag(EXECUTION_RECOMMENDER_V1)


def knowledge_store_mode() -> KnowledgeStoreMode:
    value = os.getenv(KNOWLEDGE_STORE, "off").strip().lower()
    if value not in ("off", "local", "dataweaver", "dual"):
        raise ValueError(f"{KNOWLEDGE_STORE} must be one of off|local|dataweaver|dual, got {value!r}")
    return value  # type: ignore[return-value]


def snapshot() -> dict:
    """For release_readiness.json / diagnostics."""
    return {
        SCIENCE_GATE_V1: science_gate_enabled(),
        EXECUTION_FABRIC_V1: execution_fabric_enabled(),
        EXECUTION_RECOMMENDER_V1: execution_recommender_enabled(),
        KNOWLEDGE_STORE: knowledge_store_mode(),
    }
