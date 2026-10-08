"""
CapabilityRegistry (Phase 3A) — the provider-neutral index of what the platform
can do, filtered by the active execution profile and runtime health.

Three layers stay separate (see ``shared/capability_registry.py`` and
``backend/config/execution_profile.py``):

* descriptors — static facts (``backend/config/capabilities/*.yaml``),
* profile — which providers THIS deployment enables and how to reach them,
* health — adapter reachability probed at startup.

``find()`` returns only descriptors whose provider is profile-enabled and
(optionally) healthy; disabled providers are never returned. The planner
resolves a step's ``required_capabilities`` by id or by tool-name alias, with a
deterministic fallback to ``local_compute:{tool_name}``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

from backend.config.execution_profile import ExecutionProfile
from shared.capability_registry import CapabilityDescriptor

logger = logging.getLogger(__name__)

CAPABILITIES_DIR = Path(__file__).resolve().parent.parent / "config" / "capabilities"
DEFAULTS_FILE = "_defaults.yaml"
LOCAL_CAPABILITY_PREFIX = "local_compute:"


class RegistryError(RuntimeError):
    pass


def _deep_merge_defaults(descriptor: Dict[str, Any], defaults: Dict[str, Any]) -> Dict[str, Any]:
    """Apply defaults for keys the descriptor omits. Objects merge one level deep; lists are not merged."""
    out = dict(descriptor)
    for key, default_value in (defaults or {}).items():
        if key not in out or out[key] is None:
            out[key] = default_value
        elif isinstance(out[key], dict) and isinstance(default_value, dict):
            merged = dict(default_value)
            merged.update(out[key])
            out[key] = merged
    return out


def load_descriptors(capabilities_dir: Path = CAPABILITIES_DIR) -> List[CapabilityDescriptor]:
    """Load and validate every descriptor under ``capabilities_dir`` (except ``_defaults.yaml``)."""
    if not capabilities_dir.exists():
        raise RegistryError(f"capabilities directory not found: {capabilities_dir}")
    defaults: Dict[str, Any] = {}
    defaults_path = capabilities_dir / DEFAULTS_FILE
    if defaults_path.exists():
        defaults = (yaml.safe_load(defaults_path.read_text(encoding="utf-8")) or {}).get("defaults", {}) or {}

    descriptors: List[CapabilityDescriptor] = []
    seen: Dict[str, str] = {}
    for path in sorted(capabilities_dir.glob("*.yaml")):
        if path.name == DEFAULTS_FILE:
            continue
        raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        for entry in raw.get("capabilities", []) or []:
            merged = _deep_merge_defaults(entry, defaults)
            try:
                descriptor = CapabilityDescriptor.model_validate(merged)
            except Exception as exc:
                raise RegistryError(f"invalid capability descriptor in {path.name}: {exc}") from exc
            if descriptor.capability_id in seen:
                raise RegistryError(
                    f"duplicate capability_id '{descriptor.capability_id}' in {path.name} and {seen[descriptor.capability_id]}"
                )
            seen[descriptor.capability_id] = path.name
            descriptors.append(descriptor)
    return descriptors


class CapabilityRegistry:
    def __init__(
        self,
        descriptors: List[CapabilityDescriptor],
        profile: ExecutionProfile,
        *,
        health: Optional[Dict[str, bool]] = None,
    ):
        self.profile = profile
        self._health = dict(health or {})
        self._by_id: Dict[str, CapabilityDescriptor] = {}
        self._by_alias: Dict[str, List[str]] = {}
        for descriptor in descriptors:
            if descriptor.capability_id in self._by_id:
                raise RegistryError(f"duplicate capability_id '{descriptor.capability_id}'")
            self._by_id[descriptor.capability_id] = descriptor
            for alias in descriptor.tool_name_aliases:
                self._by_alias.setdefault(alias, []).append(descriptor.capability_id)

    @classmethod
    def from_config(
        cls,
        profile: ExecutionProfile,
        *,
        capabilities_dir: Path = CAPABILITIES_DIR,
        health: Optional[Dict[str, bool]] = None,
    ) -> "CapabilityRegistry":
        return cls(load_descriptors(capabilities_dir), profile, health=health)

    # ── enablement / health ────────────────────────────────────────────────

    def is_enabled(self, descriptor: CapabilityDescriptor) -> bool:
        return self.profile.is_enabled(descriptor.provider)

    def is_healthy(self, descriptor: CapabilityDescriptor) -> bool:
        # Unknown health is treated as healthy (health probing is best-effort in P3A).
        return self._health.get(descriptor.provider, True)

    def set_health(self, provider_id: str, healthy: bool) -> None:
        self._health[provider_id] = healthy

    # ── lookup ──────────────────────────────────────────────────────────────

    def get(self, capability_id: str) -> Optional[CapabilityDescriptor]:
        return self._by_id.get(capability_id)

    def all(self, *, enabled_only: bool = False, healthy_only: bool = False) -> List[CapabilityDescriptor]:
        out = list(self._by_id.values())
        if enabled_only:
            out = [d for d in out if self.is_enabled(d)]
        if healthy_only:
            out = [d for d in out if self.is_healthy(d)]
        return sorted(out, key=lambda d: d.capability_id)

    def resolve_alias(self, tool_name: str, *, enabled_only: bool = True) -> Optional[str]:
        """Resolve a Helix tool name to a capability id. Prefers an enabled provider; deterministic."""
        candidates = self._by_alias.get(tool_name, [])
        if enabled_only:
            enabled = [cid for cid in candidates if self.is_enabled(self._by_id[cid])]
            candidates = enabled or []
        if not candidates:
            return None
        return sorted(candidates)[0]

    def _matches_constraints(self, descriptor: CapabilityDescriptor, constraints: Dict[str, Any]) -> bool:
        category = constraints.get("category")
        if category and descriptor.category != category:
            return False
        system = constraints.get("system")
        if system and system not in descriptor.supported_systems:
            return False
        data_class = constraints.get("data_class")
        if data_class and data_class not in descriptor.security.data_classes_allowed:
            return False
        return True

    def find(
        self,
        capability_id: Optional[str] = None,
        *,
        tool_name: Optional[str] = None,
        constraints: Optional[Dict[str, Any]] = None,
        enabled_only: bool = True,
        healthy_only: bool = True,
    ) -> List[CapabilityDescriptor]:
        """Profile-enabled (and, by default, healthy) descriptors matching the query.

        Query by exact ``capability_id``, by ``tool_name`` alias, or neither
        (all). ``constraints`` narrows by category/system/data_class.
        """
        if capability_id is not None:
            descriptor = self._by_id.get(capability_id)
            candidates = [descriptor] if descriptor else []
        elif tool_name is not None:
            candidates = [self._by_id[cid] for cid in self._by_alias.get(tool_name, [])]
        else:
            candidates = list(self._by_id.values())

        constraints = constraints or {}
        out = []
        for descriptor in candidates:
            if descriptor is None:
                continue
            if enabled_only and not self.is_enabled(descriptor):
                continue
            if healthy_only and not self.is_healthy(descriptor):
                continue
            if not self._matches_constraints(descriptor, constraints):
                continue
            out.append(descriptor)
        return sorted(out, key=lambda d: d.capability_id)

    def resolve_for_tool(self, tool_name: str) -> str:
        """Capability id for a tool: a registered alias when available, else the local fallback."""
        resolved = self.resolve_alias(tool_name)
        return resolved or f"{LOCAL_CAPABILITY_PREFIX}{tool_name}"


def get_registry(
    profile: Optional[ExecutionProfile] = None,
    *,
    health: Optional[Dict[str, bool]] = None,
) -> CapabilityRegistry:
    """Registry for the active (or given) execution profile."""
    if profile is None:
        from backend.config.execution_profile import load_execution_profile

        profile = load_execution_profile(check_adapters=False)
    return CapabilityRegistry.from_config(profile, health=health)
