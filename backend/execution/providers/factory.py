"""
Provider factory (Phase 3A) — instantiate the ExecutionProvider adapters that
the active execution profile enables, keyed by provider name
(``descriptor.provider`` == ``ProviderConfig.id``).

The adapter class comes from the capability descriptor's
``integration.adapter`` (``module:Class``). Each provider instance receives only
the descriptors it serves and its ``ProviderConfig.config`` — never the whole
registry, never environment variables.
"""

from __future__ import annotations

import importlib
import logging
from typing import Any, Dict, List

from backend.config.execution_profile import ExecutionProfile
from backend.execution.registry import CapabilityRegistry
from shared.capability_registry import CapabilityDescriptor

logger = logging.getLogger(__name__)


class FactoryError(RuntimeError):
    pass


def _load_adapter(adapter: str):
    module_name, _, class_name = adapter.partition(":")
    if not module_name or not class_name:
        raise FactoryError(f"invalid adapter spec (expected 'module:Class'): {adapter!r}")
    try:
        module = importlib.import_module(module_name)
    except ImportError as exc:
        raise FactoryError(f"adapter module not importable: {module_name} ({exc})") from exc
    try:
        return getattr(module, class_name)
    except AttributeError as exc:
        raise FactoryError(f"adapter class not found: {adapter}") from exc


def build_providers(
    profile: ExecutionProfile,
    registry: CapabilityRegistry,
    **adapter_overrides: Any,
) -> Dict[str, Any]:
    """Build one provider instance per enabled provider name.

    ``adapter_overrides`` maps a provider name to a pre-built instance (used by
    tests to inject fakes, e.g. a fake-backed NextflowProvider).
    """
    by_provider: Dict[str, List[CapabilityDescriptor]] = {}
    adapters: Dict[str, str] = {}
    for descriptor in registry.all():
        by_provider.setdefault(descriptor.provider, []).append(descriptor)
        adapters.setdefault(descriptor.provider, descriptor.integration.adapter)

    providers: Dict[str, Any] = {}
    for provider_cfg in profile.enabled_providers():
        name = provider_cfg.id
        if name in adapter_overrides:
            providers[name] = adapter_overrides[name]
            continue
        descriptors = by_provider.get(name)
        if not descriptors:
            logger.debug("[fabric] profile enables provider %s but no descriptor serves it; skipping", name)
            continue
        adapter_spec = provider_cfg.adapter or adapters.get(name)
        if not adapter_spec:
            raise FactoryError(f"no adapter for provider {name}")
        cls = _load_adapter(adapter_spec)
        providers[name] = cls(descriptors, provider_cfg.config)
    return providers
