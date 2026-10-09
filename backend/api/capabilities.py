"""
Capabilities endpoint (Phase 3A).

    GET /capabilities        — descriptors with profile enablement + health.
    GET /capabilities/{id}   — a single descriptor.

Read-only. Reflects the active execution profile (``HELIX_EXECUTION_PROFILE``);
``enabled`` means this deployment turns the provider on, ``healthy`` reflects
best-effort adapter health (defaults to true until probing lands).
"""

from __future__ import annotations

from typing import Optional

from fastapi import APIRouter, HTTPException

from backend.execution.registry import CapabilityRegistry, get_registry

router = APIRouter(tags=["capabilities"])


def _registry() -> CapabilityRegistry:
    return get_registry()


def _descriptor_payload(registry: CapabilityRegistry, descriptor) -> dict:
    return {
        "capability_id": descriptor.capability_id,
        "provider": descriptor.provider,
        "category": descriptor.category,
        "description": descriptor.description,
        "enabled": registry.is_enabled(descriptor),
        "healthy": registry.is_healthy(descriptor),
        "execution_mode": descriptor.integration.execution_mode,
        "idempotency_support": descriptor.integration.idempotency_support,
        "human_approval_required": descriptor.security.human_approval_required,
        "screening_required": descriptor.security.screening_required,
        "data_classes_allowed": descriptor.security.data_classes_allowed,
        "tool_name_aliases": descriptor.tool_name_aliases,
    }


@router.get("/capabilities")
async def list_capabilities(enabled: Optional[bool] = None, category: Optional[str] = None):
    registry = _registry()
    descriptors = registry.all()
    out = []
    for descriptor in descriptors:
        if category and descriptor.category != category:
            continue
        is_enabled = registry.is_enabled(descriptor)
        if enabled is not None and is_enabled != enabled:
            continue
        out.append(_descriptor_payload(registry, descriptor))
    return {
        "success": True,
        "profile": registry.profile.profile,
        "count": len(out),
        "capabilities": out,
    }


@router.get("/capabilities/{capability_id:path}")
async def get_capability(capability_id: str):
    registry = _registry()
    descriptor = registry.get(capability_id)
    if descriptor is None:
        raise HTTPException(status_code=404, detail=f"capability not found: {capability_id}")
    return {"success": True, "profile": registry.profile.profile, "capability": _descriptor_payload(registry, descriptor)}
