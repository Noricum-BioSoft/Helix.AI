"""
Phase 3A — CapabilityRegistry + execution-profile behaviour.

Disabled providers are never returned; aliases resolve deterministically and
prefer an enabled provider; constraints narrow by category.
"""

from __future__ import annotations

import pytest

from backend.config.execution_profile import load_execution_profile
from backend.execution.registry import CapabilityRegistry, RegistryError, load_descriptors


def _registry(profile_name: str) -> CapabilityRegistry:
    profile = load_execution_profile(profile_name, check_adapters=False)
    return CapabilityRegistry.from_config(profile)


def test_descriptors_load_and_are_unique():
    descriptors = load_descriptors()
    ids = [d.capability_id for d in descriptors]
    assert len(ids) == len(set(ids))
    assert "local_compute:bulk_rnaseq_analysis" in ids
    assert "mock_experimental:protein_expression" in ids


def test_local_only_enabled_set_matches_profile():
    registry = _registry("local-only")
    enabled = {d.capability_id for d in registry.all(enabled_only=True)}
    assert "local_compute:bulk_rnaseq_analysis" in enabled
    # mock is a local fake and is part of the local-only profile.
    assert "mock_experimental:protein_expression" in enabled
    # data providers are not in the local-only profile, so they stay hidden.
    assert registry.find("ncbi_entrez:fetch_sequence") == []
    # ...but they still exist as descriptors when the enabled filter is dropped.
    assert registry.find("ncbi_entrez:fetch_sequence", enabled_only=False)


def test_alias_resolution_prefers_enabled_and_is_deterministic():
    # bulk_rnaseq_analysis is aliased by both local_compute and nextflow; local sorts first.
    registry = _registry("local-only")
    assert registry.resolve_alias("bulk_rnaseq_analysis") == "local_compute:bulk_rnaseq_analysis"
    assert registry.resolve_for_tool("some_unregistered_tool") == "local_compute:some_unregistered_tool"


def test_find_respects_category_constraint_and_health():
    registry = _registry("local-only")
    comp_caps = registry.find(constraints={"category": "computational"})
    assert comp_caps and all(d.category == "computational" for d in comp_caps)
    # health gate: marking a provider unhealthy removes it from find()
    registry.set_health("local_compute", False)
    assert registry.find("local_compute:bulk_rnaseq_analysis") == []
    assert registry.find("local_compute:bulk_rnaseq_analysis", healthy_only=False)


def test_unknown_profile_capability_dir_errors(tmp_path):
    from backend.config.execution_profile import load_execution_profile

    profile = load_execution_profile("local-only", check_adapters=False)
    with pytest.raises(RegistryError):
        CapabilityRegistry.from_config(profile, capabilities_dir=tmp_path / "does_not_exist")
