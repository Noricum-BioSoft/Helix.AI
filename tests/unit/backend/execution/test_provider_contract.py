"""
Phase 3A — ExecutionProvider conformance suite.

Every provider must honour the same contract: describe_capabilities returns its
descriptors; validate/prepare/execute/retrieve round-trip; idempotency_support
is declared; and the lifecycle types validate. Run against Local, Mock and a
fake-backed Nextflow (no Nextflow binary in CI).
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Dict, List

import pytest

from backend.config.execution_profile import load_execution_profile
from backend.contracts.execution_request import ExecutionRequest, make_idempotency_key
from backend.contracts.ids import new_id
from backend.execution.providers.base import (
    ExecutionHandle,
    ExecutionResult,
    ExecutionStatus,
    PreparedExecution,
    PROVIDER_METHODS,
    ValidationResult,
)
from backend.execution.providers.legacy_broker import LegacyBrokerProvider
from backend.execution.providers.local_compute import LocalComputeProvider
from backend.execution.providers.mock_experimental import MockExperimentalProvider
from backend.execution.providers.nextflow import NextflowProvider
from backend.execution.registry import CapabilityRegistry


@pytest.fixture()
def registry() -> CapabilityRegistry:
    # local-only enables local_compute, nextflow and the mock experimental provider.
    profile = load_execution_profile("local-only", check_adapters=False)
    return CapabilityRegistry.from_config(profile)


def _descriptors(registry: CapabilityRegistry, provider: str):
    return [d for d in registry.all() if d.provider == provider]


def _request(capability_id: str, *, approval: bool = True, parameters: Dict[str, Any] | None = None) -> ExecutionRequest:
    return ExecutionRequest(
        execution_request_id=new_id(),
        trace_id="orch_" + "a" * 32,
        execution_intent_id=new_id(),
        execution_intent_hash="d" * 64,
        approval_id=(new_id() if approval else None),
        idempotency_key=make_idempotency_key("d" * 64),
        provider_id="legacy:Local",
        capability_id=capability_id,
        parameters=parameters or {},
    )


def _fake_nextflow(registry: CapabilityRegistry) -> NextflowProvider:
    submitted: Dict[str, str] = {}

    def submit(run_name, pipeline, revision, params, config):
        submitted[run_name] = f"nfhandle:{run_name}"
        return submitted[run_name]

    def status(handle):
        return "succeeded"

    def results(handle):
        return {"status": "succeeded", "outputs": [f"file:///work/{handle}/results"], "cost_actual_usd": 0.0}

    return NextflowProvider(
        _descriptors(registry, "nextflow"), {"executor": "local"},
        submit_fn=submit, status_fn=status, results_fn=results,
    )


def _providers(registry: CapabilityRegistry):
    local = LocalComputeProvider(
        _descriptors(registry, "local_compute"), {},
        tool_runner=lambda tool, args: {"status": "success", "text": f"{tool} ok"},
    )
    mock = MockExperimentalProvider(_descriptors(registry, "mock_experimental_provider"), {})
    nextflow = _fake_nextflow(registry)
    legacy = LegacyBrokerProvider(
        _descriptors(registry, "ncbi_entrez"), {},
        tool_runner=lambda tool, args: {"status": "success", "text": f"{tool} fetched"},
    )
    return {
        "local_compute:bulk_rnaseq_analysis": local,
        "mock_experimental:protein_expression": mock,
        "nextflow:nf-core/rnaseq": nextflow,
        "ncbi_entrez:fetch_sequence": legacy,
    }


@pytest.mark.parametrize("capability_id", [
    "local_compute:bulk_rnaseq_analysis",
    "mock_experimental:protein_expression",
    "nextflow:nf-core/rnaseq",
    "ncbi_entrez:fetch_sequence",
])
def test_provider_conforms_to_contract(registry, capability_id):
    provider = _providers(registry)[capability_id]

    # 1. Implements every protocol method.
    for method in PROVIDER_METHODS:
        assert callable(getattr(provider, method)), f"{provider.provider_id} missing {method}"
    assert provider.idempotency_support in ("native", "client_dedup", "none")

    # 2. describe_capabilities returns this provider's descriptors.
    described = provider.describe_capabilities()
    assert described and all(d.provider == described[0].provider for d in described)

    # 3. validate → prepare → execute → status → retrieve round-trip.
    request = _request(capability_id)
    validation = provider.validate_request(request)
    assert isinstance(validation, ValidationResult) and validation.ok, validation.errors

    estimate = provider.estimate(request)
    assert len(estimate.cost_range_usd) == 2 and estimate.cost_range_usd[0] <= estimate.cost_range_usd[1]

    prepared = provider.prepare_execution(request)
    assert isinstance(prepared, PreparedExecution) and prepared.idempotency_key == request.idempotency_key

    handle = provider.execute(prepared)
    assert isinstance(handle, ExecutionHandle) and handle.idempotency_key == request.idempotency_key

    status = provider.get_status(handle)
    assert isinstance(status, ExecutionStatus)

    result = provider.retrieve_results(handle)
    assert isinstance(result, ExecutionResult) and result.status == "succeeded"

    prov = provider.provenance(handle)
    assert prov.provider_id == provider.provider_id


def test_mock_native_idempotency_same_key_same_run(registry):
    mock = MockExperimentalProvider(_descriptors(registry, "mock_experimental_provider"), {})
    request = _request("mock_experimental:protein_expression")
    h1 = mock.execute(mock.prepare_execution(request))
    # resubmit with the same key
    h2 = mock.execute(mock.prepare_execution(request))
    assert h1.provider_handle == h2.provider_handle
    r1 = mock.retrieve_results(h1).metrics["measurements"]
    r2 = mock.retrieve_results(h2).metrics["measurements"]
    assert r1 == r2  # deterministic, single synthetic run


def test_mock_refuses_without_approval(registry):
    mock = MockExperimentalProvider(_descriptors(registry, "mock_experimental_provider"), {})
    request = _request("mock_experimental:protein_expression", approval=False)
    validation = mock.validate_request(request)
    assert not validation.ok and any("approval" in e for e in validation.errors)


def test_nextflow_without_backend_raises(registry):
    nextflow = NextflowProvider(_descriptors(registry, "nextflow"), {})
    request = _request("nextflow:nf-core/rnaseq")
    from backend.execution.providers.base import ProviderNotConfigured

    with pytest.raises(ProviderNotConfigured):
        nextflow.execute(nextflow.prepare_execution(request))


def test_nextflow_run_name_is_keyed_on_idempotency(registry):
    seen: List[str] = []
    nf = NextflowProvider(
        _descriptors(registry, "nextflow"), {},
        submit_fn=lambda name, *a: (seen.append(name) or f"h:{name}"),
        status_fn=lambda h: "succeeded",
        results_fn=lambda h: {"status": "succeeded", "outputs": []},
    )
    request = _request("nextflow:nf-core/rnaseq")
    nf.execute(nf.prepare_execution(request))
    assert seen == [f"helix_{request.idempotency_key[:12]}"]
