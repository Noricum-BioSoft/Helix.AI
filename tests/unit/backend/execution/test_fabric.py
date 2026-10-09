"""
Phase 3A — ExecutionFabric: approval binding, fail-closed behaviour, and
idempotent submission (a duplicate run with the same key produces exactly one
provider submission, even across a simulated timeout).
"""

from __future__ import annotations

import pytest

from backend.config.execution_profile import load_execution_profile
from backend.contracts.execution_request import ExecutionRequest, make_idempotency_key
from backend.contracts.human_approval import HumanApproval, Principal
from backend.contracts.ids import new_id
from backend.execution.fabric import ExecutionFabric, ProviderTimeout
from backend.execution.providers.base import ExecutionHandle
from backend.execution.providers.factory import build_providers
from backend.execution.providers.local_compute import LocalComputeProvider
from backend.execution.providers.mock_experimental import MockExperimentalProvider
from backend.execution.registry import CapabilityRegistry
from backend.orchestration.ledger import LocalLedger
from backend.orchestration.intent_builder import build_intent
from backend.orchestration.objective_builder import build_objective_deterministic
from backend.orchestration.plan_staging import build_scientific_plan

TRACE = "orch_" + "b" * 32
PRINCIPAL = Principal(subject_id="alice", identity_provider="dev_header", auth_method="header", roles=["approver"])


def _plan_and_intent(capability_id: str, tool_name: str):
    cmd = "Run the analysis"
    objective = build_objective_deterministic(cmd, {"session_id": "s"}, None, trace_id=TRACE)
    plan = build_scientific_plan(
        {"version": "v1", "steps": [{"id": "s1", "action_type": "run_analysis", "tool_name": tool_name, "arguments": {}}]},
        objective, cmd,
    )
    intent = build_intent(plan, None, capability_id=capability_id)
    return plan, intent


def _approval(intent, plan, *, approved: bool = True) -> HumanApproval:
    return HumanApproval(
        approval_id=new_id(), trace_id=intent.trace_id,
        execution_intent_id=intent.execution_intent_id, execution_intent_hash=intent.execution_intent_hash,
        plan_id=plan.plan_id, plan_hash=plan.plan_hash,
        decision="approved" if approved else "rejected", principal=PRINCIPAL,
    )


def _request(intent, approval) -> ExecutionRequest:
    return ExecutionRequest(
        execution_request_id=new_id(), trace_id=intent.trace_id,
        execution_intent_id=intent.execution_intent_id, execution_intent_hash=intent.execution_intent_hash,
        approval_id=(approval.approval_id if approval else None),
        idempotency_key=make_idempotency_key(intent.execution_intent_hash),
        provider_id=intent.provider_id, capability_id=intent.capability_id,
    )


@pytest.fixture()
def ledger(tmp_path) -> LocalLedger:
    return LocalLedger(tmp_path)


def _local_fabric(ledger, *, runner=None):
    profile = load_execution_profile("local-only", check_adapters=False)
    registry = CapabilityRegistry.from_config(profile)
    local = LocalComputeProvider(
        [d for d in registry.all() if d.provider == "local_compute"], {},
        tool_runner=runner or (lambda tool, args: {"status": "success", "text": f"{tool} ok"}),
    )
    providers = build_providers(profile, registry, local_compute=local)
    return ExecutionFabric(registry, providers, ledger), local


def _mock_fabric(ledger):
    profile = load_execution_profile("local-only", check_adapters=False)
    registry = CapabilityRegistry.from_config(profile)
    providers = build_providers(profile, registry)
    return ExecutionFabric(registry, providers, ledger)


def test_local_execution_with_approval_succeeds(ledger):
    plan, intent = _plan_and_intent("local_compute:bulk_rnaseq_analysis", "bulk_rnaseq_analysis")
    approval = _approval(intent, plan)
    fabric, _ = _local_fabric(ledger)
    run = fabric.run("s", intent, approval, _request(intent, approval))
    assert run.status == "succeeded"
    assert ledger.load_run("s", run.execution_run_id) is not None


def test_fail_closed_without_approval(ledger):
    plan, intent = _plan_and_intent("local_compute:bulk_rnaseq_analysis", "bulk_rnaseq_analysis")
    fabric, _ = _local_fabric(ledger)
    from backend.orchestration.invariants import InvariantViolation

    with pytest.raises(InvariantViolation):
        fabric.run("s", intent, None, _request(intent, None))


def test_fail_closed_on_rejected_approval(ledger):
    plan, intent = _plan_and_intent("local_compute:bulk_rnaseq_analysis", "bulk_rnaseq_analysis")
    approval = _approval(intent, plan, approved=False)
    fabric, _ = _local_fabric(ledger)
    from backend.orchestration.invariants import InvariantViolation

    with pytest.raises(InvariantViolation):
        fabric.run("s", intent, approval, _request(intent, approval))


def test_duplicate_execute_same_key_is_one_submission(ledger):
    calls = {"n": 0}

    def counting_runner(tool, args):
        calls["n"] += 1
        return {"status": "success", "text": "ok"}

    plan, intent = _plan_and_intent("local_compute:bulk_rnaseq_analysis", "bulk_rnaseq_analysis")
    approval = _approval(intent, plan)
    fabric, _ = _local_fabric(ledger, runner=counting_runner)

    req1 = _request(intent, approval)
    run1 = fabric.run("s", intent, approval, req1)
    # retry with a fresh request id but the SAME idempotency key (same intent hash)
    req2 = req1.model_copy(update={"execution_request_id": new_id()})
    run2 = fabric.run("s", intent, approval, req2)

    assert run1.execution_run_id == run2.execution_run_id  # re-attached, not re-run
    assert calls["n"] == 1  # provider executed exactly once
    runs = ledger.records_by_trace("s", intent.trace_id)["execution_runs"]
    assert len(runs) == 1


def test_timeout_after_execute_does_not_blind_retry(ledger):
    def timing_out_runner(tool, args):
        raise ProviderTimeout("network timeout")

    plan, intent = _plan_and_intent("local_compute:bulk_rnaseq_analysis", "bulk_rnaseq_analysis")
    approval = _approval(intent, plan)
    fabric, _ = _local_fabric(ledger, runner=timing_out_runner)
    req = _request(intent, approval)
    run = fabric.run("s", intent, approval, req)
    assert run.status == "unknown"  # ambiguous; recorded, not retried blindly
    # the request was persisted before execute, so a later re-query can re-attach
    assert ledger.request_by_idempotency_key("s", req.idempotency_key) is not None


def test_experimental_requires_approval(ledger):
    plan, intent = _plan_and_intent("mock_experimental:protein_expression", "protein_expression")
    fabric = _mock_fabric(ledger)
    from backend.orchestration.invariants import InvariantViolation

    with pytest.raises(InvariantViolation):
        fabric.run("s", intent, None, _request(intent, None))


def test_mock_experimental_with_approval_yields_measurements(ledger):
    plan, intent = _plan_and_intent("mock_experimental:protein_expression", "protein_expression")
    approval = _approval(intent, plan)
    fabric = _mock_fabric(ledger)
    run = fabric.run("s", intent, approval, _request(intent, approval))
    assert run.status == "succeeded"
    assert run.outputs and run.outputs[0].uri.startswith("mock://")
