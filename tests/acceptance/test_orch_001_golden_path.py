"""
ORCH-001 — Scientific Orchestration Golden Path (acceptance test).

Created in Phase 0 as a skeleton: every step is skipped with the phase that
will implement it. The skip list may only shrink; by Phase 9 nothing is
skipped and the whole loop runs on deterministic fixtures in CI.

Flow: Objective → persisted Plan → Security Assessment (envelope) →
Capability resolution → Recommendation → Provider Authorization →
ExecutionIntent → Human Approval → ExecutionProvider → Result →
persist to DataWeaver → Observation → Interpretation → Decision.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.orch001

# Steps and the phase that un-skips them. Keep this list in sync with the plan.
STEPS = [
    ("objective_created", "P1"),
    ("plan_persisted_with_hash", "P1"),
    ("execution_intent_built", "P1"),
    ("approval_bound_to_intent_hash", "P1"),
    ("security_assessment_with_envelope", "P2"),
    ("denied_state_blocks_intent", "P2"),
    ("capability_resolved_from_registry", "P3A"),
    ("provider_executes_local", "P3A"),
    ("provider_executes_mock_experimental_with_approval", "P3A"),
    ("retry_after_timeout_creates_no_duplicate", "P3A"),
    ("recommendation_with_candidate_set", "P4"),
    ("provider_authorization_against_envelope", "P4"),
    ("intent_provider_in_candidate_set", "P4"),
    ("one_real_computational_backend", "P3B/P3A-nextflow"),
    ("persisted_to_dataweaver", "P7"),
    ("ids_stable_across_stores", "P7"),
    ("trace_query_returns_full_loop", "P7"),
    ("observation_has_no_claims", "P8"),
    ("decision_references_evidence_ids", "P8"),
    ("proposed_objective_not_auto_activated", "P8"),
    ("provenance_from_decision_to_inputs", "P8"),
]


@pytest.mark.parametrize("step,phase", STEPS, ids=[s for s, _ in STEPS])
def test_orch_001_step(step: str, phase: str) -> None:
    pytest.skip(f"ORCH-001 step '{step}' is implemented in {phase}")


def test_orch_001_contracts_importable() -> None:
    """The only non-skipped assertion in Phase 0: the contract surface exists."""
    from backend.contracts import (  # noqa: F401
        evidence_assessment,
        execution_intent,
        execution_recommendation,
        execution_request,
        human_approval,
        provenance,
        provider_authorization,
        scientific_objective,
        scientific_plan,
        security_assessment,
    )
    from backend.execution.providers.base import ExecutionProvider  # noqa: F401
    from backend.orchestration import invariants  # noqa: F401
    from shared.capability_registry import CapabilityDescriptor  # noqa: F401
