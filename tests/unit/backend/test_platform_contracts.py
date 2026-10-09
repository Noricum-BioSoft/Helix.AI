"""Phase 0 — platform contract tests: round-trip, hash stability, structured rationale, schema snapshot."""

from __future__ import annotations

import json

import pytest
from pydantic import ValidationError

from backend.contracts import schema_export
from backend.contracts.evidence_assessment import EvidenceAssessment, Interpretation, NextDecision, Observation
from backend.contracts.execution_intent import ExecutionIntent
from backend.contracts.execution_recommendation import CriteriaScores, ExecutionRecommendation
from backend.contracts.execution_request import ExecutionRequest, make_idempotency_key
from backend.contracts.human_approval import HumanApproval, Principal
from backend.contracts.ids import new_id, new_trace_id, stable_hash
from backend.contracts.provenance import Actor, ProvenanceRecord
from backend.contracts.rationale import RationaleItem
from backend.contracts.scientific_objective import ScientificObjective
from backend.contracts.scientific_plan import ScientificPlan
from backend.contracts.security_assessment import (
    AppliedPolicy,
    AssessmentAudit,
    PolicyEnvelope,
    SecurityAssessment,
    combine_outcomes,
)
from backend.plan_ir import Plan, PlanStep

TRACE = new_trace_id()


def _plan(tool: str = "fastqc_quality_analysis") -> ScientificPlan:
    ir = Plan(steps=[PlanStep(id="s1", tool_name=tool, arguments={"input": "s3://b/reads.fq"})])
    return ScientificPlan.from_plan_ir(
        ir,
        plan_id=new_id(),
        trace_id=TRACE,
        objective_id=new_id(),
        rationale=[RationaleItem(statement="Paired-end reads supplied.", conclusion="QC first.", evidence_refs=["dataset:reads"])],
    )


def _intent(plan: ScientificPlan, provider: str = "legacy:Local", **over) -> ExecutionIntent:
    data = dict(
        execution_intent_id=new_id(),
        trace_id=TRACE,
        plan_id=plan.plan_id,
        plan_hash=plan.plan_hash,
        provider_id=provider,
        capability_id="local_compute:fastqc_quality_analysis",
        input_manifest_hash=stable_hash(["s3://b/reads.fq"]),
        execution_parameters_hash=stable_hash({"input": "s3://b/reads.fq"}),
    )
    data.update(over)
    return ExecutionIntent(**data)


def _principal() -> Principal:
    return Principal(subject_id="alice", identity_provider="dev_header", auth_method="header")


# --- rationale -----------------------------------------------------------------


def test_rationale_item_requires_statement_and_conclusion():
    with pytest.raises(ValidationError):
        RationaleItem(statement="  ", conclusion="x")
    item = RationaleItem(statement="a", conclusion="b", confidence=0.5)
    assert item.evidence_refs == []


def test_no_contract_has_free_text_reasoning_field():
    for name, model in schema_export.CONTRACTS.items():
        assert "reasoning" not in model.model_fields, f"{name} must use rationale: list[RationaleItem]"


# --- plan hash ----------------------------------------------------------------------


def test_plan_hash_stable_across_roundtrip_and_ignores_rationale():
    p1 = _plan()
    p2 = ScientificPlan.model_validate_json(p1.model_dump_json())
    assert p1.plan_hash == p2.plan_hash
    p3 = p1.model_copy(update={"plan_rationale": []})
    assert ScientificPlan.model_validate(p3.model_dump()).plan_hash == p1.plan_hash


def test_plan_hash_changes_with_steps_and_rejects_tampering():
    p1, p2 = _plan("fastqc_quality_analysis"), _plan("bwa_align")
    assert p1.plan_hash != p2.plan_hash
    tampered = p1.model_dump()
    tampered["plan_hash"] = p2.plan_hash
    with pytest.raises(ValidationError, match="plan_hash"):
        ScientificPlan.model_validate(tampered)


def test_from_plan_ir_defaults_capability_to_local_compute():
    assert _plan().steps[0].required_capabilities == ["local_compute:fastqc_quality_analysis"]


# --- execution intent ----------------------------------------------------------------


def test_intent_hash_stable_and_frozen():
    plan = _plan()
    i1 = _intent(plan)
    i2 = ExecutionIntent.model_validate_json(i1.model_dump_json())
    assert i1.execution_intent_hash == i2.execution_intent_hash
    with pytest.raises(ValidationError):
        i1.provider_id = "other"  # type: ignore[misc]


@pytest.mark.parametrize(
    "field,value",
    [
        ("provider_id", "carolina_cloud"),
        ("capability_id", "nextflow:nf-core/rnaseq"),
        ("input_manifest_hash", stable_hash(["s3://other"])),
        ("execution_parameters_hash", stable_hash({"x": 1})),
        ("assessment_id", "a1"),
    ],
)
def test_intent_hash_changes_when_execution_risk_changes(field, value):
    plan = _plan()
    base = _intent(plan)
    changed = _intent(plan, **{field: value})
    assert base.execution_intent_hash != changed.execution_intent_hash


def test_intent_rejects_mismatched_hash():
    plan = _plan()
    data = _intent(plan).model_dump()
    data["execution_intent_hash"] = "0" * 64
    with pytest.raises(ValidationError, match="execution_intent_hash"):
        ExecutionIntent.model_validate(data)


# --- approval -------------------------------------------------------------------------


def test_approval_binds_to_intent_hash_and_is_frozen():
    plan = _plan()
    intent = _intent(plan)
    approval = HumanApproval(
        approval_id=new_id(),
        trace_id=TRACE,
        execution_intent_id=intent.execution_intent_id,
        execution_intent_hash=intent.execution_intent_hash,
        plan_id=plan.plan_id,
        plan_hash=plan.plan_hash,
        decision="approved",
        principal=_principal(),
    )
    assert approval.approved and intent.matches(approval.execution_intent_hash)
    with pytest.raises(ValidationError):
        approval.decision = "rejected"  # type: ignore[misc]


# --- security assessment ------------------------------------------------------------------


def test_combine_outcomes_severity_order():
    assert combine_outcomes(["ALLOW", "ALLOW_WITH_APPROVAL"]) == "ALLOW_WITH_APPROVAL"
    assert combine_outcomes(["REQUIRE_REVIEW", "ALLOW_WITH_APPROVAL"]) == "REQUIRE_REVIEW"
    assert combine_outcomes(["ALLOW", "DENY", "REQUIRE_REVIEW"]) == "DENY"
    assert combine_outcomes([]) == "ALLOW"


def test_policy_envelope_intersection_is_conservative():
    a = PolicyEnvelope(allowed_data_classes=["public", "internal"], allowed_regions=["us", "eu"], max_cost_usd=10)
    b = PolicyEnvelope(allowed_data_classes=["internal"], external_execution_allowed=False, screening_required=True, max_cost_usd=5)
    c = a.intersect(b)
    assert c.allowed_data_classes == ["internal"]
    assert c.external_execution_allowed is False and c.screening_required is True
    assert c.allowed_regions == ["us", "eu"]  # empty on one side = unrestricted
    assert c.max_cost_usd == 5


def test_assessment_hash_covers_envelope_and_policies():
    plan = _plan()
    common = dict(assessment_id=new_id(), trace_id=TRACE, plan_id=plan.plan_id, plan_hash=plan.plan_hash,
                  outcome="ALLOW_WITH_APPROVAL", audit=AssessmentAudit(actor="system", gate_version="0"))
    a1 = SecurityAssessment(**common, applied_policies=[AppliedPolicy(policy_id="p", version="1", result="ALLOW")])
    a2 = SecurityAssessment(**common, applied_policies=[AppliedPolicy(policy_id="p", version="2", result="ALLOW")])
    a3 = SecurityAssessment(**common, policy_envelope=PolicyEnvelope(screening_required=True))
    assert len({a1.assessment_hash, a2.assessment_hash, a3.assessment_hash}) == 3
    assert SecurityAssessment.model_validate_json(a1.model_dump_json()).assessment_hash == a1.assessment_hash


# --- recommendation -----------------------------------------------------------------------


def _scores(**over) -> CriteriaScores:
    base = {k: 0.5 for k in CriteriaScores.model_fields}
    base.update(over)
    return CriteriaScores(**base)


def test_recommendation_provider_must_be_in_candidate_set():
    plan = _plan()
    kwargs = dict(recommendation_id=new_id(), trace_id=TRACE, plan_id=plan.plan_id, capability_id="c",
                  criteria_scores=_scores(), score=0.5, confidence=0.8)
    ExecutionRecommendation(provider_id="local_compute", candidate_set=["local_compute", "nextflow"], **kwargs)
    with pytest.raises(ValidationError, match="candidate_set"):
        ExecutionRecommendation(provider_id="aws_emr", candidate_set=["local_compute"], **kwargs)


# --- execution request / idempotency ------------------------------------------------------------


def test_idempotency_key_deterministic_for_same_intent_and_nonce():
    plan = _plan()
    intent = _intent(plan)
    k1 = make_idempotency_key(intent.execution_intent_hash, nonce="n")
    k2 = make_idempotency_key(intent.execution_intent_hash, nonce="n")
    k3 = make_idempotency_key(intent.execution_intent_hash)
    assert k1 == k2 and k1 != k3 and len(k1) == 64
    req = ExecutionRequest(execution_request_id=new_id(), trace_id=TRACE, execution_intent_id=intent.execution_intent_id,
                           execution_intent_hash=intent.execution_intent_hash, idempotency_key=k1,
                           provider_id=intent.provider_id, capability_id=intent.capability_id)
    assert req.submitted_at is None


# --- evidence / decision ----------------------------------------------------------------------------


def test_decision_requires_evidence_and_proposed_objective_stays_proposed():
    decision_id = new_id()
    proposed = ScientificObjective(objective_id=new_id(), trace_id=TRACE, objective="o", question="q",
                                   status="proposed", parent_decision_id=decision_id)
    NextDecision(decision_id=decision_id, trace_id=TRACE, kind="repeat_experiment",
                 rationale=[RationaleItem(statement="s", conclusion="c")], evidence_ids=["obs1"], proposed_objective=proposed)
    with pytest.raises(ValidationError):
        NextDecision(decision_id=decision_id, trace_id=TRACE, kind="stop", rationale=[], evidence_ids=["obs1"])
    with pytest.raises(ValidationError, match="proposed"):
        NextDecision(decision_id=decision_id, trace_id=TRACE, kind="stop",
                     rationale=[RationaleItem(statement="s", conclusion="c")], evidence_ids=["obs1"],
                     proposed_objective=proposed.model_copy(update={"status": "active"}))


def test_observation_is_a_measured_fact():
    obs = Observation(observation_id=new_id(), trace_id=TRACE, metric="expression_yield", value=12.5, unit="mg/L")
    assert "rationale" not in obs.model_fields and "conclusion" not in obs.model_fields
    EvidenceAssessment(objective_id=new_id(), trace_id=TRACE, observations=[obs], interpretation=Interpretation(),
                       decision=NextDecision(decision_id=new_id(), trace_id=TRACE, kind="stop",
                                             rationale=[RationaleItem(statement="s", conclusion="c")], evidence_ids=[obs.observation_id]))


# --- provenance / trace ---------------------------------------------------------------------------------


def test_provenance_record_versions_and_trace_format():
    rec = ProvenanceRecord(event_id=new_id(), trace_id=TRACE, event_type="plan_created", actor=Actor(kind="agent", id="planner"),
                           prompt_template_id="planner_v3", prompt_template_hash="abc", code_commit="deadbeef")
    assert rec.schema_version == 1 and "prompt_text" not in rec.model_fields
    with pytest.raises(ValidationError, match="trace_id"):
        ProvenanceRecord(event_id=new_id(), trace_id="not-a-trace", event_type="x", actor=Actor(kind="system", id="s"))


# --- schema snapshot ----------------------------------------------------------------------------------------


def test_contract_schema_snapshot_is_current():
    problems = schema_export.stale()
    assert not problems, f"run `python -m backend.contracts.schema_export` — stale: {problems}"


def test_every_contract_schema_is_valid_json_with_id():
    for name in schema_export.CONTRACTS:
        schema = json.loads(schema_export.render(name))
        assert schema["$id"].endswith(f"{name}.json")
        assert "properties" in schema
