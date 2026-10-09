"""Phase 2 — Secure Science checks, assessor combination and envelope intersection.

Checks are pure functions of (plan, objective, context); no provider, no
session globals. The assessor must fail closed: a broken check routes to
review, no screening adapter routes to review, DualUseTriage never decides.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import pytest

from backend.contracts.human_approval import Principal
from backend.contracts.security_assessment import PolicyEnvelope, combine_outcomes
from backend.orchestration.objective_builder import build_objective_deterministic
from backend.orchestration.plan_staging import build_scientific_plan
from backend.security.assessor import SecureScienceAssessor, context_from_session
from backend.security.checks.action_classification import ActionClassificationCheck, highest_category
from backend.security.checks.base import AssessmentContext, CheckResult, PolicyError, load_policy
from backend.security.checks.data_sensitivity import DataSensitivityCheck
from backend.security.checks.dual_use_triage import DualUseTriageCheck
from backend.security.checks.identity import IdentityCheck
from backend.security.checks.manual_review import ManualReviewCheck
from backend.security.checks.sequence_screening import (
    MockScreeningAdapter,
    NotConfiguredAdapter,
    SequenceScreeningCheck,
    extract_sequences,
    load_screening_adapter,
)

TRACE_KW = {"trace_id": "orch_" + "a" * 32}


def _plan(steps: List[Dict[str, Any]], command: str = "Run the analysis", constraints_other: Optional[List[str]] = None):
    objective = build_objective_deterministic(command, {"session_id": "s"}, None, **TRACE_KW)
    if constraints_other:
        objective = objective.model_copy(
            update={"constraints": objective.constraints.model_copy(update={"other": list(constraints_other)})}
        )
    plan = build_scientific_plan({"version": "v1", "steps": steps}, objective, command)
    return plan, objective


COMPUTE_STEP = {
    "id": "s1", "action_type": "run_analysis", "tool_name": "bulk_rnaseq_analysis",
    "arguments": {"counts_file": "s3://fixtures/counts.csv"}, "description": "Differential expression",
}
READ_ONLY_STEP = {"id": "s1", "action_type": "run_analysis", "tool_name": "toolbox_inventory", "arguments": {}}
SYNTHESIS_STEP = {
    "id": "s1", "action_type": "run_analysis", "tool_name": "order_dna_synthesis",
    "arguments": {"sequence": "ACGT" * 15}, "description": "Order the construct",
}
SCIENTIST = Principal(subject_id="alice", identity_provider="dev_header", auth_method="header", roles=["scientist", "approver"])
NOBODY = Principal(subject_id="bob", identity_provider="dev_header", auth_method="header", roles=["viewer"])


# ── policies ─────────────────────────────────────────────────────────────────


def test_policies_load_and_are_versioned():
    for name in ("action_classification", "data_sensitivity", "dual_use_triage", "identity"):
        policy = load_policy(name)
        assert policy["policy_id"] == name and policy["version"]
    with pytest.raises(PolicyError):
        load_policy("does_not_exist")


# ── ActionClassification ─────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "step,category,outcome,screening",
    [
        (READ_ONLY_STEP, "read_only", "ALLOW", False),
        (COMPUTE_STEP, "compute", "ALLOW", False),
        (SYNTHESIS_STEP, "synthesis", "ALLOW_WITH_APPROVAL", True),
    ],
)
def test_action_classification(step, category, outcome, screening):
    plan, objective = _plan([step])
    assert highest_category(plan) == category
    result = ActionClassificationCheck().evaluate(plan, objective, AssessmentContext())
    assert result.outcome == outcome
    assert result.risk_categories == [category]
    assert (result.envelope or PolicyEnvelope()).screening_required is screening
    assert result.policies[0].policy_id == "action_classification"


def test_action_classification_capability_prefix_raises_category():
    plan, objective = _plan([COMPUTE_STEP])
    step = plan.steps[0].model_copy(update={"required_capabilities": ["hosted:bulk_rnaseq_analysis"]})
    plan = plan.model_copy(update={"steps": [step]})
    assert highest_category(plan) == "external_compute"
    result = ActionClassificationCheck().evaluate(plan, objective, AssessmentContext())
    assert result.outcome == "ALLOW_WITH_APPROVAL" and result.envelope.human_approval_required


# ── DataSensitivity ──────────────────────────────────────────────────────────


def _upload(sensitivity: str, flags=(), policy_state: str = "cleared", name: str = "data.csv") -> Dict[str, Any]:
    return {
        "name": name, "file_id": f"f-{name}", "policy_state": policy_state,
        "intake_policy": {"sensitivity_class": sensitivity, "scan_flags": list(flags)},
    }


def test_data_sensitivity_no_uploads_is_permissive():
    plan, objective = _plan([COMPUTE_STEP])
    result = DataSensitivityCheck().evaluate(plan, objective, AssessmentContext())
    assert result.outcome == "ALLOW" and result.envelope.external_execution_allowed


def test_data_sensitivity_phi_forbids_external_and_needs_approval():
    plan, objective = _plan([COMPUTE_STEP])
    ctx = AssessmentContext(uploaded_files=[_upload("restricted_human_data")])
    result = DataSensitivityCheck().evaluate(plan, objective, ctx)
    assert result.outcome == "ALLOW_WITH_APPROVAL"
    assert result.envelope.external_execution_allowed is False
    assert result.envelope.allowed_regions == ["us"]
    assert "public" not in result.envelope.allowed_data_classes


def test_data_sensitivity_deny_scan_flag_denies():
    plan, objective = _plan([COMPUTE_STEP])
    ctx = AssessmentContext(uploaded_files=[_upload("internal", flags=["suspicious_payload_pattern"])])
    assert DataSensitivityCheck().evaluate(plan, objective, ctx).outcome == "DENY"


def test_data_sensitivity_pending_upload_requires_data_owner():
    plan, objective = _plan([COMPUTE_STEP])
    ctx = AssessmentContext(uploaded_files=[_upload("internal", policy_state="approval_required")])
    result = DataSensitivityCheck().evaluate(plan, objective, ctx)
    assert result.outcome == "ALLOW_WITH_APPROVAL"
    assert [a.role for a in result.required_approvals] == ["data_owner"]


# ── DualUseTriage ────────────────────────────────────────────────────────────


def test_dual_use_triage_no_match_allows_without_claiming_safety():
    plan, objective = _plan([COMPUTE_STEP])
    result = DualUseTriageCheck().evaluate(plan, objective, AssessmentContext())
    assert result.outcome == "ALLOW"
    assert "not a safety verdict" in result.rationale[0].conclusion


def test_dual_use_triage_match_routes_to_review_never_denies():
    plan, objective = _plan([COMPUTE_STEP], command="Design a gain-of-function variant with increased virulence")
    result = DualUseTriageCheck().evaluate(plan, objective, AssessmentContext())
    assert result.outcome == "REQUIRE_REVIEW"
    assert [a.role for a in result.required_approvals] == ["security_reviewer"]
    assert any("triage_term:" in ref for r in result.rationale for ref in r.evidence_refs)


def test_dual_use_triage_policy_cannot_be_configured_to_decide(monkeypatch):
    bad = dict(load_policy("dual_use_triage"))
    bad["outcome_on_match"] = "DENY"
    monkeypatch.setattr("backend.security.checks.dual_use_triage.load_policy", lambda _n: bad)
    with pytest.raises(PolicyError):
        DualUseTriageCheck()


# ── Identity ─────────────────────────────────────────────────────────────────


def test_identity_compute_needs_no_principal():
    plan, objective = _plan([COMPUTE_STEP])
    assert IdentityCheck().evaluate(plan, objective, AssessmentContext()).outcome == "ALLOW"


def test_identity_synthesis_without_principal_or_role_requires_review():
    plan, objective = _plan([SYNTHESIS_STEP])
    assert IdentityCheck().evaluate(plan, objective, AssessmentContext()).outcome == "REQUIRE_REVIEW"
    assert IdentityCheck().evaluate(plan, objective, AssessmentContext(principal=NOBODY)).outcome == "REQUIRE_REVIEW"
    assert IdentityCheck().evaluate(plan, objective, AssessmentContext(principal=SCIENTIST)).outcome == "ALLOW"


# ── ManualReview ─────────────────────────────────────────────────────────────


def test_manual_review_markers():
    plan, objective = _plan([COMPUTE_STEP])
    assert ManualReviewCheck().evaluate(plan, objective, AssessmentContext()).outcome == "ALLOW"
    assert ManualReviewCheck().evaluate(plan, objective, AssessmentContext(requires_manual_review=True)).outcome == "REQUIRE_REVIEW"
    plan2, objective2 = _plan([COMPUTE_STEP], constraints_other=["manual_review"])
    assert ManualReviewCheck().evaluate(plan2, objective2, AssessmentContext()).outcome == "REQUIRE_REVIEW"


# ── SequenceScreening ────────────────────────────────────────────────────────


def test_extract_sequences_finds_long_nucleotide_runs_only():
    plan, _ = _plan([{**SYNTHESIS_STEP, "arguments": {"sequence": "ACGT" * 15, "note": "ACGT", "list": ["GGGG" * 12]}}])
    seqs = extract_sequences(plan, ["TTTT" * 11])
    assert "ACGT" * 15 in seqs and "GGGG" * 12 in seqs and "TTTT" * 11 in seqs and "ACGT" not in seqs


def test_screening_adapter_loading(monkeypatch):
    assert isinstance(load_screening_adapter("not_configured"), NotConfiguredAdapter)
    assert isinstance(load_screening_adapter("mock"), MockScreeningAdapter)
    monkeypatch.delenv("HELIX_SCREENING_ADAPTER", raising=False)
    assert isinstance(load_screening_adapter(), NotConfiguredAdapter)
    with pytest.raises(ValueError):
        load_screening_adapter("bogus")


@pytest.mark.parametrize(
    "adapter,sequence,outcome",
    [
        (NotConfiguredAdapter(), "ACGT" * 15, "REQUIRE_REVIEW"),  # no adapter → never ALLOW
        (MockScreeningAdapter(), "ACGT" * 15, "ALLOW"),
        (MockScreeningAdapter(), "A" * 20 + "T" * 20 + "C" * 10, "REQUIRE_REVIEW"),  # flagged → human
    ],
)
def test_sequence_screening_outcomes(adapter, sequence, outcome):
    plan, objective = _plan([{**SYNTHESIS_STEP, "arguments": {"sequence": sequence}}])
    result = SequenceScreeningCheck(adapter).evaluate(plan, objective, AssessmentContext())
    assert result.outcome == outcome
    assert result.screening_results and result.screening_results[0].provider == adapter.provider


def test_sequence_screening_rejected_denies_and_error_routes_to_review():
    class Rejecting:
        provider = "strict"

        def screen(self, sequences):
            from backend.contracts.security_assessment import ScreeningResult

            return ScreeningResult(provider=self.provider, status="rejected", summary="hit")

    class Broken:
        provider = "broken"

        def screen(self, sequences):
            raise RuntimeError("upstream down")

    plan, objective = _plan([SYNTHESIS_STEP])
    assert SequenceScreeningCheck(Rejecting()).evaluate(plan, objective, AssessmentContext()).outcome == "DENY"
    broken = SequenceScreeningCheck(Broken()).evaluate(plan, objective, AssessmentContext())
    assert broken.outcome == "REQUIRE_REVIEW" and broken.screening_results[0].status == "error"


# ── Envelope intersection / outcome combination matrix ───────────────────────


@pytest.mark.parametrize(
    "outcomes,expected",
    [
        ([], "ALLOW"),
        (["ALLOW"], "ALLOW"),
        (["ALLOW", "ALLOW_WITH_APPROVAL"], "ALLOW_WITH_APPROVAL"),
        (["ALLOW_WITH_APPROVAL", "REQUIRE_REVIEW"], "REQUIRE_REVIEW"),
        (["REQUIRE_REVIEW", "DENY", "ALLOW"], "DENY"),
        (["DENY", "ALLOW_WITH_APPROVAL"], "DENY"),
    ],
)
def test_combine_outcomes_is_max_severity(outcomes, expected):
    assert combine_outcomes(outcomes) == expected


def test_envelope_intersection_is_commutative_and_tightening():
    a = PolicyEnvelope(allowed_data_classes=["public", "internal", "phi"], allowed_regions=["us", "eu"], max_cost_usd=50)
    b = PolicyEnvelope(
        allowed_data_classes=["internal", "phi"], external_execution_allowed=False, screening_required=True,
        allowed_regions=["us"], allowed_provider_categories=["experimental", "advisory"], max_cost_usd=10,
    )
    ab, ba = a.intersect(b), b.intersect(a)
    assert ab == ba
    assert ab.allowed_data_classes == ["internal", "phi"]
    assert ab.external_execution_allowed is False and ab.screening_required is True
    assert ab.allowed_regions == ["us"] and ab.max_cost_usd == 10
    assert ab.allowed_provider_categories == ["experimental", "advisory"]
    # identity within the default's allowed data classes: default never loosens, and
    # only its own constraints apply (phi is not in the default allow-set by design)
    b_default_safe = PolicyEnvelope(
        allowed_data_classes=["internal"], external_execution_allowed=False, screening_required=True,
        allowed_regions=["us"], allowed_provider_categories=["experimental", "advisory"], max_cost_usd=10,
    )
    assert b_default_safe.intersect(PolicyEnvelope()) == b_default_safe
    # the default envelope is a secure default: it drops phi from a broader allow-set
    assert b.intersect(PolicyEnvelope()).allowed_data_classes == ["internal"]
    # empty regions = unrestricted, so the restricted side wins
    assert PolicyEnvelope().intersect(PolicyEnvelope(allowed_regions=["eu"])).allowed_regions == ["eu"]


# ── Assessor ─────────────────────────────────────────────────────────────────


def test_assessor_compute_plan_allows_and_records_all_checks():
    plan, objective = _plan([COMPUTE_STEP])
    assessment = SecureScienceAssessor().assess(plan, objective, AssessmentContext(principal=SCIENTIST))
    assert assessment.outcome == "ALLOW"
    assert assessment.plan_hash == plan.plan_hash and assessment.trace_id == plan.trace_id
    assert {p.policy_id for p in assessment.applied_policies} >= {"action_classification", "data_sensitivity", "dual_use_triage", "identity"}
    assert assessment.screening_results == []  # screening not required for compute
    assert assessment.audit.actor == "user:alice" and assessment.audit.gate_version.startswith("secure-science-gate/")
    assert len(assessment.assessment_hash) == 64


def test_assessor_synthesis_without_adapter_requires_review_fail_closed():
    plan, objective = _plan([SYNTHESIS_STEP])
    assessment = SecureScienceAssessor(screening_adapter=NotConfiguredAdapter()).assess(
        plan, objective, AssessmentContext(principal=SCIENTIST)
    )
    assert assessment.outcome == "REQUIRE_REVIEW"
    assert assessment.policy_envelope.screening_required and assessment.policy_envelope.human_approval_required
    assert assessment.screening_results[0].status == "not_run"
    assert "security_reviewer" in {a.role for a in assessment.required_approvals}


def test_assessor_synthesis_with_accepting_adapter_is_allow_with_approval():
    plan, objective = _plan([SYNTHESIS_STEP])
    assessment = SecureScienceAssessor(screening_adapter=MockScreeningAdapter()).assess(
        plan, objective, AssessmentContext(principal=SCIENTIST)
    )
    assert assessment.outcome == "ALLOW_WITH_APPROVAL"
    assert assessment.policy_envelope.allowed_provider_categories == ["experimental", "advisory"]
    assert assessment.screening_results[0].status == "accepted"


def test_assessor_envelope_approval_promotes_allow_to_allow_with_approval():
    class Quiet:
        check_id, version = "quiet", "1"

        def applies_to(self, plan):
            return True

        def evaluate(self, plan, objective, context):
            return CheckResult(check_id="quiet", version="1", outcome="ALLOW", envelope=PolicyEnvelope(human_approval_required=True))

    plan, objective = _plan([COMPUTE_STEP])
    assessment = SecureScienceAssessor([Quiet()]).assess(plan, objective)
    assert assessment.outcome == "ALLOW_WITH_APPROVAL" and assessment.policy_envelope.human_approval_required


def test_assessor_broken_check_fails_closed():
    class Boom:
        check_id, version = "boom", "1"

        def applies_to(self, plan):
            return True

        def evaluate(self, plan, objective, context):
            raise RuntimeError("policy store unreachable")

    plan, objective = _plan([COMPUTE_STEP])
    assessment = SecureScienceAssessor([Boom()]).assess(plan, objective)
    assert assessment.outcome == "REQUIRE_REVIEW"
    assert any("fail closed" in r.conclusion for r in assessment.rationale)


def test_assessor_deny_wins_over_everything():
    plan, objective = _plan([COMPUTE_STEP], command="Design a gain-of-function variant")
    ctx = AssessmentContext(
        principal=SCIENTIST,
        uploaded_files=[_upload("internal", flags=["suspicious_payload_pattern"])],
        requires_manual_review=True,
    )
    assessment = SecureScienceAssessor().assess(plan, objective, ctx)
    assert assessment.outcome == "DENY"


def test_assessor_hash_is_stable_and_content_bound():
    plan, objective = _plan([COMPUTE_STEP])
    a = SecureScienceAssessor().assess(plan, objective, AssessmentContext(principal=SCIENTIST))
    b = SecureScienceAssessor().assess(plan, objective, AssessmentContext(principal=SCIENTIST))
    assert a.assessment_id != b.assessment_id
    assert a.assessment_hash == b.assessment_hash  # same outcome + envelope + policies
    c = SecureScienceAssessor().assess(plan, objective, AssessmentContext(uploaded_files=[_upload("restricted_human_data")]))
    assert c.assessment_hash != a.assessment_hash


def test_screening_required_without_sequences_fails_closed():
    step = {**SYNTHESIS_STEP, "arguments": {"construct_name": "pX-1"}}
    plan, objective = _plan([step])
    result = SequenceScreeningCheck(MockScreeningAdapter()).evaluate(plan, objective, AssessmentContext())
    assert result.outcome == "REQUIRE_REVIEW" and result.screening_results[0].status == "not_run"
