"""
SecureScienceAssessor — runs the checks and combines them into one SecurityAssessment.

Two passes: all non-screening checks first, their envelopes intersected;
sequence screening runs only when the combined envelope says
``screening_required``. Outcome = max severity; envelope = intersection;
``assessment_hash`` over outcome + envelope + applied policies.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Sequence

from backend.contracts.human_approval import Principal
from backend.contracts.ids import new_id
from backend.contracts.rationale import RationaleItem
from backend.contracts.scientific_objective import ScientificObjective
from backend.contracts.scientific_plan import ScientificPlan
from backend.contracts.security_assessment import (
    AppliedPolicy,
    AssessmentAudit,
    PolicyEnvelope,
    RequiredApproval,
    ScreeningResult,
    SecurityAssessment,
    combine_outcomes,
)
from backend.security.checks.action_classification import ActionClassificationCheck
from backend.security.checks.base import AssessmentContext, CheckResult, SecurityCheck
from backend.security.checks.data_sensitivity import DataSensitivityCheck
from backend.security.checks.dual_use_triage import DualUseTriageCheck
from backend.security.checks.identity import IdentityCheck
from backend.security.checks.manual_review import ManualReviewCheck
from backend.security.checks.sequence_screening import SequenceScreeningAdapter, SequenceScreeningCheck

logger = logging.getLogger(__name__)

GATE_VERSION = "secure-science-gate/2.0"


def default_checks() -> List[SecurityCheck]:
    return [ActionClassificationCheck(), DataSensitivityCheck(), DualUseTriageCheck(), IdentityCheck(), ManualReviewCheck()]


class SecureScienceAssessor:
    def __init__(
        self,
        checks: Optional[Sequence[SecurityCheck]] = None,
        *,
        screening_adapter: Optional[SequenceScreeningAdapter] = None,
        screening_check: Optional[SequenceScreeningCheck] = None,
        gate_version: str = GATE_VERSION,
        profile: Optional[str] = None,
    ):
        self.checks: List[SecurityCheck] = list(checks) if checks is not None else default_checks()
        self.screening = screening_check or SequenceScreeningCheck(screening_adapter)
        self.gate_version = gate_version
        self.profile = profile

    def assess(self, plan: ScientificPlan, objective: ScientificObjective, context: Optional[AssessmentContext] = None) -> SecurityAssessment:
        context = context or AssessmentContext()
        results: List[CheckResult] = []
        for check in self.checks:
            if not check.applies_to(plan):
                continue
            try:
                results.append(check.evaluate(plan, objective, context))
            except Exception as exc:  # a broken check must not open the gate
                logger.exception("security check %s failed", getattr(check, "check_id", check))
                results.append(
                    CheckResult(
                        check_id=str(getattr(check, "check_id", type(check).__name__)),
                        version=str(getattr(check, "version", "?")),
                        outcome="REQUIRE_REVIEW",
                        rationale=[RationaleItem(statement=f"check raised {type(exc).__name__}", conclusion="REQUIRE_REVIEW (fail closed)")],
                        required_approvals=[RequiredApproval(role="security_reviewer", reason=f"check {getattr(check, 'check_id', '?')} errored")],
                    )
                )

        envelope = PolicyEnvelope()
        for r in results:
            if r.envelope is not None:
                envelope = envelope.intersect(r.envelope)

        if envelope.screening_required:
            results.append(self.screening.evaluate(plan, objective, context))

        outcome = combine_outcomes([r.outcome for r in results])
        if outcome == "ALLOW_WITH_APPROVAL" or envelope.human_approval_required:
            envelope = envelope.model_copy(update={"human_approval_required": True})
            if outcome == "ALLOW":
                outcome = "ALLOW_WITH_APPROVAL"

        rationale: List[RationaleItem] = []
        policies: List[AppliedPolicy] = []
        screening: List[ScreeningResult] = []
        approvals: List[RequiredApproval] = []
        for r in results:
            rationale.extend(r.rationale)
            policies.extend(r.policies)
            screening.extend(r.screening_results)
            approvals.extend(r.required_approvals)
        rationale.append(
            RationaleItem(
                statement=f"{len(results)} check(s): " + ", ".join(f"{r.check_id}={r.outcome}" for r in results),
                conclusion=f"combined outcome {outcome}",
                evidence_refs=[f"check:{r.check_id}@{r.version}" for r in results],
                confidence=1.0,
            )
        )
        actor = f"user:{context.principal.subject_id}" if context.principal else "system"
        assessment = SecurityAssessment(
            assessment_id=new_id(),
            trace_id=plan.trace_id,
            plan_id=plan.plan_id,
            plan_hash=plan.plan_hash,
            outcome=outcome,
            rationale=rationale,
            applied_policies=policies,
            screening_results=screening,
            required_approvals=approvals,
            policy_envelope=envelope,
            audit=AssessmentAudit(actor=actor, gate_version=self.gate_version, profile=self.profile or context.profile),
        )
        logger.info("[security] trace=%s plan=%s outcome=%s envelope=%s", plan.trace_id, plan.plan_id, outcome, envelope.model_dump(mode="json"))
        return assessment


def context_from_session(
    session: Optional[Dict[str, Any]],
    *,
    session_id: Optional[str] = None,
    principal: Optional[Principal] = None,
    profile: Optional[str] = None,
) -> AssessmentContext:
    """Build the context from a history_manager session dict (metadata.uploaded_files)."""
    session = session or {}
    uploads: List[Dict[str, Any]] = []
    meta_uploads = (session.get("metadata") or {}).get("uploaded_files") or []
    uploads.extend(u for u in meta_uploads if isinstance(u, dict))
    return AssessmentContext(
        session_id=session_id or session.get("session_id"),
        principal=principal,
        uploaded_files=uploads,
        requires_manual_review=bool(session.get("requires_manual_review")),
        profile=profile,
    )


def default_assessor() -> SecureScienceAssessor:
    """Assessor with default checks; screening adapter from ``HELIX_SCREENING_ADAPTER`` (default not configured)."""
    import os

    return SecureScienceAssessor(profile=os.getenv("HELIX_EXECUTION_PROFILE") or None)
