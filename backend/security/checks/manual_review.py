"""Manual review check: an explicit request for review on the context, objective or plan always wins."""

from __future__ import annotations

from backend.contracts.rationale import RationaleItem
from backend.contracts.scientific_objective import ScientificObjective
from backend.contracts.scientific_plan import ScientificPlan
from backend.contracts.security_assessment import RequiredApproval
from backend.security.checks.base import AssessmentContext, CheckResult

MARKERS = {"manual_review", "requires_manual_review", "security_review"}


def _flagged(plan: ScientificPlan, objective: ScientificObjective, context: AssessmentContext) -> str:
    if context.requires_manual_review:
        return "context.requires_manual_review"
    if objective.user_context.get("requires_manual_review"):
        return "objective.user_context.requires_manual_review"
    if MARKERS & {c.strip().lower() for c in objective.constraints.other}:
        return "objective.constraints.other"
    if MARKERS & {r.strip().lower() for r in plan.security_requirements}:
        return "plan.security_requirements"
    return ""


class ManualReviewCheck:
    check_id = "manual_review"
    version = "2026.10"

    def applies_to(self, plan: ScientificPlan) -> bool:
        return True

    def evaluate(self, plan: ScientificPlan, objective: ScientificObjective, context: AssessmentContext) -> CheckResult:
        source = _flagged(plan, objective, context)
        if not source:
            return CheckResult(check_id=self.check_id, version=self.version, outcome="ALLOW",
                               rationale=[RationaleItem(statement="no explicit manual-review request", conclusion="ALLOW")])
        return CheckResult(
            check_id=self.check_id, version=self.version, outcome="REQUIRE_REVIEW",
            rationale=[RationaleItem(statement=f"manual review requested via {source}", conclusion="REQUIRE_REVIEW", evidence_refs=[source])],
            required_approvals=[RequiredApproval(role="security_reviewer", reason="explicitly requested")],
        )
