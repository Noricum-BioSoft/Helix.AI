"""
DualUseTriageCheck — conservative routing to manual review.

This is **not a biological-risk classifier and not a substitute for
specialised screening**. It matches a curated term list
(``policies/dual_use_triage.yaml``, derived from ``docs/SAFETY_POLICY.md``)
against the objective and plan text and, on any match, routes the loop to a
human reviewer (``REQUIRE_REVIEW``). It never emits ``DENY``,
``ALLOW_WITH_APPROVAL``, or any safe/unsafe label; absence of a match is not
a safety statement. Context terms (vaccine, diagnostic, …) are recorded as
evidence for the reviewer, never used to cancel a match.
"""

from __future__ import annotations

from typing import Dict, List

from backend.contracts.rationale import RationaleItem
from backend.contracts.scientific_objective import ScientificObjective
from backend.contracts.scientific_plan import ScientificPlan
from backend.contracts.security_assessment import RequiredApproval, SecurityOutcome
from backend.security.checks.base import AssessmentContext, CheckResult, PolicyError, applied, load_policy, plan_text

POLICY = "dual_use_triage"
ALLOWED_OUTCOMES = {"ALLOW", "REQUIRE_REVIEW"}


def triage_matches(text: str) -> Dict[str, List[str]]:
    policy = load_policy(POLICY)
    lowered = text.lower()
    hits: Dict[str, List[str]] = {}
    for category, spec in (policy.get("categories") or {}).items():
        matched = [t for t in (spec.get("terms") or []) if str(t).lower() in lowered]
        if matched:
            hits[category] = matched
    return hits


def context_terms_present(text: str) -> List[str]:
    policy = load_policy(POLICY)
    lowered = text.lower()
    return [t for t in (policy.get("context_terms") or []) if str(t).lower() in lowered]


class DualUseTriageCheck:
    check_id = "dual_use_triage"

    def __init__(self) -> None:
        policy = load_policy(POLICY)
        self.version = str(policy["version"])
        self.outcome_on_match: SecurityOutcome = str(policy.get("outcome_on_match") or "REQUIRE_REVIEW")  # type: ignore[assignment]
        if self.outcome_on_match not in ALLOWED_OUTCOMES:
            raise PolicyError("dual_use_triage may only route to REQUIRE_REVIEW (or ALLOW); it is not a classifier")

    def applies_to(self, plan: ScientificPlan) -> bool:
        return True

    def evaluate(self, plan: ScientificPlan, objective: ScientificObjective, context: AssessmentContext) -> CheckResult:
        policy = load_policy(POLICY)
        text = plan_text(plan, objective)
        hits = triage_matches(text)
        evidence = [f"policy:{policy['policy_id']}@{policy['version']}"]
        if not hits:
            return CheckResult(
                check_id=self.check_id,
                version=self.version,
                outcome="ALLOW",
                rationale=[RationaleItem(statement="no triage terms matched", conclusion="ALLOW (triage only; not a safety verdict)", evidence_refs=evidence)],
                policies=[applied(policy, "ALLOW")],
            )
        outcome = self.outcome_on_match
        rationale = [
            RationaleItem(
                statement=f"triage category '{category}' matched {len(terms)} term(s)",
                conclusion=f"{outcome} (routing to manual review; not a risk classification)",
                evidence_refs=evidence + [f"triage_term:{t}" for t in terms],
            )
            for category, terms in hits.items()
        ]
        ctx = context_terms_present(text)
        if ctx:
            rationale.append(
                RationaleItem(
                    statement=f"context terms present for the reviewer: {', '.join(ctx)}",
                    conclusion="recorded as evidence only; does not change the outcome",
                    evidence_refs=[f"context_term:{t}" for t in ctx],
                )
            )
        assert outcome in ALLOWED_OUTCOMES
        return CheckResult(
            check_id=self.check_id,
            version=self.version,
            outcome=outcome,
            rationale=rationale,
            policies=[applied(policy, outcome)],
            required_approvals=[RequiredApproval(role="security_reviewer", reason=f"dual-use triage matched: {', '.join(sorted(hits))}")],
        )
