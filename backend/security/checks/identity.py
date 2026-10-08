"""
Identity check: a principal is present and its roles are allowed for the plan's risk category.

Chat/dev usage without a principal is fine for ``read_only``/``compute``
(``policies/identity.yaml: principal_required_from_category``); anything
riskier without a principal routes to review.
"""

from __future__ import annotations

from typing import List

from backend.contracts.rationale import RationaleItem
from backend.contracts.scientific_objective import ScientificObjective
from backend.contracts.scientific_plan import ScientificPlan
from backend.contracts.security_assessment import RequiredApproval
from backend.security.checks.action_classification import category_rank, highest_category
from backend.security.checks.base import AssessmentContext, CheckResult, applied, load_policy

POLICY = "identity"


def review_roles() -> List[str]:
    return list(load_policy(POLICY).get("review_roles") or [])


class IdentityCheck:
    check_id = "identity"

    def __init__(self) -> None:
        self.version = str(load_policy(POLICY)["version"])

    def applies_to(self, plan: ScientificPlan) -> bool:
        return True

    def evaluate(self, plan: ScientificPlan, objective: ScientificObjective, context: AssessmentContext) -> CheckResult:
        policy = load_policy(POLICY)
        category = highest_category(plan)
        threshold = str(policy["principal_required_from_category"])
        evidence = [f"policy:{policy['policy_id']}@{policy['version']}", f"risk_category:{category}"]
        if category_rank(category) < category_rank(threshold):
            return CheckResult(
                check_id=self.check_id, version=self.version, outcome="ALLOW",
                rationale=[RationaleItem(statement=f"category {category} does not require a principal", conclusion="ALLOW", evidence_refs=evidence)],
                policies=[applied(policy, "ALLOW")],
            )
        principal = context.principal
        if principal is None:
            return CheckResult(
                check_id=self.check_id, version=self.version, outcome="REQUIRE_REVIEW",
                rationale=[RationaleItem(statement=f"category {category} requires an authenticated principal; none present", conclusion="REQUIRE_REVIEW", evidence_refs=evidence)],
                policies=[applied(policy, "REQUIRE_REVIEW")],
                required_approvals=[RequiredApproval(role="security_reviewer", reason="no authenticated principal for a controlled action")],
            )
        required = list((policy.get("roles_required") or {}).get(category) or [])
        if required and not set(required) & set(principal.roles):
            return CheckResult(
                check_id=self.check_id, version=self.version, outcome="REQUIRE_REVIEW",
                rationale=[RationaleItem(statement=f"principal '{principal.subject_id}' lacks a required role for {category} ({', '.join(required)})", conclusion="REQUIRE_REVIEW", evidence_refs=evidence + [f"principal:{principal.identity_provider}"])],
                policies=[applied(policy, "REQUIRE_REVIEW")],
                required_approvals=[RequiredApproval(role=required[0], reason=f"role required for {category}")],
            )
        return CheckResult(
            check_id=self.check_id, version=self.version, outcome="ALLOW",
            rationale=[RationaleItem(statement=f"principal '{principal.subject_id}' ({principal.identity_provider}) allowed for {category}", conclusion="ALLOW", evidence_refs=evidence)],
            policies=[applied(policy, "ALLOW")],
        )
