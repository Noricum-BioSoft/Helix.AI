"""
Action classification: step kind / capability → risk category → envelope.

Categories (``read_only < compute < external_compute < experimental < synthesis``)
come from ``policies/action_classification.yaml``. Experimental and synthesis
steps force ``human_approval_required`` and ``screening_required`` into the
envelope regardless of what later checks say.
"""

from __future__ import annotations

from typing import Dict, List

from backend.contracts.rationale import RationaleItem
from backend.contracts.scientific_objective import ScientificObjective
from backend.contracts.scientific_plan import ScientificPlan, ScientificPlanStep
from backend.contracts.security_assessment import PolicyEnvelope, combine_outcomes
from backend.security.checks.base import AssessmentContext, CheckResult, applied, load_policy

POLICY = "action_classification"


def category_rank(category: str) -> int:
    return load_policy(POLICY)["categories"].index(category)


def classify_step(step: ScientificPlanStep) -> str:
    policy = load_policy(POLICY)
    if step.tool_name and step.tool_name in set(policy.get("synthesis_tools") or []):
        return "synthesis"
    if step.tool_name and step.tool_name in set(policy.get("read_only_tools") or []):
        return "read_only"
    category = policy["kind_to_category"].get(step.kind, "compute")
    for cap in step.required_capabilities:
        for prefix, cap_category in policy["capability_prefix_to_category"].items():
            if cap.startswith(prefix):
                if category_rank(cap_category) > category_rank(category):
                    category = cap_category
                break
    return category


def classify_plan(plan: ScientificPlan) -> Dict[str, str]:
    return {step.id: classify_step(step) for step in plan.steps}


def highest_category(plan: ScientificPlan) -> str:
    cats = list(classify_plan(plan).values())
    return max(cats, key=category_rank) if cats else "read_only"


class ActionClassificationCheck:
    check_id = "action_classification"

    def __init__(self) -> None:
        self.version = str(load_policy(POLICY)["version"])

    def applies_to(self, plan: ScientificPlan) -> bool:
        return bool(plan.steps)

    def evaluate(self, plan: ScientificPlan, objective: ScientificObjective, context: AssessmentContext) -> CheckResult:
        policy = load_policy(POLICY)
        by_step = classify_plan(plan)
        envelope = PolicyEnvelope()
        outcomes = []
        rationale: List[RationaleItem] = []
        for step_id, category in by_step.items():
            cat_policy = policy["category_policy"][category]
            outcomes.append(cat_policy["outcome"])
            envelope = envelope.intersect(PolicyEnvelope(**(cat_policy.get("envelope") or {})))
            rationale.append(
                RationaleItem(
                    statement=f"step {step_id} classified as {category}",
                    conclusion=cat_policy["outcome"],
                    evidence_refs=[f"plan_step:{step_id}", f"policy:{policy['policy_id']}@{policy['version']}"],
                    confidence=1.0,
                )
            )
        outcome = combine_outcomes(outcomes)
        return CheckResult(
            check_id=self.check_id,
            version=self.version,
            outcome=outcome,
            rationale=rationale,
            policies=[applied(policy, outcome)],
            envelope=envelope,
            risk_categories=sorted(set(by_step.values()), key=category_rank),
        )
