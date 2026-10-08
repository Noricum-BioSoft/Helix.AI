"""
Data sensitivity: reuse upload-intake results to constrain where data may go.

Reads ``intake_policy.sensitivity_class`` / ``scan_flags`` / ``policy_state``
from the session's upload metadata and maps them through
``policies/data_sensitivity.yaml`` to ``allowed_data_classes``,
``external_execution_allowed`` and ``allowed_regions``.
"""

from __future__ import annotations

from typing import Any, Dict, List

from backend.contracts.rationale import RationaleItem
from backend.contracts.scientific_objective import ScientificObjective
from backend.contracts.scientific_plan import ScientificPlan
from backend.contracts.security_assessment import PolicyEnvelope, RequiredApproval, combine_outcomes
from backend.security.checks.base import AssessmentContext, CheckResult, applied, load_policy

POLICY = "data_sensitivity"


def data_class_of(upload: Dict[str, Any]) -> str:
    policy = load_policy(POLICY)
    intake = upload.get("intake_policy") or {}
    sensitivity = str(intake.get("sensitivity_class") or upload.get("sensitivity_class") or "internal")
    return policy["sensitivity_to_data_class"].get(sensitivity, "sensitive")


class DataSensitivityCheck:
    check_id = "data_sensitivity"

    def __init__(self) -> None:
        self.version = str(load_policy(POLICY)["version"])

    def applies_to(self, plan: ScientificPlan) -> bool:
        return True  # with no uploads it contributes the permissive default

    def evaluate(self, plan: ScientificPlan, objective: ScientificObjective, context: AssessmentContext) -> CheckResult:
        policy = load_policy(POLICY)
        deny_flags = set(policy.get("deny_scan_flags") or [])
        pending_states = set(policy.get("pending_policy_states") or [])
        outcomes = []
        envelope = PolicyEnvelope(allowed_data_classes=["public", "internal", "sensitive", "phi", "restricted"])
        rationale: List[RationaleItem] = []
        approvals: List[RequiredApproval] = []

        if not context.uploaded_files:
            outcomes.append("ALLOW")
            rationale.append(
                RationaleItem(statement="no session uploads", conclusion="ALLOW", evidence_refs=[f"policy:{policy['policy_id']}@{policy['version']}"])
            )

        for upload in context.uploaded_files:
            name = str(upload.get("name") or upload.get("filename") or upload.get("file_id") or "upload")
            intake = upload.get("intake_policy") or {}
            flags = set(intake.get("scan_flags") or upload.get("scan_flags") or [])
            data_class = data_class_of(upload)
            cls_policy = policy["data_class_policy"][data_class]
            outcome = cls_policy["outcome"]
            if flags & deny_flags:
                outcome = "DENY"
            elif str(upload.get("policy_state") or "") in pending_states:
                outcome = combine_outcomes([outcome, "ALLOW_WITH_APPROVAL"])
                approvals.append(RequiredApproval(role="data_owner", reason=f"upload '{name}' awaits policy approval"))
            outcomes.append(outcome)
            envelope = envelope.intersect(PolicyEnvelope(**(cls_policy.get("envelope") or {})))
            rationale.append(
                RationaleItem(
                    statement=f"upload '{name}' has data class {data_class}" + (f"; scan flags {sorted(flags)}" if flags else ""),
                    conclusion=outcome,
                    evidence_refs=[f"upload:{upload.get('file_id') or name}", f"policy:{policy['policy_id']}@{policy['version']}"],
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
            required_approvals=approvals,
            envelope=envelope,
        )
