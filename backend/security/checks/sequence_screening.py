"""
Sequence screening adapter + check.

``SequenceScreeningAdapter.screen(sequences) -> ScreeningResult`` is the
extension point for real screening services (IBBIS Common Mechanism,
SecureDNA, synthesis-vendor APIs). Two adapters ship:

* ``MockScreeningAdapter`` (``provider="mock"``) — deterministic, for tests
  and demos; accepts everything unless a sequence contains a configured
  marker, which it flags.
* ``NotConfiguredAdapter`` (``provider="not_configured"``) — returns
  ``not_run``; the check then yields ``REQUIRE_REVIEW`` whenever screening
  is required. **No adapter never means ALLOW.**

Adapter selection: ``HELIX_SCREENING_ADAPTER=not_configured|mock`` (default
``not_configured``). Real adapters register by module path
``HELIX_SCREENING_ADAPTER=module:Class``.
"""

from __future__ import annotations

import importlib
import os
import re
from typing import Iterable, List, Optional, Protocol, runtime_checkable

from backend.contracts.rationale import RationaleItem
from backend.contracts.scientific_objective import ScientificObjective
from backend.contracts.scientific_plan import ScientificPlan
from backend.contracts.security_assessment import RequiredApproval, ScreeningResult
from backend.security.checks.base import AssessmentContext, CheckResult

ADAPTER_ENV = "HELIX_SCREENING_ADAPTER"
CHECK_VERSION = "2026.10"
_SEQ_RE = re.compile(r"\b[ACGTUN]{40,}\b", re.IGNORECASE)


@runtime_checkable
class SequenceScreeningAdapter(Protocol):
    provider: str

    def screen(self, sequences: List[str]) -> ScreeningResult: ...


class NotConfiguredAdapter:
    provider = "not_configured"

    def screen(self, sequences: List[str]) -> ScreeningResult:
        return ScreeningResult(provider=self.provider, status="not_run", summary="no sequence screening adapter configured")


class MockScreeningAdapter:
    """Deterministic test adapter. Flags sequences containing ``flag_marker``."""

    provider = "mock"

    def __init__(self, flag_marker: str = "TTTTTTTTTTTTTTTTTTTT"):
        self.flag_marker = flag_marker.upper()

    def screen(self, sequences: List[str]) -> ScreeningResult:
        flagged = [s for s in sequences if self.flag_marker in s.upper()]
        if flagged:
            return ScreeningResult(provider=self.provider, status="flagged", summary=f"{len(flagged)} sequence(s) matched mock marker", reference_id="mock-flag")
        return ScreeningResult(provider=self.provider, status="accepted", summary=f"{len(sequences)} sequence(s) accepted by mock screener", reference_id="mock-ok")


def load_screening_adapter(spec: Optional[str] = None) -> SequenceScreeningAdapter:
    spec = (spec if spec is not None else os.getenv(ADAPTER_ENV, "not_configured")).strip()
    if spec in ("", "not_configured", "none"):
        return NotConfiguredAdapter()
    if spec == "mock":
        return MockScreeningAdapter()
    if ":" in spec:
        module_name, class_name = spec.split(":", 1)
        adapter = getattr(importlib.import_module(module_name), class_name)()
        if not isinstance(adapter, SequenceScreeningAdapter):
            raise TypeError(f"{spec} does not implement SequenceScreeningAdapter")
        return adapter
    raise ValueError(f"unknown screening adapter {spec!r}; use not_configured | mock | module:Class")


def extract_sequences(plan: ScientificPlan, extra: Iterable[str] = ()) -> List[str]:
    found: List[str] = []
    for step in plan.steps:
        for value in step.arguments.values():
            values = value if isinstance(value, list) else [value]
            for v in values:
                if isinstance(v, str):
                    found.extend(m.group(0) for m in _SEQ_RE.finditer(v))
    found.extend(s for s in extra if isinstance(s, str) and s)
    return sorted(set(found))


class SequenceScreeningCheck:
    """Runs only when the combined envelope says ``screening_required`` (the assessor decides)."""

    check_id = "sequence_screening"
    version = CHECK_VERSION

    def __init__(self, adapter: Optional[SequenceScreeningAdapter] = None):
        self.adapter = adapter or load_screening_adapter()

    def applies_to(self, plan: ScientificPlan) -> bool:
        return True

    def evaluate(self, plan: ScientificPlan, objective: ScientificObjective, context: AssessmentContext) -> CheckResult:
        sequences = extract_sequences(plan, context.sequences)
        if not sequences:
            # Screening is required but there is nothing concrete to screen: an adapter
            # "accepting" an empty list must not open the gate.
            return CheckResult(
                check_id=self.check_id,
                version=self.version,
                outcome="REQUIRE_REVIEW",
                rationale=[RationaleItem(statement="screening required but no sequence found in plan or context", conclusion="REQUIRE_REVIEW (nothing to screen)", evidence_refs=["screening:none:not_run"])],
                screening_results=[ScreeningResult(provider=getattr(self.adapter, "provider", "unknown"), status="not_run", summary="no sequences to screen")],
                required_approvals=[RequiredApproval(role="security_reviewer", reason="screening required but no sequence to screen")],
            )
        try:
            result = self.adapter.screen(sequences)
        except Exception as exc:  # adapter failure is never ALLOW
            result = ScreeningResult(provider=getattr(self.adapter, "provider", "unknown"), status="error", summary=str(exc)[:500])
        evidence = [f"screening:{result.provider}:{result.status}"] + ([f"screening_ref:{result.reference_id}"] if result.reference_id else [])
        if result.status == "accepted":
            outcome, conclusion, approvals = "ALLOW", "screening accepted", []
        elif result.status == "rejected":
            outcome, conclusion, approvals = "DENY", "screening rejected", []
        else:  # flagged | not_run | error → human
            outcome, conclusion = "REQUIRE_REVIEW", f"screening {result.status}; routing to manual review"
            approvals = [RequiredApproval(role="security_reviewer", reason=conclusion)]
        return CheckResult(
            check_id=self.check_id,
            version=self.version,
            outcome=outcome,
            rationale=[RationaleItem(statement=f"{len(sequences)} sequence(s) submitted to '{result.provider}'", conclusion=conclusion, evidence_refs=evidence)],
            screening_results=[result],
            required_approvals=approvals,
        )
