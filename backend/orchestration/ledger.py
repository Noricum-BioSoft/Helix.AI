"""
Platform ledger — the single write path for shared domain records (Phase 1).

Objective, plan, execution intent and approval records are written here and
nowhere else. Phase 7 routes these ``record_*`` helpers through
``KnowledgeStore`` (local | dataweaver | dual); until then the local JSON
ledger under ``sessions/{session_id}/platform/`` is the only backend.

Layout::

    sessions/{sid}/platform/objectives/{objective_id}.v{n}.json
    sessions/{sid}/platform/plans/{plan_id}.v{n}.json
    sessions/{sid}/platform/intents/{execution_intent_id}.json
    sessions/{sid}/platform/approvals/{approval_id}.json

IDs are minted by the callers at the domain layer (``backend.contracts.ids``);
the ledger never mints IDs for shared objects.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Type, TypeVar

from pydantic import BaseModel

from backend.contracts.execution_intent import ExecutionIntent
from backend.contracts.execution_request import ExecutionRequest, ExecutionRun
from backend.contracts.human_approval import HumanApproval
from backend.contracts.scientific_objective import ScientificObjective
from backend.contracts.scientific_plan import ScientificPlan
from backend.contracts.security_assessment import SecurityAssessment
from backend.contracts.security_review import SecurityReview

T = TypeVar("T", bound=BaseModel)

PLATFORM_DIR = "platform"
KIND_OBJECTIVE = "objectives"
KIND_PLAN = "plans"
KIND_INTENT = "intents"
KIND_APPROVAL = "approvals"
KIND_ASSESSMENT = "assessments"
KIND_REVIEW = "security_reviews"
KIND_REQUEST = "execution_requests"
KIND_RUN = "execution_runs"
AUDIT_FILE = "audit/events.jsonl"

RECORD_KINDS = (
    (KIND_OBJECTIVE, ScientificObjective),
    (KIND_PLAN, ScientificPlan),
    (KIND_ASSESSMENT, SecurityAssessment),
    (KIND_REVIEW, SecurityReview),
    (KIND_INTENT, ExecutionIntent),
    (KIND_APPROVAL, HumanApproval),
    (KIND_REQUEST, ExecutionRequest),
    (KIND_RUN, ExecutionRun),
)


class LedgerError(RuntimeError):
    pass


class LocalLedger:
    """File-backed ledger scoped to a sessions storage directory."""

    def __init__(self, storage_dir: Path):
        self.storage_dir = Path(storage_dir)

    # ── paths ────────────────────────────────────────────────────────────────

    def _dir(self, session_id: str, kind: str) -> Path:
        return self.storage_dir / session_id / PLATFORM_DIR / kind

    @staticmethod
    def _file_name(record_id: str, version: Optional[int]) -> str:
        return f"{record_id}.v{version}.json" if version is not None else f"{record_id}.json"

    # ── generic write/read ───────────────────────────────────────────────────

    def _write(self, session_id: str, kind: str, record_id: str, version: Optional[int], model: BaseModel) -> Path:
        target_dir = self._dir(session_id, kind)
        target_dir.mkdir(parents=True, exist_ok=True)
        target = target_dir / self._file_name(record_id, version)
        tmp = target.with_suffix(".json.tmp")
        tmp.write_text(model.model_dump_json(indent=2), encoding="utf-8")
        os.replace(tmp, target)
        return target

    def _read(self, session_id: str, kind: str, record_id: str, version: Optional[int], model: Type[T]) -> Optional[T]:
        path = self._dir(session_id, kind) / self._file_name(record_id, version)
        if not path.exists():
            return None
        return model.model_validate_json(path.read_text(encoding="utf-8"))

    def _iter(self, session_id: str, kind: str, model: Type[T]) -> Iterator[T]:
        target_dir = self._dir(session_id, kind)
        if not target_dir.exists():
            return
        for path in sorted(target_dir.glob("*.json")):
            try:
                yield model.model_validate_json(path.read_text(encoding="utf-8"))
            except Exception:  # corrupt file must not hide the rest of the ledger
                continue

    # ── objectives ───────────────────────────────────────────────────────────

    def record_objective(self, session_id: str, objective: ScientificObjective) -> Path:
        return self._write(session_id, KIND_OBJECTIVE, objective.objective_id, objective.version, objective)

    def load_objective(self, session_id: str, objective_id: str, version: int = 1) -> Optional[ScientificObjective]:
        return self._read(session_id, KIND_OBJECTIVE, objective_id, version, ScientificObjective)

    # ── plans ────────────────────────────────────────────────────────────────

    def record_plan(self, session_id: str, plan: ScientificPlan) -> Path:
        return self._write(session_id, KIND_PLAN, plan.plan_id, plan.version, plan)

    def load_plan(self, session_id: str, plan_id: str, version: int = 1) -> Optional[ScientificPlan]:
        return self._read(session_id, KIND_PLAN, plan_id, version, ScientificPlan)

    def list_plans(self, session_id: str) -> List[ScientificPlan]:
        return list(self._iter(session_id, KIND_PLAN, ScientificPlan))

    # ── execution intents ────────────────────────────────────────────────────

    def record_intent(self, session_id: str, intent: ExecutionIntent) -> Path:
        existing = self.load_intent(session_id, intent.execution_intent_id)
        if existing is not None and existing.execution_intent_hash != intent.execution_intent_hash:
            raise LedgerError(f"ExecutionIntent {intent.execution_intent_id} is immutable; refusing to overwrite")
        return self._write(session_id, KIND_INTENT, intent.execution_intent_id, None, intent)

    def load_intent(self, session_id: str, execution_intent_id: str) -> Optional[ExecutionIntent]:
        return self._read(session_id, KIND_INTENT, execution_intent_id, None, ExecutionIntent)

    # ── approvals ────────────────────────────────────────────────────────────

    def record_approval(self, session_id: str, approval: HumanApproval) -> Path:
        if self.load_approval(session_id, approval.approval_id) is not None:
            raise LedgerError(f"HumanApproval {approval.approval_id} already recorded")
        return self._write(session_id, KIND_APPROVAL, approval.approval_id, None, approval)

    def load_approval(self, session_id: str, approval_id: str) -> Optional[HumanApproval]:
        return self._read(session_id, KIND_APPROVAL, approval_id, None, HumanApproval)

    def approvals_for_intent(self, session_id: str, execution_intent_id: str) -> List[HumanApproval]:
        return [
            a for a in self._iter(session_id, KIND_APPROVAL, HumanApproval) if a.execution_intent_id == execution_intent_id
        ]

    # ── security assessments / reviews (P2) ─────────────────────────────────

    def record_assessment(self, session_id: str, assessment: SecurityAssessment) -> Path:
        existing = self.load_assessment(session_id, assessment.assessment_id)
        if existing is not None and existing.assessment_hash != assessment.assessment_hash:
            raise LedgerError(f"SecurityAssessment {assessment.assessment_id} is immutable; refusing to overwrite")
        return self._write(session_id, KIND_ASSESSMENT, assessment.assessment_id, None, assessment)

    def load_assessment(self, session_id: str, assessment_id: str) -> Optional[SecurityAssessment]:
        return self._read(session_id, KIND_ASSESSMENT, assessment_id, None, SecurityAssessment)

    def record_review(self, session_id: str, review: SecurityReview) -> Path:
        if self.load_review(session_id, review.review_id) is not None:
            raise LedgerError(f"SecurityReview {review.review_id} already recorded")
        return self._write(session_id, KIND_REVIEW, review.review_id, None, review)

    def load_review(self, session_id: str, review_id: str) -> Optional[SecurityReview]:
        return self._read(session_id, KIND_REVIEW, review_id, None, SecurityReview)

    def reviews_for_assessment(self, session_id: str, assessment_id: str) -> List[SecurityReview]:
        return [r for r in self._iter(session_id, KIND_REVIEW, SecurityReview) if r.assessment_id == assessment_id]

    # ── execution requests / runs (P3A) ──────────────────────────────────────

    def record_request(self, session_id: str, request: ExecutionRequest) -> Path:
        """Persist an ExecutionRequest. Immutable once written (idempotency anchor)."""
        existing = self.load_request(session_id, request.execution_request_id)
        if existing is not None and existing.idempotency_key != request.idempotency_key:
            raise LedgerError(f"ExecutionRequest {request.execution_request_id} is immutable; refusing to overwrite")
        return self._write(session_id, KIND_REQUEST, request.execution_request_id, None, request)

    def load_request(self, session_id: str, execution_request_id: str) -> Optional[ExecutionRequest]:
        return self._read(session_id, KIND_REQUEST, execution_request_id, None, ExecutionRequest)

    def request_by_idempotency_key(self, session_id: str, idempotency_key: str) -> Optional[ExecutionRequest]:
        """The existing submission for a key, if any — the heart of fail-safe retries."""
        for req in self._iter(session_id, KIND_REQUEST, ExecutionRequest):
            if req.idempotency_key == idempotency_key:
                return req
        return None

    def record_run(self, session_id: str, run: ExecutionRun) -> Path:
        return self._write(session_id, KIND_RUN, run.execution_run_id, None, run)

    def load_run(self, session_id: str, execution_run_id: str) -> Optional[ExecutionRun]:
        return self._read(session_id, KIND_RUN, execution_run_id, None, ExecutionRun)

    def run_for_request(self, session_id: str, execution_request_id: str) -> Optional[ExecutionRun]:
        for run in self._iter(session_id, KIND_RUN, ExecutionRun):
            if run.execution_request_id == execution_request_id:
                return run
        return None

    # ── audit events (append-only JSONL, one line per event) ─────────────────

    def record_audit_event(self, session_id: str, event: Dict) -> Path:
        target = self.storage_dir / session_id / PLATFORM_DIR / AUDIT_FILE
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(event, default=str, sort_keys=True) + "\n")
        return target

    def audit_events(self, session_id: str, trace_id: Optional[str] = None) -> List[Dict]:
        target = self.storage_dir / session_id / PLATFORM_DIR / AUDIT_FILE
        if not target.exists():
            return []
        events = [json.loads(line) for line in target.read_text(encoding="utf-8").splitlines() if line.strip()]
        return [e for e in events if trace_id is None or e.get("trace_id") == trace_id]

    # ── trace queries ────────────────────────────────────────────────────────

    def records_by_trace(self, session_id: str, trace_id: str) -> Dict[str, List[BaseModel]]:
        """All records of a loop, grouped by kind. Used by tests and (P7) the trace endpoint."""
        out: Dict[str, List[BaseModel]] = {}
        for kind, model in RECORD_KINDS:
            out[kind] = [r for r in self._iter(session_id, kind, model) if r.trace_id == trace_id]
        return out


def get_ledger() -> LocalLedger:
    """Ledger bound to the live ``history_manager`` storage dir (tests monkeypatch that dir)."""
    from backend.history_manager import history_manager

    return LocalLedger(Path(history_manager.storage_dir))
