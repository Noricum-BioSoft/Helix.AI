"""Phase 1.2 — plan staging + ledger: versions, hashes, supersedes, structured rationale, one trace_id."""

from __future__ import annotations

import copy
import json

import pytest

from backend.contracts.rationale import RationaleItem
from backend.orchestration.ledger import LedgerError, LocalLedger
from backend.orchestration.plan_staging import (
    build_scientific_plan,
    plan_ir_from_dict,
    rationale_from_plan,
    stage_plan,
)


@pytest.fixture(autouse=True)
def _mock_mode(monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")


@pytest.fixture()
def ledger(tmp_path):
    return LocalLedger(tmp_path / "sessions")


ROUTER_PLAN = {
    "version": "v1",
    "steps": [
        {
            "id": "step1",
            "action_type": "run_analysis",
            "tool_name": "fastqc_quality_analysis",
            "arguments": {
                "session_id": "s1",
                "fastq_file": "s3://bucket/reads.fastq.gz",
                "router_reasoning": {"deterministic_path": "fastqc_quality_analysis", "intent": "single-tool", "suggested_steps": []},
            },
            "description": "Run FastQC on my reads",
        }
    ],
}

TABULAR_PLAN = {
    "type": "tabular_analysis",
    "title": "Correlation",
    "goal": "Which numeric columns co-vary?",
    "steps": [
        {"id": 1, "name": "Load", "description": "Load the table", "type": "load", "operation": "read_csv"},
        {"id": 2, "name": "Correlate", "description": "Pearson", "type": "compute", "operation": "Pearson correlation matrix"},
        {"id": 3, "name": "Interpret", "description": "Summarise", "type": "interpret", "operation": "Summarise top pairs"},
    ],
    "expected_outputs": ["heatmap"],
}


def test_plan_ir_from_router_plan_and_tabular_plan():
    ir = plan_ir_from_dict(ROUTER_PLAN)
    assert ir.steps[0].tool_name == "fastqc_quality_analysis"
    ir_t = plan_ir_from_dict(TABULAR_PLAN)
    assert [s.id for s in ir_t.steps] == ["1", "2", "3"]
    assert {s.tool_name for s in ir_t.steps} == {"tabular_analysis"}
    with pytest.raises(ValueError):
        plan_ir_from_dict({"steps": []})


def test_rationale_is_structured_items_not_prose():
    items = rationale_from_plan(ROUTER_PLAN, "Run FastQC")
    assert items and all(isinstance(i, RationaleItem) for i in items)
    assert items[0].evidence_refs == ["router:fastqc_quality_analysis"]
    assert items[0].confidence == 1.0  # deterministic route
    t_items = rationale_from_plan(TABULAR_PLAN, "x")
    assert t_items[0].statement == "Which numeric columns co-vary?"
    assert [i.evidence_refs for i in t_items[1:]] == [["plan_step:1"], ["plan_step:2"], ["plan_step:3"]]


def test_stage_plan_writes_objective_plan_intent_with_one_trace(ledger):
    staged = stage_plan("s1", "Run FastQC on my reads", ROUTER_PLAN, {"session_id": "s1"}, ledger=ledger)
    assert staged.objective.trace_id == staged.plan.trace_id == staged.intent.trace_id
    assert staged.plan.objective_id == staged.objective.objective_id
    assert staged.intent.plan_hash == staged.plan.plan_hash
    assert staged.plan.version == 1 and staged.plan.status == "draft"
    assert staged.intent.provider_id == "legacy:Local"
    assert staged.intent.capability_id == "local_compute:fastqc_quality_analysis"

    by_trace = ledger.records_by_trace("s1", staged.trace_id)
    assert {k: len(v) for k, v in by_trace.items()} == {"objectives": 1, "plans": 1, "intents": 1, "approvals": 0}
    plan_file = ledger.storage_dir / "s1" / "platform" / "plans" / f"{staged.plan.plan_id}.v1.json"
    assert plan_file.exists()
    on_disk = json.loads(plan_file.read_text())
    assert on_disk["plan_hash"] == staged.plan.plan_hash
    assert on_disk["plan_rationale"][0]["statement"]
    assert "reasoning" not in on_disk

    fields = staged.checkpoint_fields()
    assert fields["pending_plan_hash"] == staged.plan.plan_hash
    assert fields["pending_execution_intent_hash"] == staged.intent.execution_intent_hash
    assert fields["approval_id"] is None


def test_plan_hash_excludes_rationale_but_covers_steps():
    from backend.orchestration.objective_builder import build_objective

    obj = build_objective("x")
    a = build_scientific_plan(ROUTER_PLAN, obj, "x")
    b = build_scientific_plan(ROUTER_PLAN, obj, "a different command phrasing")
    assert a.plan_hash == b.plan_hash  # rationale/command wording does not matter
    changed = copy.deepcopy(ROUTER_PLAN)
    changed["steps"][0]["arguments"]["fastq_file"] = "s3://bucket/other.fastq.gz"
    c = build_scientific_plan(changed, obj, "x")
    assert c.plan_hash != a.plan_hash


def test_revision_creates_new_version_with_supersedes(ledger):
    v1 = stage_plan("s1", "Run FastQC", ROUTER_PLAN, ledger=ledger)
    changed = copy.deepcopy(ROUTER_PLAN)
    changed["steps"][0]["arguments"]["fastq_file"] = "s3://bucket/other.fastq.gz"
    v2 = stage_plan("s1", "Run FastQC on the other file", changed, supersedes=v1.plan, ledger=ledger)
    assert v2.plan.plan_id == v1.plan.plan_id
    assert v2.plan.version == 2
    assert v2.plan.supersedes_plan_id == v1.plan.plan_id
    assert v2.plan.plan_hash != v1.plan.plan_hash
    assert v2.objective.objective_id == v1.objective.objective_id
    assert v2.trace_id == v1.trace_id
    assert v2.intent.execution_intent_id != v1.intent.execution_intent_id
    assert ledger.load_plan("s1", v1.plan.plan_id, 1) is not None
    assert ledger.load_plan("s1", v1.plan.plan_id, 2).version == 2


def test_tabular_plan_step_kinds_and_expected_outputs(ledger):
    staged = stage_plan("s1", "correlate", TABULAR_PLAN, ledger=ledger)
    kinds = {s.id: s.kind for s in staged.plan.steps}
    assert kinds == {"1": "computational", "2": "computational", "3": "reasoning"}
    assert staged.plan.expected_outputs == ["heatmap"]
    assert staged.intent.capability_id == "local_compute:tabular_analysis"


def test_intent_is_immutable_in_ledger(ledger):
    staged = stage_plan("s1", "Run FastQC", ROUTER_PLAN, ledger=ledger)
    other = staged.intent.model_copy(update={"provider_id": "legacy:EMR", "execution_intent_hash": ""})
    # model_copy skips validation; rebuild so the hash is recomputed for the new provider
    from backend.contracts.execution_intent import ExecutionIntent

    other = ExecutionIntent(**{**other.model_dump(), "execution_intent_hash": ""})
    assert other.execution_intent_id == staged.intent.execution_intent_id
    with pytest.raises(LedgerError):
        ledger.record_intent("s1", other)
    # idempotent re-record of the identical intent is fine
    ledger.record_intent("s1", staged.intent)


def test_intent_creation_blocked_in_security_states(ledger):
    from backend.orchestration.invariants import InvariantViolation

    with pytest.raises(InvariantViolation):
        stage_plan("s1", "x", ROUTER_PLAN, current_state="WAITING_FOR_SECURITY_REVIEW", ledger=ledger)
    with pytest.raises(InvariantViolation):
        stage_plan("s1", "x", ROUTER_PLAN, current_state="DENIED", ledger=ledger)
