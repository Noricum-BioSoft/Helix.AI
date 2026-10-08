"""
Plan staging (Phase 1.2 + 1.3).

Called from the staging branch of ``/execute`` when ``HELIX_SCIENCE_GATE_V1``
is on. Converts the legacy plan dict (plan IR, or the tabular analysis plan)
into persisted platform records sharing one ``trace_id``:

    ScientificObjective → ScientificPlan (hashed, versioned) → ExecutionIntent

and returns the identifiers the ``WorkflowCheckpoint`` carries. The legacy
``pending_plan`` dict remains authoritative for execution; these records are
what approval binds to.

Rationale is structured (``RationaleItem``), derived from the router's
structured fields or the tabular planner's ``goal``/``steps`` — never from
persisted model prose.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from backend.contracts.execution_intent import ExecutionIntent
from backend.contracts.ids import new_id
from backend.contracts.rationale import RationaleItem
from backend.contracts.scientific_objective import ScientificObjective
from backend.contracts.scientific_plan import ScientificPlan, ScientificPlanStep, StepKind
from backend.orchestration.intent_builder import build_intent
from backend.orchestration.ledger import LocalLedger, get_ledger
from backend.orchestration.objective_builder import build_objective
from backend.plan_ir import Plan as PlanIR
from backend.plan_ir import PlanStep as IRStep

logger = logging.getLogger(__name__)

TABULAR_PLAN_TYPE = "tabular_analysis"
TABULAR_TOOL_NAME = "tabular_analysis"
_TABULAR_STEP_KINDS: Dict[str, StepKind] = {
    "load": "computational",
    "filter": "computational",
    "compute": "computational",
    "visualize": "computational",
    "interpret": "reasoning",
}


@dataclass(frozen=True)
class StagedRecords:
    objective: ScientificObjective
    plan: ScientificPlan
    intent: ExecutionIntent

    @property
    def trace_id(self) -> str:
        return self.plan.trace_id

    def checkpoint_fields(self) -> Dict[str, Optional[str]]:
        return {
            "trace_id": self.trace_id,
            "objective_id": self.objective.objective_id,
            "pending_plan_id": self.plan.plan_id,
            "pending_plan_hash": self.plan.plan_hash,
            "pending_execution_intent_id": self.intent.execution_intent_id,
            "pending_execution_intent_hash": self.intent.execution_intent_hash,
            "approval_id": None,
        }


# ── plan dict → IR / steps ───────────────────────────────────────────────────


def _is_tabular_plan(plan_dict: Dict[str, Any]) -> bool:
    return plan_dict.get("type") == TABULAR_PLAN_TYPE or (
        "goal" in plan_dict and any("operation" in s for s in plan_dict.get("steps", []) if isinstance(s, dict))
    )


def plan_ir_from_dict(plan_dict: Dict[str, Any]) -> PlanIR:
    """Normalise a legacy plan dict to ``plan_ir.Plan``. Tabular plans map each step to the tabular tool."""
    steps: List[IRStep] = []
    tabular = _is_tabular_plan(plan_dict)
    for index, raw in enumerate(plan_dict.get("steps") or []):
        if not isinstance(raw, dict):
            continue
        step_id = str(raw.get("id") or f"step{index + 1}")
        if tabular:
            steps.append(
                IRStep(
                    id=step_id,
                    action_type=str(raw.get("type") or "compute"),
                    tool_name=TABULAR_TOOL_NAME,
                    arguments={k: raw[k] for k in ("name", "operation", "type") if k in raw},
                    description=raw.get("description") or raw.get("name"),
                )
            )
        else:
            steps.append(
                IRStep(
                    id=step_id,
                    action_type=str(raw.get("action_type") or "execute_tool"),
                    tool_name=raw.get("tool_name") or raw.get("tool"),
                    arguments=dict(raw.get("arguments") or raw.get("params") or {}),
                    description=raw.get("description"),
                )
            )
    if not steps:
        raise ValueError("plan has no steps")
    return PlanIR(version=str(plan_dict.get("version") or "v1"), steps=steps)


def _step_kinds(plan_dict: Dict[str, Any], ir: PlanIR) -> Dict[str, StepKind]:
    kinds: Dict[str, StepKind] = {}
    if _is_tabular_plan(plan_dict):
        for step in ir.steps:
            kinds[step.id] = _TABULAR_STEP_KINDS.get(step.action_type, "computational")
    return kinds


# ── structured rationale ─────────────────────────────────────────────────────


def rationale_from_plan(plan_dict: Dict[str, Any], command: str, router_params: Optional[Dict[str, Any]] = None) -> List[RationaleItem]:
    items: List[RationaleItem] = []
    steps = [s for s in plan_dict.get("steps") or [] if isinstance(s, dict)]
    if _is_tabular_plan(plan_dict):
        goal = str(plan_dict.get("goal") or command).strip()
        if goal:
            items.append(
                RationaleItem(
                    statement=goal[:1000],
                    conclusion=f"Tabular analysis plan with {len(steps)} step(s)",
                    evidence_refs=[f"planner:{TABULAR_PLAN_TYPE}"],
                )
            )
        for step in steps:
            op = str(step.get("operation") or step.get("description") or "").strip()
            if op:
                items.append(
                    RationaleItem(
                        statement=op[:1000],
                        conclusion=f"step {step.get('id')}: {step.get('type') or 'compute'}",
                        evidence_refs=[f"plan_step:{step.get('id')}"],
                    )
                )
        return items

    reasoning = (router_params or {}).get("router_reasoning") or {}
    if not reasoning and steps:
        reasoning = (steps[0].get("arguments") or {}).get("router_reasoning") or {}
    tools = sorted({str(s.get("tool_name") or s.get("tool")) for s in steps if s.get("tool_name") or s.get("tool")})
    deterministic = reasoning.get("deterministic_path")
    items.append(
        RationaleItem(
            statement=(
                f"Router matched the request deterministically to '{deterministic}'"
                if deterministic
                else f"Router classified the request as intent '{reasoning.get('intent') or 'execute'}'"
            ),
            conclusion=f"Plan with {len(steps)} step(s) using {', '.join(tools) or 'router placeholder'}",
            evidence_refs=[f"router:{t}" for t in tools],
            confidence=1.0 if deterministic else None,
        )
    )
    for suggested in list(reasoning.get("suggested_steps") or [])[:10]:
        text = str(suggested).strip()
        if text:
            items.append(RationaleItem(statement=text[:500], conclusion="suggested step", evidence_refs=["router:suggested_steps"]))
    return items


# ── public API ───────────────────────────────────────────────────────────────


def build_scientific_plan(
    plan_dict: Dict[str, Any],
    objective: ScientificObjective,
    command: str,
    *,
    router_params: Optional[Dict[str, Any]] = None,
    supersedes: Optional[ScientificPlan] = None,
) -> ScientificPlan:
    ir = plan_ir_from_dict(plan_dict)
    expected_outputs = [str(o) for o in plan_dict.get("expected_outputs") or [] if str(o).strip()]
    return ScientificPlan.from_plan_ir(
        ir,
        plan_id=supersedes.plan_id if supersedes else new_id(),
        trace_id=objective.trace_id,
        objective_id=objective.objective_id,
        rationale=rationale_from_plan(plan_dict, command, router_params),
        step_kinds=_step_kinds(plan_dict, ir),
        version=(supersedes.version + 1) if supersedes else 1,
        supersedes_plan_id=supersedes.plan_id if supersedes else None,
        expected_outputs=expected_outputs,
        status="draft",
    )


def stage_plan(
    session_id: str,
    command: str,
    plan_dict: Dict[str, Any],
    session_context: Optional[Dict[str, Any]] = None,
    *,
    router_params: Optional[Dict[str, Any]] = None,
    infra_decision: Optional[Any] = None,
    current_state: Optional[str] = None,
    supersedes: Optional[ScientificPlan] = None,
    objective: Optional[ScientificObjective] = None,
    ledger: Optional[LocalLedger] = None,
) -> StagedRecords:
    """Persist objective, plan and intent for a staged plan and return their identifiers.

    ``supersedes`` makes this a revision: same ``plan_id``/``trace_id``/objective,
    ``version + 1``. ``objective`` reuses an existing objective (revisions);
    otherwise a new one (and a new ``trace_id``) is minted.
    """
    ledger = ledger or get_ledger()
    if objective is None:
        if supersedes is not None:
            objective = ledger.load_objective(session_id, supersedes.objective_id)
        if objective is None:
            routed_tool = None
            steps = plan_dict.get("steps") or []
            if steps and isinstance(steps[0], dict):
                routed_tool = steps[0].get("tool_name") or steps[0].get("tool")
            objective = build_objective(command, session_context, {"tool": routed_tool} if routed_tool else None)
            ledger.record_objective(session_id, objective)

    plan = build_scientific_plan(plan_dict, objective, command, router_params=router_params, supersedes=supersedes)
    ledger.record_plan(session_id, plan)

    intent = build_intent(plan, infra_decision, current_state=current_state)
    ledger.record_intent(session_id, intent)

    logger.info(
        "[platform] staged trace=%s objective=%s plan=%s v%s intent=%s",
        plan.trace_id, objective.objective_id, plan.plan_id, plan.version, intent.execution_intent_id,
    )
    return StagedRecords(objective=objective, plan=plan, intent=intent)


def try_stage_plan(*args: Any, **kwargs: Any) -> Optional[StagedRecords]:
    """``stage_plan`` that never raises into the request path; failures are logged and yield ``None``.

    With the gate on, a staging failure means no intent exists, so the NL
    approval path refuses to execute (fail closed) rather than running an
    unrecorded plan.
    """
    try:
        return stage_plan(*args, **kwargs)
    except Exception as exc:
        logger.warning("[platform] plan staging failed: %s", exc, exc_info=True)
        return None
