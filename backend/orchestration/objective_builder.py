"""
Objective builder (Phase 1.1).

``build_objective`` turns the user's command plus session context into a
``ScientificObjective`` and mints ``objective_id`` + ``trace_id`` — the only
place a loop's ``trace_id`` originates.

Two layers:

* **Deterministic** (always runs; the only layer in ``HELIX_MOCK_MODE``):
  objective/question from the command, ``known_inputs`` from session uploads
  and registered artifacts, execution constraints from the active execution
  profile policy.
* **LLM enrichment** (``HELIX_OBJECTIVE_LLM=1``, off by default): asks the
  model for ``question``, ``biological_system``, ``desired_evidence``,
  ``expected_outputs`` and ``success_criteria`` as strict JSON. Any failure
  falls back to the deterministic objective. Enrichment never changes
  ``objective``, ``known_inputs``, ``constraints`` or the IDs.
"""

from __future__ import annotations

import json
import logging
import os
from typing import Any, Dict, List, Optional

from backend.contracts.ids import new_id, new_trace_id
from backend.contracts.scientific_objective import (
    ArtifactRef,
    ExecutionConstraints,
    ObjectiveConstraints,
    ScientificObjective,
)

logger = logging.getLogger(__name__)

OBJECTIVE_LLM_FLAG = "HELIX_OBJECTIVE_LLM"
_MAX_KNOWN_INPUTS = 50

_ENRICH_PROMPT = """You extract a structured scientific objective from a bioinformatics request.
Return ONLY a JSON object with keys:
  question (string), biological_system (string or null),
  desired_evidence (array of short strings), expected_outputs (array of short strings),
  success_criteria (array of short strings).
Do not include explanations. Do not invent inputs or data that were not mentioned."""


def _known_inputs_from_context(session_context: Optional[Dict[str, Any]]) -> List[ArtifactRef]:
    refs: List[ArtifactRef] = []
    if not isinstance(session_context, dict):
        return refs
    for entry in session_context.get("uploaded_files") or []:
        if isinstance(entry, dict):
            uri = entry.get("stored_path") or entry.get("path") or entry.get("uri") or entry.get("filename") or entry.get("name")
            if uri:
                refs.append(ArtifactRef(uri=str(uri), kind="file", content_hash=entry.get("sha256") or entry.get("content_hash")))
        elif isinstance(entry, str) and entry:
            refs.append(ArtifactRef(uri=entry, kind="file"))
    artifacts = session_context.get("artifacts")
    if isinstance(artifacts, dict):
        for artifact_id, record in artifacts.items():
            if not isinstance(record, dict):
                continue
            refs.append(ArtifactRef(artifact_id=str(artifact_id), uri=record.get("uri"), kind=record.get("type")))
    # Deterministic order + bound
    refs.sort(key=lambda r: (r.artifact_id or "", r.uri or ""))
    return refs[:_MAX_KNOWN_INPUTS]


def _execution_constraints_from_profile() -> Optional[ExecutionConstraints]:
    try:
        from backend.config.execution_profile import load_execution_profile

        profile = load_execution_profile(check_adapters=False)
    except Exception as exc:  # profile problems must not block objective creation
        logger.debug("objective_builder: execution profile unavailable: %s", exc)
        return None
    return ExecutionConstraints(
        allowed_providers=[p.id for p in profile.enabled_providers()],
        preferred_provider=profile.default_provider,
        data_residency=profile.policy.data_residency,
        max_cost_usd=profile.policy.max_cost_usd_per_run,
    )


def build_objective_deterministic(
    command: str,
    session_context: Optional[Dict[str, Any]] = None,
    intent_result: Optional[Dict[str, Any]] = None,
    *,
    trace_id: Optional[str] = None,
) -> ScientificObjective:
    text = (command or "").strip() or "Unspecified objective"
    intent_result = intent_result or {}
    user_context: Dict[str, Any] = {}
    if isinstance(session_context, dict) and session_context.get("session_id"):
        user_context["session_id"] = session_context["session_id"]
    if intent_result.get("intent"):
        user_context["intent"] = intent_result["intent"]
    if intent_result.get("tool"):
        user_context["routed_tool"] = intent_result["tool"]
    return ScientificObjective(
        objective_id=new_id(),
        trace_id=trace_id or new_trace_id(),
        objective=text[:2000],
        question=text[:2000],
        constraints=ObjectiveConstraints(execution=_execution_constraints_from_profile()),
        known_inputs=_known_inputs_from_context(session_context),
        user_context=user_context,
        status="active",
    )


def _llm_enabled() -> bool:
    if os.getenv("HELIX_MOCK_MODE") == "1":
        return False
    return os.getenv(OBJECTIVE_LLM_FLAG, "").strip().lower() in {"1", "true", "yes", "on"}


def _enrich_with_llm(objective: ScientificObjective, command: str) -> ScientificObjective:
    from backend.orchestration.approval_classifier import _get_llm  # same lazy pattern, same key handling

    llm = _get_llm()
    response = llm.invoke(
        [
            {"role": "system", "content": _ENRICH_PROMPT},
            {"role": "user", "content": command},
        ]
    )
    raw = (getattr(response, "content", "") or "").strip()
    start, end = raw.find("{"), raw.rfind("}") + 1
    if start < 0 or end <= start:
        raise ValueError("objective enrichment returned no JSON object")
    data = json.loads(raw[start:end])

    def _strs(value: Any) -> List[str]:
        return [str(v)[:500] for v in value if str(v).strip()][:20] if isinstance(value, list) else []

    return objective.model_copy(
        update={
            "question": (str(data.get("question") or "").strip() or objective.question)[:2000],
            "biological_system": (str(data["biological_system"]).strip()[:200] if data.get("biological_system") else None),
            "desired_evidence": _strs(data.get("desired_evidence")),
            "expected_outputs": _strs(data.get("expected_outputs")),
            "success_criteria": _strs(data.get("success_criteria")),
        }
    )


def build_objective(
    command: str,
    session_context: Optional[Dict[str, Any]] = None,
    intent_result: Optional[Dict[str, Any]] = None,
    *,
    trace_id: Optional[str] = None,
) -> ScientificObjective:
    """Build the loop's ScientificObjective. Deterministic in mock mode; LLM-enriched when enabled."""
    objective = build_objective_deterministic(command, session_context, intent_result, trace_id=trace_id)
    if not _llm_enabled():
        return objective
    try:
        return _enrich_with_llm(objective, command)
    except Exception as exc:
        logger.warning("objective_builder: LLM enrichment failed, using deterministic objective: %s", exc)
        return objective
