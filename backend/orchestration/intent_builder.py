"""
Execution-intent builder (Phase 1.3 shape).

Builds the immutable ``ExecutionIntent`` a human approves. Until the
recommender/authorizer exist (P2/P4) the provider comes from the legacy
``InfraDecision`` (``provider_id = legacy:<infrastructure>``) and the
capability is ``local_compute:<tool_name>``; assessment/recommendation/
authorization fields stay ``None`` while their flags are off.

Hashes:

* ``input_manifest_hash`` — stable hash over the resolved inputs: for each
  argument that looks like a file/URI binding, ``{key, uri, size, content_hash}``
  where size/hash are filled only when the file is local and readable.
* ``execution_parameters_hash`` — stable hash over every step's
  ``{id, tool_name, action_type, arguments}``.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from backend.contracts.execution_intent import ExecutionIntent
from backend.contracts.ids import new_id, stable_hash
from backend.contracts.scientific_plan import ScientificPlan
from backend.orchestration.invariants import check_intent_creation_allowed

LEGACY_PROVIDER_PREFIX = "legacy:"
LOCAL_CAPABILITY_PREFIX = "local_compute:"
_URI_SCHEMES = ("s3://", "gs://", "az://", "http://", "https://", "file://")
_INPUT_KEY_HINTS = ("file", "path", "uri", "url", "input", "fastq", "fasta", "reads", "matrix", "counts", "sheet", "bam", "vcf")
_MAX_LOCAL_HASH_BYTES = 256 * 1024 * 1024


def legacy_provider_id(infrastructure: Optional[str]) -> str:
    return f"{LEGACY_PROVIDER_PREFIX}{infrastructure or 'Local'}"


def capability_id_for_plan(plan: ScientificPlan) -> str:
    """P1: one capability per plan — the first concrete tool; multi-step plans keep step capabilities on the plan."""
    for step in plan.steps:
        if step.required_capabilities:
            return step.required_capabilities[0]
        if step.tool_name:
            return f"{LOCAL_CAPABILITY_PREFIX}{step.tool_name}"
    return f"{LOCAL_CAPABILITY_PREFIX}__plan__"


def _looks_like_input(key: str, value: Any) -> bool:
    if not isinstance(value, str) or not value:
        return False
    if value.startswith(_URI_SCHEMES):
        return True
    lowered = key.lower()
    return any(h in lowered for h in _INPUT_KEY_HINTS)


def _local_size_and_hash(uri: str) -> Tuple[Optional[int], Optional[str]]:
    path = Path(uri[len("file://"):] if uri.startswith("file://") else uri)
    try:
        if not path.is_file():
            return None, None
        size = path.stat().st_size
        if size > _MAX_LOCAL_HASH_BYTES:
            return size, None
        digest = hashlib.sha256()
        with path.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                digest.update(chunk)
        return size, digest.hexdigest()
    except OSError:
        return None, None


def collect_input_manifest(plan: ScientificPlan, extra_inputs: Optional[Iterable[Dict[str, Any]]] = None) -> List[Dict[str, Any]]:
    """Resolved inputs referenced by the plan (deterministically ordered)."""
    manifest: List[Dict[str, Any]] = []
    for step in plan.steps:
        for key, value in step.arguments.items():
            values = value if isinstance(value, list) else [value]
            for v in values:
                if _looks_like_input(key, v):
                    size, content_hash = _local_size_and_hash(v)
                    manifest.append({"step_id": step.id, "key": key, "uri": v, "size": size, "content_hash": content_hash})
    for item in extra_inputs or []:
        if isinstance(item, dict) and item.get("uri"):
            manifest.append(
                {
                    "step_id": item.get("step_id"),
                    "key": item.get("key"),
                    "uri": item["uri"],
                    "size": item.get("size"),
                    "content_hash": item.get("content_hash"),
                }
            )
    manifest.sort(key=lambda m: (str(m.get("step_id")), str(m.get("key")), str(m.get("uri"))))
    return manifest


def compute_input_manifest_hash(manifest: List[Dict[str, Any]]) -> str:
    return stable_hash(manifest)


def compute_execution_parameters_hash(plan: ScientificPlan) -> str:
    return stable_hash(
        [{"id": s.id, "tool_name": s.tool_name, "action_type": s.action_type, "arguments": s.arguments} for s in plan.steps]
    )


def build_intent(
    plan: ScientificPlan,
    infra_decision: Optional[Any] = None,
    inputs: Optional[Iterable[Dict[str, Any]]] = None,
    *,
    current_state: Optional[str] = None,
    provider_id: Optional[str] = None,
    capability_id: Optional[str] = None,
    assessment_id: Optional[str] = None,
    assessment_hash: Optional[str] = None,
    recommendation_id: Optional[str] = None,
    authorization_id: Optional[str] = None,
) -> ExecutionIntent:
    """Build the P1 ExecutionIntent for a staged plan.

    ``infra_decision`` may be an ``InfraDecision`` or anything with an
    ``infrastructure`` attribute; ``None`` means Local.
    """
    if current_state is not None:
        check_intent_creation_allowed(current_state)
    infrastructure = getattr(infra_decision, "infrastructure", None) if infra_decision is not None else None
    manifest = collect_input_manifest(plan, inputs)
    return ExecutionIntent(
        execution_intent_id=new_id(),
        trace_id=plan.trace_id,
        plan_id=plan.plan_id,
        plan_hash=plan.plan_hash,
        assessment_id=assessment_id,
        assessment_hash=assessment_hash,
        recommendation_id=recommendation_id,
        authorization_id=authorization_id,
        provider_id=provider_id or legacy_provider_id(infrastructure),
        capability_id=capability_id or capability_id_for_plan(plan),
        input_manifest_hash=compute_input_manifest_hash(manifest),
        execution_parameters_hash=compute_execution_parameters_hash(plan),
    )
