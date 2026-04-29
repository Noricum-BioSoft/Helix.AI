"""
Canonical advisory response schema and normalizer.

The LLM may return advisory/explanation/planning responses in several ad-hoc
shapes.  This module defines a single canonical ``HelixAdvisory`` Pydantic
model and a ``normalize_advisory_text`` function that converts any recognised
shape into it, adding a top-level ``"helix_type": "advisory"`` sentinel that
the frontend uses for deterministic rendering.

Supported input shapes:
  1. ``{ classification, requirements, recommended_workflow, ... }``
     (ChIP-seq / planning advisor from the BioAgent prompt)
  2. ``{ type: "answer", title, summary, sections, next_steps }``
     (explanation / "Explain results" response)
  3. ``{ domain, task_type, user_friendly_summary, details_markdown, ... }``
     (legacy agent.md §11 schema — NOT rewritten; returned as-is)
"""
from __future__ import annotations

import json
import re
from typing import Any, Dict, List, Optional, Union

from pydantic import BaseModel, Field


# ---------------------------------------------------------------------------
# Canonical schema
# ---------------------------------------------------------------------------

class AdvisoryClassification(BaseModel):
    domain: Optional[str] = None
    task_type: Optional[str] = None
    feasible: Optional[bool] = None


class AdvisoryItem(BaseModel):
    """Generic label+description item used in requirements, questions, etc."""
    label: str
    description: Optional[str] = None
    examples: Optional[List[str]] = None
    tools: Optional[List[str]] = None


class AdvisorySection(BaseModel):
    """A narrative section with optional list items."""
    heading: str
    content: Optional[str] = None
    items: Optional[List[Union[AdvisoryItem, str]]] = None


class AdvisoryWorkflowStep(BaseModel):
    step: int
    name: str
    description: Optional[str] = None


class HelixAdvisory(BaseModel):
    """Canonical advisory response envelope returned by the backend."""
    helix_type: str = Field(default="advisory", frozen=True)
    title: str
    summary: str
    classification: Optional[AdvisoryClassification] = None
    sections: List[AdvisorySection] = Field(default_factory=list)
    workflow_steps: List[AdvisoryWorkflowStep] = Field(default_factory=list)
    requirements: List[AdvisoryItem] = Field(default_factory=list)
    questions_for_user: List[AdvisoryItem] = Field(default_factory=list)
    next_steps: List[str] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _extract_json(text: str) -> Optional[Dict[str, Any]]:
    """Return the first JSON object parsed from *text* (may be fenced)."""
    stripped = text.strip()
    # Strip fenced code blocks
    if stripped.startswith("```"):
        m = re.match(r"```(?:json)?\s*([\s\S]*?)\s*```", stripped)
        if m:
            stripped = m.group(1).strip()
    if not stripped.startswith("{"):
        return None
    try:
        return json.loads(stripped)
    except json.JSONDecodeError:
        return None


def _is_advisory(obj: Dict[str, Any]) -> bool:
    """Return True if *obj* looks like any known advisory shape."""
    keys = set(obj.keys())
    # Already canonical
    if obj.get("helix_type") == "advisory":
        return True
    # Shape 1 – planning advisor (classification may be dict or string)
    if keys & {"classification", "requirements", "recommended_workflow",
                "minimum_information_needed_from_you"}:
        return True
    # Shape 3 – gpt-4.1 plan+answer shape
    if "plan" in keys and "answer" in keys:
        return True
    # Shape 2 – explanation advisor
    if (obj.get("type") in ("answer", "advisory")
            or (isinstance(obj.get("sections"), list) and "summary" in keys)):
        return True
    return False


def _item_from_dict(d: Any) -> AdvisoryItem:
    """Convert an arbitrary dict to AdvisoryItem (best-effort)."""
    if isinstance(d, str):
        return AdvisoryItem(label=d)
    label = (d.get("item") or d.get("step") or d.get("metric")
             or d.get("term") or d.get("label") or d.get("name") or "")
    description = (d.get("details") or d.get("description") or d.get("purpose")
                   or d.get("content") or d.get("meaning") or d.get("interpretation") or "")
    examples = d.get("examples") or []
    tools = d.get("tools") or []
    return AdvisoryItem(
        label=str(label),
        description=str(description) if description else None,
        examples=examples if isinstance(examples, list) else [str(examples)],
        tools=tools if isinstance(tools, list) else [],
    )


def _section_from_shape2(s: Dict[str, Any]) -> AdvisorySection:
    """Convert a Shape-2 section dict to AdvisorySection."""
    heading = s.get("heading", "")
    content = s.get("content")
    raw_items = s.get("items")
    items: Optional[List[Union[AdvisoryItem, str]]] = None
    if raw_items:
        items = [_item_from_dict(i) if isinstance(i, dict) else str(i)
                 for i in raw_items]
    return AdvisorySection(heading=heading, content=content, items=items)


def _normalise_shape1(obj: Dict[str, Any]) -> HelixAdvisory:
    """Normalise the planning advisory (classification / requirements / workflow) shape."""
    # Classification
    raw_cls = obj.get("classification")
    classification = None
    if isinstance(raw_cls, dict):
        classification = AdvisoryClassification(
            domain=raw_cls.get("domain"),
            task_type=raw_cls.get("task_type"),
            feasible=raw_cls.get("feasible"),
        )

    # Summary
    raw_ans = obj.get("answer", {})
    summary = (
        raw_ans if isinstance(raw_ans, str)
        else (raw_ans.get("summary") if isinstance(raw_ans, dict) else None)
        or obj.get("summary", "")
    )

    # Title – derive from title field, then task_type, then classification
    raw_title = obj.get("title", "")
    task_type = (raw_cls or {}).get("task_type", "") if isinstance(raw_cls, dict) else ""
    title = raw_title or task_type or "Analysis Plan"

    # Requirements — may be a flat list OR a nested dict of sub-sections
    req = obj.get("requirements", {})
    sections: List[AdvisorySection] = []
    requirements: List[AdvisoryItem] = []
    if isinstance(req, list):
        # Flat list → populate requirements directly
        requirements = [_item_from_dict(i) for i in req]
    elif isinstance(req, dict):
        for section_key, heading in [
            ("input_data", "Required Input Data"),
            ("reference_files", "Reference Files"),
            ("software_or_pipeline_components", "Pipeline Components"),
        ]:
            items_raw = req.get(section_key, [])
            if items_raw:
                sections.append(AdvisorySection(
                    heading=heading,
                    items=[_item_from_dict(i) for i in items_raw],
                ))

    # Workflow steps
    wf_raw = obj.get("recommended_workflow", [])
    workflow_steps: List[AdvisoryWorkflowStep] = []
    for w in (wf_raw if isinstance(wf_raw, list) else []):
        try:
            workflow_steps.append(AdvisoryWorkflowStep(
                step=int(w.get("step", 0)),
                name=w.get("name", ""),
                description=w.get("description"),
            ))
        except Exception:
            pass

    # Questions for user
    q_raw = obj.get("minimum_information_needed_from_you", [])
    questions = [_item_from_dict({"label": q.get("question", ""), **q})
                 for q in (q_raw if isinstance(q_raw, list) else [])]

    # Next steps — may be "next_steps" (list) or legacy "next_step" (str/dict)
    next_steps_raw = obj.get("next_steps") or []
    next_steps: List[str] = [ns if isinstance(ns, str) else str(ns)
                              for ns in (next_steps_raw if isinstance(next_steps_raw, list)
                                         else [next_steps_raw])]
    if not next_steps:
        next_step = obj.get("next_step", {})
        if isinstance(next_step, dict) and next_step.get("message"):
            next_steps.append(next_step["message"])
        elif isinstance(next_step, str):
            next_steps.append(next_step)

    return HelixAdvisory(
        title=title,
        summary=str(summary),
        classification=classification,
        sections=sections,
        workflow_steps=workflow_steps,
        requirements=requirements,
        questions_for_user=questions,
        next_steps=next_steps,
    )


def _normalise_shape3(obj: Dict[str, Any]) -> HelixAdvisory:
    """Normalise the gpt-4.1 plan+answer shape.

    Handles two sub-variants:
      a) { plan: [str, ...], answer: { section_key: { ... }, ... } }
      b) { plan: { title, summary, workflow_type }, answer: { sections: [...], next_steps: [...] } }
    """
    raw_cls = obj.get("classification")
    task_type = ""
    if isinstance(raw_cls, dict):
        task_type = raw_cls.get("task_type", "")
    elif isinstance(raw_cls, str):
        task_type = raw_cls

    # plan may be a list of step strings OR a dict with title/summary/workflow_type
    plan_raw = obj.get("plan", [])
    if isinstance(plan_raw, dict):
        title = (plan_raw.get("title")
                 or (plan_raw.get("workflow_type") or task_type or "").replace("_", " ").title()
                 or "Analysis Plan")
        summary = plan_raw.get("summary", "")
        # Prefer workflow_type from plan dict if top-level classification was absent
        if not task_type:
            task_type = plan_raw.get("workflow_type", "")
    else:
        title = task_type.replace("_", " ").title() if task_type else "Analysis Plan"
        summary = " ".join(plan_raw) if isinstance(plan_raw, list) else ""

    # Walk the answer dict
    answer = obj.get("answer", {})
    sections: List[AdvisorySection] = []
    workflow_steps: List[AdvisoryWorkflowStep] = []
    requirements: List[AdvisoryItem] = []
    questions: List[AdvisoryItem] = []

    if isinstance(answer, dict):
        # Sub-variant b: answer has explicit "sections" / "next_steps" keys
        if "sections" in answer or "next_steps" in answer:
            for s in (answer.get("sections") or []):
                if isinstance(s, dict):
                    sections.append(_section_from_shape2(s))
        else:
            # Sub-variant a: freeform section_key → section_value mapping
            step_counter = 1
            for section_key, section_val in answer.items():
                heading = section_key.replace("_", " ").title()
                if isinstance(section_val, dict):
                    items: List[Union[AdvisoryItem, str]] = []
                    for sub_key, sub_val in section_val.items():
                        sub_heading = sub_key.replace("_", " ").title()
                        if isinstance(sub_val, list):
                            for entry in sub_val:
                                if isinstance(entry, str):
                                    items.append(AdvisoryItem(label=entry))
                                elif isinstance(entry, dict):
                                    items.append(_item_from_dict(entry))
                        elif isinstance(sub_val, str):
                            items.append(AdvisoryItem(label=sub_heading, description=sub_val))
                    sections.append(AdvisorySection(heading=heading, items=items if items else None))
                elif isinstance(section_val, list):
                    items_list: List[Union[AdvisoryItem, str]] = []
                    for entry in section_val:
                        if isinstance(entry, str):
                            workflow_steps.append(AdvisoryWorkflowStep(
                                step=step_counter, name=entry
                            ))
                            step_counter += 1
                        elif isinstance(entry, dict):
                            items_list.append(_item_from_dict(entry))
                    if items_list:
                        sections.append(AdvisorySection(heading=heading, items=items_list))
                elif isinstance(section_val, str):
                    sections.append(AdvisorySection(heading=heading, content=section_val))
    elif isinstance(answer, str):
        summary = answer

    # next_steps from answer.next_steps (sub-variant b) or obj.next_steps
    next_steps_raw = (
        (answer.get("next_steps") if isinstance(answer, dict) else None)
        or obj.get("next_steps")
        or obj.get("next_step")
        or []
    )
    next_steps: List[str] = [ns if isinstance(ns, str) else str(ns)
                              for ns in (next_steps_raw if isinstance(next_steps_raw, list)
                                         else [next_steps_raw])]

    return HelixAdvisory(
        title=title,
        summary=summary,
        sections=sections,
        workflow_steps=workflow_steps,
        requirements=requirements,
        questions_for_user=questions,
        next_steps=next_steps,
    )


def _normalise_shape2(obj: Dict[str, Any]) -> HelixAdvisory:
    """Normalise the explanation advisory (type/answer, sections, next_steps) shape."""
    title = obj.get("title", "Analysis Results")
    summary = obj.get("summary", "")
    sections = [_section_from_shape2(s)
                for s in (obj.get("sections") or [])
                if isinstance(s, dict)]
    next_steps_raw = obj.get("next_steps") or obj.get("next_step") or []
    next_steps: List[str] = []
    for ns in (next_steps_raw if isinstance(next_steps_raw, list) else [next_steps_raw]):
        if isinstance(ns, str):
            next_steps.append(ns)
        elif isinstance(ns, dict):
            msg = ns.get("message") or ns.get("text") or str(ns)
            next_steps.append(msg)

    return HelixAdvisory(
        title=str(title),
        summary=str(summary),
        sections=sections,
        next_steps=next_steps,
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def normalize_advisory_text(text: str) -> str:
    """
    If *text* is an advisory JSON response, normalise it to the canonical
    ``HelixAdvisory`` JSON shape and return the serialised result.
    Otherwise return *text* unchanged.
    """
    if not text or not text.strip():
        return text
    obj = _extract_json(text)
    if obj is None or not isinstance(obj, dict):
        return text
    if not _is_advisory(obj):
        return text

    try:
        # Already canonical — validate fields exist and pass through unchanged
        if obj.get("helix_type") == "advisory":
            # Ensure required fields are present; fill defaults if LLM omitted them
            obj.setdefault("title", "Analysis Plan")
            obj.setdefault("summary", "")
            obj.setdefault("sections", [])
            obj.setdefault("workflow_steps", [])
            obj.setdefault("requirements", [])
            obj.setdefault("questions_for_user", [])
            obj.setdefault("next_steps", [])
            return json.dumps(obj)

        # Determine which legacy shape we have
        keys = set(obj.keys())
        if "plan" in keys and "answer" in keys:
            advisory = _normalise_shape3(obj)
        elif keys & {"classification", "requirements", "recommended_workflow",
                     "minimum_information_needed_from_you"}:
            advisory = _normalise_shape1(obj)
        else:
            advisory = _normalise_shape2(obj)

        return advisory.model_dump_json(exclude_none=True)
    except Exception:
        # Safety: if normalization fails, return original text unchanged
        return text
