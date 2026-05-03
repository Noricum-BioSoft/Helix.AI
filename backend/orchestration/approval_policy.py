"""Approval-staging policy.

Determines whether a command should be staged for user approval before
execution.  All intent classification is LLM-based — there are no keyword
or regex fallbacks.  If the LLM is unavailable, ``should_stage_for_approval``
raises ``StagingClassificationError`` rather than silently guessing.
"""

from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Configuration (policy-level constants, NOT keyword matchers)
# ---------------------------------------------------------------------------

READ_ONLY_ROUTER_TOOLS = {
    "toolbox_inventory",
    "session_run_io_summary",
    "visualize_job_results",
    "fetch_ncbi_sequence",
    "query_uniprot",
    "lookup_go_term",
    "unsupported_tool",
    "bio_diff_runs",
    # Q&A tools — read-only, never need approval
    "tabular_qa",
}

HIGH_IMPACT_ACTION_TYPES = {
    "correct_metadata",
}


# ---------------------------------------------------------------------------
# LLM-based approval re-export (drop-in for legacy callers)
# ---------------------------------------------------------------------------

from backend.orchestration.approval_classifier import (  # noqa: E402
    is_approval_command,
)


# ---------------------------------------------------------------------------
# LLM staging classifier
# ---------------------------------------------------------------------------

class StagingClassificationError(RuntimeError):
    """Raised when the LLM is unavailable and staging cannot be determined."""


@dataclass(frozen=True)
class StagingDecision:
    requires_approval: bool
    has_execute_intent: bool
    is_planning_request: bool
    method: str
    reason: str


_STAGING_SYSTEM_PROMPT = """\
You are a staging-intent classifier for a bioinformatics AI assistant.

Given a user command, decide three things:

1. "requires_approval": Is this a high-impact correction or metadata change that
   should be reviewed before execution? Examples: relabeling samples, correcting
   design formulas, changing key experimental parameters.

2. "has_execute_intent": Does the user explicitly want to execute something
   immediately? Examples: "run it", "execute the analysis", "rerun with these
   parameters", "go ahead and run", "now do X".

3. "is_planning_request": Is this a request to design or propose a concrete
   analytical workflow that should be staged for the user to review before
   anything runs? This is ONLY true when the user asks the system to plan a
   specific pipeline or workflow. Examples that ARE planning requests:
   "design a workflow for my RNA-seq data", "propose a plan for this analysis",
   "what would a full scRNA-seq pipeline look like for this dataset".

   IMPORTANT — these are NOT planning requests (set false):
   - Advisory or meta questions: "what should I do next?", "what are my options?",
     "what do you recommend?", "can you help me?", "where do I go from here?"
   - Conceptual/educational questions: "what is BLAST?", "explain this result"
   - Simple follow-up questions after a completed analysis step

Respond with ONLY a JSON object on a single line:
{"requires_approval": bool, "has_execute_intent": bool, "is_planning_request": bool}

No explanation. No markdown. Just the JSON object.
"""


def _get_llm():
    """Lazy LLM initializer for staging classification."""
    if os.getenv("HELIX_MOCK_MODE") == "1":
        raise StagingClassificationError("LLM disabled in HELIX_MOCK_MODE")

    openai_key = os.getenv("OPENAI_API_KEY", "").strip()
    deepseek_key = os.getenv("DEEPSEEK_API_KEY", "").strip()

    openai_enabled = openai_key not in ("", "disabled", "your_openai_api_key_here", "none")
    deepseek_enabled = deepseek_key not in ("", "disabled", "your_deepseek_api_key_here", "none")

    if openai_enabled:
        from langchain.chat_models import init_chat_model
        model = os.getenv("HELIX_STAGING_OPENAI_MODEL", "openai:gpt-4.1-mini").strip()
        if ":" not in model:
            model = f"openai:{model}"
        return init_chat_model(model, temperature=0, max_tokens=30)

    if deepseek_enabled:
        from langchain_deepseek import ChatDeepSeek
        return ChatDeepSeek(
            model="deepseek-chat", temperature=0, max_tokens=30, timeout=8.0, max_retries=1
        )

    raise StagingClassificationError("No LLM API key configured for staging classification.")


def _classify_staging_intent(command: str) -> StagingDecision:
    """Ask the LLM for a staging decision on *command*.

    Raises ``StagingClassificationError`` if the LLM is unavailable.
    """
    messages = [
        {"role": "system", "content": _STAGING_SYSTEM_PROMPT},
        {"role": "user", "content": f"User command: {command!r}"},
    ]
    try:
        llm = _get_llm()
        response = llm.invoke(messages)
        raw = (response.content or "").strip()
    except StagingClassificationError:
        raise
    except Exception as exc:
        raise StagingClassificationError(
            f"Staging LLM call failed: {exc}"
        ) from exc

    start = raw.find("{")
    end = raw.rfind("}") + 1
    if start >= 0 and end > start:
        obj = json.loads(raw[start:end])
        return StagingDecision(
            requires_approval=bool(obj.get("requires_approval", False)),
            has_execute_intent=bool(obj.get("has_execute_intent", False)),
            is_planning_request=bool(obj.get("is_planning_request", False)),
            method="llm",
            reason="llm_classified",
        )
    raise StagingClassificationError(f"LLM returned unparseable response: {raw!r}")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def has_explicit_execute_intent(command: str) -> bool:
    """Return True if the user's command explicitly signals they want to execute now."""
    decision = _classify_staging_intent(command)
    return decision.has_execute_intent


def requires_approval_semantics(command: str, action_type: Optional[str] = None) -> bool:
    """Return True if this command requires user approval before execution.

    Checks the action_type label first (fast, free).  If action_type is not
    available or is not a high-impact type, asks the LLM.
    """
    if (action_type or "").lower() in HIGH_IMPACT_ACTION_TYPES:
        return True
    decision = _classify_staging_intent(command)
    return decision.requires_approval


def should_stage_for_approval(
    tool_name: str,
    command: str,
    params: Optional[Dict[str, Any]] = None,
    *,
    action_type: Optional[str] = None,
    has_pending_plan: bool = False,
) -> bool:
    """Decide whether to stage a plan and request user approval before execution.

    Returns True when the command should be reviewed before the system runs it.

    Raises ``StagingClassificationError`` if the LLM is unavailable and the
    decision cannot be made by structural checks alone.
    """
    if not tool_name or tool_name in READ_ONLY_ROUTER_TOOLS:
        return False
    if tool_name == "__plan__":
        return False
    if isinstance(params, dict) and params.get("session_resolution_error"):
        return False

    # High-impact action types always require approval (fast path, no LLM).
    # This check runs BEFORE is_approval_command to avoid an unnecessary LLM
    # call when action_type already tells us the answer.
    if (action_type or "").lower() in HIGH_IMPACT_ACTION_TYPES:
        return True

    # If this IS an approval of a pending plan, don't re-stage.
    if is_approval_command(command, has_pending_plan=has_pending_plan):
        return False

    # Concrete routed tools execute directly unless the LLM says otherwise.
    decision = _classify_staging_intent(command)
    logger.debug(
        "staging_classifier: tool=%r → requires_approval=%s has_execute_intent=%s "
        "is_planning_request=%s",
        tool_name,
        decision.requires_approval,
        decision.has_execute_intent,
        decision.is_planning_request,
    )

    if decision.has_execute_intent:
        return False
    if decision.requires_approval:
        return True
    # Router-level fallbacks always need a plan card before execution.
    # ``multi_step_workflow`` always represents a composite request that needs
    # the approval-then-agent re-route path; ``handle_natural_command`` only
    # stages when the LLM has flagged it as planning.
    if tool_name == "multi_step_workflow":
        return True
    if tool_name == "handle_natural_command" and decision.is_planning_request:
        return True
    return False
