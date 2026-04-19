from __future__ import annotations

import re
from typing import Any, Dict, Optional


APPROVAL_COMMANDS = {
    "approve",
    "approve.",
    "approved",
    "approved.",
    "yes, approve",
    "yes approve",
    "yes, run it",
    "yes run it",
    "proceed",
    "proceed.",
    "run it",
    "run it.",
}

READ_ONLY_ROUTER_TOOLS = {
    "toolbox_inventory",
    "session_run_io_summary",
    "visualize_job_results",
    "fetch_ncbi_sequence",
    "query_uniprot",
    "lookup_go_term",
    "unsupported_tool",
    "bio_diff_runs",
}

HIGH_IMPACT_ACTION_TYPES = {
    "correct_metadata",
}


def is_approval_command(command: str) -> bool:
    normalized = " ".join((command or "").strip().lower().split())
    return normalized in APPROVAL_COMMANDS or normalized.startswith("approve ")


def has_explicit_execute_intent(command: str) -> bool:
    c = (command or "").lower()
    if any(tok in c for tok in (" and then ", " then ", "->", "→", ";")):
        return True
    return bool(
        re.search(
            r"\b(run|execute|rerun|re-run|regenerate|recreate|reconstruct|start|launch|fix|update|exclude|highlight|color|make|generate|plot|heatmap|show)\b",
            c,
        )
    )


def requires_approval_semantics(command: str, action_type: Optional[str] = None) -> bool:
    c = (command or "").lower()
    if (action_type or "").lower() in HIGH_IMPACT_ACTION_TYPES:
        return True
    triggers = (
        "mislabeled",
        "should be",
        "looks wrong",
        "reversed",
        "focus only",
        "only female",
        "only male",
        "adjusting for",
        "design formula",
        "correction",
    )
    return any(t in c for t in triggers)


def _looks_like_planning_analysis_request(command: str, params: Optional[Dict[str, Any]] = None) -> bool:
    """
    Detect high-level analysis requests that should be staged as a plan before execution.
    Keep this intentionally narrow so concrete/explicit tool commands still run directly.
    """
    c = (command or "").lower()
    planning_cues = (
        "analyze this",
        "analyse this",
        "what is going on",
        "tell me what is going on",
        "design the workflow",
        "before execution",
        "propose expected",
    )
    if not any(cue in c for cue in planning_cues):
        return False
    if has_explicit_execute_intent(command):
        return False
    # NOTE: keep this command-text driven. Router parameters may be auto-filled
    # from demos/defaults, but planning-style user phrasing should still stage.
    return True


def should_stage_for_approval(
    tool_name: str,
    command: str,
    params: Optional[Dict[str, Any]] = None,
    *,
    action_type: Optional[str] = None,
) -> bool:
    """
    Decide whether to stage a plan and request approval before execution.
    Policy is action- and impact-based rather than workflow-name-based.
    """
    if not tool_name or tool_name in READ_ONLY_ROUTER_TOOLS:
        return False
    if tool_name == "__plan__":
        return False
    if isinstance(params, dict) and params.get("session_resolution_error"):
        return False
    if is_approval_command(command):
        return False
    # Concrete routed tools should execute directly. The approval gate is primarily
    # for ambiguous natural-language execution requests (handle_natural_command)
    # and explicit high-impact correction semantics.
    if tool_name != "handle_natural_command" and not requires_approval_semantics(command, action_type):
        return _looks_like_planning_analysis_request(command, params)
    if tool_name == "handle_natural_command" and not requires_approval_semantics(command, action_type):
        return False
    if requires_approval_semantics(command, action_type):
        return True
    if has_explicit_execute_intent(command):
        return False
    return True

