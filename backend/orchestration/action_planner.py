from __future__ import annotations

from typing import Any, Dict, Optional

from backend.action_plan import infer_action_type, map_action_to_tool


def build_single_step_plan(command: str, tool_name: str, params: Dict[str, Any]) -> Dict[str, Any]:
    """Build a single-step action plan.

    When ``tool_name`` is ``handle_natural_command`` and the router returned
    ``suggested_steps``, expand those into labelled steps so the plan card
    shows a meaningful workflow instead of a bare "handle_natural_command" line.
    The execution path (approval → agent re-route) handles actual dispatch.
    """
    params = params or {}
    router_reasoning = params.get("router_reasoning") or {}
    suggested_steps: list = router_reasoning.get("suggested_steps") or []
    action_type = infer_action_type(command, tool_name)

    if tool_name == "handle_natural_command" and len(suggested_steps) > 1:
        steps = [
            {
                "id": f"step{i + 1}",
                "action_type": action_type,
                "tool_name": "handle_natural_command",
                "arguments": params,
                "description": step_desc,
            }
            for i, step_desc in enumerate(suggested_steps)
        ]
    else:
        steps = [
            {
                "id": "step1",
                "action_type": action_type,
                "tool_name": tool_name,
                "arguments": params,
                "description": (command or "").strip(),
            }
        ]

    return {"version": "v1", "steps": steps}


def resolve_tool_for_action(
    *,
    action_type: str,
    fallback_tool: Optional[str],
    arguments: Dict[str, Any],
) -> Optional[str]:
    return map_action_to_tool(action_type, fallback_tool, arguments)


def infer_action(command: str, tool_name: Optional[str] = None) -> str:
    return infer_action_type(command, tool_name)

