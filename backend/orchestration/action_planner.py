from __future__ import annotations

from typing import Any, Dict, Optional

from backend.action_plan import infer_action_type, map_action_to_tool


def build_single_step_plan(command: str, tool_name: str, params: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "version": "v1",
        "steps": [
            {
                "id": "step1",
                "action_type": infer_action_type(command, tool_name),
                "tool_name": tool_name,
                "arguments": params or {},
                "description": (command or "").strip(),
            }
        ],
    }


def resolve_tool_for_action(
    *,
    action_type: str,
    fallback_tool: Optional[str],
    arguments: Dict[str, Any],
) -> Optional[str]:
    return map_action_to_tool(action_type, fallback_tool, arguments)


def infer_action(command: str, tool_name: Optional[str] = None) -> str:
    return infer_action_type(command, tool_name)

