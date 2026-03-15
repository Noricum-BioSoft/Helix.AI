from __future__ import annotations

from typing import Any, Dict, List

from backend.tool_schemas import get_tool_schema


def validate_tool_bindings(tool_name: str, arguments: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Validate required input bindings for a concrete tool invocation.
    Returns a list of structured issues. Empty list means validation passed.
    """
    issues: List[Dict[str, Any]] = []
    if not tool_name:
        return [{"type": "invalid_tool", "message": "Missing tool_name"}]
    if not isinstance(arguments, dict):
        return [{"type": "invalid_arguments", "message": "Arguments must be an object"}]

    schema = get_tool_schema(tool_name) or {}
    required = (
        ((schema.get("inputs") or {}).get("required"))
        if isinstance(schema, dict)
        else None
    ) or []
    required = [k for k in required if isinstance(k, str)]

    for key in required:
        val = arguments.get(key)
        if val in (None, "", [], {}):
            issues.append(
                {
                    "type": "missing_artifact",
                    "tool_name": tool_name,
                    "input": key,
                    "message": f"Required input '{key}' is missing.",
                }
            )
            continue
        if isinstance(val, str) and val.strip().lower() in {"unresolved", "__unresolved__", "unknown"}:
            issues.append(
                {
                    "type": "unresolved_historical_selector",
                    "tool_name": tool_name,
                    "input": key,
                    "message": f"Input '{key}' resolved to an unresolved selector token.",
                }
            )

    return issues

