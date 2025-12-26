from __future__ import annotations

"""
Phase 3 — Plan IR (minimal v1)

This schema is intentionally small:
- A plan is an ordered list of steps.
- Each step is a tool invocation with arguments.
- Arguments may contain references to previous step results using:
    {"$ref": "steps.<step_id>.result.<path>"}
where <path> is a dot-separated traversal into dicts/lists.
"""

from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field


PLAN_IR_VERSION = "v1"


class PlanStep(BaseModel):
    id: str = Field(..., description="Unique step id (stable within plan)")
    tool_name: str = Field(..., description="Tool name to execute")
    arguments: Dict[str, Any] = Field(default_factory=dict, description="Tool arguments")
    description: Optional[str] = Field(default=None, description="Human-readable step description")


class Plan(BaseModel):
    version: str = Field(default=PLAN_IR_VERSION)
    steps: List[PlanStep]





