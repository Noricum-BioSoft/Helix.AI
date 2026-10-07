"""
ProvenanceRecord — enough to answer "why did the platform behave differently
for this workflow than for an earlier run?"

Includes tool, model, agent, prompt-template (id/version/hash — never the
prompt text), code commit, schema version, parameters, execution backend,
inputs and outputs. No hidden reasoning is stored.
"""

from __future__ import annotations

from datetime import datetime
from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field

from backend.contracts.base import TracedModel, utcnow


class Actor(BaseModel):
    model_config = ConfigDict(extra="forbid")

    kind: Literal["user", "agent", "system"]
    id: str


class ProvenanceRecord(TracedModel):
    event_id: str
    event_type: str = Field(..., description="plan_created | assessed | recommended | authorized | approved | submitted | completed | observed | decided | …")
    actor: Actor
    tool: Optional[str] = None
    tool_version: Optional[str] = None
    model_version: Optional[str] = None
    agent_id: Optional[str] = None
    agent_version: Optional[str] = None
    prompt_template_id: Optional[str] = None
    prompt_template_version: Optional[str] = None
    prompt_template_hash: Optional[str] = None
    code_commit: Optional[str] = None
    parameters: Dict[str, Any] = Field(default_factory=dict)
    execution_backend: Optional[str] = None
    provider_id: Optional[str] = None
    inputs: List[str] = Field(default_factory=list, description="Ids/URIs consumed.")
    outputs: List[str] = Field(default_factory=list, description="Ids/URIs produced.")
    timestamp: datetime = Field(default_factory=utcnow)
