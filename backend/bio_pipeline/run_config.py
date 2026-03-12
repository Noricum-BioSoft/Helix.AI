"""BioRunConfig — immutable configuration for a single bioinformatics run."""
from __future__ import annotations

import uuid
from dataclasses import dataclass, field
from typing import Any, Dict, Optional


@dataclass
class BioRunConfig:
    """Captures everything needed to reproduce and compare one bio tool run."""

    tool_name: str
    params: Dict[str, Any]
    run_id: str = field(default_factory=lambda: str(uuid.uuid4()))
    parent_run_id: Optional[str] = None
    session_id: Optional[str] = None
    objective: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "tool_name": self.tool_name,
            "params": self.params,
            "run_id": self.run_id,
            "parent_run_id": self.parent_run_id,
            "session_id": self.session_id,
            "objective": self.objective,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "BioRunConfig":
        return cls(
            tool_name=d["tool_name"],
            params=d["params"],
            run_id=d.get("run_id", str(uuid.uuid4())),
            parent_run_id=d.get("parent_run_id"),
            session_id=d.get("session_id"),
            objective=d.get("objective", ""),
        )
