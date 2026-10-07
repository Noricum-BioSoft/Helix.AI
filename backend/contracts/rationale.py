"""
Structured, user-safe rationale.

Platform rule (plan rev 2, principle 4): persist decisions, evidence refs,
assumptions, confidence and concise explanations — never opaque model
reasoning traces. Every place a free-text ``reasoning``/``explanation``
would otherwise be stored uses ``list[RationaleItem]`` instead.
"""

from __future__ import annotations

from typing import List, Optional

from pydantic import BaseModel, Field, field_validator


class RationaleItem(BaseModel):
    """One statement → conclusion pair, optionally backed by evidence references."""

    evidence_refs: List[str] = Field(
        default_factory=list,
        description="Ids/URIs of evidence this item relies on (dataset:…, artifact:…, criterion:…, policy:…).",
    )
    statement: str = Field(..., min_length=1, description="Observed fact or assumption, stated plainly.")
    conclusion: str = Field(..., min_length=1, description="What follows from the statement for this decision.")
    confidence: Optional[float] = Field(default=None, ge=0.0, le=1.0)

    @field_validator("statement", "conclusion")
    @classmethod
    def _strip_non_empty(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("must not be blank")
        return v


Rationale = List[RationaleItem]
