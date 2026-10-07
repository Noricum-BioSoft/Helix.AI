"""
Common base for platform contracts: trace id, schema version, timestamps.

Every shared domain object carries ``trace_id`` (one per scientific loop) and
``schema_version`` so provenance can explain behavioural differences between
runs. Ids are minted via ``backend.contracts.ids`` only.
"""

from __future__ import annotations

from datetime import datetime, timezone

from pydantic import BaseModel, ConfigDict, Field, field_validator

from backend.contracts.ids import is_trace_id

CONTRACT_SCHEMA_VERSION = 1


def utcnow() -> datetime:
    return datetime.now(timezone.utc)


class TracedModel(BaseModel):
    """Base for all persisted platform records."""

    model_config = ConfigDict(extra="forbid")

    trace_id: str = Field(..., description="Shared by every record of one scientific loop (orch_<uuid7>).")
    schema_version: int = Field(default=CONTRACT_SCHEMA_VERSION, ge=1)
    created_at: datetime = Field(default_factory=utcnow)

    @field_validator("trace_id")
    @classmethod
    def _trace_format(cls, v: str) -> str:
        if not is_trace_id(v):
            raise ValueError("trace_id must look like orch_<32 hex chars>")
        return v
