"""
Scientific Capability Registry — descriptor schema.

A ``CapabilityDescriptor`` states what a provider *can* do (static, vendor-
published facts). Whether this deployment *uses* it is an execution-profile
concern, so there is deliberately no ``available``/``enabled`` field here.
Descriptors are loaded from ``backend/config/capabilities/*.yaml`` (Phase 3A).
"""

from __future__ import annotations

from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, field_validator

CapabilityCategory = Literal["computational", "experimental", "data", "advisory"]
ExecutionMode = Literal["sync", "async"]
IdempotencySupport = Literal["native", "client_dedup", "none"]
CostModel = Literal["free", "per_vcpu_hour", "per_instance_hour", "per_sample", "per_request", "per_gb", "quoted"]
DescriptorSource = Literal["manual", "generated", "provider"]


class CapabilityIO(BaseModel):
    model_config = ConfigDict(extra="forbid")

    name: str
    type: str = Field(..., description="file-path | directory-path | sequence | table | string | number | …")
    format: Optional[str] = Field(default=None, description="fastq | bam | csv | fasta | …")
    required: bool = True
    description: Optional[str] = None


class CapabilityCost(BaseModel):
    model_config = ConfigDict(extra="forbid")

    model: CostModel
    unit_usd_range: Optional[List[float]] = Field(default=None, min_length=2, max_length=2)
    notes: Optional[str] = None

    @field_validator("unit_usd_range")
    @classmethod
    def _ordered(cls, v: Optional[List[float]]) -> Optional[List[float]]:
        if v is not None and (v[0] < 0 or v[0] > v[1]):
            raise ValueError("unit_usd_range must be [min, max] with 0 <= min <= max")
        return v


class CapabilityTurnaround(BaseModel):
    model_config = ConfigDict(extra="forbid")

    estimated_minutes: Optional[float] = Field(default=None, ge=0)
    estimated_days: Optional[float] = Field(default=None, ge=0)


class CapabilityIntegration(BaseModel):
    model_config = ConfigDict(extra="forbid")

    api_available: bool = True
    adapter: str = Field(..., description="Import path 'module:Class' of the ExecutionProvider implementation.")
    execution_mode: ExecutionMode = "sync"
    idempotency_support: IdempotencySupport = "none"


class CapabilitySecurity(BaseModel):
    model_config = ConfigDict(extra="forbid")

    screening_required: bool = False
    human_approval_required: bool = False
    data_classes_allowed: List[str] = Field(default_factory=lambda: ["public", "internal"])
    regions: List[str] = Field(default_factory=list, description="Where the provider runs/stores data; empty = unspecified.")


class CapabilityProvenance(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: DescriptorSource = "manual"
    catalog_sha: Optional[str] = None
    catalog_fetched_at: Optional[str] = None
    pinned_release: Optional[Dict[str, str]] = None


class CapabilityDescriptor(BaseModel):
    model_config = ConfigDict(extra="forbid")

    provider: str = Field(..., description="Provider id, e.g. local_compute, nextflow, mock_experimental_provider.")
    capability_id: str = Field(..., description="Globally unique, e.g. 'nextflow:nf-core/rnaseq', 'mock_experimental:protein_expression'.")
    category: CapabilityCategory
    description: Optional[str] = None
    inputs: List[CapabilityIO] = Field(default_factory=list)
    outputs: List[CapabilityIO] = Field(default_factory=list)
    supported_systems: List[str] = Field(default_factory=list)
    constraints: Dict[str, Any] = Field(default_factory=dict)
    cost: CapabilityCost
    turnaround: CapabilityTurnaround = Field(default_factory=CapabilityTurnaround)
    integration: CapabilityIntegration
    security: CapabilitySecurity = Field(default_factory=CapabilitySecurity)
    prior_performance: Dict[str, Any] = Field(default_factory=dict)
    provenance: CapabilityProvenance = Field(default_factory=CapabilityProvenance)
    tool_name_aliases: List[str] = Field(default_factory=list, description="Helix tool names that resolve to this capability.")

    @field_validator("capability_id")
    @classmethod
    def _namespaced(cls, v: str) -> str:
        if ":" not in v:
            raise ValueError("capability_id must be namespaced as '<provider>:<name>'")
        return v
