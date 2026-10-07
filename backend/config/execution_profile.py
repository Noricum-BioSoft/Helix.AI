"""
Execution Profile — deployment-level provider configuration.

Three layers are kept separate:
- capability descriptors (static: what a provider can do),
- the execution profile (this module: which providers THIS deployment enables,
  how to reach them, and policy limits),
- runtime health (adapter ping at startup; Phase 3A).

Selection: ``HELIX_EXECUTION_PROFILE=<name>`` → ``backend/config/profiles/<name>.yaml``
(default ``local-only``). ``${ENV_VAR}`` placeholders are interpolated; a
missing variable leaves the value ``None`` (never the literal placeholder) so
adapters can fail with a clear ProviderNotConfigured. No secrets live in YAML.
"""

from __future__ import annotations

import importlib
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator

PROFILES_DIR = Path(__file__).resolve().parent / "profiles"
DEFAULT_PROFILE = "local-only"
ENV_VAR_NAME = "HELIX_EXECUTION_PROFILE"

CredentialSource = Literal["env", "profile", "role", "api_key", "none"]

_PLACEHOLDER = re.compile(r"^\$\{([A-Z0-9_]+)\}$")


class ProviderCredentials(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: CredentialSource = "none"
    env_var: Optional[str] = Field(default=None, description="For source=api_key: env var holding the key.")
    profile_name: Optional[str] = Field(default=None, description="For source=profile: named SDK profile.")
    role_arn: Optional[str] = None


class ProviderConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: str
    enabled: bool = True
    adapter: Optional[str] = Field(default=None, description="'module:Class'; None = resolved from the capability descriptor.")
    config: Dict[str, Any] = Field(default_factory=dict)
    credentials: ProviderCredentials = Field(default_factory=ProviderCredentials)

    @field_validator("adapter")
    @classmethod
    def _adapter_shape(cls, v: Optional[str]) -> Optional[str]:
        if v is not None and (":" not in v or v.startswith(":") or v.endswith(":")):
            raise ValueError("adapter must be 'module.path:ClassName'")
        return v


class ExecutionPolicy(BaseModel):
    model_config = ConfigDict(extra="forbid")

    allow_cloud: bool = False
    allowed_regions: List[str] = Field(default_factory=list)
    data_residency: Optional[str] = None
    max_cost_usd_per_run: Optional[float] = Field(default=None, ge=0)


class ExecutionProfile(BaseModel):
    model_config = ConfigDict(extra="forbid")

    profile: str
    default_provider: str
    policy: ExecutionPolicy = Field(default_factory=ExecutionPolicy)
    providers: List[ProviderConfig]
    storage: Dict[str, Any] = Field(default_factory=dict, description="Reserved for the artifact-storage abstraction (LATER).")

    @field_validator("providers")
    @classmethod
    def _unique_ids(cls, v: List[ProviderConfig]) -> List[ProviderConfig]:
        ids = [p.id for p in v]
        if len(ids) != len(set(ids)):
            raise ValueError("provider ids must be unique")
        return v

    def enabled_providers(self) -> List[ProviderConfig]:
        return [p for p in self.providers if p.enabled]

    def get(self, provider_id: str) -> Optional[ProviderConfig]:
        return next((p for p in self.providers if p.id == provider_id), None)

    def is_enabled(self, provider_id: str) -> bool:
        p = self.get(provider_id)
        return bool(p and p.enabled)


class ProfileError(ValueError):
    pass


def _interpolate(value: Any, env: Dict[str, str]) -> Any:
    if isinstance(value, str):
        m = _PLACEHOLDER.match(value)
        if m:
            return env.get(m.group(1))
        return value
    if isinstance(value, dict):
        return {k: _interpolate(v, env) for k, v in value.items()}
    if isinstance(value, list):
        return [_interpolate(v, env) for v in value]
    return value


def _check_adapter_importable(adapter: str) -> None:
    module_name, _, class_name = adapter.partition(":")
    try:
        module = importlib.import_module(module_name)
    except ImportError as exc:
        raise ProfileError(f"adapter module not importable: {module_name} ({exc})") from exc
    if not hasattr(module, class_name):
        raise ProfileError(f"adapter class not found: {adapter}")


def load_execution_profile(
    name: Optional[str] = None,
    *,
    profiles_dir: Path = PROFILES_DIR,
    env: Optional[Dict[str, str]] = None,
    check_adapters: bool = True,
) -> ExecutionProfile:
    env = dict(os.environ if env is None else env)
    name = name or env.get(ENV_VAR_NAME) or DEFAULT_PROFILE
    path = profiles_dir / f"{name}.yaml"
    if not path.exists():
        raise ProfileError(f"execution profile not found: {path}")
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    data = _interpolate(raw, env)
    profile = ExecutionProfile.model_validate(data)
    if profile.get(profile.default_provider) is None:
        raise ProfileError(f"default_provider '{profile.default_provider}' is not declared in profile '{name}'")
    if check_adapters:
        for p in profile.enabled_providers():
            if p.adapter:
                _check_adapter_importable(p.adapter)
    return profile
