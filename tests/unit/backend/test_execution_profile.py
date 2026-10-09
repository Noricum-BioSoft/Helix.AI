from __future__ import annotations

from pathlib import Path

import pytest

from backend.config.execution_profile import (
    DEFAULT_PROFILE,
    PROFILES_DIR,
    ENV_VAR_NAME,
    ProfileError,
    load_execution_profile,
)


def test_local_only_loads_with_empty_env_and_has_no_cloud():
    p = load_execution_profile(env={})
    assert p.profile == DEFAULT_PROFILE
    assert p.policy.allow_cloud is False
    assert {x.id for x in p.enabled_providers()} == {"local_compute", "nextflow", "mock_experimental_provider"}
    assert p.is_enabled("aws_emr") is False


def test_profile_selected_by_env_var_and_placeholders_interpolated():
    p = load_execution_profile(env={ENV_VAR_NAME: "aws-dev", "AWS_REGION": "us-east-1"})
    emr = p.get("aws_emr")
    assert emr is not None and emr.enabled
    assert emr.config["region"] == "us-east-1"
    assert emr.config["workdir_uri"] is None  # missing env var -> None, never the literal placeholder
    assert emr.credentials.source == "env"


def test_unknown_profile_and_unknown_adapter_rejected(tmp_path: Path):
    with pytest.raises(ProfileError, match="not found"):
        load_execution_profile("nope", env={})
    (tmp_path / "bad.yaml").write_text(
        "profile: bad\ndefault_provider: x\nproviders:\n  - id: x\n    adapter: backend.no_such_module:Thing\n"
    )
    with pytest.raises(ProfileError, match="not importable"):
        load_execution_profile("bad", profiles_dir=tmp_path, env={})
    (tmp_path / "bad2.yaml").write_text(
        "profile: bad2\ndefault_provider: x\nproviders:\n  - id: x\n    adapter: backend.contracts.ids:NoSuchClass\n"
    )
    with pytest.raises(ProfileError, match="class not found"):
        load_execution_profile("bad2", profiles_dir=tmp_path, env={})


def test_default_provider_must_be_declared_and_ids_unique(tmp_path: Path):
    (tmp_path / "p.yaml").write_text("profile: p\ndefault_provider: missing\nproviders:\n  - id: a\n")
    with pytest.raises(ProfileError, match="default_provider"):
        load_execution_profile("p", profiles_dir=tmp_path, env={})
    (tmp_path / "q.yaml").write_text("profile: q\ndefault_provider: a\nproviders:\n  - id: a\n  - id: a\n")
    with pytest.raises(Exception, match="unique"):
        load_execution_profile("q", profiles_dir=tmp_path, env={})


def test_shipped_profiles_all_parse():
    for path in PROFILES_DIR.glob("*.yaml"):
        load_execution_profile(path.stem, env={})
