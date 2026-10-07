from __future__ import annotations

import pytest

from backend.config import auth_mode as am


def test_defaults_are_development_dev_header(monkeypatch):
    monkeypatch.delenv(am.AUTH_MODE_VAR, raising=False)
    monkeypatch.delenv(am.ENV_VAR, raising=False)
    assert am.environment() == "development" and am.auth_mode() == "dev_header"
    am.assert_auth_mode_allowed()


def test_production_with_dev_header_refuses_to_start(monkeypatch):
    monkeypatch.setenv(am.ENV_VAR, "production")
    monkeypatch.setenv(am.AUTH_MODE_VAR, "dev_header")
    with pytest.raises(am.AuthConfigurationError, match="not permitted"):
        am.assert_auth_mode_allowed()


def test_production_with_oidc_is_allowed_by_guard(monkeypatch):
    monkeypatch.setenv(am.ENV_VAR, "production")
    monkeypatch.setenv(am.AUTH_MODE_VAR, "oidc")
    am.assert_auth_mode_allowed()


def test_invalid_values_rejected(monkeypatch):
    monkeypatch.setenv(am.AUTH_MODE_VAR, "basic")
    with pytest.raises(am.AuthConfigurationError):
        am.auth_mode()
    monkeypatch.setenv(am.AUTH_MODE_VAR, "dev_header")
    monkeypatch.setenv(am.ENV_VAR, "prod")
    with pytest.raises(am.AuthConfigurationError):
        am.environment()


def test_dev_header_principal_resolution(monkeypatch):
    monkeypatch.delenv(am.ENV_VAR, raising=False)
    monkeypatch.setenv(am.AUTH_MODE_VAR, "dev_header")
    assert am.principal_from_headers({}) is None
    p = am.principal_from_headers({am.DEV_HEADER: "alice", am.DEV_ROLES_HEADER: "scientist, reviewer"})
    assert p is not None and p.subject_id == "alice" and p.identity_provider == "dev_header"
    assert p.roles == ["scientist", "reviewer"]


def test_oidc_mode_is_explicit_extension_point(monkeypatch):
    monkeypatch.delenv(am.ENV_VAR, raising=False)
    monkeypatch.setenv(am.AUTH_MODE_VAR, "oidc")
    with pytest.raises(am.AuthConfigurationError, match="not implemented"):
        am.principal_from_headers({"Authorization": "Bearer x"})
