"""
Identity mode for approval principals.

``HELIX_AUTH_MODE``:
- ``dev_header`` (default): trust the ``X-Helix-User`` header. Development only.
- ``oidc``: authenticated principal from an OIDC bearer token (pilot+; Phase 1.3 extension point).

``HELIX_ENV``: development | staging | production.

Guard: a production deployment MUST NOT run with ``dev_header``. ``main.py``
calls ``assert_auth_mode_allowed()`` at startup and refuses to start otherwise,
so the temporary mechanism cannot silently survive into production.
"""

from __future__ import annotations

import os
from typing import Literal, Optional

from backend.contracts.human_approval import Principal

AuthMode = Literal["dev_header", "oidc"]
Environment = Literal["development", "staging", "production"]

AUTH_MODE_VAR = "HELIX_AUTH_MODE"
ENV_VAR = "HELIX_ENV"
DEV_HEADER = "X-Helix-User"
DEV_ROLES_HEADER = "X-Helix-Roles"


class AuthConfigurationError(RuntimeError):
    pass


def auth_mode() -> AuthMode:
    value = os.getenv(AUTH_MODE_VAR, "dev_header").strip().lower()
    if value not in ("dev_header", "oidc"):
        raise AuthConfigurationError(f"{AUTH_MODE_VAR} must be dev_header|oidc, got {value!r}")
    return value  # type: ignore[return-value]


def environment() -> Environment:
    value = os.getenv(ENV_VAR, "development").strip().lower()
    if value not in ("development", "staging", "production"):
        raise AuthConfigurationError(f"{ENV_VAR} must be development|staging|production, got {value!r}")
    return value  # type: ignore[return-value]


def assert_auth_mode_allowed() -> None:
    """Fail closed: production + dev_header is a misconfiguration, not a warning."""
    env = environment()
    mode = auth_mode()
    if env == "production" and mode == "dev_header":
        raise AuthConfigurationError(
            f"{AUTH_MODE_VAR}=dev_header is not permitted when {ENV_VAR}=production; configure {AUTH_MODE_VAR}=oidc"
        )


def principal_from_headers(headers) -> Optional[Principal]:
    """Resolve a principal from request headers according to the active mode.

    Returns None when no principal is present (callers decide whether that is
    an error). OIDC resolution is a Phase 1.3 extension point and currently
    raises so it cannot be mistaken for working auth.
    """
    assert_auth_mode_allowed()
    mode = auth_mode()
    if mode == "dev_header":
        subject = headers.get(DEV_HEADER)
        if not subject:
            return None
        roles = [r.strip() for r in (headers.get(DEV_ROLES_HEADER) or "").split(",") if r.strip()]
        return Principal(
            subject_id=subject.strip(),
            display_name=subject.strip(),
            identity_provider="dev_header",
            auth_method="header",
            roles=roles,
        )
    raise AuthConfigurationError("HELIX_AUTH_MODE=oidc is not implemented yet (Phase 1.3 extension point)")
