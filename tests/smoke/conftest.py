"""
Smoke test configuration.

Smoke tests require a live LLM (OPENAI_API_KEY or DEEPSEEK_API_KEY).
They are automatically skipped when neither key is set or when
HELIX_MOCK_MODE=1 is active.
"""
from __future__ import annotations

import os

import pytest


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "smoke: live-LLM smoke tests (skipped without API key)",
    )


@pytest.fixture(autouse=True, scope="session")
def _require_llm():
    """Skip the entire smoke suite if no LLM key is available."""
    has_key = bool(
        os.getenv("OPENAI_API_KEY") or os.getenv("DEEPSEEK_API_KEY")
    )
    mock_mode = os.getenv("HELIX_MOCK_MODE", "0") == "1"
    if not has_key or mock_mode:
        pytest.skip(
            "Smoke tests require OPENAI_API_KEY or DEEPSEEK_API_KEY "
            "and HELIX_MOCK_MODE != 1",
            allow_module_level=True,
        )


@pytest.fixture(autouse=True, scope="session")
def _disable_mock_mode():
    """Ensure mock mode is off for the smoke suite."""
    os.environ["HELIX_MOCK_MODE"] = "0"
    os.environ["HELIX_AGENT_DISABLED"] = "0"
    yield
    # Restore after session (best-effort; pytest isolation handles the rest)
    os.environ.pop("HELIX_MOCK_MODE", None)
    os.environ.pop("HELIX_AGENT_DISABLED", None)
