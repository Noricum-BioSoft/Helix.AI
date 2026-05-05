"""
Eval-suite conftest.

Evals exercise the *live* LLM router/intent classifier and validate that the
real model produces the expected tool routing and intent decisions on a fixed
set of curated cases. Running them under ``HELIX_MOCK_MODE=1`` (the default for
the rest of the test suite, see ``tests/conftest.py``) defeats the purpose:
``backend/command_router.py`` and ``backend/intent_classifier.py`` raise
``RoutingError`` / ``RuntimeError("LLM is disabled in HELIX_MOCK_MODE")`` rather
than silently falling back to a heuristic.

This conftest:

1. Disables ``HELIX_MOCK_MODE`` for the eval suite only (top-level
   ``tests/conftest.py`` setdefaults it to "1"; we override to "0" here so
   each subprocess-imported module sees the live-LLM setting).
2. Loads ``.env`` so ``OPENAI_API_KEY`` (or ``DEEPSEEK_API_KEY``) is available
   to the LLM clients.
3. Skips the suite cleanly if no LLM credentials are present, instead of
   producing 70+ identical "LLM unavailable" failures that drown out real
   regressions.
"""

from __future__ import annotations

import os

import pytest


def _has_llm_credentials() -> bool:
    return bool(
        os.environ.get("OPENAI_API_KEY")
        or os.environ.get("DEEPSEEK_API_KEY")
        or os.environ.get("AZURE_OPENAI_API_KEY")
    )


def pytest_configure(config):  # noqa: D401 — pytest hook signature
    """Switch to live-LLM mode and load .env credentials before collection."""
    try:
        from dotenv import load_dotenv

        load_dotenv()
    except Exception:
        pass

    os.environ["HELIX_MOCK_MODE"] = "0"


def pytest_collection_modifyitems(config, items):
    """Skip the entire eval suite when no LLM credentials are configured.

    This keeps CI green on PRs that don't carry secrets while still failing
    loudly when the eval suite is intended to run.
    """
    if _has_llm_credentials():
        return

    skip_marker = pytest.mark.skip(
        reason=(
            "Eval suite requires a live LLM; set OPENAI_API_KEY (or "
            "DEEPSEEK_API_KEY / AZURE_OPENAI_API_KEY) to enable."
        )
    )
    for item in items:
        item.add_marker(skip_marker)
