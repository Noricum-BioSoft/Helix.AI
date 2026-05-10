"""Configurable codegen ↔ sandbox retry limits for tabular analysis."""

from __future__ import annotations

import os
from typing import Optional

_ENV_KEY = "HELIX_TABULAR_CODEGEN_MAX_ATTEMPTS"
_DEFAULT = 4
_MIN = 1
_MAX = 16
_CODE_PREVIEW_LIMIT = 2000


def get_tabular_codegen_max_attempts() -> int:
    """Total LLM codegen + sandbox attempts (initial try + repairs)."""
    raw = os.getenv(_ENV_KEY, str(_DEFAULT)).strip()
    try:
        n = int(raw)
    except ValueError:
        n = _DEFAULT
    return max(_MIN, min(n, _MAX))


def truncate_code_preview(code: str, limit: int = _CODE_PREVIEW_LIMIT) -> str:
    if not code:
        return ""
    text = code.strip()
    if len(text) <= limit:
        return text
    return text[: limit - 24] + "\n…(truncated)"


def classify_tabular_execution_error(message: Optional[str]) -> str:
    if not message:
        return "unknown"
    m = message.lower()
    if "sandboximporterror" in m or "import statements are not allowed" in m:
        return "sandbox_import"
    if "timeout" in m:
        return "timeout"
    if "syntaxerror" in m or "syntax error" in m:
        return "syntax"
    return "execution"
