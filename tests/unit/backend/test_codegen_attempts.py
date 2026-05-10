"""Tests for tabular codegen attempt limits and error classification."""

from __future__ import annotations

import pytest

from backend.tabular_qa.codegen_attempts import (
    classify_tabular_execution_error,
    get_tabular_codegen_max_attempts,
    truncate_code_preview,
)


def test_get_tabular_codegen_max_attempts_default(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("HELIX_TABULAR_CODEGEN_MAX_ATTEMPTS", raising=False)
    assert get_tabular_codegen_max_attempts() == 4


def test_get_tabular_codegen_max_attempts_clamp_high(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("HELIX_TABULAR_CODEGEN_MAX_ATTEMPTS", "99")
    assert get_tabular_codegen_max_attempts() == 16


def test_get_tabular_codegen_max_attempts_clamp_low(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("HELIX_TABULAR_CODEGEN_MAX_ATTEMPTS", "0")
    assert get_tabular_codegen_max_attempts() == 1


def test_classify_sandbox_import() -> None:
    assert (
        classify_tabular_execution_error("SandboxImportError: import statements are not allowed")
        == "sandbox_import"
    )


def test_truncate_code_preview() -> None:
    long = "x" * 5000
    out = truncate_code_preview(long, limit=100)
    assert len(out) <= 100
    assert "truncated" in out.lower()


def test_truncate_short_unchanged() -> None:
    assert truncate_code_preview("print(1)") == "print(1)"
