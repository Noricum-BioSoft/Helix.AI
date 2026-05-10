"""
Unit tests for tabular analysis plan execution retry / re-codegen behaviour.

Mocks LLM codegen + sandbox ``execute_code`` so we do not call real models or
run unsafe exec — while still verifying re-generation and re-execution paths.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple
from unittest.mock import MagicMock

import pytest

from backend.tabular_qa import analysis_executor as ae


def _minimal_csv(tmp_path: Path) -> Path:
    p = tmp_path / "session_data.csv"
    p.write_text("gene,count\nA,1\nB,2\n", encoding="utf-8")
    return p


def _session_context(csv_path: Path) -> Dict[str, Any]:
    return {
        "uploaded_files": [
            {
                "filename": "session_data.csv",
                "local_path": str(csv_path),
                "schema_preview": {
                    "family": "tabular",
                    "schema": {
                        "columns": [
                            {"name": "gene", "dtype": "object"},
                            {"name": "count", "dtype": "int64"},
                        ]
                    },
                    "summary": {},
                },
            }
        ]
    }


def _plan() -> Dict[str, Any]:
    return {
        "type": "tabular_analysis",
        "title": "Test plan",
        "goal": "Summarise counts",
        "steps": [
            {
                "id": 1,
                "name": "aggregate",
                "type": "compute",
                "operation": "group by gene",
            }
        ],
    }


def test_execute_analysis_plan_regenerates_after_sandbox_failures_then_succeeds(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """First sandbox runs fail; executor must call LLM + execute again until success."""
    monkeypatch.setenv("HELIX_TABULAR_CODEGEN_MAX_ATTEMPTS", "5")

    csv_path = _minimal_csv(tmp_path)
    plan = _plan()
    session_id = "test-session-retry-ok"

    codegen_llm_calls: List[Tuple[str, str]] = []
    exec_calls: List[str] = []

    def fake_execute_code(code: str, df: Any) -> Dict[str, Any]:
        exec_calls.append(code)
        n = len(exec_calls)
        if n == 1:
            return {
                "success": False,
                "error": (
                    "SandboxImportError: import statements are not allowed in this sandbox. "
                    "Use df, pd, np directly."
                ),
            }
        if n == 2:
            return {"success": False, "error": "RuntimeError: simulated numeric failure"}
        return {
            "success": True,
            "result": {"type": "dict", "value": {"attempts_needed": n, "rows": int(len(df))}},
        }

    def fake_llm_call(llm: Any, system: str, user: str) -> str:
        if "biological interpretation" in system:
            return "Interpretation: synthetic success after retries."
        assert "Python data analyst" in system, "expected codegen system prompt"
        codegen_llm_calls.append((system, user))
        # Return distinct bodies so execute mock can branch (here execute uses call count).
        return f"# gen_{len(codegen_llm_calls)}\nresult = {{'gen': {len(codegen_llm_calls)}}}\n"

    with monkeypatch.context() as m:
        m.setattr(ae, "_get_llm", lambda: object())
        m.setattr(ae, "_llm_call", fake_llm_call)
        m.setattr("backend.tabular_qa.executor.execute_code", fake_execute_code)

        out = ae.execute_analysis_plan(plan, session_id, _session_context(csv_path))

    assert out["status"] == "success"
    assert out["attempts_used"] == 3
    assert out["attempts_max"] == 5
    assert "Interpretation" in (out.get("text") or "")
    assert len(codegen_llm_calls) == 3
    assert len(exec_calls) == 3

    # Second and third codegen prompts must include correction context after first failure.
    assert "CORRECTION REQUIRED" in codegen_llm_calls[1][1]
    assert "SandboxImportError" in codegen_llm_calls[1][1] or "import statements" in codegen_llm_calls[1][1]
    assert "CORRECTION REQUIRED" in codegen_llm_calls[2][1]
    assert "RuntimeError" in codegen_llm_calls[2][1]


def test_execute_analysis_plan_exhausts_attempts_without_success(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setenv("HELIX_TABULAR_CODEGEN_MAX_ATTEMPTS", "3")

    csv_path = _minimal_csv(tmp_path)
    plan = _plan()

    def always_fail_execute(code: str, df: Any) -> Dict[str, Any]:
        return {"success": False, "error": "SandboxImportError: import statements are not allowed"}

    def fake_llm_call(llm: Any, system: str, user: str) -> str:
        if "biological interpretation" in system:
            pytest.fail("interpretation LLM must not run when codegen never succeeds")
        return "import os\nresult = {}\n"

    with monkeypatch.context() as m:
        m.setattr(ae, "_get_llm", lambda: MagicMock())
        m.setattr(ae, "_llm_call", fake_llm_call)
        m.setattr("backend.tabular_qa.executor.execute_code", always_fail_execute)

        out = ae.execute_analysis_plan(plan, "test-session-all-fail", _session_context(csv_path))

    assert out["status"] == "error"
    assert out["attempts_used"] == 3
    assert out["attempts_max"] == 3
    assert out.get("tabular_retry_available") is True
    assert out.get("error_class") == "sandbox_import"
    assert "last_code_preview" in out
    assert "import os" in (out.get("last_code_preview") or "")
    assert "after 3 attempts" in (out.get("error") or "")
