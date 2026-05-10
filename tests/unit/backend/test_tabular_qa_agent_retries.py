"""
Unit tests for ``run_tabular_qa`` retry / re-codegen behaviour (tabular_qa/agent.py).

Uses a real CSV ingest + mocked LLM and sandbox execution so retries are fully
deterministic without calling external models.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List
from unittest.mock import MagicMock

import pytest

from backend.tabular_qa import agent as tab_agent


def _csv_and_profile(tmp_path: Path) -> tuple[str, Dict[str, Any]]:
    p = tmp_path / "qa_sample.csv"
    p.write_text("gene,count\nA,10\nB,20\n", encoding="utf-8")
    profile: Dict[str, Any] = {
        "schema": {
            "columns": [
                {"name": "gene", "dtype": "object", "n_missing": 0},
                {"name": "count", "dtype": "int64", "n_missing": 0},
            ]
        },
        "summary": {},
        "sample_rows": [{"gene": "A", "count": 10}],
    }
    return str(p), profile


def test_run_tabular_qa_retries_after_exec_failures_then_narrates(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setenv("HELIX_TABULAR_CODEGEN_MAX_ATTEMPTS", "6")

    file_path, profile = _csv_and_profile(tmp_path)
    prompts: List[str] = []
    exec_calls: List[str] = []

    def fake_execute_code(code: str, df: Any, *, timeout_s: int = 30) -> Dict[str, Any]:
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
            return {"success": False, "error": "ValueError: simulated failure on second run"}
        return {
            "success": True,
            "result": {"type": "dict", "value": {"rows": len(df), "exec_n": n}},
        }

    def fake_llm_invoke(llm: Any, prompt: str) -> str:
        prompts.append(prompt)
        if "Python code was executed and produced this result" in prompt:
            return "Narration: counts look consistent after retries."
        assert "You are a pandas data analysis assistant" in prompt
        gen_idx = len([x for x in prompts if "You are a pandas data analysis assistant" in x])
        return f"result = {{'gen': {gen_idx}}}\n"

    with monkeypatch.context() as m:
        m.setattr(tab_agent, "_get_llm", lambda: MagicMock(name="llm"))
        m.setattr(tab_agent, "_llm_invoke", fake_llm_invoke)
        m.setattr("backend.tabular_qa.executor.execute_code", fake_execute_code)

        out = tab_agent.run_tabular_qa(
            question="What are the counts?",
            session_id="qa-retry-session",
            file_path=file_path,
            profile=profile,
        )

    assert out["success"] is True
    assert out["attempts"] == 3
    assert out["attempts_max"] == 6
    assert "Narration" in (out.get("answer") or "")
    assert len(exec_calls) == 3
    assert len([p for p in prompts if "You are a pandas data analysis assistant" in p]) == 3
    narration_prompts = [p for p in prompts if "Python code was executed and produced this result" in p]
    assert len(narration_prompts) == 1

    # Retry prompts must surface the prior sandbox / runtime error to the LLM.
    second_codegen = [p for p in prompts if "You are a pandas data analysis assistant" in p][1]
    assert "The previous attempt raised this error" in second_codegen
    assert "SandboxImportError" in second_codegen or "import statements" in second_codegen
    third_codegen = [p for p in prompts if "You are a pandas data analysis assistant" in p][2]
    assert "ValueError" in third_codegen


def test_run_tabular_qa_exhausts_attempts_without_success(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setenv("HELIX_TABULAR_CODEGEN_MAX_ATTEMPTS", "3")

    file_path, profile = _csv_and_profile(tmp_path)

    def always_fail_execute(code: str, df: Any, *, timeout_s: int = 30) -> Dict[str, Any]:
        return {
            "success": False,
            "error": "SandboxImportError: import statements are not allowed in this sandbox.",
        }

    def fake_llm_invoke(llm: Any, prompt: str) -> str:
        if "Python code was executed and produced this result" in prompt:
            narration_hit.append(True)
            pytest.fail("narration must not run when execution never succeeds")
        return "import pandas as pd\nresult = {}\n"

    with monkeypatch.context() as m:
        m.setattr(tab_agent, "_get_llm", lambda: MagicMock())
        m.setattr(tab_agent, "_llm_invoke", fake_llm_invoke)
        m.setattr("backend.tabular_qa.executor.execute_code", always_fail_execute)

        out = tab_agent.run_tabular_qa(
            question="Summarise the table.",
            session_id="qa-all-fail",
            file_path=file_path,
            profile=profile,
        )

    assert out["success"] is False
    assert out["attempts"] == 3
    assert out["attempts_max"] == 3
    assert out.get("error_class") == "sandbox_import"
    assert out.get("tabular_retry_available") is False
    assert "import pandas" in (out.get("last_code_preview") or "")
