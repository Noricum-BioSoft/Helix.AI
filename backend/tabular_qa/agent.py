"""
Conversational data analysis agent (Code Interpreter pattern).

Flow
----
1. Load the DataFrame from the session-resident uploaded file.
2. Build an LLM prompt containing the file's schema + sample rows.
3. Ask the LLM to write minimal pandas code that answers the user's question.
4. Execute in sandbox (executor.py).
5. If the code errors, feed the traceback back to the LLM and retry (≤ 2 times).
6. Narrate the result back to the user in plain language.
"""
from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

MAX_RETRIES = 2
_QA_MODEL_ENV = "HELIX_TABULAR_QA_MODEL"
_DEFAULT_MODEL = "openai:gpt-4o"

# ---------------------------------------------------------------------------
# Prompt helpers
# ---------------------------------------------------------------------------

def _schema_to_text(profile: Dict[str, Any]) -> str:
    lines: List[str] = []
    cols = (profile.get("schema") or {}).get("columns") or []
    if cols:
        lines.append("Columns:")
        for c in cols:
            if isinstance(c, dict):
                name = c.get("name", "?")
                dtype = c.get("dtype", "?")
                n_missing = c.get("n_missing", 0)
                stats = ""
                if "mean" in c:
                    stats = f"  mean={c['mean']}, min={c['min']}, max={c['max']}"
                elif "top_values" in c:
                    top = list(c["top_values"].keys())[:5]
                    stats = f"  top values: {top}"
                lines.append(f"  - {name} ({dtype}){', has_nulls' if n_missing else ''}{stats}")
            else:
                lines.append(f"  - {c}")
    summary = profile.get("summary") or {}
    n_rows = summary.get("n_rows") or profile.get("n_records")
    if n_rows:
        lines.append(f"Shape: {n_rows} rows × {len(cols)} cols")
    available = profile.get("available_sheets")
    if available:
        lines.append(f"Excel sheets available: {available}")
    return "\n".join(lines)


def _sample_to_text(profile: Dict[str, Any]) -> str:
    sample = profile.get("sample") or []
    if not sample:
        return ""
    import pandas as pd
    try:
        df = pd.DataFrame(sample)
        return df.to_string(index=False, max_cols=12, max_colwidth=20)
    except Exception:
        return json.dumps(sample[:3], default=str)


def _build_code_prompt(
    question: str,
    profile: Dict[str, Any],
    error_context: Optional[str] = None,
) -> str:
    schema_text = _schema_to_text(profile)
    sample_text = _sample_to_text(profile)
    sheet = (profile.get("summary") or {}).get("source_sheet")
    sheet_note = f"\nThe data was loaded from sheet: **{sheet}**." if sheet else ""

    error_block = ""
    if error_context:
        error_block = f"""
The previous attempt raised this error — fix it:
```
{error_context}
```
"""

    return f"""You are a pandas data analysis assistant.{sheet_note}

The user has a pandas DataFrame called `df` with the following schema:
{schema_text}

First few rows:
{sample_text}

{error_block}
Write **only** Python/pandas code (no markdown fences, no explanation) that answers:
"{question}"

Rules:
- Use only `df`, `pd`, `np` — no other imports.
- Store the final answer in a variable called `result`.
- `result` must be a DataFrame, Series, scalar, or dict.
- Keep code concise and correct.
"""


def _build_narration_prompt(
    question: str,
    exec_result: Dict[str, Any],
    profile: Dict[str, Any],
) -> str:
    result_json = json.dumps(exec_result.get("result"), default=str, indent=2)
    if len(result_json) > 4000:
        result_json = result_json[:4000] + "\n…(truncated)"
    schema_text = _schema_to_text(profile)

    return f"""You are a scientific data analyst. The user asked:
"{question}"

The data has this schema:
{schema_text}

Python code was executed and produced this result:
{result_json}

Write a clear, concise plain-language answer to the user's question based on this result.
Include specific numbers, gene names, or values from the result where relevant.
If the result is a table, summarise the key findings.
Do not repeat the code. Do not use markdown headers.
"""


# ---------------------------------------------------------------------------
# LLM helper
# ---------------------------------------------------------------------------

def _get_llm():
    if os.getenv("HELIX_MOCK_MODE") == "1":
        raise RuntimeError("LLM disabled in HELIX_MOCK_MODE")

    from langchain.chat_models import init_chat_model

    openai_key = os.getenv("OPENAI_API_KEY", "").strip()
    deepseek_key = os.getenv("DEEPSEEK_API_KEY", "").strip()

    openai_bad = {"", "disabled", "your_openai_api_key_here", "none"}
    deepseek_bad = {"", "disabled", "your_deepseek_api_key_here", "none"}

    if openai_key and openai_key not in openai_bad:
        model_str = os.getenv(_QA_MODEL_ENV, _DEFAULT_MODEL)
        return init_chat_model(model_str, temperature=0)

    if deepseek_key and deepseek_key not in deepseek_bad:
        from langchain_deepseek import ChatDeepSeek
        return ChatDeepSeek(model="deepseek-chat", temperature=0, max_tokens=2048)

    raise ValueError("No LLM API key configured.")


def _llm_invoke(llm: Any, prompt: str) -> str:
    from langchain_core.messages import HumanMessage
    response = llm.invoke([HumanMessage(content=prompt)])
    return response.content.strip()


def _strip_code_fences(text: str) -> str:
    """Remove ```python ... ``` fencing if the LLM adds it despite instructions."""
    import re
    text = re.sub(r"^```[a-z]*\n?", "", text.strip(), flags=re.MULTILINE)
    text = re.sub(r"\n?```$", "", text.strip(), flags=re.MULTILINE)
    return text.strip()


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_tabular_qa(
    *,
    question: str,
    session_id: str,
    file_path: str,
    profile: Dict[str, Any],
    sheet: Optional[str] = None,
    timeout_s: int = 30,
) -> Dict[str, Any]:
    """
    Answer *question* about the file at *file_path* using the Code Interpreter loop.

    Parameters
    ----------
    question :   The user's natural-language question.
    session_id : For logging / audit.
    file_path :  Absolute path to the uploaded file on disk.
    profile :    FileProfile returned by ``profile_file()`` at upload time.
    sheet :      For Excel, the sheet to load (None → first sheet).
    timeout_s :  Per-execution sandbox timeout.

    Returns
    -------
    dict with keys:
        success      bool
        answer       str  — narrated plain-language answer
        result       dict — raw serialised execution result
        code         str  — the generated + executed code
        attempts     int  — number of LLM/exec iterations used
        error        str  — set only on final failure
    """
    from backend.ds_pipeline.pipelines.ingest import ingest_tabular
    from backend.tabular_qa.executor import execute_code

    # --- Load DataFrame ---
    try:
        _conn, df, _meta = ingest_tabular(file_path, sheet=sheet)
    except Exception as exc:
        return {
            "success": False,
            "answer": f"Could not load the file: {exc}",
            "result": None,
            "code": "",
            "attempts": 0,
            "error": str(exc),
        }

    # --- Code generation + execution loop ---
    try:
        llm = _get_llm()
    except RuntimeError:
        # Mock mode: return a canned response for testability
        return {
            "success": True,
            "answer": "LLM unavailable (mock mode). Code execution skipped.",
            "result": {"type": "null", "value": None},
            "code": "result = None  # mock mode",
            "attempts": 0,
            "error": None,
        }

    last_error: Optional[str] = None
    last_exec: Optional[Dict[str, Any]] = None
    code = ""
    attempt = 0  # guard: keeps attempt defined if the loop range is ever empty

    for attempt in range(1, MAX_RETRIES + 2):  # 1 initial try + MAX_RETRIES retries
        prompt = _build_code_prompt(question, profile, error_context=last_error)
        try:
            raw = _llm_invoke(llm, prompt)
        except Exception as exc:
            logger.warning("tabular_qa: LLM call failed on attempt %d: %s", attempt, exc)
            return {
                "success": False,
                "answer": f"LLM error: {exc}",
                "result": None,
                "code": code,
                "attempts": attempt,
                "error": str(exc),
            }

        code = _strip_code_fences(raw)
        exec_result = execute_code(code, df, timeout_s=timeout_s)
        last_exec = exec_result

        if exec_result["success"]:
            break

        last_error = exec_result["error"]
        logger.info("tabular_qa: execution failed (attempt %d): %s", attempt, last_error)

    if not (last_exec and last_exec.get("success")):
        return {
            "success": False,
            "answer": f"Could not produce a valid result after {attempt} attempts. Last error: {last_error}",
            "result": last_exec.get("result") if last_exec else None,
            "code": code,
            "attempts": attempt,
            "error": last_error,
        }

    # --- Narration ---
    narration_prompt = _build_narration_prompt(question, last_exec, profile)
    try:
        answer = _llm_invoke(llm, narration_prompt)
    except Exception as exc:
        # Narration failure is non-fatal — fall back to repr of result
        answer = f"Result computed successfully. (Narration failed: {exc})"

    return {
        "success": True,
        "answer": answer,
        "result": last_exec.get("result"),
        "code": code,
        "attempts": attempt,
        "error": None,
    }
