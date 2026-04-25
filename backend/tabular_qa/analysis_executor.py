"""
Execute a tabular analysis plan produced by analysis_planner.py.

On approval the executor:
1. Generates a comprehensive pandas/numpy/matplotlib script from the plan
2. Runs it in the existing sandbox (executor.py)
3. Generates a plain-language biological interpretation of the results
"""
from __future__ import annotations

import json
import logging
import os
import re
import shutil
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)

_EXECUTOR_MODEL_ENV = "HELIX_ANALYSIS_EXECUTOR_MODEL"
_DEFAULT_MODEL = "openai:gpt-4o"

# ---------------------------------------------------------------------------
# Prompts
# ---------------------------------------------------------------------------

CODEGEN_SYSTEM_PROMPT = """\
You are an expert Python data analyst writing code for a sandboxed execution environment.

Rules — MUST follow:
- The variable `df` (a pandas DataFrame) is already defined. Do NOT load any file.
- Pre-bound names (use directly, NO import statements):
    df       — pandas DataFrame with the data
    pd       — pandas
    np       — numpy
    plt      — matplotlib.pyplot
    sns      — seaborn
    stats    — scipy.stats  (e.g. stats.pearsonr, stats.ttest_ind, stats.spearmanr)
    scipy    — scipy root   (e.g. scipy.stats, scipy.cluster)
- Do NOT write any import/from statements — they will raise ImportError in the sandbox.
- Do NOT use: open(), exec(), eval(), os, subprocess, sys, __import__
- Store all final results in a dict named `result` with descriptive string keys.
  Every value must be JSON-serialisable: use .to_dict(), .tolist(), float(), str() as needed.
- For every plot: create the figure, then call
      plt.savefig('/tmp/helix_plot.png', bbox_inches='tight', dpi=150); plt.close()
  Only one plot per run.
- Handle missing values: use .dropna() where scientifically appropriate.
- Keep code concise and correct.

Return ONLY Python code — no markdown fences, no comments, no explanation.
"""

INTERPRETATION_SYSTEM_PROMPT = """\
You are an expert bioinformatics analyst writing for a scientist audience.

Given an analysis plan, the dataset schema, and the computed results, write a concise
biological interpretation of the findings.

Requirements:
- 2–4 paragraphs of plain prose (no bullet points, no markdown headers)
- Lead with the key scientific finding
- Include specific numbers or names from the results
- Explain biological significance where relevant
- Note any important caveats (sample size, missing data, assumptions)
- Do NOT repeat raw code or JSON
"""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _get_llm():
    bad = {"", "disabled", "your_openai_api_key_here", "none"}
    openai_key = os.getenv("OPENAI_API_KEY", "").strip()
    deepseek_key = os.getenv("DEEPSEEK_API_KEY", "").strip()

    if openai_key and openai_key not in bad:
        from langchain.chat_models import init_chat_model
        return init_chat_model(os.getenv(_EXECUTOR_MODEL_ENV, _DEFAULT_MODEL), temperature=0)

    if deepseek_key and deepseek_key not in bad:
        from langchain_deepseek import ChatDeepSeek
        return ChatDeepSeek(model="deepseek-chat", temperature=0, max_tokens=4096)

    raise ValueError("No LLM API key configured for analysis executor.")


def _llm_call(llm, system: str, user: str) -> str:
    from langchain_core.messages import HumanMessage, SystemMessage
    response = llm.invoke([SystemMessage(content=system), HumanMessage(content=user)])
    return response.content.strip()


def _strip_fences(text: str) -> str:
    text = re.sub(r"^```[a-z]*\n?", "", text.strip(), flags=re.MULTILINE)
    return re.sub(r"\n?```$", "", text.strip(), flags=re.MULTILINE).strip()


def _schema_text(profile: Dict[str, Any]) -> str:
    from backend.tabular_qa.agent import _schema_to_text
    return _schema_to_text(profile)


def _load_df(local_path: str, profile: Dict[str, Any]):
    """Load the DataFrame from the uploaded file path."""
    import pandas as pd
    sheet = (profile.get("summary") or {}).get("source_sheet")
    if local_path.endswith((".xlsx", ".xls")):
        return pd.read_excel(local_path, sheet_name=sheet or 0)
    return pd.read_csv(local_path, sep=None, engine="python")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def execute_analysis_plan(
    plan: Dict[str, Any],
    session_id: str,
    session_context: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Execute a tabular analysis plan and return interpreted findings.

    Returns
    -------
    {
        "status": "success" | "error",
        "text":          str  — biological interpretation (on success),
        "result_data":   Any  — raw computed results,
        "plot_path":     str | None — absolute path to a generated PNG,
        "plan_title":    str,
        "plan_goal":     str,
        "error":         str  — (on error only)
    }
    """
    from backend.tabular_qa.executor import execute_code

    # Resolve file(s) — prefer metadata embedded in the plan
    plan_files = plan.get("_files") or []
    if not plan_files:
        uploaded = session_context.get("uploaded_files", [])
        plan_files = [
            {
                "name": f.get("name") or f.get("original_filename"),
                "local_path": f.get("local_path"),
                "schema_preview": f.get("schema_preview"),
            }
            for f in uploaded
            if (f.get("schema_preview") or {}).get("family") == "tabular"
        ]

    if not plan_files:
        return {"status": "error", "error": "No tabular files available for execution."}

    file_info = plan_files[0]
    local_path: Optional[str] = file_info.get("local_path")
    profile: Dict[str, Any] = file_info.get("schema_preview") or {}
    filename: str = file_info.get("name") or "data"

    if not local_path or not Path(local_path).exists():
        return {"status": "error", "error": f"File not found on server: {local_path}"}

    try:
        df = _load_df(local_path, profile)
    except Exception as exc:
        return {"status": "error", "error": f"Could not load file: {exc}"}

    schema_text = _schema_text(profile)

    # Build step descriptions for the code prompt (exclude interpret steps)
    steps_text = "\n".join(
        f"{s['id']}. {s['name']}: {s.get('operation') or s.get('description', '')}"
        for s in (plan.get("steps") or [])
        if s.get("type") != "interpret"
    )

    code_prompt = (
        f"Analysis goal: {plan.get('goal', '')}\n\n"
        f"Dataset ({filename}):\n{schema_text}\n\n"
        f"Steps to implement:\n{steps_text}\n\n"
        f"Write the complete Python analysis script."
    )

    llm = _get_llm()

    # Remove any stale plot from a previous run
    _plot_tmp = "/tmp/helix_plot.png"
    if Path(_plot_tmp).exists():
        try:
            os.remove(_plot_tmp)
        except OSError:
            pass

    # Code generation → execution (one retry on failure)
    last_error = ""
    exec_result: Dict[str, Any] = {}
    for attempt in range(2):
        retry_suffix = (
            f"\n\nThe previous attempt raised:\n{last_error}\nFix the error."
            if attempt > 0 else ""
        )
        code = _strip_fences(_llm_call(llm, CODEGEN_SYSTEM_PROMPT, code_prompt + retry_suffix))
        exec_result = execute_code(code, df)

        if exec_result.get("error"):
            last_error = exec_result["error"]
            logger.warning(
                "[analysis_executor] code attempt %d failed: %s", attempt + 1, last_error
            )
        else:
            break
    else:
        return {
            "status": "error",
            "error": f"Analysis execution failed after 2 attempts: {last_error}",
        }

    # Persist plot and encode as base64 for inline display in the frontend
    plot_path: Optional[str] = None
    plot_base64: Optional[str] = None
    if Path(_plot_tmp).exists():
        dest_dir = Path(local_path).parent.parent / "plots"
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest = dest_dir / f"analysis_{session_id[:8]}.png"
        shutil.copy2(_plot_tmp, dest)
        plot_path = str(dest)
        try:
            import base64 as _b64
            with open(dest, "rb") as _pf:
                plot_base64 = "data:image/png;base64," + _b64.b64encode(_pf.read()).decode()
        except Exception:
            plot_base64 = None

    # Generate biological interpretation
    result_json = json.dumps(exec_result.get("result"), default=str, indent=2)
    if len(result_json) > 5000:
        result_json = result_json[:5000] + "\n…(truncated)"

    interp_prompt = (
        f"Analysis goal: {plan.get('goal', '')}\n\n"
        f"Dataset schema:\n{schema_text}\n\n"
        f"Computed results:\n{result_json}\n\n"
        f"Write the biological interpretation."
    )

    try:
        interpretation = _llm_call(llm, INTERPRETATION_SYSTEM_PROMPT, interp_prompt)
    except Exception as exc:
        logger.warning("[analysis_executor] interpretation failed: %s", exc)
        interpretation = f"Analysis completed. Results: {result_json[:500]}"

    return {
        "status": "success",
        "text": interpretation,
        "result_data": exec_result.get("result"),
        "plot_path": plot_path,
        "plot_base64": plot_base64,
        "plan_title": plan.get("title"),
        "plan_goal": plan.get("goal"),
    }
