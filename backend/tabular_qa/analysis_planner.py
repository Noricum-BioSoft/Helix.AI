"""
Analysis planner for tabular data.

Given a user command and uploaded file profiles, uses an LLM to generate
a structured analysis plan that is presented to the user for approval before
any code is executed.
"""
from __future__ import annotations

import json
import logging
import os
import re
from typing import Any, Dict, List

logger = logging.getLogger(__name__)

_PLANNER_MODEL_ENV = "HELIX_ANALYSIS_PLANNER_MODEL"
_DEFAULT_MODEL = "openai:gpt-4o"

PLANNER_SYSTEM_PROMPT = """\
You are an expert bioinformatics analyst for Helix.AI, an autonomous data analysis platform.

Given a user's analysis request and the profile of their uploaded tabular data, generate a
concise, executable analysis plan.

The plan MUST:
1. Have a clear biological/scientific goal stated in one sentence
2. Break the analysis into 3–6 concrete, ordered steps
3. Each step must specify a distinct computation or visualisation
4. Use only pandas, numpy, scipy, matplotlib, seaborn — no external databases
5. Be achievable with the available data columns

Return ONLY valid JSON — no markdown, no preamble, no explanation:
{
  "title": "Brief title (max 10 words)",
  "goal": "One sentence: the scientific question being answered",
  "steps": [
    {
      "id": 1,
      "name": "Step name (3–6 words)",
      "description": "What this step does and why it matters",
      "type": "load | filter | compute | visualize | interpret",
      "operation": "Technical detail e.g. 'Pearson correlation matrix across all numeric columns'"
    }
  ],
  "expected_outputs": ["output 1", "output 2"]
}

Step type meanings:
- load     : data loading and quality validation
- filter   : subsetting rows or columns
- compute  : statistics (correlation, clustering, DE, PCA, etc.)
- visualize: plot generation (heatmap, scatter, bar, etc.)
- interpret: LLM-generated biological summary of findings
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
        return init_chat_model(os.getenv(_PLANNER_MODEL_ENV, _DEFAULT_MODEL), temperature=0.1)

    if deepseek_key and deepseek_key not in bad:
        from langchain_deepseek import ChatDeepSeek
        return ChatDeepSeek(model="deepseek-chat", temperature=0.1, max_tokens=1024)

    raise ValueError("No LLM API key configured for analysis planner.")


def _build_file_context(uploaded_files: List[Dict[str, Any]]) -> str:
    """Return a concise text description of all uploaded tabular files."""
    parts: List[str] = []
    for f in uploaded_files:
        sp = f.get("schema_preview") or {}
        if sp.get("family") != "tabular":
            continue
        name = f.get("name") or f.get("original_filename") or "unknown"
        summary = sp.get("summary") or {}
        n_rows = summary.get("n_rows") or sp.get("n_records") or "?"
        cols_data: List[Dict] = (sp.get("schema") or {}).get("columns") or []
        n_cols = len(cols_data)
        col_names = [c["name"] for c in cols_data if isinstance(c, dict)][:20]
        sheet = summary.get("source_sheet")
        sheet_info = f" [sheet: {sheet}]" if sheet else ""

        interesting: List[str] = []
        for c in cols_data[:30]:
            if not isinstance(c, dict):
                continue
            tv = c.get("top_values")
            if tv and len(tv) <= 10:
                vals = list(tv.keys())[:5]
                interesting.append(f"  - {c['name']} (categorical: {vals})")
            elif c.get("dtype") in ("float64", "float32", "int64", "int32"):
                interesting.append(
                    f"  - {c['name']} (numeric, mean={c.get('mean','?')}, "
                    f"min={c.get('min','?')}, max={c.get('max','?')})"
                )

        part = (
            f"File: {name}{sheet_info}\n"
            f"Shape: {n_rows} rows × {n_cols} columns\n"
            f"Columns: {col_names}"
        )
        if interesting:
            part += "\nKey columns:\n" + "\n".join(interesting[:10])
        parts.append(part)

    return "\n\n".join(parts) if parts else "No tabular files available."


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------


def is_analytical_request(command: str, session_context: Dict[str, Any]) -> bool:
    """
    Return True when the command is a complex analysis request that needs a plan
    (vs a simple question that the code interpreter can answer directly).
    """
    uploaded_files = session_context.get("uploaded_files", [])
    has_tabular = any(
        (f.get("schema_preview") or {}).get("family") == "tabular"
        for f in uploaded_files
    )
    if not has_tabular:
        return False

    c = command.lower().strip()

    # Always pass short lookup-style questions straight to code interpreter
    simple_patterns = [
        r"^how many (rows|columns|records|samples|patients|genes|entries)",
        r"^what (are|is) the (columns?|column names?|fields?|schema|headers?|data types?|dtypes?)",
        r"^(show|display|print|list) (me )?(the )?(first|top|last|sample|head|tail|all columns?)",
        r"^(count|sum|mean|average|max|min|median|std) (of |the )?",
        r"^(describe|summarize|summary) (the )?(data(set|frame)?|file|table)?$",
        r"^(is |are )?there (any )?(missing|null|nan|empty)",
    ]
    for pat in simple_patterns:
        if re.search(pat, c):
            return False

    # Analytical trigger keywords
    analytical_cues = (
        "correlat", "cluster", "differenti", "compare", "comparison",
        "find gene", "which gene", "identify", "analyze", "analyse",
        "predict", "classif", "enrichment", "pathway",
        "pca", "principal component",
        "heatmap", "visualize", "visualise",
        "association", "regression", "survival", "cohort",
        "treatment response", "biomarker",
        "expression pattern", "top genes", "ranked", "rank",
        "distribution of", "over time", "trend",
        "group by", "stratif", "segment",
        "feature importan", "variable importan",
    )
    for cue in analytical_cues:
        if cue in c:
            return True

    # Long commands on tabular data almost always need planning
    return len(c.split()) >= 8


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


async def plan_analysis(
    command: str,
    session_context: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Generate a structured analysis plan from a user's natural-language request
    and the profile of their uploaded tabular files.

    Returns
    -------
    {
        "status": "success" | "error" | "no_tabular_files",
        "plan":  { title, goal, steps, expected_outputs, type, _files }  # on success
        "error": str                                                       # on error
    }
    """
    uploaded_files = session_context.get("uploaded_files", [])
    tabular_files = [
        f for f in uploaded_files
        if (f.get("schema_preview") or {}).get("family") == "tabular"
    ]

    if not tabular_files:
        return {"status": "no_tabular_files", "error": "No tabular files uploaded."}

    file_context = _build_file_context(uploaded_files)

    user_message = (
        f'User request: "{command}"\n\n'
        f"Uploaded data:\n{file_context}\n\n"
        f"Generate the analysis plan as JSON."
    )

    try:
        llm = _get_llm()
        from langchain_core.messages import HumanMessage, SystemMessage

        response = llm.invoke([
            SystemMessage(content=PLANNER_SYSTEM_PROMPT),
            HumanMessage(content=user_message),
        ])

        content = response.content.strip()
        # Strip markdown fences if LLM adds them
        content = re.sub(r"^```[a-z]*\n?", "", content, flags=re.MULTILINE)
        content = re.sub(r"\n?```$", "", content.strip(), flags=re.MULTILINE).strip()

        plan = json.loads(content)

        # Normalise required fields
        plan["type"] = "tabular_analysis"
        plan.setdefault("expected_outputs", [])

        # Attach file metadata so the executor doesn't need to re-resolve paths
        plan["_files"] = [
            {
                "name": f.get("name") or f.get("original_filename"),
                "local_path": f.get("local_path"),
                "schema_preview": f.get("schema_preview"),
            }
            for f in tabular_files
        ]

        logger.info(
            "[analysis_planner] plan generated: title=%r steps=%d",
            plan.get("title"), len(plan.get("steps", [])),
        )
        return {"status": "success", "plan": plan}

    except json.JSONDecodeError as exc:
        logger.error("[analysis_planner] LLM returned non-JSON: %s", exc)
        return {"status": "error", "error": "Plan generation failed: LLM returned non-JSON."}
    except Exception as exc:
        logger.error("[analysis_planner] unexpected error: %s", exc, exc_info=True)
        return {"status": "error", "error": str(exc)}
