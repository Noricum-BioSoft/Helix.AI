from __future__ import annotations

import json
import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, List

logger = logging.getLogger(__name__)

Intent = Literal["execute", "qa"]


@dataclass(frozen=True)
class IntentDecision:
    intent: Intent
    reason: str


# Load agent prompt from markdown file
PROJECT_ROOT = Path(__file__).resolve().parent.parent
INTENT_DETECTOR_PROMPT_PATH = PROJECT_ROOT / "agents" / "intent-detector-agent.md"

try:
    if not INTENT_DETECTOR_PROMPT_PATH.exists():
        raise FileNotFoundError(f"Intent detector prompt file not found at: {INTENT_DETECTOR_PROMPT_PATH}")
    
    INTENT_DETECTOR_SYSTEM_PROMPT = INTENT_DETECTOR_PROMPT_PATH.read_text(encoding='utf-8')
    
    if not INTENT_DETECTOR_SYSTEM_PROMPT or not INTENT_DETECTOR_SYSTEM_PROMPT.strip():
        raise ValueError(f"Intent detector prompt file is empty: {INTENT_DETECTOR_PROMPT_PATH}")
    
    logger.info(f"Loaded intent detector prompt from {INTENT_DETECTOR_PROMPT_PATH} ({len(INTENT_DETECTOR_SYSTEM_PROMPT)} chars)")
except Exception as e:
    logger.warning(f"Could not load intent-detector-agent.md from {INTENT_DETECTOR_PROMPT_PATH}: {e}")
    logger.warning("Using fallback prompt")
    INTENT_DETECTOR_SYSTEM_PROMPT = """You are an Intent Detector. Given a single user prompt, classify the user's intent into a small set of multi-label categories. Return only a JSON object with "prompt" and "intent" (array of labels: "question", "action", "data", "workflow").
"""


# Heuristic classifier constants (used as fallback)
_QUESTION_STARTERS = (
    "what ",
    "why ",
    "how ",
    "explain ",
    "tell me ",
    "describe ",
    "compare ",
    "difference ",
    "help me understand ",
    "can you explain",
    "is it ",
    "are there ",
)

_EXEC_VERBS = (
    "run ",
    "execute ",
    "perform ",
    "analyze ",
    "analyse ",
    "align ",
    "mutate ",
    "trim ",
    "merge ",
    "fetch ",
    "download ",
    "visualize ",
    "create ",
    "generate ",
    "submit ",
    "build ",
    "calculate ",
    "compute ",
)


def _get_llm():
    """
    Lazily initialize and return the LLM for intent classification.
    
    Uses the same pattern as backend.agent._get_llm() for consistency.
    Handles HELIX_MOCK_MODE and supports both OpenAI and DeepSeek.
    """
    if os.getenv("HELIX_MOCK_MODE") == "1":
        raise RuntimeError("LLM is disabled in HELIX_MOCK_MODE")

    openai_key = os.getenv("OPENAI_API_KEY", "").strip()
    deepseek_key = os.getenv("DEEPSEEK_API_KEY", "").strip()

    openai_enabled = openai_key and openai_key not in ["", "disabled", "your_openai_api_key_here", "none"]
    deepseek_enabled = deepseek_key and deepseek_key not in ["", "disabled", "your_deepseek_api_key_here", "none"]

    if openai_enabled:
        # Local import to avoid loading heavy deps at module import time
        from langchain.chat_models import init_chat_model
        # Load model name from .env, default to gpt-5.2
        # LangChain init_chat_model expects provider-prefixed names like "openai:<model>"
        openai_model = os.getenv("HELIX_INTENT_OPENAI_MODEL", "openai:gpt-5.2").strip()
        # Be forgiving: allow HELIX_INTENT_OPENAI_MODEL="gpt-5.2" as well (add prefix if missing)
        if ":" not in openai_model:
            openai_model = f"openai:{openai_model}"
        return init_chat_model(openai_model, temperature=0)

    if deepseek_enabled:
        # Local import to avoid loading heavy deps at module import time
        from langchain_deepseek import ChatDeepSeek
        return ChatDeepSeek(
            model="deepseek-chat",
            temperature=0,
            max_tokens=None,
            timeout=30.0,  # Shorter timeout for classification task
            max_retries=1,
        )

    raise ValueError(
        "No API keys found for intent classification. "
        "Please set either OPENAI_API_KEY or DEEPSEEK_API_KEY "
        "in your environment variables, or the classifier will fall back to heuristics."
    )


def _map_intent_labels_to_binary(intent_labels: List[str]) -> tuple[Intent, str]:
    """
    Map multi-label intent classification to binary decision (execute or qa).
    
    Args:
        intent_labels: List of intent labels from LLM (e.g., ["question", "action", "data"])
        
    Returns:
        Tuple of (intent, reason) where:
        - intent: "execute" or "qa"
        - reason: String describing the classification reason
    """
    has_action = "action" in intent_labels
    has_question = "question" in intent_labels
    
    if has_action:
        # If action is present, it's an execute request
        intent: Intent = "execute"
        reason = f"llm_classified_action_{','.join(intent_labels)}"
    elif has_question and not has_action:
        # Only question, no action → Q&A
        intent = "qa"
        reason = f"llm_classified_question_{','.join(intent_labels)}"
    else:
        # Default to execute for ambiguous cases (safer for command-style UI)
        intent = "execute"
        reason = f"llm_classified_default_{','.join(intent_labels) if intent_labels else 'no_labels'}"
    
    return intent, reason


def _classify_intent_with_llm(text: str) -> IntentDecision:
    """
    Classify intent using LLM agent based on intent-detector-agent.md prompt.
    
    Returns IntentDecision with intent ("execute" or "qa") and reason.
    Maps multi-label output from agent to binary decision:
    - If "action" is in labels → "execute"
    - If only "question" (no "action") → "qa"
    - Otherwise → "execute" (safer default for command-style UI)
    """
    try:
        llm = _get_llm()
        
        # Build the classification request
        # The system prompt defines the task and format. Pass the user prompt directly
        # as the agent description specifies - it should return JSON with the prompt included
        user_message = text.strip()
        
        # Invoke LLM with system prompt and user prompt
        # LangChain ChatModel accepts messages as dicts with "role" and "content"
        # The system prompt (INTENT_DETECTOR_SYSTEM_PROMPT) instructs the LLM to return JSON
        # with "prompt" (the original user message) and "intent" (array of labels)
        
        # Debug: Verify prompt is loaded
        logger.debug(f"Loaded system prompt from {INTENT_DETECTOR_PROMPT_PATH}: {INTENT_DETECTOR_SYSTEM_PROMPT}")

        if not INTENT_DETECTOR_SYSTEM_PROMPT or not INTENT_DETECTOR_SYSTEM_PROMPT.strip():
            logger.error(f"INTENT_DETECTOR_SYSTEM_PROMPT is empty! Path was: {INTENT_DETECTOR_PROMPT_PATH}")
            raise ValueError("System prompt is empty - cannot classify intent")
        
        logger.debug(f"Using system prompt ({len(INTENT_DETECTOR_SYSTEM_PROMPT)} chars) for intent classification")
        
        messages = [
            {"role": "system", "content": INTENT_DETECTOR_SYSTEM_PROMPT},
            {"role": "user", "content": user_message}
        ]
        
        response = llm.invoke(messages)
        response_content = response.content.strip()
        
        # Parse JSON response
        # Try to extract JSON from response (handle cases where LLM adds extra text)
        json_start = response_content.find("{")
        json_end = response_content.rfind("}") + 1
        
        if json_start >= 0 and json_end > json_start:
            json_str = response_content[json_start:json_end]
            result = json.loads(json_str)
        else:
            # Fallback: try parsing the whole response
            result = json.loads(response_content)
        
        # Extract intent labels
        intent_labels: List[str] = result.get("intent", [])
        if not isinstance(intent_labels, list):
            intent_labels = []
        
        # Map multi-label to binary decision
        intent, reason = _map_intent_labels_to_binary(intent_labels)
        
        logger.info(f"Intent classified by LLM: {intent} (labels: {intent_labels}, reason: {reason})")
        return IntentDecision(intent=intent, reason=reason)
        
    except Exception as e:
        logger.warning(f"LLM intent classification failed: {e}, falling back to heuristics")
        raise  # Re-raise to trigger fallback


def _classify_intent_with_heuristics(text: str) -> IntentDecision:
    """
    Heuristic intent classifier (fallback when LLM is unavailable).
    
    Rules of thumb:
    - If the user provides obvious inputs (S3 URIs / file paths / FASTA headers), assume execute.
    - If the user asks a question without any execution cues, assume Q&A.
    - Otherwise, default to execute (safer for command-style UI).
    """
    t = (text or "").strip()
    if not t:
        return IntentDecision(intent="qa", reason="heuristic_empty")

    tl = t.lower()

    # Strong execute cues: explicit "execute/run" / data references / common bio file extensions.
    if "s3://" in tl:
        return IntentDecision(intent="execute", reason="heuristic_s3_uri")
    if re.search(r"(^|\s)(/[^\s]+)", t):
        return IntentDecision(intent="execute", reason="heuristic_local_path")
    if re.search(r"\.(fastq|fq|fasta|fa|bam|sam|vcf|csv|tsv|json)(\b|$)", tl):
        return IntentDecision(intent="execute", reason="heuristic_file_extension")
    if re.search(r"(^|\n)\s*>\S+", t):
        return IntentDecision(intent="execute", reason="heuristic_fasta_header")

    # Question cues.
    if "?" in t:
        # If this is an explanatory question ("what/why/how..."), treat as Q&A unless there are
        # strong execute cues (paths/URIs/FASTA already handled above).
        if tl.startswith(_QUESTION_STARTERS):
            # Allow "Can you ...?" style questions to be treated as execute requests.
            if any(v in tl for v in _EXEC_VERBS) and any(
                p in tl for p in ("can you ", "could you ", "please ")
            ):
                return IntentDecision(intent="execute", reason="heuristic_question_starter_polite_exec_request")
            return IntentDecision(intent="qa", reason="heuristic_question_mark_explanatory")

        # Non-explanatory questions like "Align these?" are often commands phrased as a question.
        if any(v in tl for v in _EXEC_VERBS):
            return IntentDecision(intent="execute", reason="heuristic_question_with_exec_verbs")
        return IntentDecision(intent="qa", reason="heuristic_question_mark")

    if tl.startswith(_QUESTION_STARTERS):
        # If they start with question words but also contain clear imperative verbs, treat as execute.
        if any(v in tl for v in _EXEC_VERBS):
            return IntentDecision(intent="execute", reason="heuristic_question_starter_with_exec_verbs")
        return IntentDecision(intent="qa", reason="heuristic_question_starter")

    # Imperative cues.
    if any(v in tl for v in _EXEC_VERBS):
        return IntentDecision(intent="execute", reason="heuristic_exec_verbs")

    # Default: treat as Q&A if it's short and looks like plain conversation; otherwise execute.
    if len(tl.split()) <= 5:
        return IntentDecision(intent="qa", reason="heuristic_short_ambiguous")
    return IntentDecision(intent="execute", reason="heuristic_default")


def classify_intent(
    text: str,
    *,
    session_context: Optional[dict] = None,
    workflow_state: Optional[str] = None,
) -> IntentDecision:
    """
    Classify user intent into "execute" or "qa" using session-aware rules first,
    then LLM, then heuristics.

    Session-aware shortcuts (applied before any LLM call):
    - WAITING_FOR_APPROVAL + short affirmative-looking text → "execute" (approval)
    - WAITING_FOR_CLARIFICATION → "execute" (user is answering a question)
    - WAITING_FOR_INPUTS → "execute" (user is supplying inputs)
    - Prior runs in session + "again / rerun / redo" → "execute"
    """
    from typing import Optional  # noqa: PLC0415 – local to avoid circular at module level

    # ── Checkpoint-state shortcuts ────────────────────────────────────────────
    if workflow_state:
        _ws = (workflow_state or "").upper()
        if _ws in {"WAITING_FOR_APPROVAL", "WAITING_FOR_CLARIFICATION", "WAITING_FOR_INPUTS"}:
            return IntentDecision(intent="execute", reason=f"checkpoint_state_{_ws.lower()}")

    # ── Session-history shortcuts ─────────────────────────────────────────────
    if session_context:
        _runs = session_context.get("runs") or []
        _history = session_context.get("history") or []
        _has_prior_runs = bool(_runs or _history)
        _tl = (text or "").lower().strip()
        _rerun_cues = ("again", "rerun", "re-run", "redo", "try again", "retry",
                       "same but", "this time", "now run", "update it")
        if _has_prior_runs and any(c in _tl for c in _rerun_cues):
            return IntentDecision(intent="execute", reason="session_rerun_cue")

    # ── LLM then heuristic fallback ───────────────────────────────────────────
    try:
        return _classify_intent_with_llm(text)
    except Exception as e:
        logger.debug(f"Falling back to heuristic classification: {e}")
        return _classify_intent_with_heuristics(text)


# Keep Optional importable at module level for type annotations
from typing import Optional  # noqa: E402
