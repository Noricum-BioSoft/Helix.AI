"""LLM-based approval intent classifier.

Classifies whether a user message is an approval of a pending analysis plan,
using a layered approach with no keyword fallback.

Architecture
------------
1. **Keyword fast-path** – Common, unambiguous approval phrases are detected
   instantly without an LLM call (e.g. "yes", "proceed", "go ahead").

2. **Early-rejection fast-path** – Clearly analytical commands (long sentences,
   bioinformatics terms, file references) are rejected without an LLM call.

3. **LLM path** – Short / ambiguous commands are sent to a lightweight
   classification call.  If the LLM is unavailable the exception propagates —
   there is intentionally no keyword fallback, because a silent wrong answer
   is worse than a visible error.

Public API
----------
``classify_approval(command, has_pending_plan=False) -> ApprovalDecision``
``is_approval_command(command, has_pending_plan=False) -> bool``

The optional ``has_pending_plan`` hint biases the classifier: when there IS a
pending plan any short affirmative phrase is much more likely to be an approval.
"""

from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass
from typing import Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Legacy keyword set – kept as fast-path / fallback
# ---------------------------------------------------------------------------

_KEYWORD_APPROVALS: frozenset[str] = frozenset({
    "approve",
    "approve.",
    "approved",
    "approved.",
    "yes",
    "yes.",
    "yes, approve",
    "yes approve",
    "yes, run it",
    "yes run it",
    "yes, proceed",
    "yes proceed",
    "yes, go ahead",
    "yes go ahead",
    "yes, looks good",
    "yes looks good",
    "go ahead",
    "go ahead.",
    "go for it",
    "sounds good",
    "sounds good.",
    "looks good",
    "looks good.",
    "ok",
    "ok.",
    "okay",
    "okay.",
    "sure",
    "sure.",
    "proceed",
    "proceed.",
    "run it",
    "run it.",
    "execute it",
    "execute.",
    "do it",
    "do it.",
    "let's do it",
    "let's go",
    "yep",
    "yep.",
    "yup",
    "yup.",
    "affirmative",
})


def _keyword_match(command: str) -> bool:
    """Return True if the command is an exact keyword approval."""
    normalized = " ".join((command or "").strip().lower().split())
    return normalized in _KEYWORD_APPROVALS or normalized.startswith("approve ")


# ---------------------------------------------------------------------------
# Early-rejection heuristics (no LLM needed)
# ---------------------------------------------------------------------------

_ANALYTICAL_PATTERNS = re.compile(
    r"""
    s3://                                     # S3 URIs
    | \.(?:fastq|fq|fasta|fa|bam|sam|vcf      # bioinformatics file extensions
          |csv|tsv|xlsx?|h5ad?)\b
    | (?:^|\s)>                               # FASTA headers
    | \b(?:analyze|analyse|compute|calculate  # analytical / execution verbs
          |compare|align|normalize            # Note: the affirmative-prefix guard runs first,
          |differential|heatmap|violin        # so "yes, generate a plot" still reaches the LLM
          |annotate|quantify|correlate        # but plain "Generate a volcano plot" is rejected.
          |cluster|generate|visualize|pca     # Action-command verbs (not approval verbs):
          |plots?|volcano|enrichment          # "trim these reads", "run FastQC", "fetch seq..."
          |trim|merge|fetch|download          # Modification/iteration verbs:
          |create|build|submit|perform        # "change the scale", "update the plot", etc.
          |run|execute                        # These are all clearly NOT approvals.
          |change|update|modify|adjust|fix
          |switch|replace|rename|remove|delete
          |scale|set|revert|rerun|patch)\b
    """,
    re.VERBOSE | re.IGNORECASE,
)

# Words that signal "I'm approving what was just shown" — if the command starts
# with one of these, we skip early rejection and let the LLM decide.
_AFFIRMATIVE_PREFIXES = re.compile(
    r"^(yes\b|yeah\b|yep\b|yup\b|ok\b|okay\b|sure\b|great\b|perfect\b|"
    r"approved?\b|go\s+ahead\b|go\s+for\b|sounds\b|looks\b|that\b|let'?s\b|"
    r"affirmative\b|proceed\b|do\s+it\b|run\s+it\b|execute\s+it\b)",
    re.IGNORECASE,
)

_WORD_COUNT_THRESHOLD = 15  # commands longer than this are unlikely approvals

# Interrogative words that open a question — questions are never approvals and
# never need an LLM call to classify.
_INTERROGATIVE_PREFIXES = re.compile(
    r"^(?:what|how|why|where|when|who|which|whose|is|are|am|was|were|"
    r"can|could|may|might|will|would|shall|should|do|does|did|has|have|had)\b",
    re.IGNORECASE,
)


def _is_clearly_not_approval(command: str) -> bool:
    """Return True if the command is clearly an analytical request, not an approval.

    Skips early rejection when the command starts with an affirmative prefix,
    so phrases like "yes, execute the plan" or "approved, run the analysis"
    still reach the LLM path instead of being rejected as analytical commands.

    Questions (interrogative opener or trailing `?`) are ALWAYS rejected here
    — they never require LLM confirmation.
    """
    text = (command or "").strip()
    if not text:
        return False
    # Interrogative questions are never approvals.
    if _INTERROGATIVE_PREFIXES.match(text) or text.endswith("?"):
        return True
    # Don't early-reject short affirmative phrases even if they contain
    # action verbs — the user is consenting to an existing plan, not issuing a
    # new analytical command.
    if _AFFIRMATIVE_PREFIXES.match(text):
        return False
    if len(text.split()) > _WORD_COUNT_THRESHOLD:
        return True
    if _ANALYTICAL_PATTERNS.search(text):
        return True
    return False


# ---------------------------------------------------------------------------
# LLM classification
# ---------------------------------------------------------------------------

_SYSTEM_PROMPT = """\
You are an approval-intent classifier for a bioinformatics AI assistant.

The assistant has just shown the user a proposed analysis plan and is waiting
for the user to approve it before executing.

Your job: decide whether the user's message is an APPROVAL of that plan.

An APPROVAL is any message whose clear intent is to say "yes, go ahead and
execute the plan as proposed".  This includes:
- Explicit confirmations: "yes", "approve", "proceed", "go ahead", "ok", "sure"
- Natural affirmations: "looks good", "that works", "sounds right", "let's do it"
- Implicit consent: "run it", "execute", "start", "do it"

NOT an approval:
- A new analytical request ("analyze my data for X instead")
- A question about the plan ("what does step 3 do?")
- A rejection or revision ("no, use a different method")
- Unrelated commands

Respond with ONLY a JSON object on a single line:
{"approval": true}   or   {"approval": false}

No explanation. No markdown. Just the JSON object.
"""


def _get_llm():
    """Lazy LLM initializer — same pattern as intent_classifier._get_llm()."""
    if os.getenv("HELIX_MOCK_MODE") == "1":
        raise RuntimeError("LLM disabled in HELIX_MOCK_MODE")

    openai_key = os.getenv("OPENAI_API_KEY", "").strip()
    deepseek_key = os.getenv("DEEPSEEK_API_KEY", "").strip()

    openai_enabled = openai_key not in ("", "disabled", "your_openai_api_key_here", "none")
    deepseek_enabled = deepseek_key not in ("", "disabled", "your_deepseek_api_key_here", "none")

    if openai_enabled:
        from langchain.chat_models import init_chat_model
        model = os.getenv("HELIX_APPROVAL_OPENAI_MODEL", "openai:gpt-4.1-mini").strip()
        if ":" not in model:
            model = f"openai:{model}"
        return init_chat_model(model, temperature=0, max_tokens=10)

    if deepseek_enabled:
        from langchain_deepseek import ChatDeepSeek
        return ChatDeepSeek(model="deepseek-chat", temperature=0, max_tokens=10,
                            timeout=8.0, max_retries=1)

    raise ValueError("No LLM API key configured for approval classification.")


def _classify_with_llm(command: str, has_pending_plan: bool) -> bool:
    """Ask the LLM whether the command is an approval. Returns True/False."""
    context = (
        "[Context: the user has a pending analysis plan waiting for approval.]\n\n"
        if has_pending_plan
        else ""
    )
    messages = [
        {"role": "system", "content": _SYSTEM_PROMPT},
        {"role": "user", "content": f"{context}User message: {command!r}"},
    ]
    llm = _get_llm()
    response = llm.invoke(messages)
    raw = (response.content or "").strip()

    # Parse {"approval": true/false} — be defensive
    import json
    start = raw.find("{")
    end = raw.rfind("}") + 1
    if start >= 0 and end > start:
        obj = json.loads(raw[start:end])
        return bool(obj.get("approval", False))
    # Tolerate bare "true"/"false"
    return raw.lower().startswith("true")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ApprovalDecision:
    is_approval: bool
    method: str        # "keyword_fast_path" | "early_rejection" | "llm"
    reason: str


def classify_approval(
    command: str,
    *,
    has_pending_plan: bool = False,
) -> ApprovalDecision:
    """Classify whether *command* is an approval of a pending plan.

    Parameters
    ----------
    command:
        The raw user command string.
    has_pending_plan:
        True when the session currently has a plan waiting for approval.
        Used to bias the classifier: any short affirmative is much more
        likely to be an approval in this context.
    """
    # 1. Keyword fast-path (free, instant)
    if _keyword_match(command):
        return ApprovalDecision(
            is_approval=True,
            method="keyword_fast_path",
            reason=f"matched_keyword: {command!r}",
        )

    # 2. Early rejection (free, instant)
    if _is_clearly_not_approval(command):
        return ApprovalDecision(
            is_approval=False,
            method="early_rejection",
            reason="analytical_pattern_or_too_long",
        )

    # 3. LLM classification — raises on failure; no keyword fallback.
    result = _classify_with_llm(command, has_pending_plan)
    return ApprovalDecision(
        is_approval=result,
        method="llm",
        reason="llm_classified",
    )


def is_approval_command(
    command: str,
    *,
    has_pending_plan: bool = False,
) -> bool:
    """Drop-in replacement for the legacy ``is_approval_command`` in approval_policy.py.

    Returns True if the command's intent is to approve a pending plan.
    """
    decision = classify_approval(command, has_pending_plan=has_pending_plan)
    logger.debug(
        "approval_classifier: command=%r → is_approval=%s method=%s reason=%s",
        command[:80],
        decision.is_approval,
        decision.method,
        decision.reason,
    )
    return decision.is_approval
