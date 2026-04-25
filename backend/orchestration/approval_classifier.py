"""LLM-based approval intent classifier.

Replaces the hard-coded ``APPROVAL_COMMANDS`` keyword set with an LLM call
that understands any natural approval phrasing ("yes, looks good", "let's do
it", "go ahead", "execute that", "sounds right — run it", etc.).

Architecture
------------
1. **Keyword fast-path** – If the command exactly matches the legacy keyword
   set it is classified as an approval instantly, with no LLM call.  This
   covers the most common, unambiguous cases at zero cost.

2. **Early-rejection fast-path** – Commands that are clearly analytical
   (long sentences, bioinformatics terms, file references) are rejected
   without an LLM call.

3. **LLM path** – Short / ambiguous commands are sent to a fast, cheap
   classification call (≤ 5 tokens of output expected).  The system prompt is
   intentionally minimal.

4. **Fallback** – If the LLM is unavailable (``HELIX_MOCK_MODE=1``, no API
   key, timeout) the keyword set alone is used, preserving the prior behavior.

Public API
----------
``classify_approval(command, has_pending_plan=False) -> ApprovalDecision``
``is_approval_command(command, has_pending_plan=False) -> bool``  ← drop-in

The optional ``has_pending_plan`` hint biases the classifier: when there IS a
pending plan any short affirmative phrase ("sure", "ok", "yes") is much more
likely to be an approval than an analytical command.
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
    | \b(?:analyze|analyse|compute|calculate  # analytical verbs — safe to reject because
          |compare|align|normalize            # the affirmative-prefix guard runs first, so
          |differential|heatmap|violin        # "yes, generate..." still gets through
          |annotate|quantify|correlate        # but "Generate a volcano plot..." is rejected
          |cluster|generate|visualize|pca
          |plot|volcano|enrichment)\b
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


def _is_clearly_not_approval(command: str) -> bool:
    """Return True if the command is clearly an analytical request, not an approval.

    Skips early rejection when the command starts with an affirmative prefix,
    so phrases like "yes, execute the plan" or "approved, run the analysis"
    still reach the LLM path instead of being rejected as analytical commands.
    """
    text = (command or "").strip()
    if not text:
        return False
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
    method: str        # "keyword_fast_path" | "early_rejection" | "llm" | "keyword_fallback"
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

    # 3. LLM classification
    try:
        result = _classify_with_llm(command, has_pending_plan)
        return ApprovalDecision(
            is_approval=result,
            method="llm",
            reason="llm_classified",
        )
    except Exception as exc:
        logger.debug("Approval LLM classification failed (%s); falling back to keyword set.", exc)

    # 4. Keyword-only fallback (same behavior as before this module existed)
    return ApprovalDecision(
        is_approval=_keyword_match(command),
        method="keyword_fallback",
        reason="llm_unavailable",
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
