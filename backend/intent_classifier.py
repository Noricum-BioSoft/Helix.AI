from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Literal, Tuple

Intent = Literal["execute", "qa"]


@dataclass(frozen=True)
class IntentDecision:
    intent: Intent
    reason: str


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
    "fastqc",
    "fetch ",
    "download ",
    "visualize ",
    "create ",
    "generate ",
    "submit ",
    "build ",
)


def classify_intent(text: str) -> IntentDecision:
    """
    Heuristic intent classifier to prevent unintended tool generation.

    Rules of thumb:
    - If the user provides obvious inputs (S3 URIs / file paths / FASTA headers), assume execute.
    - If the user asks a question without any execution cues, assume Q&A.
    - Otherwise, default to execute (safer for command-style UI).
    """
    t = (text or "").strip()
    if not t:
        return IntentDecision(intent="qa", reason="empty")

    tl = t.lower()

    # Strong execute cues: explicit "execute/run" / data references / common bio file extensions.
    if "s3://" in tl:
        return IntentDecision(intent="execute", reason="s3_uri")
    if re.search(r"(^|\\s)(/[^\\s]+)", t):
        return IntentDecision(intent="execute", reason="local_path")
    if re.search(r"\\.(fastq|fq|fasta|fa|bam|sam|vcf|csv|tsv|json)(\\b|$)", tl):
        return IntentDecision(intent="execute", reason="file_extension")
    if re.search(r"(^|\\n)\\s*>\\S+", t):
        return IntentDecision(intent="execute", reason="fasta_header")

    # Question cues.
    if "?" in t:
        # Unless it's clearly a command with execute verbs.
        if any(v in tl for v in _EXEC_VERBS):
            return IntentDecision(intent="execute", reason="question_with_exec_verbs")
        return IntentDecision(intent="qa", reason="question_mark")

    if tl.startswith(_QUESTION_STARTERS):
        # If they start with question words but also contain clear imperative verbs, treat as execute.
        if any(v in tl for v in _EXEC_VERBS):
            return IntentDecision(intent="execute", reason="question_starter_with_exec_verbs")
        return IntentDecision(intent="qa", reason="question_starter")

    # Imperative cues.
    if any(v in tl for v in _EXEC_VERBS):
        return IntentDecision(intent="execute", reason="exec_verbs")

    # Default: treat as Q&A if it's short and looks like plain conversation; otherwise execute.
    if len(tl.split()) <= 5:
        return IntentDecision(intent="qa", reason="short_ambiguous")
    return IntentDecision(intent="execute", reason="default")


