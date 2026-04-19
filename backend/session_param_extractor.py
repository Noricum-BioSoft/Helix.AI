"""
Session-aware parameter extraction for user requests.

Extracts parameters from session history (previous runs, produced artifacts,
pipeline outputs) when the user's command references them implicitly
(e.g. "merged reads", "the output", "previous results").
"""

from __future__ import annotations

import logging
import re
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


def _get_history(session_context: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Get history entries from session, normalizing structure."""
    if not session_context:
        return []
    history = session_context.get("history") or []
    if isinstance(history, list):
        return history
    return []


def _get_artifacts(session_context: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Get artifact list from session (from runs or artifacts registry)."""
    if not session_context:
        return []
    arts = []
    # From artifacts registry
    reg = session_context.get("artifacts") or {}
    if isinstance(reg, dict):
        arts.extend(reg.values())
    # From runs' produced_artifacts
    for run in session_context.get("runs") or []:
        if isinstance(run, dict):
            arts.extend(run.get("produced_artifacts") or [])
    # From history entries' result.produced_artifacts / result.links
    for entry in _get_history(session_context):
        result = entry.get("result") or {}
        if isinstance(result, dict):
            arts.extend(result.get("produced_artifacts") or [])
            for link in result.get("links") or []:
                if isinstance(link, dict) and link.get("url") or link.get("uri"):
                    arts.append({
                        "uri": link.get("url") or link.get("uri"),
                        "type": link.get("type", "file"),
                        "title": link.get("label", link.get("title", "")),
                        "format": link.get("format"),
                    })
    return arts


def _uri_looks_like_fastq(uri: str) -> bool:
    """Check if URI suggests FASTQ format."""
    if not uri or not isinstance(uri, str):
        return False
    u = uri.lower()
    return bool(re.search(r"\.(?:fastq|fq)(?:\.gz)?(?:$|[?#])", u))


def _uri_looks_like_fasta(uri: str) -> bool:
    """Check if URI suggests FASTA format."""
    if not uri or not isinstance(uri, str):
        return False
    u = uri.lower()
    return bool(re.search(r"\.(?:fasta|fa|fas|fna)(?:\.gz)?(?:$|[?#])", u))


def _extract_uris_from_result(result: Any) -> List[str]:
    """Recursively extract S3/local URIs from a result dict."""
    uris = []
    if isinstance(result, str) and (result.startswith("s3://") or result.startswith("/")):
        uris.append(result)
    elif isinstance(result, dict):
        for k, v in result.items():
            if k in ("uri", "url", "output", "output_uri", "merged_output", "path"):
                if isinstance(v, str) and (v.startswith("s3://") or v.startswith("/")):
                    uris.append(v)
            uris.extend(_extract_uris_from_result(v))
    elif isinstance(result, list):
        for item in result:
            uris.extend(_extract_uris_from_result(item))
    return uris


def get_fastqc_inputs_from_session(
    session_context: Dict[str, Any],
    command: str,
) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
    """
    Extract FastQC input_r1, input_r2, output from session when user references
    "merged reads", "previous output", "pipeline results", etc.

    Returns:
        (r1_uri, r2_uri, output_uri, explanation_or_error)
        - If found: (r1, r2, output, None)
        - If merged/FASTA detected and incompatible: (None, None, None, error_message)
        - If not found: (None, None, None, None)
    """
    cmd_lower = (command or "").lower()
    history = _get_history(session_context)
    artifacts = _get_artifacts(session_context)

    # User explicitly asked for "merged reads" or "merged output" or gave merged.fasta
    merged_ref = bool(
        re.search(r"\bmerged\s*(?:reads?|output|file|sequences?)\b", cmd_lower)
        or "merged.fasta" in cmd_lower
        or "merged.fa" in cmd_lower
    )
    single_uri_in_cmd = None
    for m in re.finditer(r"s3://[^\s]+|/[^\s]+\.(?:fasta?|fq|fastq)(?:\.gz)?", command):
        single_uri_in_cmd = m.group(0).strip()
        break

    if merged_ref or (single_uri_in_cmd and _uri_looks_like_fasta(single_uri_in_cmd)):
        # User wants FastQC on merged output - incompatible
        merged_uri = single_uri_in_cmd
        if not merged_uri:
            for a in artifacts:
                u = a.get("uri") or ""
                if _uri_looks_like_fasta(u) and "merged" in u.lower():
                    merged_uri = u
                    break
        return (
            None,
            None,
            None,
            "FastQC cannot run on merged reads. Merged output is in **FASTA** format, which has no "
            "per-base quality scores. FastQC requires **paired FASTQ** files (R1 and R2) with "
            "Phred quality values.\n\n"
            "**What you can do:**\n"
            "1. Run FastQC on the **original raw reads** (R1 and R2) before merging — "
            "that will show base quality, adapter content, and per-sequence metrics.\n"
            "2. For post-merge QC, use tools designed for FASTA (e.g. sequence length distribution, "
            "chimera checks) — Helix does not yet offer MultiQC or merged-FASTA QC."
        )

    # Look for R1/R2 from pipeline steps and command text
    r1_uri = None
    r2_uri = None
    output_uri = None

    # First: extract from command text in history (e.g. "R1: s3://... R2: s3://...")
    for entry in reversed(history):
        cmd = entry.get("command") or ""
        if "s3://" in cmd:
            m1 = re.search(r"(?:r1|forward|read\s*1)[^:]*:\s*(s3://[^\s]+)", cmd, re.IGNORECASE)
            m2 = re.search(r"(?:r2|reverse|read\s*2)[^:]*:\s*(s3://[^\s]+)", cmd, re.IGNORECASE)
            if m1 and _uri_looks_like_fastq(m1.group(1)):
                r1_uri = r1_uri or m1.group(1).strip()
            if m2 and _uri_looks_like_fastq(m2.group(1)):
                r2_uri = r2_uri or m2.group(1).strip()
            if not r1_uri:
                m = re.search(r"s3://[^\s]*(?:R1|r1|_1\.fq|mate_R1)[^\s]*", cmd, re.IGNORECASE)
                if m and _uri_looks_like_fastq(m.group(0)):
                    r1_uri = r1_uri or m.group(0).strip()
            if not r2_uri:
                m = re.search(r"s3://[^\s]*(?:R2|r2|_2\.fq|mate_R2)[^\s]*", cmd, re.IGNORECASE)
                if m and _uri_looks_like_fastq(m.group(0)):
                    r2_uri = r2_uri or m.group(0).strip()
        if r1_uri and r2_uri:
            break

    # Second: from result structure (steps, args, nested result)
    if not r1_uri or not r2_uri:
        for entry in reversed(history):
            result = entry.get("result") or {}
            if not isinstance(result, dict):
                continue
            meta = entry.get("metadata") or {}
            args = meta.get("tool_args") or result.get("arguments") or {}
            steps = result.get("steps") or []
            for step in steps:
                if isinstance(step, dict):
                    a = step.get("arguments") or step.get("tool_args") or {}
                    if isinstance(a, dict):
                        u1 = a.get("forward_reads") or a.get("input_r1")
                        u2 = a.get("reverse_reads") or a.get("input_r2")
                        if isinstance(u1, str) and u1.startswith("s3://") and _uri_looks_like_fastq(u1):
                            r1_uri = r1_uri or u1
                        if isinstance(u2, str) and u2.startswith("s3://") and _uri_looks_like_fastq(u2):
                            r2_uri = r2_uri or u2
            if isinstance(args, dict):
                u1 = args.get("forward_reads") or args.get("input_r1")
                u2 = args.get("reverse_reads") or args.get("input_r2")
                if isinstance(u1, str) and u1.startswith("s3://") and _uri_looks_like_fastq(u1):
                    r1_uri = r1_uri or u1
                if isinstance(u2, str) and u2.startswith("s3://") and _uri_looks_like_fastq(u2):
                    r2_uri = r2_uri or u2
            inner = result.get("result") or result
            if isinstance(inner, dict):
                u1 = inner.get("forward_reads") or inner.get("input_r1")
                u2 = inner.get("reverse_reads") or inner.get("input_r2")
                if isinstance(u1, str) and u1.startswith("s3://") and _uri_looks_like_fastq(u1):
                    r1_uri = r1_uri or u1
                if isinstance(u2, str) and u2.startswith("s3://") and _uri_looks_like_fastq(u2):
                    r2_uri = r2_uri or u2
            for u in _extract_uris_from_result(result):
                if _uri_looks_like_fastq(u):
                    if not r1_uri and ("r1" in u.lower() or "_1" in u or "mate_r1" in u.lower()):
                        r1_uri = u
                    elif not r2_uri and ("r2" in u.lower() or "_2" in u or "mate_r2" in u.lower()):
                        r2_uri = u
            if r1_uri and r2_uri:
                break

    # Also scan command for R1/R2 in the original prompt (e.g. from pipeline description)
    if not r1_uri or not r2_uri:
        for a in artifacts:
            u = (a.get("uri") or "").strip()
            if not u or not _uri_looks_like_fastq(u):
                continue
            if not r1_uri and ("r1" in u.lower() or "_1." in u.lower() or "mate_r1" in u.lower()):
                r1_uri = u
            elif not r2_uri and ("r2" in u.lower() or "_2." in u.lower() or "mate_r2" in u.lower()):
                r2_uri = u
            if r1_uri and r2_uri:
                break

    if r1_uri and r2_uri:
        return (r1_uri, r2_uri, output_uri, None)

    return (None, None, None, None)
