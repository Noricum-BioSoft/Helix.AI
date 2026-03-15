"""
Format and structure validation for tool inputs.

Validates file formats (FASTQ vs FASTA), required structure (paired vs single),
and returns clear, actionable error messages.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import List, Optional

FASTQ_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
FASTA_EXTS = (".fasta", ".fa", ".fas", ".fna", ".fa.gz", ".fasta.gz")


@dataclass
class ValidationResult:
    """Result of input validation."""
    valid: bool
    message: Optional[str] = None
    suggestion: Optional[str] = None


def _ext_from_uri(uri: str) -> str:
    """Get lowercase extension from URI (handles .gz)."""
    if not uri or not isinstance(uri, str):
        return ""
    name = uri.split("/")[-1].split("?")[0].lower()
    if name.endswith(".fastq.gz") or name.endswith(".fq.gz"):
        return ".fastq.gz" if "fastq" in name else ".fq.gz"
    if name.endswith(".fasta.gz") or name.endswith(".fa.gz"):
        return ".fasta.gz" if "fasta" in name else ".fa.gz"
    if "." in name:
        return "." + name.rsplit(".", 1)[-1]
    return ""


def _is_fastq(uri: str) -> bool:
    if not uri:
        return False
    ext = _ext_from_uri(uri)
    return any(uri.lower().endswith(e) or ext == e for e in FASTQ_EXTS)


def _is_fasta(uri: str) -> bool:
    if not uri:
        return False
    ext = _ext_from_uri(uri)
    return any(uri.lower().endswith(e) or ext == e for e in FASTA_EXTS)


def validate_fastqc_inputs(
    input_r1: str,
    input_r2: str,
    command_hint: Optional[str] = None,
) -> ValidationResult:
    """
    Validate FastQC inputs: must be paired FASTQ, not FASTA.

    Returns ValidationResult with clear message when invalid.
    """
    r1 = (input_r1 or "").strip()
    r2 = (input_r2 or "").strip()

    if not r1 and not r2:
        return ValidationResult(
            valid=False,
            message=(
                "FastQC requires **two** input files: forward (R1) and reverse (R2) reads."
            ),
            suggestion=(
                "Provide S3 URIs for both files, e.g.:\n"
                "`input_r1: s3://bucket/sample_R1.fastq.gz`\n"
                "`input_r2: s3://bucket/sample_R2.fastq.gz`\n\n"
                "If you ran a pipeline in this session, you can say "
                "\"Run FastQC on the raw reads\" and Helix will use the R1/R2 from the pipeline."
            ),
        )

    if not r1 or not r2:
        provided = r1 or r2
        return ValidationResult(
            valid=False,
            message=(
                f"FastQC requires **both** R1 and R2. You provided: `{provided}`"
            ),
            suggestion=(
                "FastQC analyzes paired-end reads. Provide both:\n"
                "- **input_r1**: forward/read 1 FASTQ\n"
                "- **input_r2**: reverse/read 2 FASTQ"
            ),
        )

    if r1 == r2:
        return ValidationResult(
            valid=False,
            message="R1 and R2 must be **different** files. You provided the same file for both.",
            suggestion="Use the forward reads for R1 and the reverse reads for R2.",
        )

    # Format validation
    if _is_fasta(r1) or _is_fasta(r2):
        which = []
        if _is_fasta(r1):
            which.append(f"R1 ({r1})")
        if _is_fasta(r2):
            which.append(f"R2 ({r2})")
        return ValidationResult(
            valid=False,
            message=(
                f"FastQC requires **FASTQ** files, not FASTA. "
                f"The following appear to be FASTA: {', '.join(which)}"
            ),
            suggestion=(
                "FASTA files do not contain per-base quality scores. FastQC needs FASTQ "
                "(e.g. .fastq, .fq, .fastq.gz). Use the raw sequencing reads (R1/R2) before "
                "any conversion to FASTA."
            ),
        )

    if not _is_fastq(r1) or not _is_fastq(r2):
        return ValidationResult(
            valid=False,
            message=(
                "FastQC expects FASTQ format (e.g. .fastq, .fq, .fastq.gz). "
                f"R1: `{r1}` | R2: `{r2}`"
            ),
            suggestion="Ensure both files are FASTQ. If you have merged FASTA, FastQC cannot run on it.",
        )

    return ValidationResult(valid=True)
