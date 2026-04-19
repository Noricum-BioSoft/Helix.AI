from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List


BLOCKED_CONTENT_TYPES = {
    "application/x-dosexec",
    "application/x-msdownload",
    "application/x-sh",
    "application/x-mach-binary",
}

SUSPICIOUS_TEXT_PATTERNS = (
    "<script",
    "<?php",
    "powershell -",
    "rm -rf",
    "os.system(",
    "subprocess.popen(",
    "#!/bin/bash",
)

RESTRICTED_NAME_MARKERS = (
    "patient",
    "clinical",
    "subject",
    "cohort_id",
    "mrn",
)


@dataclass
class UploadPolicyDecision:
    profile: str
    decision: str
    reason: str
    sensitivity_class: str
    scan_flags: List[str]
    approval_required: bool

    def to_dict(self) -> Dict[str, Any]:
        return {
            "profile": self.profile,
            "decision": self.decision,
            "reason": self.reason,
            "sensitivity_class": self.sensitivity_class,
            "scan_flags": list(self.scan_flags),
            "approval_required": bool(self.approval_required),
        }


def get_policy_profile() -> str:
    return (os.getenv("HELIX_POLICY_PROFILE", "p0") or "p0").strip().lower()


def policy_profile_defaults(profile: str) -> Dict[str, Any]:
    normalized = (profile or "p0").strip().lower()
    if normalized != "p0":
        normalized = "p0"
    return {
        "profile": normalized,
        "upload_intake_enabled": True,
        "block_suspicious_payloads": True,
        "approval_for_restricted_human_data": True,
        "audit_enabled": True,
    }


def evaluate_upload_intake(
    *,
    filename: str,
    content_type: str,
    first_chunk: bytes,
) -> UploadPolicyDecision:
    profile = get_policy_profile()
    defaults = policy_profile_defaults(profile)

    scan_flags: List[str] = []
    lower_name = (filename or "").lower()
    lower_content_type = (content_type or "").lower()
    suffix = Path(lower_name).suffix

    sensitivity_class = "internal"
    if any(marker in lower_name for marker in RESTRICTED_NAME_MARKERS):
        sensitivity_class = "restricted_human_data"
        scan_flags.append("restricted_name_marker")

    if lower_content_type in BLOCKED_CONTENT_TYPES:
        scan_flags.append("blocked_content_type")

    # Binary bioinformatics formats legitimately contain null bytes in their
    # magic bytes / headers — exclude them from the binary null-byte scan.
    _BINARY_BIO_SUFFIXES = {
        ".xlsx", ".xls",   # ZIP-based Office
        ".bam", ".cram",   # compressed alignment
        ".bcf",            # compressed VCF
        ".h5ad", ".loom", ".h5",  # HDF5-based
    }
    _is_known_binary = suffix in _BINARY_BIO_SUFFIXES
    if not _is_known_binary and first_chunk and b"\x00" in first_chunk[:4096]:
        scan_flags.append("binary_null_bytes")

    preview = (first_chunk or b"")[:4096].decode("utf-8", errors="ignore").lower()
    if preview:
        for pattern in SUSPICIOUS_TEXT_PATTERNS:
            if pattern in preview:
                scan_flags.append("suspicious_payload_pattern")
                break

    if defaults["block_suspicious_payloads"] and any(
        flag in scan_flags for flag in ("blocked_content_type", "suspicious_payload_pattern")
    ):
        return UploadPolicyDecision(
            profile=profile,
            decision="block",
            reason="Upload blocked by policy scan.",
            sensitivity_class=sensitivity_class,
            scan_flags=scan_flags,
            approval_required=False,
        )

    approval_required = bool(
        defaults["approval_for_restricted_human_data"]
        and sensitivity_class == "restricted_human_data"
    )

    if approval_required:
        return UploadPolicyDecision(
            profile=profile,
            decision="allow_with_approval",
            reason="Restricted data marker detected; approval required before high-risk execution.",
            sensitivity_class=sensitivity_class,
            scan_flags=scan_flags,
            approval_required=True,
        )

    return UploadPolicyDecision(
        profile=profile,
        decision="allow",
        reason="Upload passed intake policy checks.",
        sensitivity_class=sensitivity_class,
        scan_flags=scan_flags,
        approval_required=False,
    )
