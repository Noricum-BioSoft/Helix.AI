"""
Known unsupported tools with clear messages and alternatives.

When a user requests a tool that Helix does not support, return a helpful
response instead of failing with a generic error.
"""

from __future__ import annotations

from typing import Dict, Optional

# Tool name (lowercase) -> {reason, alternatives, docs}
UNSUPPORTED_TOOLS: Dict[str, Dict[str, str]] = {
    "multiqc": {
        "reason": (
            "MultiQC aggregates QC reports from multiple tools (FastQC, STAR, etc.) "
            "into a single HTML report. Helix does not yet implement MultiQC."
        ),
        "alternatives": (
            "1. **Run FastQC** on your raw paired-end reads (R1 and R2) — Helix supports this. "
            "Each sample produces an HTML report you can view.\n"
            "2. **View individual FastQC reports** — After running FastQC, use the download links "
            "to get each sample's report.\n"
            "3. **Request MultiQC** — If you need aggregated reports, we can add MultiQC as a "
            "future feature."
        ),
        "docs": "https://multiqc.info/",
    },
    "kraken": {
        "reason": "Kraken is a taxonomic classification tool. Helix does not yet implement Kraken.",
        "alternatives": (
            "For 16S amplicon analysis, Helix supports: read trimming, merging, and quality reports. "
            "Taxonomic classification (Kraken, QIIME2) may be added in a future release."
        ),
        "docs": "",
    },
    "qiime": {
        "reason": "QIIME2 is a microbiome analysis platform. Helix does not yet implement QIIME2.",
        "alternatives": (
            "Helix supports 16S amplicon preprocessing (trim, merge, QC). "
            "For full QIIME2-style diversity analysis, use QIIME2 locally or in a separate pipeline."
        ),
        "docs": "",
    },
}


def get_unsupported_response(tool_name: str) -> Optional[Dict[str, str]]:
    """
    If the requested tool is known to be unsupported, return a dict with
    keys: reason, alternatives, (optional) docs.
    Otherwise return None.
    """
    if not tool_name:
        return None
    key = tool_name.lower().strip().replace(" ", "_")
    return UNSUPPORTED_TOOLS.get(key)
