"""
test_demo_followup_flows.py

End-to-end smoke tests for each demo's "Load & run" followUpPrompt.

These are the exact strings submitted when a user clicks the
"▶ Load & run" button in the demo panel, tested against the live backend.

Run against localhost (default):
    pytest tests/integration/test_demo_followup_flows.py -v -m integration

Run against the beta server:
    BACKEND_URL=https://helix-beta.noricum-biosoft.com \
      pytest tests/integration/test_demo_followup_flows.py -v -m integration

A test PASSES when:
  - HTTP response is 200
  - Response is not a bare error string (no "Error:" prefix, no 4xx/5xx bubbled up)
  - The response contains meaningful scientific content (not an empty body)
"""
from __future__ import annotations

import json
import os
import time
from typing import Any, Dict, Optional

import pytest
import requests

BACKEND_URL = os.getenv("BACKEND_URL", "http://localhost:8001").rstrip("/")

# ---------------------------------------------------------------------------
# Demo followUpPrompts — kept in sync with frontend/src/data/demoScenarios.ts
# ---------------------------------------------------------------------------

DEMOS = [
    {
        "id": "bulk-rnaseq-factorial",
        "tool": "bulk_rnaseq_analysis",
        "followUpPrompt": (
            "Run bulk RNA-seq differential expression analysis with these inputs:\n"
            "count_matrix: s3://noricum-ngs-data/demo/rnaseq/tgondii_counts.csv\n"
            "sample_metadata: s3://noricum-ngs-data/demo/rnaseq/tgondii_metadata.csv\n"
            "design_formula: ~infection_status + time_point + infection_status:time_point"
        ),
        "expect_no_error": True,
        "expect_content_keywords": ["rna", "differential", "expression", "deseq", "gene"],
    },
    {
        "id": "scrna-sle-pbmc",
        "tool": "single_cell_analysis",
        "followUpPrompt": (
            "data_file: s3://noricum-ngs-data/demo/scrna/sle_pbmc_filtered_feature_bc_matrix.h5\n"
            "data_format: 10x\n"
            "resolution: 0.5\n"
            "steps: all"
        ),
        "expect_no_error": True,
        "expect_content_keywords": ["cell", "cluster", "umap", "single"],
    },
    {
        "id": "phylogenetics-sarscov2",
        "tool": "phylogenetic_tree",
        "followUpPrompt": (
            "Build a phylogenetic tree and pairwise amino acid identity matrix for the "
            "SARS-CoV-2 spike protein across 8 variants of concern using these pre-aligned sequences:\n"
            "sequences: s3://noricum-ngs-data/demo/phylo/sars_cov2_spike.fasta\n"
            "variants: Wuhan-Hu-1, Alpha B.1.1.7, Beta B.1.351, Gamma P.1, "
            "Delta B.1.617.2, Omicron BA.1, Omicron BA.4/5, XBB.1.5"
        ),
        "expect_no_error": True,
        "expect_content_keywords": ["tree", "phylogen", "spike", "variant", "identity"],
    },
]

# Known error strings that should NEVER appear in a successful demo response.
_ERROR_PATTERNS = [
    "Error: No sequences provided",
    "Error: At least 2 sequences",
    "at least 2 sequences are required",
    "No sequences provided",
    "HTTP 404",
    "HTTP 500",
    "internal server error",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _create_session() -> str:
    r = requests.post(f"{BACKEND_URL}/create_session", timeout=15)
    r.raise_for_status()
    return r.json()["session_id"]


def _execute_stream(command: str, session_id: str, timeout: int = 90) -> Dict[str, Any]:
    """POST to /execute/stream and return the final 'result' event data."""
    resp = requests.post(
        f"{BACKEND_URL}/execute/stream",
        json={"command": command, "session_id": session_id},
        stream=True,
        timeout=timeout,
    )
    resp.raise_for_status()

    result: Optional[Dict[str, Any]] = None
    buffer = ""
    for chunk in resp.iter_content(chunk_size=None):
        buffer += chunk.decode("utf-8", errors="replace")
        while "\n\n" in buffer:
            event_text, buffer = buffer.split("\n\n", 1)
            for line in event_text.splitlines():
                if not line.startswith("data: "):
                    continue
                try:
                    event = json.loads(line[6:])
                except json.JSONDecodeError:
                    continue
                if event.get("type") == "result":
                    result = event.get("data", {})
                    return result
                if event.get("type") == "error":
                    pytest.fail(
                        f"SSE error event received: {event.get('detail', event)}"
                    )
    return result or {}


def _flatten_text(result: dict) -> str:
    """Pull all text content out of a result dict for assertion scanning."""
    parts = []
    for key in ("text", "output", "result", "advisory", "message", "explanation"):
        val = result.get(key)
        if isinstance(val, str):
            parts.append(val)
        elif isinstance(val, dict):
            parts.append(json.dumps(val))
    # Also include nested result.result
    inner = result.get("result") or result.get("data") or {}
    if isinstance(inner, dict):
        for key in ("text", "output", "advisory"):
            if isinstance(inner.get(key), str):
                parts.append(inner[key])
    return " ".join(parts).lower()


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def session():
    return _create_session()


# ---------------------------------------------------------------------------
# Parametrised demo tests
# ---------------------------------------------------------------------------

@pytest.mark.integration
@pytest.mark.demo
@pytest.mark.parametrize("demo", DEMOS, ids=[d["id"] for d in DEMOS])
def test_demo_followup_no_error(demo, session):
    """
    Each demo's followUpPrompt ('Load & run') must not return a known error
    string and must produce meaningful content.
    """
    result = _execute_stream(demo["followUpPrompt"], session)

    assert result, (
        f"Demo '{demo['id']}': got an empty result from /execute/stream"
    )

    text = _flatten_text(result)

    # Hard error patterns must not appear anywhere in the response
    if demo["expect_no_error"]:
        for pattern in _ERROR_PATTERNS:
            assert pattern.lower() not in text, (
                f"Demo '{demo['id']}': response contained forbidden error string "
                f"'{pattern}'.\nFull text (truncated):\n{text[:800]}"
            )

    # At least one expected scientific keyword must appear
    keywords = demo.get("expect_content_keywords", [])
    if keywords:
        found = any(kw in text for kw in keywords)
        assert found, (
            f"Demo '{demo['id']}': none of the expected keywords "
            f"{keywords} found in response.\nFull text (truncated):\n{text[:800]}"
        )


@pytest.mark.integration
@pytest.mark.demo
@pytest.mark.parametrize("demo", DEMOS, ids=[d["id"] for d in DEMOS])
def test_demo_followup_routes_to_expected_tool(demo, session):
    """The result must reference the expected tool name (not a fallback tool)."""
    result = _execute_stream(demo["followUpPrompt"], session)

    # Tool can appear at top level or nested
    tool = (
        result.get("tool")
        or (result.get("result") or {}).get("tool")
        or (result.get("data") or {}).get("tool")
        or ""
    )

    expected = demo["tool"]
    # Allow partial match (e.g. "single_cell_analysis" matches "single_cell")
    assert expected in tool or expected.split("_")[0] in tool or expected.split("_")[0] in _flatten_text(result), (
        f"Demo '{demo['id']}': expected tool '{expected}' but got '{tool}'.\n"
        f"Full result keys: {list(result.keys())}"
    )
