"""
Pytest fixtures for the USER_SCENARIOS testbed.

Uses TestClient against the FastAPI app so tests run without a live server.
Default: HELIX_MOCK_MODE=1 for deterministic, offline-friendly runs.
"""

import os
from pathlib import Path

import pytest
from fastapi.testclient import TestClient


# Ensure project root and tools are on path (same as top-level conftest)
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in __import__("sys").path:
    __import__("sys").path.insert(0, str(PROJECT_ROOT))
TOOLS_PATH = PROJECT_ROOT / "tools"
if str(TOOLS_PATH) not in __import__("sys").path:
    __import__("sys").path.append(str(TOOLS_PATH))


@pytest.fixture
def client(tmp_path, monkeypatch):
    """TestClient for main FastAPI app; isolated session storage in tmp_path."""
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")
    monkeypatch.setenv("HELIX_DEMO_MODE", "0")
    monkeypatch.setenv("HELIX_SANDBOX_HOST_FALLBACK", "1")

    from backend.history_manager import history_manager
    session_dir = tmp_path / "sessions"
    session_dir.mkdir(parents=True, exist_ok=True)
    history_manager.storage_dir = session_dir
    history_manager.sessions = {}
    history_manager._sessions_loaded = True

    from backend.main import app
    return TestClient(app)


@pytest.fixture
def session_id(client):
    """Create a session and return its session_id."""
    r = client.post("/create_session", json={})
    assert r.status_code == 200, r.text
    data = r.json()
    sid = data.get("session_id")
    assert sid, data
    return sid


# Demo payloads from USER_SCENARIOS / demoScenarios (abbreviated for tests)
BULK_RNASEQ_PROMPT = (
    "I have bulk RNA-seq count data and sample metadata. "
    "Run differential expression with design ~infection_status + time_point. "
    "Use count_matrix: s3://noricum-ngs-data/demo/rnaseq/tgondii_counts.csv "
    "and sample_metadata: s3://noricum-ngs-data/demo/rnaseq/tgondii_metadata.csv "
    "and design_formula: ~infection_status + time_point + infection_status:time_point"
)

BULK_RNASEQ_NEEDS_INPUTS_PROMPT = (
    "I have bulk RNA-seq count data and sample metadata. What can you do with it?"
)

SINGLECELL_NEEDS_INPUTS_PROMPT = "What do I need to run single-cell RNA-seq analysis?"

FASTQC_NEEDS_INPUTS_PROMPT = "How do I run FastQC on my FASTQ files?"

QA_PROMPT = "What's the difference between bulk and single-cell RNA-seq?"

FASTA_THREE = ">s1\nATGCGATCG\n>s2\nATGCGATC\n>s3\nATGCGATCGATC"
