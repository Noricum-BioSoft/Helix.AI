"""
Phase 3A — GET /capabilities catalog endpoint.

The endpoint is a read-only catalog of every descriptor, annotated with whether
the active execution profile enables it. ``?enabled=true`` narrows to the
providers this deployment turns on; ``?category=`` narrows by category.
"""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient


@pytest.fixture()
def client(monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")
    from backend.main import app

    return TestClient(app)


def test_lists_full_catalog_with_enablement(client):
    body = client.get("/capabilities").json()
    assert body["success"] and body["profile"] == "local-only"
    ids = {c["capability_id"] for c in body["capabilities"]}
    assert "local_compute:bulk_rnaseq_analysis" in ids
    assert "mock_experimental:protein_expression" in ids
    assert "ncbi_entrez:fetch_sequence" in ids  # catalog shows it even when disabled
    by_id = {c["capability_id"]: c for c in body["capabilities"]}
    assert by_id["local_compute:bulk_rnaseq_analysis"]["enabled"] is True
    assert by_id["ncbi_entrez:fetch_sequence"]["enabled"] is False  # not in local-only profile


def test_enabled_filter_reflects_profile(client):
    enabled_ids = {c["capability_id"] for c in client.get("/capabilities?enabled=true").json()["capabilities"]}
    assert "local_compute:bulk_rnaseq_analysis" in enabled_ids
    assert "ncbi_entrez:fetch_sequence" not in enabled_ids


def test_category_filter(client):
    body = client.get("/capabilities?category=experimental").json()
    assert body["capabilities"] and all(c["category"] == "experimental" for c in body["capabilities"])


def test_single_capability_and_404(client):
    ok = client.get("/capabilities/local_compute:bulk_rnaseq_analysis")
    assert ok.status_code == 200 and ok.json()["capability"]["provider"] == "local_compute"
    assert client.get("/capabilities/does_not:exist").status_code == 404
