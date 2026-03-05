import os
import sys
from pathlib import Path
from fastapi.testclient import TestClient

PROJECT_ROOT = Path(__file__).resolve().parents[2]

# Ensure mock mode to avoid network/tooling dependencies
os.environ["HELIX_MOCK_MODE"] = "1"

# Ensure repo root is importable even when pytest is invoked from a subdir
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from backend.main_with_mcp import app


client = TestClient(app)


def _assert_standard_envelope(payload: dict):
    # Required top-level keys
    for key in [
        "version",
        "success",
        "session_id",
        "prompt",
        "tool",
        "status",
        "text",
        "data",
        "logs",
        "errors",
        "mcp",
        "raw_result",
        "timestamp",
    ]:
        assert key in payload, f"Missing key in response: {key}"

    assert isinstance(payload["data"], dict)
    assert isinstance(payload["logs"], list)
    assert isinstance(payload["errors"], list)
    assert isinstance(payload["mcp"], dict)


def test_execute_returns_standard_envelope():
    cmd = "fetch sequence NC_000001.11"
    resp = client.post("/execute", json={"command": cmd})
    assert resp.status_code == 200
    data = resp.json()

    _assert_standard_envelope(data)
    assert data["success"] is True
    assert data["tool"] == "fetch_ncbi_sequence"
    assert data["prompt"] == cmd
    assert data["status"] != "error"

    # Check that mock data appears in raw_result
    raw = data.get("raw_result", {})
    assert isinstance(raw, dict)
    # With the ExecutionBroker envelope, the original tool output is under raw_result["result"].
    accession = (
        raw.get("accession")
        or raw.get("result", {}).get("accession")
        or raw.get("result", {}).get("result", {}).get("accession")
    )
    status = (
        raw.get("status")
        or raw.get("result", {}).get("status")
        or raw.get("result", {}).get("result", {}).get("status")
    )
    assert accession == "NC_000001.11"
    assert status == "success"


def test_health():
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json().get("status") == "healthy"
