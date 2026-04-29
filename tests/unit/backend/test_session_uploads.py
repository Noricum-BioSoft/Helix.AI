from pathlib import Path

import pytest
from fastapi.testclient import TestClient


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")
    monkeypatch.setenv("HELIX_LOCAL_SESSIONS_ONLY", "1")
    monkeypatch.setenv("HELIX_MAX_UPLOAD_MB", "20")
    # Keyword routing for tests without a live LLM; production uses LLM-first.

    from backend.history_manager import history_manager

    session_dir = tmp_path / "sessions"
    session_dir.mkdir(parents=True, exist_ok=True)
    history_manager.storage_dir = session_dir
    history_manager.sessions = {}
    history_manager._sessions_loaded = True

    from backend.main import app

    return TestClient(app)


def _create_session(client: TestClient) -> str:
    response = client.post("/create_session", json={})
    assert response.status_code == 200, response.text
    session_id = response.json().get("session_id")
    assert session_id
    return session_id


def test_upload_session_file_creates_directories_and_metadata(client: TestClient):
    session_id = _create_session(client)

    files = [
        ("files", ("sample.fastq", b"@r1\nACGT\n+\nIIII\n", "text/plain")),
    ]
    response = client.post(f"/session/{session_id}/uploads", files=files)
    assert response.status_code == 200, response.text
    body = response.json()
    assert body["success"] is True
    assert body["uploaded_count"] == 1

    uploaded = body["files"][0]
    assert uploaded["filename"].endswith(".fastq")
    assert uploaded["size"] > 0
    assert uploaded["local_path"]

    from backend.history_manager import history_manager

    session = history_manager.get_session(session_id)
    assert session is not None
    metadata = session.get("metadata", {})
    upload_paths = metadata.get("upload_paths", {})
    assert upload_paths.get("raw")
    assert upload_paths.get("processed")
    assert upload_paths.get("meta")
    assert Path(uploaded["local_path"]).exists()

    uploaded_files = metadata.get("uploaded_files", [])
    assert len(uploaded_files) == 1
    assert uploaded_files[0]["file_id"] == uploaded["file_id"]

    list_response = client.get(f"/session/{session_id}/files")
    assert list_response.status_code == 200, list_response.text
    list_body = list_response.json()
    assert list_body["uploaded_count"] == 1
    assert list_body["files"][0]["filename"] == uploaded["filename"]


def test_upload_rejects_file_over_limit(client: TestClient, monkeypatch):
    monkeypatch.setenv("HELIX_MAX_UPLOAD_MB", "0.0005")
    session_id = _create_session(client)

    # ~1KB file with 0.0005MB (~524B) server-side limit should be rejected.
    payload = b"A" * 1024
    files = [
        ("files", ("too_large.csv", payload, "text/csv")),
    ]
    response = client.post(f"/session/{session_id}/uploads", files=files)
    assert response.status_code == 413, response.text
    assert "exceeds the allowed upload size" in response.json()["detail"].lower()


def test_upload_rejects_unsupported_extension(client: TestClient):
    session_id = _create_session(client)

    files = [
        ("files", ("malware.exe", b"binary", "application/octet-stream")),
    ]
    response = client.post(f"/session/{session_id}/uploads", files=files)
    assert response.status_code == 415, response.text
    assert "unsupported file type" in response.json()["detail"].lower()


def test_upload_blocks_suspicious_payload(client: TestClient):
    session_id = _create_session(client)
    files = [
        ("files", ("notes.txt", b"<script>alert('xss')</script>", "text/plain")),
    ]
    response = client.post(f"/session/{session_id}/uploads", files=files)
    assert response.status_code == 400, response.text
    assert "upload blocked by intake policy" in response.json()["detail"].lower()


def test_restricted_upload_requires_approval_and_can_be_approved(client: TestClient):
    session_id = _create_session(client)
    files = [
        ("files", ("patient_cohort.csv", b"gene,score\nTP53,0.9\n", "text/csv")),
    ]
    response = client.post(f"/session/{session_id}/uploads", files=files)
    assert response.status_code == 200, response.text
    uploaded = response.json()["files"][0]
    assert uploaded["requires_policy_approval"] is True
    assert uploaded["policy_state"] == "approval_required"
    assert uploaded["intake_policy"]["approval_required"] is True

    execute_response = client.post(
        "/execute",
        # Use a non-approval command so the approval-bypass path is not triggered
        json={"session_id": session_id, "command": "run analysis", "execute_plan": True},
    )
    # The endpoint raises HTTPException(409) when execute_plan=True with pending policy uploads
    assert execute_response.status_code == 409, execute_response.text
    detail = execute_response.json().get("detail", "")
    assert "policy approval" in detail.lower()

    approve_response = client.post(f"/session/{session_id}/uploads/approve")
    assert approve_response.status_code == 200, approve_response.text
    assert approve_response.json()["approved_count"] == 1

    from backend.history_manager import history_manager

    session = history_manager.get_session(session_id)
    uploaded_files = (session or {}).get("metadata", {}).get("uploaded_files", [])
    assert uploaded_files[0]["policy_state"] == "approved_for_execution"
