"""
test_session_schema_version.py — Slice E: session schema version contract.

Verifies that:
 - Every new session carries a `schema_version` field equal to
   CURRENT_SESSION_SCHEMA_VERSION.
 - Loading an existing session whose stored JSON lacks `schema_version`
   (i.e. an "old" session) does not raise an exception — the field is
   simply absent, and callers should treat its absence as version 0.
 - CURRENT_SESSION_SCHEMA_VERSION is accessible and is an int >= 1.
"""
from __future__ import annotations

import json
import uuid

import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_history_manager(tmp_path):
    """Return a HistoryManager backed by *tmp_path* (no S3, no network)."""
    from backend.history_manager import HistoryManager

    hm = HistoryManager.__new__(HistoryManager)
    hm.sessions = {}
    hm.storage_dir = tmp_path
    hm.s3_client = None
    hm.s3_bucket_name = None
    hm._session_locks = {}
    hm._session_locks_lock = __import__("threading").Lock()
    return hm


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestCurrentSchemaVersion:
    def test_constant_exists_and_is_int(self):
        from backend.history_manager import CURRENT_SESSION_SCHEMA_VERSION
        assert isinstance(CURRENT_SESSION_SCHEMA_VERSION, int)
        assert CURRENT_SESSION_SCHEMA_VERSION >= 1, (
            "Schema version must start at 1 or higher"
        )


class TestNewSessionHasSchemaVersion:
    def test_schema_version_present_in_session_data(self, tmp_path, monkeypatch):
        from backend.history_manager import HistoryManager, CURRENT_SESSION_SCHEMA_VERSION

        monkeypatch.setenv("HELIX_LOCAL_SESSIONS_ONLY", "1")
        hm = HistoryManager(storage_dir=str(tmp_path))
        session_id = hm.create_session()
        session = hm.get_session(session_id)

        assert "schema_version" in session, (
            "New sessions must include a 'schema_version' field"
        )
        assert session["schema_version"] == CURRENT_SESSION_SCHEMA_VERSION, (
            f"Expected schema_version={CURRENT_SESSION_SCHEMA_VERSION}, "
            f"got {session['schema_version']}"
        )

    def test_schema_version_persisted_to_disk(self, tmp_path, monkeypatch):
        """Verify schema_version survives the JSON serialization round-trip."""
        from backend.history_manager import HistoryManager, CURRENT_SESSION_SCHEMA_VERSION

        monkeypatch.setenv("HELIX_LOCAL_SESSIONS_ONLY", "1")
        hm = HistoryManager(storage_dir=str(tmp_path))
        session_id = hm.create_session()

        session_file = tmp_path / f"{session_id}.json"
        assert session_file.exists(), "Session file must be created on disk"

        with session_file.open() as f:
            on_disk = json.load(f)

        assert on_disk.get("schema_version") == CURRENT_SESSION_SCHEMA_VERSION, (
            "schema_version must be persisted to session.json"
        )


class TestOldSessionCompatibility:
    """Ensure sessions created before schema_version was introduced load cleanly."""

    def test_loading_session_without_schema_version_does_not_raise(self, tmp_path, monkeypatch):
        """Write a session.json without schema_version and confirm get_session works."""
        from backend.history_manager import HistoryManager

        monkeypatch.setenv("HELIX_LOCAL_SESSIONS_ONLY", "1")
        hm = HistoryManager(storage_dir=str(tmp_path))

        # Manually write a legacy session file (no schema_version field).
        # HistoryManager persists sessions as {storage_dir}/{session_id}.json.
        session_id = str(uuid.uuid4())
        legacy_data = {
            "session_id": session_id,
            "history": [],
            "runs": [],
            "artifacts": {},
            "metadata": {},
        }
        (tmp_path / f"{session_id}.json").write_text(json.dumps(legacy_data))

        # Loading should not raise
        session = hm.get_session(session_id)
        assert session is not None
        # schema_version absent — callers should treat this as version 0
        assert session.get("schema_version", 0) == 0, (
            "Absence of schema_version should be treated as version 0 by callers"
        )
