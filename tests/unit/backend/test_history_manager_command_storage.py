import json
from pathlib import Path

from backend.history_manager import HistoryManager, sanitize_command_for_storage


def test_sanitize_command_for_storage_truncates_large_payload(monkeypatch):
    monkeypatch.setenv("HELIX_MAX_STORED_COMMAND_CHARS", "120")
    command = "x" * 2000

    out = sanitize_command_for_storage(command)

    assert out["meta"]["truncated"] is True
    assert out["meta"]["original_chars"] == 2000
    assert out["meta"]["stored_chars"] > 120
    assert "[helix-storage-note]" in out["command"]


def test_add_history_entry_stores_trimmed_command_and_metadata(tmp_path, monkeypatch):
    monkeypatch.setenv("HELIX_MAX_STORED_COMMAND_CHARS", "150")
    hm = HistoryManager(storage_dir=str(tmp_path / "sessions"))
    session_id = hm.create_session()

    # Simulate accidental binary-ish attachment text embedded into a command.
    binary_blob = "".join(chr(i % 32) for i in range(3000))
    command = f"analyze this file payload\n{binary_blob}\nend"

    hm.add_history_entry(
        session_id=session_id,
        command=command,
        tool="handle_natural_command",
        result={"status": "success"},
    )

    session = hm.get_session(session_id)
    assert session is not None
    assert len(session["history"]) == 1

    entry = session["history"][0]
    assert len(entry["command"]) < len(command)
    assert "[helix-storage-note]" in entry["command"]
    assert entry["metadata"]["command_storage"]["truncated"] is True
    assert entry["metadata"]["command_storage"]["probable_binary"] is True

    run = session["runs"][0]
    assert run["metadata"]["command_storage"]["truncated"] is True

    # Ensure persisted session JSON also contains the compact command form.
    persisted = json.loads(Path(tmp_path / "sessions" / f"{session_id}.json").read_text())
    assert "[helix-storage-note]" in persisted["history"][0]["command"]
