"""Bundle ZIP contents for tabular_analysis (inputs, script, plots, tables)."""

from __future__ import annotations

import io
import zipfile
from pathlib import Path
from typing import Any, Dict, List

import pytest

from backend import bundle_generator
from backend.history_manager import history_manager


def _zip_names(buf: io.BytesIO) -> set[str]:
    buf.seek(0)
    with zipfile.ZipFile(buf, "r") as zf:
        return set(zf.namelist())


@pytest.fixture()
def tabular_session_on_disk(tmp_path: Path) -> tuple[str, str, List[Dict[str, Any]]]:
    session_id = "sess-bundle-tabular-01"
    run_id = "run_testbundle1"
    session_dir = tmp_path / session_id
    raw = session_dir / "uploads" / "raw"
    raw.mkdir(parents=True)
    (raw / "sheet.xlsx").write_bytes(b"FAKEXLSX")

    run_root = session_dir / "runs" / run_id
    plots = run_root / "plots"
    tables = run_root / "tables"
    plots.mkdir(parents=True)
    tables.mkdir(parents=True)
    script_path = run_root / "analysis.py"
    script_path.write_text("# sandbox tabular script\nresult = {}\n", encoding="utf-8")
    plot_path = plots / "analysis_sess-bu.png"
    plot_path.write_bytes(b"\x89PNG\r\n\x1a\n\x00")
    table_path = tables / "tabular_result.json"
    table_path.write_text('{"k": 1}', encoding="utf-8")

    runs: List[Dict[str, Any]] = [
        {
            "run_id": run_id,
            "tool": "tabular_analysis",
            "command": "Approve.",
            "tool_args": {"plan_title": "Demo", "plan_goal": "Test"},
            "produced_artifacts": [
                {
                    "type": "script",
                    "title": "analysis.py",
                    "uri": str(script_path.resolve()),
                    "format": "python",
                },
                {
                    "type": "plot",
                    "title": plot_path.name,
                    "uri": str(plot_path.resolve()),
                    "format": "png",
                },
                {
                    "type": "table",
                    "title": "tabular_result.json",
                    "uri": str(table_path.resolve()),
                    "format": "json",
                },
            ],
        }
    ]
    return session_id, run_id, runs


def test_tabular_bundle_includes_inputs_script_plots_tables(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    tabular_session_on_disk: tuple[str, str, List[Dict[str, Any]]],
) -> None:
    session_id, run_id, runs = tabular_session_on_disk

    monkeypatch.setattr(history_manager, "storage_dir", tmp_path)

    def _list_runs(sid: str) -> List[Dict[str, Any]]:
        return runs if sid == session_id else []

    monkeypatch.setattr(history_manager, "list_runs", _list_runs)

    buf, filename = bundle_generator.build_bundle(session_id, run_id)
    assert filename.startswith("tabular_analysis_bundle_")

    names = _zip_names(buf)
    prefix = f"tabular_analysis_bundle_{run_id[:8]}"
    assert f"{prefix}/analysis.py" in names
    assert f"{prefix}/plots/{runs[0]['produced_artifacts'][1]['title']}" in names
    assert f"{prefix}/tables/tabular_result.json" in names
    assert f"{prefix}/inputs/raw/sheet.xlsx" in names
    assert f"{prefix}/run_manifest.json" in names
    assert f"{prefix}/README.md" in names
    # inputs/raw comes from history_manager.storage_dir, not repo ./sessions
    assert history_manager.storage_dir == tmp_path
