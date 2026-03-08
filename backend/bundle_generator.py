"""
bundle_generator.py — Create a self-contained reproducibility ZIP for a run.

Bundle layout
-------------
<tool>_bundle_<run_id_short>/
├── README.md               Human-readable walkthrough (prompts, steps, I/O)
├── run_manifest.json       Machine-readable metadata for every run in the chain
├── analysis.py             Re-runnable script for the target run
├── plots/                  PNG plots (skipped if > LARGE_FILE_BYTES)
├── tables/                 JSON / CSV tables (skipped if > LARGE_FILE_BYTES)
├── iteration_history/      One sub-folder per ancestor run in the chain
│   └── <run_id_short>/
│       ├── analysis.py
│       └── params.json
└── large_files.txt         Paths / URIs of files that exceeded the size limit

The bundle is assembled in-memory as a BytesIO zip so the API endpoint can
stream it directly without touching the filesystem.
"""

from __future__ import annotations

import io
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

# Files larger than this are referenced in large_files.txt rather than embedded.
LARGE_FILE_BYTES = 5 * 1024 * 1024  # 5 MB

_PROJECT_ROOT = Path(__file__).parent.parent.resolve()
_SESSIONS_ROOT = _PROJECT_ROOT / "sessions"


# ── public entry point ────────────────────────────────────────────────────────

def build_bundle(
    session_id: str,
    run_id: Optional[str] = None,
) -> tuple[io.BytesIO, str]:
    """Build a reproducibility ZIP and return (buffer, suggested_filename).

    Parameters
    ----------
    session_id : str
        Active session ID.
    run_id : str | None
        Target run.  When None, the most-recent scriptable run is used.

    Returns
    -------
    (BytesIO, str)
        In-memory ZIP buffer (position reset to 0) and a suggested filename
        like ``bulk_rnaseq_bundle_a1b2c3d4.zip``.
    """
    from backend.history_manager import history_manager

    all_runs: List[Dict[str, Any]] = history_manager.list_runs(session_id) or []

    # Resolve target run
    target = _find_run(all_runs, run_id)
    if target is None:
        raise ValueError(
            f"No run found for session={session_id!r} run_id={run_id!r}"
        )

    target_run_id = target.get("run_id", run_id or "unknown")
    tool_name = target.get("tool", "analysis")

    # Build the ancestor chain (oldest → newest)
    chain = _build_chain(all_runs, target_run_id)

    short_id = target_run_id[:8]
    dir_name = f"{tool_name}_bundle_{short_id}"
    suggested_filename = f"{dir_name}.zip"

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        large_refs: List[str] = []

        # ── analysis.py for the target run ───────────────────────────────────
        _add_script(zf, target, dir_name, "analysis.py", large_refs)

        # ── plots & tables from the target run ───────────────────────────────
        _add_artifacts(zf, target, dir_name, large_refs, session_id, target_run_id)

        # ── iteration history ─────────────────────────────────────────────────
        if len(chain) > 1:
            for ancestor in chain[:-1]:  # all but the last (= target)
                a_run_id = ancestor.get("run_id", "")
                a_short = a_run_id[:8]
                sub = f"{dir_name}/iteration_history/{a_short}"
                _add_script(zf, ancestor, sub, "analysis.py", large_refs)
                # params.json
                params_bytes = json.dumps(
                    ancestor.get("tool_args") or {}, indent=2
                ).encode()
                zf.writestr(f"{sub}/params.json", params_bytes)

        # ── run_manifest.json ─────────────────────────────────────────────────
        manifest = _build_manifest(session_id, target_run_id, chain)
        zf.writestr(
            f"{dir_name}/run_manifest.json",
            json.dumps(manifest, indent=2).encode(),
        )

        # ── README.md ─────────────────────────────────────────────────────────
        readme = _build_readme(session_id, tool_name, chain, large_refs)
        zf.writestr(f"{dir_name}/README.md", readme.encode())

        # ── large_files.txt ───────────────────────────────────────────────────
        if large_refs:
            zf.writestr(
                f"{dir_name}/large_files.txt",
                "\n".join(large_refs) + "\n",
            )

    buf.seek(0)
    return buf, suggested_filename


# ── chain helpers ─────────────────────────────────────────────────────────────

def _find_run(runs: List[Dict], run_id: Optional[str]) -> Optional[Dict]:
    """Return the run dict matching run_id, or the most-recent scriptable run."""
    if run_id:
        for r in reversed(runs):
            if isinstance(r, dict) and r.get("run_id") == run_id:
                return r
        return None
    # Fall back to most-recent run that has a saved analysis.py
    for r in reversed(runs):
        if not isinstance(r, dict):
            continue
        arts = r.get("produced_artifacts") or []
        if any(
            isinstance(a, dict) and a.get("type") == "script"
            for a in arts
        ):
            return r
    # Last resort: any run
    return runs[-1] if runs else None


def _build_chain(runs: List[Dict], target_run_id: str) -> List[Dict]:
    """Walk parent_run_id links upward and return the list oldest → target."""
    by_id = {r.get("run_id"): r for r in runs if isinstance(r, dict)}

    chain: List[Dict] = []
    current_id: Optional[str] = target_run_id
    seen: set = set()

    while current_id and current_id not in seen:
        node = by_id.get(current_id)
        if node is None:
            break
        chain.append(node)
        seen.add(current_id)
        current_id = node.get("parent_run_id") or node.get("metadata", {}).get("parent_run_id")

    chain.reverse()  # oldest first
    return chain


# ── zip helpers ───────────────────────────────────────────────────────────────

def _add_script(
    zf: zipfile.ZipFile,
    run: Dict,
    dest_dir: str,
    filename: str,
    large_refs: List[str],
) -> None:
    arts = run.get("produced_artifacts") or []
    for a in arts:
        if not isinstance(a, dict):
            continue
        if a.get("type") == "script" and str(a.get("uri", "")).endswith(".py"):
            p = Path(str(a["uri"]))
            if p.exists():
                size = p.stat().st_size
                if size > LARGE_FILE_BYTES:
                    large_refs.append(f"[script] {p}")
                else:
                    zf.writestr(f"{dest_dir}/{filename}", p.read_bytes())
            return

    # Fallback: try to locate analysis.py from sessions/ tree
    run_id = run.get("run_id", "")
    if run_id:
        sid = run.get("session_id", "")
        for candidate in _SESSIONS_ROOT.rglob(f"*{run_id[:8]}*/analysis.py"):
            zf.writestr(f"{dest_dir}/{filename}", candidate.read_bytes())
            return


def _add_artifacts(
    zf: zipfile.ZipFile,
    run: Dict,
    dest_dir: str,
    large_refs: List[str],
    session_id: str,
    run_id: str,
) -> None:
    """Add plots and tables from the run's artifact files."""
    # Check produced_artifacts list
    arts = run.get("produced_artifacts") or []
    for a in arts:
        if not isinstance(a, dict):
            continue
        atype = a.get("type", "")
        if atype not in ("plot", "table"):
            continue
        p = Path(str(a.get("uri", "")))
        if not p.exists():
            continue
        size = p.stat().st_size
        sub = "plots" if atype == "plot" else "tables"
        if size > LARGE_FILE_BYTES:
            large_refs.append(f"[{atype}] {p}")
        else:
            zf.writestr(f"{dest_dir}/{sub}/{p.name}", p.read_bytes())

    # Also sweep run_dir/plots/ and run_dir/tables/ on disk
    run_dir = _SESSIONS_ROOT / session_id / "runs" / run_id
    for sub in ("plots", "tables"):
        sub_dir = run_dir / sub
        if not sub_dir.exists():
            continue
        for f in sorted(sub_dir.iterdir()):
            if not f.is_file():
                continue
            size = f.stat().st_size
            if size > LARGE_FILE_BYTES:
                large_refs.append(f"[{sub[:-1]}] {f}")
            else:
                zf.writestr(f"{dest_dir}/{sub}/{f.name}", f.read_bytes())


# ── manifest ──────────────────────────────────────────────────────────────────

def _build_manifest(
    session_id: str,
    target_run_id: str,
    chain: List[Dict],
) -> Dict[str, Any]:
    runs_meta = []
    for r in chain:
        runs_meta.append({
            "run_id":        r.get("run_id"),
            "parent_run_id": r.get("parent_run_id") or r.get("metadata", {}).get("parent_run_id"),
            "tool":          r.get("tool"),
            "prompt":        r.get("command") or r.get("prompt"),
            "params":        r.get("tool_args") or {},
            "artifacts": [
                {
                    "type":  a.get("type"),
                    "title": a.get("title") or a.get("name"),
                    "uri":   a.get("uri"),
                }
                for a in (r.get("produced_artifacts") or [])
                if isinstance(a, dict)
            ],
            "timestamp": r.get("timestamp"),
        })

    return {
        "schema_version": "1.0",
        "created_at":     datetime.now(timezone.utc).isoformat(),
        "session_id":     session_id,
        "target_run_id":  target_run_id,
        "n_iterations":   len(chain),
        "runs":           runs_meta,
    }


# ── README ────────────────────────────────────────────────────────────────────

_TOOL_TITLES = {
    "bulk_rnaseq_analysis":   "Bulk RNA-seq Differential Expression Analysis",
    "single_cell_analysis":   "Single-Cell RNA-seq Analysis",
    "phylogenetic_tree":      "Phylogenetic Tree Analysis",
    "fastqc_quality_analysis": "FastQC / Amplicon Quality Assessment",
}


def _build_readme(
    session_id: str,
    tool_name: str,
    chain: List[Dict],
    large_refs: List[str],
) -> str:
    title = _TOOL_TITLES.get(tool_name, tool_name.replace("_", " ").title())
    lines: List[str] = []

    lines += [
        f"# {title}",
        "",
        f"**Session:** `{session_id}`  ",
        f"**Generated:** {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}  ",
        f"**Iterations:** {len(chain)}",
        "",
        "---",
        "",
        "## How to reproduce",
        "",
        "1. Make sure the Helix.AI Python environment is active:",
        "   ```bash",
        "   cd /path/to/Helix.AI",
        "   source .venv/bin/activate",
        "   ```",
        "2. Run the analysis script for the final state:",
        "   ```bash",
        "   python analysis.py",
        "   ```",
        "3. To reproduce an earlier iteration, run the script from",
        "   `iteration_history/<run_id_short>/analysis.py`.",
        "",
        "Each script is self-contained — just edit the `# ── Parameters ──` block",
        "at the top to adjust settings, then re-run.",
        "",
        "---",
        "",
        "## Analysis steps",
        "",
    ]

    for i, run in enumerate(chain, start=1):
        run_id = run.get("run_id", "?")
        prompt = run.get("command") or run.get("prompt") or "(no prompt recorded)"
        # Strip internal prefix like "[bulk_rnaseq_analysis] ..."
        if prompt.startswith("[") and "]" in prompt:
            prompt = prompt.split("]", 1)[-1].strip()

        params = run.get("tool_args") or {}
        arts = [
            a for a in (run.get("produced_artifacts") or [])
            if isinstance(a, dict) and a.get("type") in ("plot", "table")
        ]

        lines.append(f"### Step {i} — `{run_id[:12]}{'…' if len(run_id) > 12 else ''}`")
        lines.append("")
        if i == 1:
            lines.append("**Type:** Initial run  ")
        else:
            parent = run.get("parent_run_id") or run.get("metadata", {}).get("parent_run_id") or "?"
            lines.append(f"**Type:** Iterative update (parent: `{str(parent)[:12]}…`)  ")

        if prompt:
            lines.append(f"**Prompt:** {prompt}  ")

        if params:
            lines.append("")
            lines.append("**Parameters:**")
            lines.append("")
            lines.append("| Parameter | Value |")
            lines.append("|-----------|-------|")
            for k, v in params.items():
                lines.append(f"| `{k}` | `{v}` |")

        if arts:
            lines.append("")
            lines.append("**Outputs:**")
            lines.append("")
            for a in arts:
                atype = a.get("type", "file")
                name = a.get("title") or a.get("name") or Path(a.get("uri", "unknown")).name
                uri = a.get("uri", "")
                lines.append(f"- [{atype}] `{name}` → `{uri}`")

        lines.append("")

    if large_refs:
        lines += [
            "---",
            "",
            "## Large file references",
            "",
            "The following files exceeded the 5 MB size limit and were not",
            "embedded in this bundle. Their paths/URIs are listed in `large_files.txt`.",
            "",
        ]
        for ref in large_refs:
            lines.append(f"- `{ref}`")
        lines.append("")

    lines += [
        "---",
        "",
        "*Generated by [Helix.AI](https://github.com/your-org/Helix.AI) — "
        "Iterative bioinformatics workflows.*",
        "",
    ]

    return "\n".join(lines)
