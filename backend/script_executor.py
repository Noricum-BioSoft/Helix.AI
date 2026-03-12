"""
script_executor.py — Execute a generated analysis script and collect artifacts.

Responsibilities
----------------
1. Run a saved ``analysis.py`` in a subprocess (isolated from the server process).
2. Parse the JSON summary printed to stdout.
3. Discover all files written under ``run_dir/`` and return them as artifact records.
4. Compute a human-readable diff between two run directories (for iterative edits).
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import logging

logger = logging.getLogger(__name__)

# Maximum time (seconds) a script may run before it is killed.
DEFAULT_TIMEOUT = 300


# ── execution ─────────────────────────────────────────────────────────────────

def execute(
    script_path: str | Path,
    timeout: int = DEFAULT_TIMEOUT,
    env: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """Run *script_path* and return a structured result dict.

    The script is expected to print a single JSON object to stdout as its last
    line.  Any preceding output is captured in ``logs``.

    Returns
    -------
    dict with keys:
        status      "success" | "error"
        summary     parsed JSON object from stdout (or {})
        artifacts   list of artifact dicts (type, uri, name, size_bytes)
        logs        combined stdout/stderr text
        elapsed_s   wall-clock seconds
        error       error message (only present on failure)
    """
    script_path = Path(script_path)
    if not script_path.exists():
        return {"status": "error", "error": f"Script not found: {script_path}"}

    run_env = {**os.environ, **(env or {})}

    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            [sys.executable, str(script_path)],
            capture_output=True,
            text=True,
            timeout=timeout,
            env=run_env,
        )
    except subprocess.TimeoutExpired:
        return {"status": "error", "error": f"Script timed out after {timeout}s"}
    except Exception as exc:
        return {"status": "error", "error": str(exc)}

    elapsed = round(time.perf_counter() - t0, 2)
    combined_logs = (proc.stdout or "") + (proc.stderr or "")

    if proc.returncode != 0:
        return {
            "status":   "error",
            "error":    f"Script exited with code {proc.returncode}",
            "logs":     combined_logs[-4000:],
            "elapsed_s": elapsed,
        }

    # Parse the last non-empty line as JSON summary
    summary: Dict[str, Any] = {}
    for line in reversed((proc.stdout or "").splitlines()):
        line = line.strip()
        if line.startswith("{"):
            try:
                summary = json.loads(line)
                break
            except json.JSONDecodeError:
                continue

    run_dir_str = summary.get("run_dir") or str(script_path.parent)
    artifacts = _collect_artifacts(Path(run_dir_str))

    return {
        "status":    summary.get("status", "success"),
        "summary":   summary,
        "artifacts": artifacts,
        "logs":      combined_logs[-2000:],
        "elapsed_s": elapsed,
    }


# ── artifact collection ───────────────────────────────────────────────────────

_ARTIFACT_TYPES: Dict[str, str] = {
    ".png":  "plot",
    ".svg":  "plot",
    ".json": "table",
    ".csv":  "table",
    ".tsv":  "table",
    ".nwk":  "newick",
    ".txt":  "text",
    ".py":   "script",
}


def _collect_artifacts(run_dir: Path) -> List[Dict[str, Any]]:
    if not run_dir.exists():
        return []
    artifacts = []
    for f in sorted(run_dir.rglob("*")):
        if not f.is_file():
            continue
        suffix = f.suffix.lower()
        atype = _ARTIFACT_TYPES.get(suffix, "file")
        artifacts.append({
            "type":       atype,
            "uri":        str(f),
            "name":       f.stem,
            "size_bytes": f.stat().st_size,
            "relative":   str(f.relative_to(run_dir)),
        })
    return artifacts


# ── script patching ───────────────────────────────────────────────────────────

def patch_parameters(script_text: str, patch: Dict[str, Any]) -> str:
    """Apply *patch* to the ``# ── Parameters ──`` block of *script_text*.

    Handles simple ``KEY = value`` assignments.  Each value is rendered via
    ``repr()`` so strings stay quoted, numbers stay bare, and None becomes None.

    Parameters
    ----------
    script_text : str
        Full source of the analysis script.
    patch : dict
        Mapping of parameter name → new value, e.g. ``{"ALPHA": 0.01}``.

    Returns
    -------
    str — modified script text.
    """
    lines = script_text.splitlines()
    new_lines: List[str] = []
    in_params = False

    for line in lines:
        stripped = line.strip()

        # Track entry/exit of the Parameters block
        if "# ── Parameters" in line or "# ── parameters" in line:
            in_params = True
        elif in_params and stripped.startswith("# ──") and "Parameters" not in line:
            in_params = False

        if in_params and "=" in stripped and not stripped.startswith("#"):
            # Extract the variable name (left of the first "=")
            var_name = stripped.split("=", 1)[0].strip()
            if var_name in patch:
                new_val = repr(patch[var_name])
                # Preserve the original indentation
                indent = len(line) - len(line.lstrip())
                comment = ""
                if "#" in stripped.split("=", 1)[-1]:
                    comment = "  " + stripped.split("=", 1)[-1].split("#", 1)[-1].strip()
                    comment = f"  # {comment}"
                new_line = " " * indent + f"{var_name} = {new_val}{comment}"
                new_lines.append(new_line)
                continue

        new_lines.append(line)

    return "\n".join(new_lines)


# ── diff helpers ──────────────────────────────────────────────────────────────

def parameter_diff(old_script: str, new_script: str) -> List[str]:
    """Return a list of human-readable change strings between two scripts.

    Only compares lines inside the ``# ── Parameters ──`` block.

    Example output::

        ["ALPHA: 0.05 → 0.01", "X_SCALE: 'log2' → 'linear'"]
    """
    def _extract_params(src: str) -> Dict[str, str]:
        params: Dict[str, str] = {}
        in_block = False
        for line in src.splitlines():
            if "# ── Parameters" in line or "# ── parameters" in line:
                in_block = True
                continue
            if in_block and line.strip().startswith("# ──") and "Parameters" not in line:
                break
            if in_block and "=" in line and not line.strip().startswith("#"):
                parts = line.split("=", 1)
                key = parts[0].strip()
                val = parts[1].split("#")[0].strip()
                params[key] = val
        return params

    old_p = _extract_params(old_script)
    new_p = _extract_params(new_script)
    changes = []
    for k in new_p:
        if old_p.get(k) != new_p[k]:
            changes.append(f"{k}: {old_p.get(k, '?')} → {new_p[k]}")
    return changes


def output_diff(
    old_run_dir: str | Path,
    new_run_dir: str | Path,
) -> Dict[str, Any]:
    """Compare the artifacts in two run directories.

    Returns a dict::

        {
            "new_files":     ["plots/volcano_....png"],
            "removed_files": [],
            "changed_files": ["tables/de_....json"],
        }
    """
    def _rel_set(d: Path) -> Dict[str, int]:
        if not d.exists():
            return {}
        return {
            str(f.relative_to(d)): f.stat().st_size
            for f in d.rglob("*") if f.is_file()
        }

    old_files = _rel_set(Path(old_run_dir))
    new_files = _rel_set(Path(new_run_dir))

    added   = [f for f in new_files if f not in old_files]
    removed = [f for f in old_files if f not in new_files]
    changed = [
        f for f in new_files
        if f in old_files and new_files[f] != old_files[f]
    ]

    return {
        "new_files":     sorted(added),
        "removed_files": sorted(removed),
        "changed_files": sorted(changed),
    }
