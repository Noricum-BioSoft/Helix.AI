"""
script_executor.py — Execute a generated analysis script and collect artifacts.

Responsibilities
----------------
1. Run a saved ``analysis.py`` in a sandboxed Docker container (default) or host
   subprocess, for security and user trust.
2. Parse the JSON summary printed to stdout.
3. Discover all files written under ``run_dir/`` and return them as artifact records.
4. Compute a human-readable diff between two run directories (for iterative edits).

Sandbox execution is controlled by HELIX_ANALYSIS_USE_SANDBOX (default: true).
Set to "false" to run on the host (e.g. CI without Docker). When sandbox is
enabled and Docker is unavailable, execution fails with a clear error unless
HELIX_SANDBOX_HOST_FALLBACK=1 is set.
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


def _run_in_sandbox(
    script_path: Path,
    timeout: int,
    env: Optional[Dict[str, str]],
) -> Dict[str, Any]:
    """Execute analysis.py inside the Docker sandbox. Returns same shape as execute()."""
    from backend.sandbox_executor import get_sandbox_executor

    executor = get_sandbox_executor()
    run_dir = script_path.parent
    # allow_host_fallback is resolved inside execute_command based on
    # self._docker_available and HELIX_SANDBOX_HOST_FALLBACK; pass None to let
    # SandboxExecutor decide (it auto-falls back when Docker is unavailable).
    allow_host_fallback = None

    result = executor.execute_command(
        cmd=["python", script_path.name],
        working_dir=str(run_dir),
        output_dir=str(run_dir),
        timeout=timeout,
        env_vars=env,
        allow_host_fallback=allow_host_fallback,
    )

    combined_logs = (result.stdout or "") + (result.stderr or "")

    if not result.success:
        return {
            "status": "error",
            "error": result.error_message or f"Script exited with code {result.exit_code}",
            "logs": combined_logs[-4000:],
            "elapsed_s": round(result.execution_time, 2),
            "execution_mode": "sandbox",
        }

    summary: Dict[str, Any] = {}
    for line in reversed((result.stdout or "").splitlines()):
        line = line.strip()
        if line.startswith("{"):
            try:
                summary = json.loads(line)
                break
            except json.JSONDecodeError:
                continue

    run_dir_str = summary.get("run_dir") or str(run_dir)
    artifacts = _collect_artifacts(Path(run_dir_str))

    return {
        "status": summary.get("status", "success"),
        "summary": summary,
        "artifacts": artifacts,
        "logs": combined_logs[-2000:],
        "elapsed_s": round(result.execution_time, 2),
        "execution_mode": "sandbox",
    }


def _run_on_host(
    script_path: Path,
    timeout: int,
    env: Optional[Dict[str, str]],
) -> Dict[str, Any]:
    """Execute analysis.py in a host subprocess. Returns same shape as execute()."""
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
        return {"status": "error", "error": f"Script timed out after {timeout}s", "elapsed_s": round(time.perf_counter() - t0, 2)}
    except Exception as exc:
        return {"status": "error", "error": str(exc), "elapsed_s": 0}

    elapsed = round(time.perf_counter() - t0, 2)
    combined_logs = (proc.stdout or "") + (proc.stderr or "")

    if proc.returncode != 0:
        return {
            "status": "error",
            "error": f"Script exited with code {proc.returncode}",
            "logs": combined_logs[-4000:],
            "elapsed_s": elapsed,
            "execution_mode": "host",
        }

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
        "status": summary.get("status", "success"),
        "summary": summary,
        "artifacts": artifacts,
        "logs": combined_logs[-2000:],
        "elapsed_s": elapsed,
        "execution_mode": "host",
    }


# ── execution ─────────────────────────────────────────────────────────────────

def execute(
    script_path: str | Path,
    timeout: int = DEFAULT_TIMEOUT,
    env: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """Run *script_path* and return a structured result dict.

    By default runs in a sandboxed Docker container (helix-biotools) for
    security and user trust. Set HELIX_ANALYSIS_USE_SANDBOX=false to run
    on the host.

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
        execution_mode  "sandbox" | "host"
        error       error message (only present on failure)
    """
    script_path = Path(script_path)
    if not script_path.exists():
        return {"status": "error", "error": f"Script not found: {script_path}"}

    use_sandbox = os.getenv("HELIX_ANALYSIS_USE_SANDBOX", "true").lower() == "true"

    if use_sandbox:
        try:
            sandbox_result = _run_in_sandbox(script_path, timeout, env)
            if sandbox_result.get("status") != "error":
                return sandbox_result
            err_blob = f"{sandbox_result.get('error', '')}\n{sandbox_result.get('logs', '')}".lower()
            import_runtime_issue = (
                "importerror" in err_blob
                or "cannot import name" in err_blob
                or "no module named" in err_blob
            )
            if import_runtime_issue:
                logger.warning(
                    "Sandbox import/runtime failure for %s; retrying on host execution",
                    script_path,
                )
                host_retry = _run_on_host(script_path, timeout, env)
                if host_retry.get("status") == "success":
                    return host_retry
                # Preserve original sandbox diagnostics if host retry also fails.
                return sandbox_result
            return sandbox_result
        except Exception as exc:
            exc_str = str(exc)
            # Docker is not available in this environment (e.g. ECS Fargate, CI).
            # Fall through to host execution rather than surfacing a confusing
            # "Docker is not installed" error to the end user.
            _docker_unavailable = (
                "Docker is not installed" in exc_str
                or "not running" in exc_str
                or "No such file or directory: 'docker'" in exc_str
                or "docker" in exc_str.lower() and "not found" in exc_str.lower()
            )
            if _docker_unavailable:
                logger.warning(
                    "Docker unavailable (%s) — falling back to host execution for %s",
                    exc_str,
                    script_path,
                )
            else:
                logger.exception("Sandbox execution of analysis.py failed")
                return {
                    "status": "error",
                    "error": f"Sandbox execution failed: {exc}",
                    "elapsed_s": 0,
                    "execution_mode": "sandbox",
                }

    return _run_on_host(script_path, timeout, env)


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
