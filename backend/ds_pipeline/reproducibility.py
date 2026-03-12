"""
Reproducibility utilities: capture environment metadata and hash data/configs.

Used to make every run uniquely identifiable and re-runnable.
"""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path
from typing import Any, Dict, Optional


def capture_env() -> Dict[str, Any]:
    """Capture Python version, platform, key package versions, and git SHA."""
    env: Dict[str, Any] = {
        "python_version": sys.version,
        "platform": platform.platform(),
        "packages": _capture_packages(),
        "git_sha": _git_sha(),
    }
    return env


def _capture_packages() -> Dict[str, Optional[str]]:
    packages = {}
    _TRACKED = ["pandas", "numpy", "scikit-learn", "duckdb", "matplotlib"]
    for pkg in _TRACKED:
        try:
            import importlib.metadata
            packages[pkg] = importlib.metadata.version(pkg)
        except Exception:
            packages[pkg] = None
    return packages


def _git_sha() -> Optional[str]:
    try:
        import subprocess
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=5,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except Exception:
        pass
    return None


def hash_file(path: str | Path) -> str:
    """Compute SHA-256 of a file's raw bytes."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_dataframe(df: Any) -> str:
    """Compute a stable SHA-256 of a pandas DataFrame (via its CSV representation)."""
    csv_bytes = df.to_csv(index=False).encode("utf-8")
    return hashlib.sha256(csv_bytes).hexdigest()


def hash_config(config: Dict[str, Any]) -> str:
    """Compute a stable SHA-256 of a config dict (sorted keys for stability)."""
    canonical = json.dumps(config, sort_keys=True, default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()
