from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict

import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
POLICY_PATH = REPO_ROOT / "benchmarks" / "compliance" / "SECURITY_POLICY_PROFILE_P0.yaml"
SANDBOX_BASELINE_PATH = REPO_ROOT / "benchmarks" / "compliance" / "SANDBOX_BASELINE.md"


def _load_yaml(path: Path) -> Dict[str, Any]:
    return yaml.safe_load(path.read_text()) or {}


def validate() -> int:
    if not POLICY_PATH.exists():
        print(f"Missing policy profile file: {POLICY_PATH}")
        return 1
    if not SANDBOX_BASELINE_PATH.exists():
        print(f"Missing sandbox baseline file: {SANDBOX_BASELINE_PATH}")
        return 1

    data = _load_yaml(POLICY_PATH)
    if str(data.get("profile_id") or "").strip().lower() != "p0":
        print("SECURITY_POLICY_PROFILE_P0.yaml must set profile_id: p0")
        return 1

    defaults = data.get("defaults") or {}
    required_defaults = {
        "upload_intake_enabled",
        "block_suspicious_payloads",
        "approval_for_restricted_human_data",
        "audit_enabled",
        "sandbox_required_for_untrusted_execution",
    }
    missing = sorted(required_defaults.difference(defaults.keys()))
    if missing:
        print("Missing required defaults:", ", ".join(missing))
        return 1

    enforcement_points = data.get("enforcement_points") or []
    if len(enforcement_points) < 3:
        print("Expected at least 3 policy enforcement points.")
        return 1

    print("Security policy validation passed.")
    return 0


if __name__ == "__main__":
    sys.exit(validate())
