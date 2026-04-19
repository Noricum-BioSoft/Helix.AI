from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, List, Set

import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
CASES_DIR = REPO_ROOT / "benchmarks" / "cases"
REGISTRY_PATH = REPO_ROOT / "benchmarks" / "compliance" / "fixture_approval_registry.yaml"


def _load_yaml(path: Path) -> Dict[str, Any]:
    return yaml.safe_load(path.read_text()) or {}


def _collect_case_fixture_paths() -> Set[str]:
    paths: Set[str] = set()
    for scenario_path in CASES_DIR.rglob("*.yaml"):
        data = _load_yaml(scenario_path)
        for fixture in data.get("input_fixtures", []) or []:
            fixture_path = str(fixture.get("path") or "").strip()
            if fixture_path:
                paths.add(fixture_path)
    return paths


def _index_registry_entries(registry: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    index: Dict[str, Dict[str, Any]] = {}
    duplicates: List[str] = []
    for entry in registry.get("fixtures", []) or []:
        fixture_path = str(entry.get("fixture_path") or "").strip()
        if not fixture_path:
            continue
        if fixture_path in index:
            duplicates.append(fixture_path)
        index[fixture_path] = entry
    if duplicates:
        raise ValueError(f"Duplicate fixture_path entries in registry: {sorted(set(duplicates))}")
    return index


def validate() -> int:
    case_paths = _collect_case_fixture_paths()
    registry = _load_yaml(REGISTRY_PATH)
    index = _index_registry_entries(registry)

    missing_registry = sorted(case_paths.difference(index.keys()))
    if missing_registry:
        print("Missing fixture registry entries:")
        for p in missing_registry:
            print(f"  - {p}")
        return 1

    invalid_status = []
    allowed_status = {"approved", "approved_with_restrictions", "pending", "rejected"}
    for fixture_path in sorted(case_paths):
        status = str(index[fixture_path].get("approval_status") or "")
        if status not in allowed_status:
            invalid_status.append((fixture_path, status))
    if invalid_status:
        print("Invalid approval_status values:")
        for p, s in invalid_status:
            print(f"  - {p}: {s}")
        return 1

    print(f"Fixture compliance validation passed for {len(case_paths)} fixture paths.")
    pending = [
        p
        for p in sorted(case_paths)
        if str(index[p].get("approval_status") or "") == "pending"
    ]
    if pending:
        print("Pending fixtures (allowed for planning, blocked for release):")
        for p in pending:
            print(f"  - {p}")
    return 0


if __name__ == "__main__":
    sys.exit(validate())
