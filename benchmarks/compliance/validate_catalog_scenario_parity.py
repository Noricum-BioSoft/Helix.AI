"""Verify every use_cases[] entry in use_case_catalog.yaml has a scenario file."""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict

import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
CATALOG_PATH = REPO_ROOT / "benchmarks" / "use_case_catalog.yaml"
CASES_DIR = REPO_ROOT / "benchmarks" / "cases"


def _load_yaml(path: Path) -> Dict[str, Any]:
    return yaml.safe_load(path.read_text()) or {}


def _source_use_case_ids_from_scenario(data: Dict[str, Any]) -> set[str]:
    """Extract all source_use_case_id values from source_refs."""
    ids: set[str] = set()
    for ref in data.get("source_refs") or []:
        uid = str(ref.get("source_use_case_id") or "").strip()
        if uid:
            ids.add(uid)
    return ids


def validate() -> int:
    catalog = _load_yaml(CATALOG_PATH)
    catalog_ids = {
        str(uc.get("id") or "").strip()
        for uc in (catalog.get("use_cases") or [])
        if uc.get("id")
    }

    covered_catalog_ids: set[str] = set()
    for path in CASES_DIR.rglob("*.yaml"):
        data = _load_yaml(path)
        covered_catalog_ids.update(_source_use_case_ids_from_scenario(data))

    missing = sorted(catalog_ids.difference(covered_catalog_ids))

    if missing:
        print("Catalog entries missing scenario files (no scenario references source_use_case_id):")
        for m in missing:
            print(f"  - {m}")
        return 1

    print(
        f"Catalog-scenario parity validated: "
        f"{len(catalog_ids)} catalog entries covered by scenario files."
    )
    return 0


if __name__ == "__main__":
    sys.exit(validate())
