from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


@dataclass(frozen=True)
class EvalCase:
    id: str
    input: Dict[str, Any]
    expect: Dict[str, Any]
    tags: List[str]


def load_jsonl(path: str | Path) -> List[dict]:
    p = Path(path)
    out: List[dict] = []
    for idx, line in enumerate(p.read_text(encoding="utf-8").splitlines(), start=1):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        try:
            out.append(json.loads(s))
        except Exception as e:
            raise ValueError(f"Invalid JSON on {p} line {idx}: {e}") from e
    return out


def load_cases(path: str | Path) -> List[EvalCase]:
    raw = load_jsonl(path)
    cases: List[EvalCase] = []
    for obj in raw:
        cases.append(
            EvalCase(
                id=str(obj["id"]),
                input=dict(obj.get("input") or {}),
                expect=dict(obj.get("expect") or {}),
                tags=list(obj.get("tags") or []),
            )
        )
    # Stable ordering
    cases.sort(key=lambda c: c.id)
    return cases


def assert_subset(expected_subset: Dict[str, Any], actual: Dict[str, Any]) -> None:
    """
    Assert that all keys/values in expected_subset are present in actual.
    This is intentionally strict for scalar values, and supports nested dicts.
    """
    for k, v in expected_subset.items():
        if k not in actual:
            raise AssertionError(f"Missing key '{k}' in actual params. expected_subset={expected_subset} actual={actual}")
        av = actual[k]
        if isinstance(v, dict):
            if not isinstance(av, dict):
                raise AssertionError(f"Key '{k}' expected dict but got {type(av)}. actual={actual}")
            assert_subset(v, av)
        else:
            if av != v:
                raise AssertionError(f"Mismatch for key '{k}': expected {v!r} got {av!r}. actual={actual}")


def filter_cases(cases: Iterable[EvalCase], required_tags: Optional[List[str]] = None) -> List[EvalCase]:
    if not required_tags:
        return list(cases)
    required = set(required_tags)
    return [c for c in cases if required.issubset(set(c.tags))]





