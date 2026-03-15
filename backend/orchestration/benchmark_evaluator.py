from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

from benchmarks.bio_benchmark_scorer import score_benchmark_run, score_run_file


def evaluate_benchmark_run(run_data: Dict[str, Any], benchmark_spec: Dict[str, Any]) -> Dict[str, Any]:
    return score_benchmark_run(run_data, benchmark_spec)


def evaluate_benchmark_run_file(run_file: Path, benchmark_yaml: Path) -> Dict[str, Any]:
    return score_run_file(run_file, benchmark_yaml)

