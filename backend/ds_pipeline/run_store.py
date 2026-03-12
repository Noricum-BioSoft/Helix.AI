"""
Run store: writes and reads per-run artifacts and the append-only experiment log.

Artifacts layout:
    {base_dir}/artifacts/{run_id}/
        run.json        <- full run metadata (all fields from the iteration contract)
        report.md       <- narrative markdown report (written by report pipeline)
        feedback.txt    <- optional user feedback
        figures/        <- matplotlib PNGs
    {base_dir}/experiments/
        experiment_log.csv   <- one row per completed run (append-only)
        decisions.md         <- human-readable decision + next-steps log
"""
from __future__ import annotations

import csv
import json
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional


def new_run_id() -> str:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    short = uuid.uuid4().hex[:6]
    return f"run_{ts}_{short}"


_BASE_LOG_FIELDS = [
    "run_id",
    "timestamp",
    "objective",
    "hypothesis",
    "changes",
    "data_hash",
    "config_hash",
    "git_sha",
    "steps_run",
    "decision",
]


class RunStore:
    """Persists run metadata, experiment log, and artifact files to disk."""

    def __init__(
        self,
        base_dir: str | Path = ".",
        artifacts_subdir: str = "artifacts",
        experiments_subdir: str = "experiments",
    ) -> None:
        self.base_dir = Path(base_dir)
        self.artifacts_dir = self.base_dir / artifacts_subdir
        self.experiments_dir = self.base_dir / experiments_subdir
        self.artifacts_dir.mkdir(parents=True, exist_ok=True)
        self.experiments_dir.mkdir(parents=True, exist_ok=True)
        self._ensure_experiment_log()

    # ── Experiment log ────────────────────────────────────────────────────────

    def _ensure_experiment_log(self) -> None:
        log_path = self.experiments_dir / "experiment_log.csv"
        if not log_path.exists():
            with open(log_path, "w", newline="") as f:
                csv.DictWriter(f, fieldnames=_BASE_LOG_FIELDS, extrasaction="ignore").writeheader()

    def _log_fields(self) -> List[str]:
        log_path = self.experiments_dir / "experiment_log.csv"
        if log_path.exists():
            with open(log_path, "r", newline="") as f:
                reader = csv.reader(f)
                try:
                    return next(reader)
                except StopIteration:
                    pass
        return _BASE_LOG_FIELDS[:]

    def append_experiment_log(self, run_data: Dict[str, Any]) -> None:
        """Append a row to experiment_log.csv, adding metric columns as needed."""
        metrics = run_data.get("metrics") or {}
        row: Dict[str, Any] = {
            "run_id": run_data.get("run_id", ""),
            "timestamp": run_data.get("timestamp", ""),
            "objective": run_data.get("objective", ""),
            "hypothesis": run_data.get("hypothesis", ""),
            "changes": run_data.get("changes", ""),
            "data_hash": run_data.get("data_hash", ""),
            "config_hash": run_data.get("config_hash", ""),
            "git_sha": (run_data.get("env_info") or {}).get("git_sha", ""),
            "steps_run": ",".join(run_data.get("steps_run") or []),
            "decision": run_data.get("decision", ""),
        }
        for k, v in metrics.items():
            row[f"metric_{k}"] = v

        existing_fields = self._log_fields()
        new_fields = [k for k in row if k not in existing_fields]
        all_fields = existing_fields + new_fields

        log_path = self.experiments_dir / "experiment_log.csv"
        if new_fields:
            existing_rows: List[Dict[str, Any]] = []
            if log_path.exists():
                with open(log_path, "r", newline="") as f:
                    existing_rows = list(csv.DictReader(f))
            with open(log_path, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=all_fields, extrasaction="ignore")
                w.writeheader()
                w.writerows(existing_rows)

        with open(log_path, "a", newline="") as f:
            csv.DictWriter(f, fieldnames=all_fields, extrasaction="ignore").writerow(row)

    def read_experiment_log(self) -> List[Dict[str, Any]]:
        log_path = self.experiments_dir / "experiment_log.csv"
        if not log_path.exists():
            return []
        with open(log_path, "r", newline="") as f:
            return list(csv.DictReader(f))

    # ── Run directory helpers ─────────────────────────────────────────────────

    def run_dir(self, run_id: str) -> Path:
        d = self.artifacts_dir / run_id
        d.mkdir(parents=True, exist_ok=True)
        return d

    def figures_dir(self, run_id: str) -> Path:
        d = self.run_dir(run_id) / "figures"
        d.mkdir(parents=True, exist_ok=True)
        return d

    # ── run.json ──────────────────────────────────────────────────────────────

    def save_run(self, run_data: Dict[str, Any]) -> Path:
        run_id = run_data["run_id"]
        path = self.run_dir(run_id) / "run.json"
        with open(path, "w") as f:
            json.dump(run_data, f, indent=2, default=str)
        return path

    def load_run(self, run_id: str) -> Optional[Dict[str, Any]]:
        path = self.artifacts_dir / run_id / "run.json"
        if not path.exists():
            return None
        with open(path, "r") as f:
            return json.load(f)

    def list_runs(self) -> List[str]:
        return sorted(d.name for d in self.artifacts_dir.iterdir() if d.is_dir())

    def save_feedback(self, run_id: str, feedback: str) -> Path:
        path = self.run_dir(run_id) / "feedback.txt"
        path.write_text(feedback)
        run_data = self.load_run(run_id)
        if run_data:
            run_data["user_feedback"] = feedback
            self.save_run(run_data)
        return path

    # ── decisions.md ─────────────────────────────────────────────────────────

    def append_decision(self, run_id: str, decision: str, next_steps: List[str]) -> None:
        path = self.experiments_dir / "decisions.md"
        ts = datetime.now().isoformat(timespec="seconds")
        lines = [
            f"\n## {ts} — {run_id}\n\n",
            f"**Decision:** {decision}\n\n",
            "**Next steps:**\n",
        ]
        for i, step in enumerate(next_steps, 1):
            lines.append(f"{i}. {step}\n")
        with open(path, "a") as f:
            f.writelines(lines)
