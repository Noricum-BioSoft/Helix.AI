"""
DataScienceOrchestrator: the iterative plan → execute → evaluate → review loop.

State machine:
  INGEST → VALIDATE → CLEAN → EDA → [TRAIN → EVALUATE] → REPORT → PLAN → DONE

Each step writes structured results into run_data.
The final run_data is persisted to artifacts/{run_id}/run.json and
a row is appended to experiments/experiment_log.csv.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


@dataclass
class DSRunConfig:
    """Configuration for one analysis iteration."""
    data_path: str
    target_col: Optional[str] = None
    task_type: str = "auto"          # "classification" | "regression" | "auto"
    objective: str = "Analyze dataset"
    hypothesis: str = ""
    changes: str = "Initial run"
    random_seed: int = 42
    test_size: float = 0.2
    missingness_threshold: float = 0.5
    error_threshold: float = 0.9
    time_col: Optional[str] = None
    entity_col: Optional[str] = None
    split_strategy: str = "random"   # "random" | "temporal"
    model_kwargs: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "data_path": self.data_path,
            "target_col": self.target_col,
            "task_type": self.task_type,
            "objective": self.objective,
            "hypothesis": self.hypothesis,
            "changes": self.changes,
            "random_seed": self.random_seed,
            "test_size": self.test_size,
            "missingness_threshold": self.missingness_threshold,
            "error_threshold": self.error_threshold,
            "time_col": self.time_col,
            "entity_col": self.entity_col,
            "split_strategy": self.split_strategy,
            "model_kwargs": self.model_kwargs,
        }


class DataScienceOrchestrator:
    """
    Runs the full iterative data science workflow for a given configuration.

    Parameters
    ----------
    base_dir : str | Path
        Root directory for artifacts and experiment logs.
        Defaults to the current working directory.
    session_id : str | None
        Optional Helix.AI session ID for integrating with the run ledger.
    """

    def __init__(
        self,
        base_dir: str | Path = ".",
        session_id: Optional[str] = None,
    ) -> None:
        from backend.ds_pipeline.run_store import RunStore
        self.base_dir = Path(base_dir)
        self.session_id = session_id
        self.store = RunStore(base_dir=self.base_dir)

    # ── Public API ────────────────────────────────────────────────────────────

    def run(
        self,
        config: DSRunConfig,
        feedback: Optional[str] = None,
        parent_run_id: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Execute one full iteration of the data science loop.

        Returns run_data — a dict conforming to the iteration contract.
        """
        from backend.ds_pipeline.run_store import new_run_id
        from backend.ds_pipeline.reproducibility import capture_env, hash_file, hash_config

        run_id = new_run_id()
        ts = datetime.now().isoformat(timespec="seconds")
        logger.info(f"[DS] Starting run {run_id}")

        env_info = capture_env()
        data_hash = hash_file(config.data_path)
        config_dict = config.to_dict()
        config_hash = hash_config(config_dict)

        run_data: Dict[str, Any] = {
            "run_id": run_id,
            "timestamp": ts,
            "parent_run_id": parent_run_id,
            "env_info": env_info,
            "seeds": {"random": config.random_seed},
            "data_hash": data_hash,
            "config_hash": config_hash,
            "config": config_dict,
            "objective": config.objective,
            "hypothesis": config.hypothesis,
            "changes": config.changes,
            "steps_run": [],
            "metrics": {},
            "slice_metrics": {},
            "validation_results": {},
            "comparison": {},
            "decision": "",
            "next_steps": [],
            "user_feedback": feedback,
            "artifacts": {},
        }

        steps = [
            ("ingest",   self._step_ingest),
            ("validate", self._step_validate),
            ("clean",    self._step_clean),
            ("eda",      self._step_eda),
        ]
        if config.target_col:
            steps += [
                ("train",    self._step_train),
                ("evaluate", self._step_evaluate),
            ]
        steps += [
            ("report",   self._step_report),
            ("plan",     self._step_plan),
        ]

        state: Dict[str, Any] = {"config": config, "run_data": run_data}
        for step_name, step_fn in steps:
            try:
                logger.info(f"[DS]   step: {step_name}")
                step_fn(state)
                run_data["steps_run"].append(step_name)
            except _HardStop as e:
                logger.warning(f"[DS] Hard stop at {step_name}: {e}")
                run_data["steps_run"].append(f"{step_name}:STOPPED")
                run_data["decision"] = f"hard_stop:{step_name}"
                break
            except Exception as e:
                logger.error(f"[DS] Step {step_name} failed: {e}", exc_info=True)
                run_data["steps_run"].append(f"{step_name}:ERROR")
                run_data[f"{step_name}_error"] = str(e)

        # Persist run.json + experiment log
        run_dir = self.store.run_dir(run_id)
        self.store.save_run(run_data)
        self.store.append_experiment_log(run_data)
        decision = run_data.get("decision", "")
        self.store.append_decision(
            run_id,
            decision,
            [
                s.get("name", str(s)) if isinstance(s, dict) else str(s)
                for s in (run_data.get("next_steps") or [])
            ],
        )

        if feedback:
            self.store.save_feedback(run_id, feedback)

        # Register with Helix.AI history_manager (best-effort)
        self._register_with_helix(run_data, run_dir)

        logger.info(f"[DS] Run {run_id} complete. Decision: {decision}")
        return run_data

    def reproduce(self, run_id: str) -> Dict[str, Any]:
        """Re-run using the recorded config for the given run_id."""
        run_data = self.store.load_run(run_id)
        if not run_data:
            raise ValueError(f"Run not found: {run_id}")
        config_dict = run_data.get("config", {})
        config = DSRunConfig(**{k: v for k, v in config_dict.items() if k in DSRunConfig.__dataclass_fields__})
        config.changes = f"Reproduce of {run_id}"
        return self.run(config, parent_run_id=run_id)

    def diff(self, run_id_a: str, run_id_b: str) -> Dict[str, Any]:
        """Compare two runs: configs, metrics, next_steps."""
        a = self.store.load_run(run_id_a)
        b = self.store.load_run(run_id_b)
        if not a:
            raise ValueError(f"Run not found: {run_id_a}")
        if not b:
            raise ValueError(f"Run not found: {run_id_b}")

        def _diff_dicts(x, y):
            keys = set(list(x.keys()) + list(y.keys()))
            return {k: {"a": x.get(k), "b": y.get(k)} for k in keys if x.get(k) != y.get(k)}

        config_diff = _diff_dicts(a.get("config", {}), b.get("config", {}))
        metric_diff = _diff_dicts(a.get("metrics", {}), b.get("metrics", {}))

        return {
            "run_id_a": run_id_a,
            "run_id_b": run_id_b,
            "config_diff": config_diff,
            "metric_diff": metric_diff,
            "steps_a": a.get("steps_run", []),
            "steps_b": b.get("steps_run", []),
            "decision_a": a.get("decision", ""),
            "decision_b": b.get("decision", ""),
        }

    # ── Pipeline steps ────────────────────────────────────────────────────────

    def _step_ingest(self, state: Dict[str, Any]) -> None:
        from backend.ds_pipeline.pipelines.ingest import ingest_csv
        config: DSRunConfig = state["config"]
        conn, df, metadata = ingest_csv(config.data_path)
        state["conn"] = conn
        state["df"] = df
        state["run_data"]["ingestion_metadata"] = metadata

    def _step_validate(self, state: Dict[str, Any]) -> None:
        from backend.ds_pipeline.pipelines.validate import validate, ValidationError
        config: DSRunConfig = state["config"]
        result = validate(
            state["df"],
            target_col=config.target_col,
            time_col=config.time_col,
            entity_col=config.entity_col,
            split_strategy=config.split_strategy,
            missingness_threshold=config.missingness_threshold,
            error_threshold=config.error_threshold,
            raise_on_error=False,  # we handle hard stops ourselves below
        )
        state["run_data"]["validation_results"] = result
        if not result["passed"]:
            raise _HardStop(result["summary"])

    def _step_clean(self, state: Dict[str, Any]) -> None:
        from backend.ds_pipeline.pipelines.clean import clean
        config: DSRunConfig = state["config"]
        conn, df, report = clean(
            state["df"],
            state["conn"],
            drop_high_missingness=config.error_threshold,
        )
        state["conn"] = conn
        state["df"] = df
        state["run_data"]["cleaning_report"] = report

    def _step_eda(self, state: Dict[str, Any]) -> None:
        from backend.ds_pipeline.pipelines.eda import run_eda
        config: DSRunConfig = state["config"]
        eda_summary = run_eda(
            state["conn"],
            state["df"],
            target_col=config.target_col,
        )
        state["run_data"]["eda_summary"] = eda_summary

    def _step_train(self, state: Dict[str, Any]) -> None:
        from backend.ds_pipeline.pipelines.train import train_baseline
        config: DSRunConfig = state["config"]
        model, X_train, X_test, train_info = train_baseline(
            state["df"],
            target_col=config.target_col,
            task_type=config.task_type,
            test_size=config.test_size,
            random_seed=config.random_seed,
            model_kwargs=config.model_kwargs,
        )
        state["model"] = model
        # Store non-array fields only in run_data (arrays stay in state)
        serialisable_info = {
            k: v for k, v in train_info.items()
            if k not in {"y_train", "y_test", "X_train", "X_test", "label_encoder"}
        }
        state["train_info"] = train_info
        state["run_data"]["train_info"] = serialisable_info

    def _step_evaluate(self, state: Dict[str, Any]) -> None:
        from backend.ds_pipeline.pipelines.evaluate import evaluate_model
        from backend.ds_pipeline.evaluator import evaluate_run

        metrics = evaluate_model(state["model"], state["train_info"])
        state["run_data"]["metrics"] = metrics

        experiment_log = self.store.read_experiment_log()
        eval_result = evaluate_run(
            metrics,
            experiment_log,
            primary_metric=None,  # auto-detect
        )
        state["run_data"]["comparison"] = eval_result.get("comparison", {})
        state["run_data"]["decision"] = eval_result.get("decision", "first_run")

    def _step_report(self, state: Dict[str, Any]) -> None:
        from backend.ds_pipeline.pipelines.report import generate_report
        run_data = state["run_data"]
        run_id = run_data["run_id"]
        run_dir = self.store.run_dir(run_id)
        figures_dir = self.store.figures_dir(run_id)
        report_path = generate_report(run_data, run_dir, figures_dir)
        run_data["artifacts"]["report_md"] = report_path

        # Collect figure paths
        fig_paths = run_data.get("figure_paths") or []
        run_data["artifacts"]["figures"] = fig_paths

    def _step_plan(self, state: Dict[str, Any]) -> None:
        from backend.ds_pipeline.planner import propose_next_steps
        experiment_log = self.store.read_experiment_log()
        next_steps = propose_next_steps(state["run_data"], experiment_log)
        state["run_data"]["next_steps"] = next_steps

        # If decision wasn't set by evaluate (EDA-only mode), set it now
        if not state["run_data"].get("decision"):
            state["run_data"]["decision"] = "first_run"

    # ── Helix.AI integration ─────────────────────────────────────────────────

    def _register_with_helix(self, run_data: Dict[str, Any], run_dir: Path) -> None:
        """Best-effort: register the run and its artifacts with Helix.AI history_manager."""
        if not self.session_id:
            return
        try:
            from backend.history_manager import history_manager
            from backend.ds_pipeline.reviewer import review as make_review

            summary_text = make_review(run_data)

            artifacts = []
            report_path = run_data.get("artifacts", {}).get("report_md")
            if report_path and Path(report_path).exists():
                artifacts.append({
                    "type": "report",
                    "title": "ds_report",
                    "uri": report_path,
                    "format": "md",
                })
            for fig_path in run_data.get("figure_paths") or []:
                if Path(fig_path).exists():
                    artifacts.append({
                        "type": "plot",
                        "title": Path(fig_path).stem,
                        "uri": fig_path,
                        "format": "png",
                    })
            run_json_path = str(run_dir / "run.json")
            artifacts.append({
                "type": "run_json",
                "title": "run_metadata",
                "uri": run_json_path,
                "format": "json",
            })

            history_manager.add_history_entry(
                self.session_id,
                command=run_data.get("objective", "ds_run"),
                tool="ds_run_analysis",
                result={
                    "status": "success",
                    "text": summary_text,
                    "run_id": run_data["run_id"],
                    "metrics": run_data.get("metrics"),
                    "decision": run_data.get("decision"),
                    "next_steps": run_data.get("next_steps"),
                    "visualization_type": "ds_report",
                    "report_md": run_data.get("artifacts", {}).get("report_md"),
                },
                metadata={
                    "produced_artifacts": artifacts,
                    "parent_run_id": run_data.get("parent_run_id"),
                    "run_id": run_data["run_id"],
                },
            )
        except Exception as e:
            logger.warning(f"[DS] Could not register run with history_manager: {e}")


class _HardStop(Exception):
    """Internal: raised when a pipeline step hits a hard-stop policy violation."""
