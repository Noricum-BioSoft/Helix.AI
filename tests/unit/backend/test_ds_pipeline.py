"""
Unit tests for the data science pipeline.

Tests run fully offline (no LLM, no Docker, no AWS).
Each test uses a small synthetic CSV so we don't depend on real data files.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def csv_classification(tmp_path: Path) -> Path:
    """Tiny classification CSV: 60 rows, 3 numeric features, binary target."""
    import random
    random.seed(42)
    path = tmp_path / "clf_data.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["feat_a", "feat_b", "feat_c", "label"])
        for i in range(60):
            a = random.uniform(0, 10)
            b = random.uniform(-5, 5)
            c = random.uniform(0, 1)
            label = "pos" if (a + b) > 5 else "neg"
            w.writerow([a, b, c, label])
    return path


@pytest.fixture
def csv_regression(tmp_path: Path) -> Path:
    """Tiny regression CSV: 60 rows, 3 numeric features, continuous target."""
    import random
    random.seed(7)
    path = tmp_path / "reg_data.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x1", "x2", "x3", "y"])
        for i in range(60):
            x1 = random.uniform(0, 1)
            x2 = random.uniform(0, 1)
            x3 = random.uniform(0, 1)
            y = 2 * x1 + 3 * x2 + random.gauss(0, 0.1)
            w.writerow([x1, x2, x3, y])
    return path


@pytest.fixture
def csv_with_missing(tmp_path: Path) -> Path:
    """CSV with ~30% missing values in one column."""
    import random
    random.seed(1)
    path = tmp_path / "missing_data.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["a", "b_missing", "target"])
        for i in range(40):
            b = "" if random.random() < 0.3 else random.uniform(0, 10)
            w.writerow([random.uniform(0, 5), b, random.randint(0, 1)])
    return path


# ── Core modules ──────────────────────────────────────────────────────────────

class TestReproducibility:
    def test_hash_file_consistent(self, csv_classification: Path):
        from backend.ds_pipeline.reproducibility import hash_file
        h1 = hash_file(csv_classification)
        h2 = hash_file(csv_classification)
        assert h1 == h2
        assert len(h1) == 64  # SHA-256 hex

    def test_hash_config_stable(self):
        from backend.ds_pipeline.reproducibility import hash_config
        cfg = {"a": 1, "b": "hello", "c": [1, 2, 3]}
        h1 = hash_config(cfg)
        h2 = hash_config(cfg)
        assert h1 == h2

    def test_capture_env_returns_dict(self):
        from backend.ds_pipeline.reproducibility import capture_env
        env = capture_env()
        assert "python_version" in env
        assert "packages" in env
        assert isinstance(env["packages"], dict)


class TestPolicies:
    def test_missingness_warning(self):
        import pandas as pd
        from backend.ds_pipeline.policies import check_missingness, PolicyLevel
        df = pd.DataFrame({"a": [1, None, None, None], "b": [1, 2, 3, 4]})
        results = check_missingness(df, threshold=0.5)
        levels = [r.level for r in results]
        assert PolicyLevel.WARNING in levels or PolicyLevel.ERROR in levels

    def test_duplicates_detected(self):
        import pandas as pd
        from backend.ds_pipeline.policies import check_duplicates
        df = pd.DataFrame({"a": [1, 1, 2], "b": [3, 3, 4]})
        results = check_duplicates(df)
        assert len(results) == 1
        assert results[0].details["n_duplicates"] == 1

    def test_no_duplicates(self):
        import pandas as pd
        from backend.ds_pipeline.policies import check_duplicates
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
        assert check_duplicates(df) == []

    def test_pii_detected(self):
        import pandas as pd
        from backend.ds_pipeline.policies import check_pii
        df = pd.DataFrame({"email": ["a@b.com"], "value": [1]})
        results = check_pii(df)
        assert any(r.check == "pii" for r in results)

    def test_leakage_target_in_feature_name(self):
        import pandas as pd
        from backend.ds_pipeline.policies import check_leakage
        df = pd.DataFrame({"label": [0, 1], "label_score": [0.1, 0.9]})
        results = check_leakage(df, target_col="label")
        assert any("leakage" in r.check for r in results)

    def test_temporal_split_warning(self):
        import pandas as pd
        from backend.ds_pipeline.policies import check_leakage
        df = pd.DataFrame({"time": [1, 2], "feat": [1, 2], "y": [0, 1]})
        results = check_leakage(df, target_col="y", time_col="time", split_strategy="random")
        assert any(r.check == "leakage" for r in results)


class TestRunStore:
    def test_save_and_load_run(self, tmp_path: Path):
        from backend.ds_pipeline.run_store import RunStore
        store = RunStore(base_dir=tmp_path)
        run = {"run_id": "run_test_001", "timestamp": "2026-01-01", "metrics": {"accuracy": 0.9}}
        store.save_run(run)
        loaded = store.load_run("run_test_001")
        assert loaded is not None
        assert loaded["metrics"]["accuracy"] == 0.9

    def test_experiment_log_append(self, tmp_path: Path):
        from backend.ds_pipeline.run_store import RunStore
        store = RunStore(base_dir=tmp_path)
        run = {
            "run_id": "run_001", "timestamp": "2026-01-01T00:00:00",
            "objective": "test", "hypothesis": "", "changes": "",
            "data_hash": "abc", "config_hash": "def",
            "env_info": {"git_sha": "ghi"},
            "steps_run": ["ingest", "eda"],
            "decision": "first_run",
            "metrics": {"accuracy": 0.85},
        }
        store.append_experiment_log(run)
        log = store.read_experiment_log()
        assert len(log) == 1
        assert log[0]["run_id"] == "run_001"
        assert log[0].get("metric_accuracy") == "0.85"

    def test_list_runs(self, tmp_path: Path):
        from backend.ds_pipeline.run_store import RunStore
        store = RunStore(base_dir=tmp_path)
        for rid in ["run_aaa", "run_bbb"]:
            store.save_run({"run_id": rid})
        runs = store.list_runs()
        assert "run_aaa" in runs
        assert "run_bbb" in runs


# ── Pipeline steps ────────────────────────────────────────────────────────────

class TestIngest:
    def test_ingest_csv(self, csv_classification: Path):
        from backend.ds_pipeline.pipelines.ingest import ingest_csv
        conn, df, meta = ingest_csv(csv_classification)
        assert meta["n_rows"] == 60
        assert meta["n_cols"] == 4
        assert "feat_a" in [c["name"] for c in meta["columns"]]
        # DuckDB table accessible
        result = conn.execute("SELECT COUNT(*) FROM data").fetchone()
        assert result[0] == 60

    def test_ingest_missing_file(self, tmp_path: Path):
        from backend.ds_pipeline.pipelines.ingest import ingest_csv
        with pytest.raises(FileNotFoundError):
            ingest_csv(tmp_path / "nonexistent.csv")


class TestValidate:
    def test_clean_data_passes(self, csv_classification: Path):
        from backend.ds_pipeline.pipelines.ingest import ingest_csv
        from backend.ds_pipeline.pipelines.validate import validate
        _, df, _ = ingest_csv(csv_classification)
        result = validate(df, raise_on_error=False)
        assert result["n_errors"] == 0

    def test_missing_data_warns(self, csv_with_missing: Path):
        from backend.ds_pipeline.pipelines.ingest import ingest_csv
        from backend.ds_pipeline.pipelines.validate import validate
        _, df, _ = ingest_csv(csv_with_missing)
        # Use a low threshold so the ~30% missingness in b_missing triggers a warning
        result = validate(df, missingness_threshold=0.2, raise_on_error=False)
        assert result["n_warnings"] > 0


class TestClean:
    def test_removes_duplicates(self, tmp_path: Path):
        import pandas as pd
        import duckdb
        from backend.ds_pipeline.pipelines.clean import clean
        df = pd.DataFrame({"a": [1, 1, 2], "b": [3, 3, 4]})
        conn = duckdb.connect(":memory:")
        conn.register("data", df)
        _, cleaned_df, report = clean(df, conn)
        assert report["dropped_duplicate_rows"] == 1
        assert len(cleaned_df) == 2

    def test_imputes_numeric_missing(self, tmp_path: Path):
        import pandas as pd
        import duckdb
        from backend.ds_pipeline.pipelines.clean import clean
        df = pd.DataFrame({"a": [1.0, None, 3.0], "b": [4.0, 5.0, 6.0]})
        conn = duckdb.connect(":memory:")
        conn.register("data", df)
        _, cleaned_df, report = clean(df, conn)
        assert cleaned_df["a"].isnull().sum() == 0
        assert "a" in report["imputed_columns"]


class TestEDA:
    def test_eda_summary_keys(self, csv_classification: Path):
        from backend.ds_pipeline.pipelines.ingest import ingest_csv
        from backend.ds_pipeline.pipelines.eda import run_eda
        conn, df, _ = ingest_csv(csv_classification)
        summary = run_eda(conn, df, target_col="label")
        assert "shape" in summary
        assert "missingness" in summary
        assert "target_distribution" in summary
        assert summary["shape"] == [60, 4]

    def test_eda_numeric_target_stats(self, csv_regression: Path):
        from backend.ds_pipeline.pipelines.ingest import ingest_csv
        from backend.ds_pipeline.pipelines.eda import run_eda
        conn, df, _ = ingest_csv(csv_regression)
        summary = run_eda(conn, df, target_col="y")
        dist = summary["target_distribution"]
        assert dist["type"] == "numeric"
        assert "mean" in dist


class TestTrainEvaluate:
    def test_classification_baseline(self, csv_classification: Path):
        from backend.ds_pipeline.pipelines.ingest import ingest_csv
        from backend.ds_pipeline.pipelines.train import train_baseline
        from backend.ds_pipeline.pipelines.evaluate import evaluate_model
        _, df, _ = ingest_csv(csv_classification)
        model, _, _, train_info = train_baseline(df, target_col="label", random_seed=42)
        assert train_info["task_type"] == "classification"
        metrics = evaluate_model(model, train_info)
        assert "accuracy" in metrics
        assert 0 <= metrics["accuracy"] <= 1

    def test_regression_baseline(self, csv_regression: Path):
        from backend.ds_pipeline.pipelines.ingest import ingest_csv
        from backend.ds_pipeline.pipelines.train import train_baseline
        from backend.ds_pipeline.pipelines.evaluate import evaluate_model
        _, df, _ = ingest_csv(csv_regression)
        model, _, _, train_info = train_baseline(df, target_col="y", random_seed=42)
        assert train_info["task_type"] == "regression"
        metrics = evaluate_model(model, train_info)
        assert "r2" in metrics
        assert "rmse" in metrics


class TestPlanner:
    def test_propose_three_steps(self):
        from backend.ds_pipeline.planner import propose_next_steps
        run_data = {
            "eda_summary": {"shape": [100, 5], "n_categorical": 2, "n_numeric": 3},
            "validation_results": {"results": []},
            "config": {},
        }
        steps = propose_next_steps(run_data)
        assert len(steps) == 3
        for step in steps:
            assert "name" in step
            assert "score" in step
            assert step["score"] > 0

    def test_triggered_imbalance_boosted(self):
        from backend.ds_pipeline.planner import propose_next_steps
        run_data = {
            "eda_summary": {
                "shape": [200, 5], "n_categorical": 0, "n_numeric": 5,
                "target_distribution": {"value_counts": {"pos": 180, "neg": 20}},
            },
            "validation_results": {"results": []},
            "config": {},
            "train_info": {"task_type": "classification"},
        }
        steps = propose_next_steps(run_data)
        names = [s["name"] for s in steps]
        assert "handle_class_imbalance" in names


# ── Full orchestrator integration ─────────────────────────────────────────────

class TestOrchestrator:
    def test_eda_only_run(self, csv_classification: Path, tmp_path: Path):
        from backend.ds_pipeline.orchestrator import DataScienceOrchestrator, DSRunConfig
        config = DSRunConfig(
            data_path=str(csv_classification),
            target_col=None,  # EDA only
            objective="EDA only test",
        )
        orch = DataScienceOrchestrator(base_dir=tmp_path)
        run_data = orch.run(config)

        assert "run_id" in run_data
        assert "eda" in run_data["steps_run"]
        assert "train" not in run_data["steps_run"]
        assert (tmp_path / "artifacts" / run_data["run_id"] / "run.json").exists()

    def test_classification_full_run(self, csv_classification: Path, tmp_path: Path):
        from backend.ds_pipeline.orchestrator import DataScienceOrchestrator, DSRunConfig
        config = DSRunConfig(
            data_path=str(csv_classification),
            target_col="label",
            objective="Test classification",
        )
        orch = DataScienceOrchestrator(base_dir=tmp_path)
        run_data = orch.run(config)

        assert run_data["decision"] in {"first_run", "new_best", "no_change", "regressed", ""} or "hard_stop" in run_data["decision"]
        assert run_data.get("metrics") or "train:ERROR" in str(run_data["steps_run"])
        assert (tmp_path / "experiments" / "experiment_log.csv").exists()
        log = (tmp_path / "experiments" / "experiment_log.csv").read_text()
        assert run_data["run_id"] in log

    def test_report_generated(self, csv_regression: Path, tmp_path: Path):
        from backend.ds_pipeline.orchestrator import DataScienceOrchestrator, DSRunConfig
        config = DSRunConfig(
            data_path=str(csv_regression),
            target_col="y",
            objective="Report test",
        )
        orch = DataScienceOrchestrator(base_dir=tmp_path)
        run_data = orch.run(config)

        report_path = run_data.get("artifacts", {}).get("report_md")
        assert report_path and Path(report_path).exists()
        report_text = Path(report_path).read_text()
        assert "## Dataset Summary" in report_text
        assert "## Baseline Model Metrics" in report_text

    def test_reproduce_run(self, csv_classification: Path, tmp_path: Path):
        from backend.ds_pipeline.orchestrator import DataScienceOrchestrator, DSRunConfig
        config = DSRunConfig(
            data_path=str(csv_classification),
            target_col="label",
        )
        orch = DataScienceOrchestrator(base_dir=tmp_path)
        run1 = orch.run(config)
        run2 = orch.reproduce(run1["run_id"])

        assert run2["run_id"] != run1["run_id"]
        assert run2["parent_run_id"] == run1["run_id"]
        assert run2["data_hash"] == run1["data_hash"]

    def test_diff_runs(self, csv_classification: Path, csv_regression: Path, tmp_path: Path):
        from backend.ds_pipeline.orchestrator import DataScienceOrchestrator, DSRunConfig
        orch = DataScienceOrchestrator(base_dir=tmp_path)
        r1 = orch.run(DSRunConfig(data_path=str(csv_classification), target_col="label", changes="run A"))
        r2 = orch.run(DSRunConfig(data_path=str(csv_classification), target_col="label", changes="run B"))
        diff = orch.diff(r1["run_id"], r2["run_id"])
        assert "config_diff" in diff
        assert diff["config_diff"].get("changes") is not None

    def test_reviewer_text(self, csv_classification: Path, tmp_path: Path):
        from backend.ds_pipeline.orchestrator import DataScienceOrchestrator, DSRunConfig
        from backend.ds_pipeline.reviewer import review as make_review
        config = DSRunConfig(data_path=str(csv_classification), target_col="label")
        orch = DataScienceOrchestrator(base_dir=tmp_path)
        run_data = orch.run(config)
        text = make_review(run_data)
        assert "Run Summary" in text
        assert run_data["run_id"] in text
