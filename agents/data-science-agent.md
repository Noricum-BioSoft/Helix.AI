You are an AI coding assistant. Generate a complete starter repository for a “Local-First Iterative Data Analysis Agent” that helps a user analyze a dataset through an explicit plan→execute→evaluate→review loop. The repo must be runnable locally (Python + DuckDB + pandas + scikit-learn), and it must be designed for later migration to production tooling (e.g., MLflow registry, orchestration, Spark/warehouse, feature store) without rewriting core logic.

PRIMARY GOALS
1) Iterative workflow: implement a state-machine style orchestrator that runs:
   - data audit/validation
   - cleaning
   - EDA summary
   - baseline model
   - evaluation + comparison to prior best
   - report generation
   - user feedback capture
   - next-steps planning (top 3 options)
   Each iteration MUST write structured artifacts and append to an experiment log.

2) Reproducibility: make it easy to reproduce the same results given the same input data.
   - Pin dependencies (pyproject + lock file).
   - Provide a Dockerfile for deterministic environment.
   - Capture run metadata: git SHA, package versions, python version, random seeds.
   - Compute and store hashes: data hash + config hash.
   - Provide “reproduce run” command by run_id.

3) Monitoring/traceability: every agent action and pipeline step must be logged.
   - JSONL logs with trace_id/run_id.
   - Step duration + success/failure.
   - Store per-run artifacts under artifacts/<run_id>/...

4) Extensibility: keep core logic cleanly separated from infrastructure choices.
   - Define minimal “ports/adapters” interfaces for: data access, tracking, execution.
   - Default adapters: local files + DuckDB, local artifact store, lightweight tracking.
   - Provide stubs/adapters placeholders for later: MLflow tracking, Spark/Databricks, warehouse, orchestration (Prefect/Airflow).

CONSTRAINTS
- Do NOT overengineer. Provide a minimal, usable MVP with a clear upgrade path.
- Use Python 3.11+.
- Use DuckDB for local SQL and parquet/csv handling.
- Use pandas for dataframes where needed, but prefer DuckDB for large-ish scans.
- Use scikit-learn for baseline modeling.
- Use matplotlib for any plots (no seaborn).
- No external web calls.
- Provide sensible defaults that work on a small CSV.

DELIVERABLES (MUST CREATE ALL)
A) Repository structure:
   analysis-agent/
     README.md
     AGENTS.md
     SKILLS.md
     pyproject.toml (+ lock file using uv or poetry; pick ONE and be consistent)
     .gitignore
     configs/
       project.yaml
       data.yaml
       run.yaml
       eval.yaml
     data/
       raw/ (empty, with README)
       interim/
       processed/
       README.md
     src/analysis_agent/
       __init__.py
       cli.py
       orchestrator.py
       planner.py
       executor.py
       evaluator.py
       reviewer.py
       policies.py
       pipelines/
         ingest.py
         validate.py
         clean.py
         eda.py
         features.py
         train.py
         evaluate.py
         report.py
       tracking/
         logging.py
         run_store.py
         mlflow_adapter.py (stub, optional dependency)
       adapters/
         data_local.py
         tracking_local.py
         execution_local.py
         registry_stub.py
       utils/
         hashing.py
         env_capture.py
         io.py
         schema.py
     experiments/
       experiment_log.csv (initialize with headers)
       decisions.md (starter)
     artifacts/ (empty; created on run)
     scripts/
       run_local.sh
       repro_run.sh
     infra/
       Dockerfile
       Makefile
     tests/
       test_hashing.py
       test_policies_leakage.py

B) CLI behavior (must work):
   - `python -m analysis_agent run --config configs/run.yaml`
       Runs one iteration end-to-end; creates a new run_id; writes artifacts.
   - `python -m analysis_agent reproduce --run-id <RUN_ID>`
       Re-runs using recorded configs and hashes; verifies match where feasible.
   - `python -m analysis_agent report --run-id <RUN_ID>`
       Regenerates the markdown report from stored artifacts.
   - `python -m analysis_agent diff --a <RUN_ID> --b <RUN_ID>`
       Shows what changed: configs, features, metrics summary.

C) Iteration Contract (MUST IMPLEMENT):
   - Write `artifacts/<run_id>/run.json` capturing:
     run_id, timestamp, git_sha, env info, seeds, data_hash, config_hash,
     objective, hypothesis, changes, steps_run, metrics, slice_metrics,
     validation_results, decision, next_steps, user_feedback (if any).
   - Append a row to experiments/experiment_log.csv with key fields.

D) Minimal Planning + Review:
   - Planner must propose top 3 next experiments using a heuristic:
     score = (expected_impact * confidence) / cost
   - Reviewer must write a concise summary report:
     what ran, what changed, key results, what to do next
   - Feedback capture:
     store user feedback in `artifacts/<run_id>/feedback.txt` and `run.json`.
     In local mode, accept feedback via CLI flag `--feedback "..."` and also
     allow interactive prompt if flag is absent (safe basic input()).

E) Data Audit + Policy Guardrails:
   Implement policies with hard stops or warnings:
   - missingness thresholds
   - duplicate row detection
   - leakage heuristics:
     * if time column exists, warn if random split chosen
     * if entity id exists, warn if entity leakage likely
     * warn if target appears in feature names
   - basic PII detection (email/phone/name-like columns) => warn and list columns.

F) Reporting:
   - Generate `artifacts/<run_id>/report.md`
   - Include: dataset summary, data quality, EDA highlights, baseline metrics,
     comparisons vs prior best (if exists), next steps.
   - Generate a small set of plots into `artifacts/<run_id>/figures/` (matplotlib).
     Keep it minimal (e.g., target distribution, missingness bar chart).

DEFAULTS / ASSUMPTIONS
- Assume input data is a CSV placed in data/raw/ and referenced in configs/data.yaml.
- Support both supervised (if target specified) and unsupervised EDA-only mode.
- Baseline model:
  - classification: logistic regression
  - regression: ridge regression
  Auto-detect based on target type in configs (or basic heuristics).

QUALITY BAR
- Code must be clean, typed where reasonable, and documented.
- Provide clear README: setup, quickstart, how to add data, how to run.
- Provide minimal tests that actually pass.
- Ensure the project runs without modification after `pip install -e .` (or chosen tool).
- Prefer small, understandable modules over cleverness.

OUTPUT FORMAT
- Provide the full repository tree with file contents.
- Ensure all referenced modules and imports exist.
- Do not leave TODO placeholders except for clearly marked adapters (mlflow/registry stubs).