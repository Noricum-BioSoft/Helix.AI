# What writes to the `sessions` directory

## Where is `sessions`?

- **HistoryManager** (session state, history, results): uses `storage_dir = "sessions"` by default — a **relative path**. So the actual directory is **current working directory (CWD) + "sessions"**.
  - If you start the server or run pytest from the **repo root** → `Helix.AI/sessions/`.
  - If you start from **backend/** → `backend/sessions/`.
- **Tool generator** (debug copy of generated code): writes to `PROJECT_ROOT / "sessions" / "generated_code"`, i.e. **repo root** `sessions/generated_code/`, regardless of CWD.
- **main** (e.g. script download, _save_analysis_script): uses `Path(__file__).parent.parent / "sessions"` → **repo root** `sessions/`.

So you can end up with session data under **both** `backend/sessions` (if CWD is backend) and **repo** `sessions/` (tool generator and main), depending on how you run the app or tests.

## What creates files and dirs

| Source | What it creates |
|--------|------------------|
| **HistoryManager** | `{storage_dir}/{session_id}.json`, `{storage_dir}/{session_id}.json.tmp`, and directory `{storage_dir}/{session_id}/`. |
| **Tool generator** | One `.py` file per execution under `sessions/generated_code/` (repo root). |
| **main _save_analysis_script** | `sessions/{session_id}/runs/{run_id}/analysis.py` when BioOrchestrator/code_generator path runs. |
| **agent_tools** (`_get_session_local_dir`) | `sessions/{session_id}/` and, when tools run, `sessions/{session_id}/runs/{run_id}/workspace/`, `artifacts/`, etc. |
| **RunStore** (when used with session_dir) | `{session_dir}/artifacts/`, `{session_dir}/experiments/`, `experiment_log.csv`, and per run `artifacts/{run_id}/`. |
| **JobManager** | `sessions/jobs.json` and per-job dirs under session dir when FastQC jobs are used. |

A single pipeline run (e.g. integration test) normally creates only a **small number** of items (one session JSON, one session dir, maybe one `generated_code` file, and optionally one run dir). It should **not** create 100+ files by itself.

## Why you might see 100+ files

1. **CWD and two “sessions” trees**  
   If the process CWD is `backend/`, HistoryManager writes to `backend/sessions/`, while the tool generator and main still use repo-root `sessions/`. You may be clearing one tree and then seeing the other, or mixing both.

2. **Many runs or sessions over time**  
   Old session JSONs, session dirs, and `generated_code` files are never auto-deleted. Repeated test runs or use of the app will keep adding files.

3. **External tools**  
   LangChain/LangSmith tracing (if enabled) or other tools may write under the project; they usually use their own dirs (e.g. `.langgraph_api`), not `sessions/`, but it’s worth checking.

4. **Run the test from repo root**  
   To avoid splitting session data between `backend/sessions` and repo `sessions/`, run pytest from the **repository root**:  
   `pytest tests/integration/test_handle_natural_command_pipeline.py ...`

## Limiting test pollution

The integration test can use a **dedicated sessions subdir** so only that test’s session data lives there. See the test file: it can set `history_manager.storage_dir` to a path like `sessions/integration-test-runs` (or a temp dir) for the duration of the test so the default `sessions/` (or `backend/sessions`) is not filled by the test.
