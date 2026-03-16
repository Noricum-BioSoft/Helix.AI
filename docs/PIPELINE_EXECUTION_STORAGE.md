# Pipeline execution: where it runs and what we store

## Where the pipeline is executed

- **Entry point:** `backend.execution_broker.ExecutionBroker._execute_plan_sync`
- **Invoked from:** `backend.main_with_mcp` when `POST /execute` receives a command that `_looks_like_workflow(req.command)` (e.g. "align … then calculate the consensus").
- **Flow:** `execute()` → `broker.execute_tool(ExecutionRequest(tool_name="__plan__", arguments={"plan": plan, "session_id": ...}))` → `_execute_plan_sync(plan)` → for each step `_tool_executor(step.tool_name, resolved_args)` (which is `call_mcp_tool`).

So the "code" that runs the entire pipeline is **ExecutionBroker._execute_plan_sync**; individual steps are run by **call_mcp_tool** (main_with_mcp).

## Where step directories and artifacts live

- **Session storage:** `history_manager.get_session_storage_paths(session_id)` gives:
  - **storage_root:** base directory (e.g. `sessions/`)
  - **session_file:** `{storage_root}/{session_id}.json` — session state (history, results, pipeline_execution).
  - **session_dir:** `{storage_root}/{session_id}/` — per-session directory for artifacts and runs.

- **Per-step directories:** There is no dedicated "step 1 dir", "step 2 dir" today. Generated code from the tool generator is written under `sessions/generated_code/` for debugging; execution may use temp dirs or the sandbox. The stored **pipeline_execution.steps[].step_dir** is reserved for when we add per-step working directories (currently `null`).

## What we store for pipeline runs (session JSON)

We only persist **pipeline execution information**, not full system prompts or generated code blocks:

- **history[].result** for `tool: "__plan__"` is built by `_build_pipeline_execution_storage_result()` and contains:
  - **result.pipeline_execution:**
    - **entry_point:** string describing where the pipeline runs.
    - **storage:** `storage_root`, `session_file`, `session_dir`.
    - **steps:** for each step: `step_index`, `id`, `tool_name`, `inputs`, `outputs`, `status`, `step_dir` (null for now).
  - **result.steps:** trimmed step list: `id`, `tool_name`, `arguments` (without `previous_plan_steps` blob), `result` (without `full_response`, `code_preview`).
  - **result.result:** trimmed last step result (no `full_response`/`code_preview`).

So in `sessions/{session_id}.json` you get:
- **Where the pipeline runs:** `result.pipeline_execution.entry_point`
- **Where session data lives:** `result.pipeline_execution.storage` (session_file, session_dir)
- **Per-step inputs/outputs:** `result.pipeline_execution.steps[].inputs` and `result.pipeline_execution.steps[].outputs`
- **Per-step directory:** `result.pipeline_execution.steps[].step_dir` (reserved; currently null)

## Summary

| Question | Answer |
|----------|--------|
| Code that executes the entire pipeline | `backend.execution_broker.ExecutionBroker._execute_plan_sync` (invoked from `POST /execute` in main_with_mcp when the command looks like a workflow). |
| Directory of step 1, 2, …, N | No dedicated step dirs yet; `pipeline_execution.steps[].step_dir` is reserved (null). Session-level dir: `session_dir` from `pipeline_execution.storage`. |
| Inputs of each step | `pipeline_execution.steps[].inputs` (summary of arguments; no `previous_plan_steps` blob). Full trimmed args in `result.steps[].arguments`. |
| Outputs of each step | `pipeline_execution.steps[].outputs` (short summary); full trimmed result in `result.steps[].result`. |
