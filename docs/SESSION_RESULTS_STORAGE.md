# Where execute and plan results are stored (session handling)

Results from `/execute` (single-tool and multi-step plans) are stored the same way as other session operations.

## Storage path

1. **History entry**  
   Every execute response goes through `_dispatch_result()`, which calls  
   `history_manager.add_history_entry(session_id, command, tool, result, metadata)` when `record_history=True` (default).

2. **Where it lives**
   - **In memory:** `history_manager.sessions[session_id]["history"]` — list of `{ "timestamp", "command", "tool", "result", "metadata", "run_id", "iteration_index" }`.
   - **In memory:** `history_manager.sessions[session_id]["results"]` — map of keys like `"sequence_alignment_1"`, `"__plan___2"` to the serialized result.
   - **On disk:** Session is persisted under `history_manager.storage_dir` (default `sessions/`) per session ID; see `HistoryManager._save_session()`.

3. **Plan execution (`tool="__plan__"`)**
   - The **full plan result** (broker output: `type: "plan_result"`, `steps: [...]`, last step result) is stored as one history entry with `tool="__plan__"` and one entry in `results` (e.g. `__plan___N`).
   - **Session context** is updated from each plan step via `_apply_session_context_side_effects(session_id, "__plan__", result)`: for every step (e.g. `sequence_alignment`, `mutate_sequence`), the same logic as single-tool runs (e.g. `aligned_sequences`, `mutated_sequences` on the session) so downstream tools and the UI see the same session state as after a single alignment or mutation.

So: **plan results are stored like single-tool results** (history + results dict + optional run ledger), and **per-step outputs are written into the session** so later commands in that session can use them.

## Where files live on disk

- **Paths for a session:** `history_manager.get_session_storage_paths(session_id)` returns:
  - `storage_root`: base directory (e.g. `sessions/`)
  - `session_file`: path to the session JSON (history, results, prompts)
  - `session_dir`: per-session directory (artifacts, runs)
- Session state (inputs, prompts, results) is in **session_file**; run artifacts live under **session_dir**.

## How to read them back

- **Latest result for a tool:** `history_manager.get_latest_result(session_id, "sequence_alignment")` or `"__plan__"`.
- **All results for a tool:** `history_manager.get_all_results(session_id, "__plan__")`.
- **Session summary (counts, result keys):** `history_manager.get_session_summary(session_id)` → `available_results`, `tool_usage`, etc.
- **Aligned sequences / mutated sequences for downstream:** `session_context["aligned_sequences"]`, `session_context["mutated_sequences"]` (set by `_apply_session_context_side_effects` after alignment/mutation or after a plan that includes those steps).
