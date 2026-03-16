# Routing Flow

This page documents the request routing in `backend/main_with_mcp.py` and `backend/command_router.py`.

## Current Flow

```text
User request
  -> /execute
  -> S3 browse fast-path        (deterministic keyword match)
  -> Approval command check      (_is_approval_command → executes pending checkpoint)
  -> Universal approval gate     (CommandRouter + _should_stage_for_approval → WorkflowCheckpoint)
  -> Agent path                  (CommandProcessor, LangGraph ReAct)
  -> Fallback router path        (CommandRouter, if agent fails/disabled)
```

### Current Path Roles

| Path | Purpose | When Used |
|---|---|---|
| S3 browse fast-path | Read/list result artifacts quickly | S3 browse/list/display prompts |
| Approval command check | Detect "Approve / I approve / Yes, proceed" and execute the staged `WorkflowCheckpoint` | User approves a previously staged plan |
| Universal approval gate | Route to CommandRouter, check if response requires user approval, stage `WorkflowCheckpoint`, return `workflow_planned` | Any analytical command with clear intent (bulk RNA-seq, metadata correction, multi-step workflows) |
| Agent path | General intent/tool orchestration via LangGraph ReAct | When not short-circuited by the above paths |
| Fallback router path | Recovery when agent fails or is disabled | Agent errors / timeouts / `HELIX_AGENT_DISABLED=1` |

---

## Approval Gate Detail

The universal approval gate is the central decision point for **all analysis requests that carry scientific risk** (data modification, multi-step pipelines, parameter changes).

```text
CommandRouter.route_command(command, session_context)
    |
    +--> _should_stage_for_approval(tool, command, params)?
         |
         Yes --> _validate_plan_bindings(pending_plan)
                 |
                 +--> Always: stage WorkflowCheckpoint(WAITING_FOR_APPROVAL)
                              return status=workflow_planned
                              (binding_diagnostics added to plan text if inputs missing)
         |
         No  --> pass through to agent / fallback router
```

Key invariant: **`workflow_planned` is always returned when intent is clear**, even if required file bindings are absent. This decouples *scientific intent clarity* from *physical file availability*.

---

## Session-Aware State Machine

Every session carries a `WorkflowCheckpoint` (defined in `backend/workflow_checkpoint.py`) that tracks the current lifecycle state:

| `WorkflowState` | Meaning |
|---|---|
| `IDLE` | No active workflow |
| `WAITING_FOR_CLARIFICATION` | System asked a clarifying question |
| `WAITING_FOR_INPUTS` | Specific required inputs are missing |
| `WAITING_FOR_APPROVAL` | Plan staged, awaiting user approval |
| `PLANNING` | Plan is being constructed |
| `READY_TO_EXECUTE` | All inputs bound, awaiting execution trigger |
| `EXECUTING` | Execution in progress |
| `FAILED` | Execution failed |
| `FAILED_WAITING_FOR_USER` | Failed but recoverable with user input |
| `COMPLETED` | Execution completed successfully |

The checkpoint is saved by `history_manager.save_checkpoint()` and loaded at the start of every `/execute` request to restore interrupted workflows across turns.

---

## CommandRouter Fast-Paths

`CommandRouter.route_command()` applies deterministic keyword rules **before** falling back to LLM-based routing. Current fast-paths (in order):

| Fast-path | Trigger | Routes to |
|---|---|---|
| `_use_historical` | `"use … before/prior/original … dataset/data/version"` | `handle_natural_command` |
| `_recreate_historical` | `"recreate/reconstruct … figure/state/version/before"` | `bio_diff_runs` |
| bio diff / compare | `"compare … run/version/result"` | `bio_diff_runs` |
| S3 browse | explicit S3 list/display patterns | `s3_browse` |
| rerun | `"rerun/re-run … excluding/without"` | `bio_rerun` |
| patch and rerun | `"fix … and … rerun/regenerate/plot"` | `patch_and_rerun` |
| GO enrichment | `"go enrichment/pathway enrichment"` | `go_enrichment_analysis` |
| single-cell explicit | explicit `"single.cell/scrna"` | `single_cell_analysis` |
| bulk RNA-seq explicit | explicit `"bulk.*rnaseq/deseq"` | `bulk_rnaseq_analysis` |

The `_use_historical` fast-path is specifically important: it routes "use X from before Y" commands directly to `handle_natural_command`, **bypassing the approval gate** and preventing the agent from treating them as new pipeline requests. Without it, these commands would time out waiting for LLM-based routing.

The LLM-based `_route_with_llm` fallback is wrapped in `asyncio.wait_for` with a configurable timeout (`HELIX_GATE_ROUTE_TIMEOUT_S`, default 20 s) to prevent the `/execute` endpoint from hanging indefinitely.

---

## Frontend Workflow State Integration

The frontend (`frontend/src/App.tsx`) renders buttons and sections based on the `workflow_state` field in the backend response — not heuristics:

| `workflow_state` | UI shown |
|---|---|
| `WAITING_FOR_APPROVAL` | "I approve" button |
| `EXECUTING` / `READY_TO_EXECUTE` | "Execute Pipeline" button |
| `COMPLETED` (with download links) | "Download Results" section |
| anything else | no action button |

---

## Acceptance Checks

- Analytical commands with clear intent must return `workflow_planned` + `WAITING_FOR_APPROVAL`.
- Missing file bindings must not downgrade a plan to `needs_inputs` at the `workflow_state` level.
- "Use X from before Y" commands must route to `handle_natural_command` (not `bio_rerun` or `bulk_rnaseq_analysis`).
- LLM routing timeout must not block the `/execute` endpoint beyond `HELIX_GATE_ROUTE_TIMEOUT_S`.
- Approval commands must execute the staged checkpoint, not re-stage another plan.
