# Routing Flow

This page documents the request routing in `backend/main.py` and `backend/command_router.py`.

## Current Flow

```text
User request
  -> /execute
  -> S3 browse fast-path        (deterministic path rule, not LLM)
  -> Approval command check      (_is_approval_command → LLM classifier → executes pending checkpoint)
  -> Universal approval gate     (LLM staging classifier → WorkflowCheckpoint)
  -> Agent path                  (CommandProcessor, LangGraph ReAct)
  -> Fallback router path        (CommandRouter, if agent fails/disabled)
```

### Path Roles

| Path | Purpose | Classification |
|---|---|---|
| S3 browse fast-path | Read/list result artifacts quickly | Deterministic path rule |
| Approval command check | Detect user approval ("yes, proceed", "looks good") and execute the staged `WorkflowCheckpoint` | LLM `ApprovalClassifier` (keyword fast-path + early-rejection + LLM) |
| Universal approval gate | LLM staging classifier decides if the command needs a plan before execution | LLM `StagingClassifier` via `should_stage_for_approval` |
| Agent path | General intent/tool orchestration via LangGraph ReAct | LLM-driven |
| Fallback router path | Recovery when agent fails or is disabled | LLM `CommandRouter` (`HELIX_LLM_ROUTER_FIRST=1`) |

---

## LLM-only classification — no keyword fallbacks

Every classification decision in the routing path is LLM-driven. There are **no silent keyword fallbacks**: if an LLM is unavailable, the relevant classifier raises an explicit exception (`RoutingError`, `ApprovalClassificationError`, `IntentClassificationError`, `StagingClassificationError`) rather than guessing.

The only allowable non-LLM shortcuts are:
- **Keyword fast-path** in `ApprovalClassifier`: exact approval phrases ("yes", "proceed", "go ahead") are matched directly without an LLM call as a performance optimisation — not a fallback.
- **Early-rejection** in `ApprovalClassifier`: questions (`_INTERROGATIVE_PREFIXES`) and analytical commands (`_ANALYTICAL_PATTERNS`) are rejected without LLM as they can never be approvals.
- **High-impact action type fast-path** in `should_stage_for_approval`: known high-risk action types always require staging regardless of LLM.

---

## Approval Classifier Detail

`classify_approval(command, has_pending_plan)` in `backend/orchestration/approval_classifier.py`:

```text
1. Keyword fast-path      Is the command an exact approval phrase? → True (no LLM)
2. Early rejection        Is it clearly not approval (question / analytical verb)? → False (no LLM)
3. LLM classification     Ambiguous short phrase → LLM decides; raises on LLM failure
```

---

## Approval Gate Detail

The universal approval gate is the central decision point for all analysis requests.

```text
CommandRouter.route_command(command, session_context)   ← LLM-first (HELIX_LLM_ROUTER_FIRST=1)
    |
    +--> _should_stage_for_approval(tool, command, params)?
         |
         Yes --> _validate_plan_bindings(pending_plan)
                 |
                 +--> Always: stage WorkflowCheckpoint(WAITING_FOR_APPROVAL)
                              return status=workflow_planned
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

## CommandRouter

`CommandRouter.route_command()` uses LLM routing by default (`HELIX_LLM_ROUTER_FIRST=1`). If the LLM returns `None` (unavailable), it raises `RoutingError`.

Set `HELIX_LLM_ROUTER_FIRST=0` only for local debugging of specific keyword routing branches. Keyword routing is preserved in the code as a test/debug aid; it is not used in production.

---

## Frontend Workflow State Integration

The frontend (`frontend/src/App.tsx`) renders buttons and sections based on the `workflow_state` field in the backend response:

| `workflow_state` | UI shown |
|---|---|
| `WAITING_FOR_APPROVAL` | Plan display + approval action chips |
| `EXECUTING` / `READY_TO_EXECUTE` | "Execute Pipeline" button |
| `COMPLETED` (with download links) | "Download Results" section |
| anything else | no action button |

---

## Acceptance Checks

- Analytical commands with clear intent must return `workflow_planned` + `WAITING_FOR_APPROVAL`.
- Approval classifier must use keyword fast-path for exact approval phrases (no LLM call).
- Approval classifier must raise (not silently fall back) when LLM is unavailable and phrase is ambiguous.
- Missing file bindings must not downgrade a plan to `needs_inputs` at the `workflow_state` level.
- LLM routing timeout must not block the `/execute` endpoint beyond `HELIX_GATE_ROUTE_TIMEOUT_S`.
- Approval commands must execute the staged checkpoint, not re-stage another plan.
