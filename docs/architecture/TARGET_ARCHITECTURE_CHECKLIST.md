# Target Architecture Rollout Checklist

This checklist tracks the implementation progress for evolving Helix.AI from tool-specific routing into a more general bioinformatics agent with universal sync/async execution.

Legend:
- [ ] not started
- [x] done
- [~] in progress

---

## Phase 0 — Prep and documentation
- [~] Create and maintain this checklist (`docs/architecture/TARGET_ARCHITECTURE_CHECKLIST.md`)
- [x] Add architecture diagrams (`docs/architecture/*.mmd`)

---

## Phase 1 — Generic execution (central broker + policy + generalized jobs)

### 1.1 Execution broker (single entrypoint for running tools)
- [x] Add `ExecutionBroker` / `ExecutionManager` module
- [x] Route `/execute` tool runs through the broker
- [x] Ensure broker supports:
  - [x] direct sync execution (local/EC2)
  - [x] async execution returning `job_id` (EMR path)
  - [x] consistent result envelope (success/error + artifacts)

### 1.2 Routing policy v1 (sync vs async)
- [x] Implement policy evaluation order:
  - [x] (A) bytes threshold (default 100MB)
  - [x] (B) tool-class overrides (small curated list)
  - [x] (C) timeout promotion hook (design + stub; implement later if needed)
- [x] Make threshold configurable via env var (e.g. `HELIX_ASYNC_BYTES_THRESHOLD`)

### 1.3 Input discovery + size estimation
- [x] Discover inputs from tool arguments:
  - [x] S3 URIs (`s3://bucket/key`)
  - [x] local paths
  - [x] session uploads / dataset references
- [x] Implement estimators:
  - [x] S3 size via `head_object` (ContentLength)
  - [x] local size via `stat`
- [x] Document limitations (e.g. unknown sizes for inline sequences)

### 1.4 Generalize JobManager (FastQC becomes “just a job type”)
- [x] Expand job record fields (job_type/tool_name, args, infra, status, timestamps, results)
- [x] Keep existing FastQC EMR flow working unchanged
- [x] Keep existing endpoints working:
  - [x] `GET /jobs/{job_id}`
  - [x] `GET /jobs/{job_id}/results`
- [x] Ensure job outputs land in session/job S3 paths when possible

---

## Phase 2 — General-purpose async runner on EMR
- [x] Add a universal EMR runner entrypoint that:
  - [x] reads a payload (plan/tool + args) from S3 or step args
  - [x] fetches inputs from S3
  - [x] runs steps (initially single-step v1)
  - [x] writes `results.json` + logs + artifacts to S3
- [x] Add “EMR submitter” path for non-FastQC jobs

---

## Phase 3 — Multi-step workflows (Plan IR)
- [x] Define minimal Plan IR schema (steps with inputs/outputs and references)
- [x] Update planner to emit Plan IR for workflows
- [x] Broker executes steps sequentially:
  - [x] sync sequential execution
  - [x] async execution (single EMR step running the plan v1)

---

## Phase 4 — Reduce keyword fragility and unintended tool generation
- [x] Demote `CommandRouter` to last resort
- [x] Improve tool schemas for semantic mapping (clear inputs/outputs)
- [x] Prevent tool-generator from running on pure Q&A:
  - [x] explicit intent classification
  - [x] tool-gen only when user intent is “execute”

---

## Notes / decisions
- Sync target: EC2 vs local behavior should be configurable (env-driven).
- Async target: EMR should be used for large/long tasks; initial policy is bytes-based + overrides.


