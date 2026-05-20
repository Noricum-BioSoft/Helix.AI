# Helix.AI Beta Tester Overview

This document gives a high-level view of how Helix works in beta, what to expect when you submit prompts, and where to look for deeper technical details.

---

## What Helix.AI Is (In Beta)

Helix is a **domain-specific AI system for biology, biotech, biopharma, and bioinformatics**.
In practice, it behaves as a plan-gated, multi-agent bioinformatics assistant.

At a high level, each request follows this loop:

1. **Understand request + route intent**
2. **Propose a plan or select a specialist tool**
3. **Ask for approval when required**
4. **Execute and stream progress**
5. **Store outputs, artifacts, and state in session memory**

The product is designed to support iterative scientific workflows, not just one-off commands.

---

## What Happens When You Submit a Prompt

### 1) API entry
Prompts usually go to:

- `POST /execute` (single final response)
- `POST /execute/stream` (streaming response / progress events)

Session creation:

- `POST /create_session` or `POST /session/create`

Primary implementation: `backend/main.py`.

### 2) Router decision (deterministic first, then LLM)
Helix routes your prompt to the right tool/workflow in `backend/command_router.py`.

Routing strategy:

- **Deterministic pre-router** for high-confidence patterns (explicit tool directives, demo follow-up input blocks, specific sequence/data patterns)
- **LLM router** for open-ended or ambiguous prompts
- Explicit fallback outcomes:
  - `multi_step_workflow` (composite request needing planner/agent path)
  - `unknown_intent` (clarification needed)

### 3) Plan and approval gate
If a request should be staged before execution, Helix returns a plan (`__plan__`) and waits for user approval.

Important states are tracked per session:

- `WAITING_FOR_APPROVAL`
- `WAITING_FOR_INPUTS`
- `EXECUTING`
- `COMPLETED`
- `FAILED` / `FAILED_WAITING_FOR_USER`

State model: `backend/workflow_checkpoint.py`.

### 4) Execution path
Execution is brokered through `backend/execution_broker.py`.

Depending on tool, inputs, and routing policy, Helix may:

- Execute synchronously via local tool handlers
- Execute asynchronously as a job
- Execute multi-step plans (`__plan__`)
- Route into Nextflow-backed execution for specific workflow tools

### 5) Results, jobs, and streaming
You can monitor execution using job endpoints:

- `GET /jobs/{job_id}`
- `GET /jobs/{job_id}/stream`
- `GET /jobs/{job_id}/results`

And session/run/artifact endpoints:

- `GET /session/{session_id}`
- `GET /session/{session_id}/runs`
- `GET /session/{session_id}/artifacts`

---

## How Tools and Nextflow Fit Together

Helix has four execution modes, selected by router + broker policy:

1. **Direct built-in tool path (fastest path)**  
   Typical for single-tool requests like FastQC, read trimming, sequence alignment, RNA-seq analysis.  
   - Router picks a concrete tool
   - Execution broker calls the tool executor directly
   - Result is returned sync (or promoted to async if needed for timeout safety)  
   See: `backend/orchestration/tool_registry.py`, `backend/execution_broker.py`

2. **Planned multi-step path (`__plan__`)**  
   Typical for requests like "fetch -> align -> tree" or "trim -> merge -> QC".  
   - Router/agent produces a staged plan
   - User approves
   - Broker executes plan steps in sequence, carrying step outputs forward  
   See: `backend/main.py`, `backend/execution_broker.py`, `backend/plan_ir.py`

3. **Dynamic tool-generation fallback**  
   Used when no existing tool can satisfy the request safely.  
   - Tool Generator builds/executes a purpose-built script
   - Isolated temp working directory is used to avoid polluting the repo
   - Outputs are still returned through normal Helix response contracts  
   See: `backend/tool_generator_agent.py` and `agents/tool-generator-agent.md`

4. **Nextflow-backed workflow execution**  
   Used for supported pipeline families where Nextflow is the best execution engine.
   Current tool->pipeline mapping includes:
   - `chip_seq_analysis` -> `workflows/chip_seq.nf`
   - `atac_seq_analysis` -> `nf-core/atacseq`
   - `genome_assembly` -> `nf-core/mag`
   - `variant_calling` -> `nf-core/sarek`
   - `metagenomics_16s` -> `nf-core/ampliseq`
   - `metagenomics_shotgun` -> `nf-core/taxprofiler`
   - `rna_splicing_isoform` -> `nf-core/rnasplice`
   - `crispr_screen_analysis` -> `nf-core/crisprseq`  
   See: `backend/nextflow_executor.py` (`HELIX_TO_PIPELINE`)

### How Helix decides sync vs async

The broker estimates input size and chooses execution mode:

- **sync** for smaller/faster tasks
- **async job** for heavier tasks or timeout-protected tools

In beta, this is especially relevant for long-running QC/pipeline tasks where Helix should return a `job_id` quickly and stream progress via job endpoints.

### How Nextflow progress reaches the UI

When Nextflow is used, Helix launches it with weblog hooks (`-with-weblog`), and Nextflow sends lifecycle events back into Helix:

- `POST /internal/nextflow/events`

Those events are then fanned out to SSE/job streams so users can see run progress and final status without polling raw process logs.

---

## Sessions and Memory (What Helix Remembers)

Session persistence is handled by `backend/history_manager.py`.

Helix stores:

- conversation history (sanitized for size/safety)
- run metadata
- produced artifacts and lineage
- checkpoint state (approval/input/execution/failure lifecycle)

Session-aware parameter extraction allows follow-up prompts like "use previous output" by resolving context from prior runs/artifacts:

- `backend/session_param_extractor.py`

---

## What Beta Testers Should Expect

- **Plan-first behavior** for many non-trivial workflows
- **Approval gating** before sensitive or multi-step execution
- **Needs-inputs responses** when required data is missing
- **Streaming progress** for long-running operations
- **Session continuity** across iterative refinement

Helix is optimized for iterative scientific work: propose, review, run, refine.

---

## Demo Prompts and Reference Cases

Good places to see realistic prompts and expected behavior:

- Demo scenarios: `frontend/src/data/demoScenarios.ts`
- Router eval corpus: `tests/evals/cases/bioinformatics_router_tool_mapping.jsonl`

These are useful for validating user-facing behavior and regression testing prompt handling.

---

## When To Use Helix vs General Coding Agents

Helix and general coding agents (for example Cursor or Claude Code) are complementary, but they are optimized for different jobs.

### Use Helix when you want domain execution with workflow state

Use Helix when your request is primarily **bioinformatics/data-analysis execution** and you want:

- tool/workflow routing from natural language
- plan -> approve -> execute safety gates
- session-aware follow-ups ("use previous outputs", "rerun with changed params")
- built-in handling of common analysis artifacts and job states
- optional Nextflow-backed execution for registered workflow families

Typical Helix requests:

- "Run bulk RNA-seq DE analysis with these counts + metadata and compare interaction effects."
- "Fetch these accessions, align sequences, and build a phylogenetic tree."
- "Rerun the previous analysis with a stricter FDR threshold."

### Use general coding agents when you need custom software engineering

Use a general coding agent when the task is mostly **software development**, not domain workflow execution:

- writing/refactoring app code
- building custom scripts/libraries from scratch
- debugging frontend/backend architecture issues
- adding new APIs, tests, CI/CD, infrastructure glue
- generating novel pipeline code that is not already represented in Helix tools

Typical coding-agent requests:

- "Refactor this React component and add tests."
- "Implement a FastAPI endpoint and Pydantic schema."
- "Write a custom parser and benchmark it."

### Practical rule of thumb

- If the main question is **"run/iterate this analysis safely with memory and approvals"**, use **Helix**.
- If the main question is **"write or change software code"**, use a **general coding agent**.

### Hand-off pattern that works well

Many teams use both:

1. Use **general coding agents** to add/improve Helix capabilities (new tools, better UI, better tests).
2. Use **Helix** to run and iterate real analysis tasks with end users.

This keeps engineering velocity high while preserving domain-specific execution safety for analysis users.

---

## Technical Map (For Engineers)

Core backend flow:

- `backend/main.py`
- `backend/command_router.py`
- `backend/execution_broker.py`
- `backend/agent.py`

State and memory:

- `backend/history_manager.py`
- `backend/workflow_checkpoint.py`
- `backend/session_param_extractor.py`

Orchestration and policy docs:

- `agents/agent-responsibilities.md`
- `agents/handoff-policy.md`
- `agents/intent-detector-agent.md`
- `agents/workflow-planner-agent.md`
- `agents/infrastructure-decision-agent.md`
- `agents/implementation-agent.md`
- `agents/tool-generator-agent.md`

---

## Current Beta Positioning (One-Liner)

Helix beta is a robust plan-gated orchestration system for bioinformatics workflows with deterministic+LLM routing, explicit approval/input checkpoints, and persistent session memory for iterative analysis.
