# Helix.AI Agent Responsibilities & Handoff Contracts

This document defines **clear, non-overlapping responsibilities** for each agent and the **contracts** they exchange.

## System-wide Rules

### 1) Single-Owner Rule
Every “decision type” has exactly one owner. Other agents may provide suggestions, but only the owner may decide.

| Decision Type | Owner |
|---|---|
| Intent (ask vs execute) | Intent Detector |
| Workflow / task decomposition (what to do) | Bioinformatics Executor (Planner) |
| Infrastructure selection (where to do it) | Infrastructure Expert |
| Packaging / filling tool gaps (how to run when tools are missing) | Code Generator (Tooling/Packaging Agent) |
| Visualization / reporting | Data Visualizer |
| Side effects (job submission, I/O, long-running tasks) | Execution Broker (non-LLM service) |

### 2) Contract Rule
Agents exchange **typed contracts**. No agent should rely on parsing another agent’s prose.

Core contracts:
- IntentResult
- WorkflowPlan (Plan IR + metadata)
- InfraDecision
- ExecutionSpec
- ExecutionResult
- VisualizationArtifacts (or VisualizationSpec)

### 3) No-Mutation Rule
Agents do not mutate upstream contracts. They may:
- add annotations (warnings, assumptions, constraints)
- produce a new downstream contract

If an upstream artifact must change, create a new version with explicit provenance.

### 4) Clarifying Questions Policy
- Agents may ask clarifying questions only within their scope.
- If missing info can be obtained via read-only tools (e.g., file sizing), prefer tool lookup over asking the user.

---

## Agent Specifications

## 1) Intent Detector
**Owner of:** intent classification

**Input:**
- User prompt (and optionally recent conversation context)

**Output:**
- `IntentResult { intent, confidence, clarifying_questions?, rationale? }`

**Must not:**
- answer bioinformatics questions
- propose workflows
- choose infrastructure
- trigger execution

**Invocation:**
- Always first step for every user message.

---

## 2) Bioinformatics Guru (Ask Agent)
**Owner of:** answering questions (education, interpretation, guidance)

**Input:**
- User question + optional context

**Output:**
- `Answer { text, citations?, followups?, suggested_next_actions? }`

**May:**
- suggest high-level approaches and options

**Must not:**
- output a `WorkflowPlan`
- choose infrastructure
- generate runnable job specs
- execute tools or submit jobs

**Invocation:**
- Only when `IntentResult.intent == "ask"`
- Or when Planner explicitly requests domain clarification.

---

## 3) Bioinformatics Executor (Workflow Planner)
**Owner of:** creating a workflow plan for “execute” requests

**Input:**
- User request to execute a task/workflow
- Optional: dataset hints, format constraints, desired outputs

**Output:**
- `WorkflowPlan` (Plan IR + metadata)
  - steps (tool invocations + args)
  - declared inputs/outputs
  - plan metadata: scale, parallelism shape, reproducibility
  - risk flags / missing info

**May:**
- reference the tool registry to choose known tools
- propose alternative workflows as `alternatives[]` (still within `WorkflowPlan`)

**Must not:**
- choose infrastructure (no EC2/EMR/Batch/Lambda decisions)
- generate new tools/code (unless explicitly delegated to Code Generator)
- execute jobs, upload data, or perform side effects

**Invocation:**
- When `IntentResult.intent == "execute"`

---

## 4) Infrastructure Expert
**Owner of:** selecting infrastructure for a given plan

**Input:**
- `WorkflowPlan`
- `DatasetSpec` (or file metadata tool outputs)

**Output:**
- `InfraDecision`
  - recommended environment
  - confidence score
  - assumptions + warnings
  - cost range + assumptions
  - constraints (e.g., “container required”, “prefer data locality”)
  - alternatives + tradeoffs

**May:**
- suggest non-invasive plan constraints (e.g., “split by sample”) as constraints/advice

**Must not:**
- rewrite workflow steps
- generate runnable scripts / job specs
- submit jobs or execute anything

**Invocation:**
- After a `WorkflowPlan` exists
- May be re-run if dataset metadata changes.

---

## 5) Code Generator (Tooling/Packaging Agent)
**Owner of:** filling tool gaps and producing a runnable execution specification

**Entry conditions (strict):**
Invoke ONLY if:
- Planner requires a capability not in the existing toolbox/registry, OR
- Infrastructure constraints require packaging/containerization not available

**Input:**
- `WorkflowPlan`
- `InfraDecision`
- Tool registry + environment capabilities (read-only)

**Output:**
- `ExecutionSpec`
  - target environment
  - container image / dependencies / job definition / EMR steps
  - commands/entrypoints
  - declared inputs/outputs
  - retry/checkpoint strategy

**Must not:**
- change scientific intent of the workflow
- choose infrastructure
- submit jobs or execute anything

---

## 6) Data Visualizer
**Owner of:** visualization/reporting of results

**Input:**
- `ExecutionResult` (+ artifact URIs, schemas, summaries)

**Output:**
- `VisualizationArtifacts` (plots/tables/report spec/links)

**May:**
- propose a visualization plan if artifacts aren’t available yet

**Must not:**
- redesign workflow
- choose infrastructure
- generate execution code
- execute jobs

---

## Non-Agent Component: Execution Broker (Service)
**Owner of:** all side effects and job lifecycle

**Input:**
- `ExecutionSpec`

**Output:**
- `ExecutionResult` (status, logs, artifact locations)

Responsibilities include:
- job submission and monitoring
- retries / timeouts
- log collection
- artifact registration
- returning consistent execution state
