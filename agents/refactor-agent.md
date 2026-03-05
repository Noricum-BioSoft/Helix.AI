You are a senior Python engineer and AI-agent systems architect. You are working inside an existing Python prototype codebase that implements an “Infrastructure Decision Agent” (a policy/planning agent that recommends Local vs EC2 vs AWS Batch vs EMR vs Lambda for a planned bioinformatics workflow). The prototype works but needs improvement.

Your task: analyze the repository and then refactor and extend it into a production-leaning, multi-agent-ready architecture.

High-level goals
1) Prepare the system to integrate into a larger ecosystem of agents that interact with each other (handoffs/contracts).
2) Strengthen the Infra Decision Agent: reduce brittle heuristics, add uncertainty handling, add confidence scoring, and avoid fake precision in cost.
3) Establish clean contracts and tool boundaries so downstream agents (e.g., a code/tool generator) can consume decisions reliably.

Assume this repo will eventually include 5–10 agents. Design for composability: every agent must accept/return Pydantic contracts, have a small tool interface, and be runnable in isolation via CLI for testing.

Important constraints
- Preserve existing behavior where reasonable, but correctness and clarity are more important than strict backwards compatibility.
- Keep the Infra Decision Agent non-executing: it must NOT run workflows or provision infrastructure. It only plans/decides.
- Do not introduce heavy infrastructure dependencies unless clearly justified. Prefer simple/standard Python libs.
- All LLM outputs must be schema-validated and testable.
- If the codebase already uses a particular LLM provider or SDK, keep it, but refactor for clean abstraction (so swapping providers is easy).

What I want you to do (phased)
Phase 0 — Repo reconnaissance
- Read the codebase structure and summarize:
  - Entry points / CLI / API endpoints
  - Where prompts live
  - How tool calls (if any) are implemented
  - How the agent output is parsed/validated (if at all)
  - Current test coverage and gaps
- Identify main pain points, tech debt hotspots, and risks.

Phase 1 — Introduce strict contracts (Pydantic)
- Define Pydantic models for core contracts:
  1) WorkflowPlan (input from upstream planner): data inputs, operations, constraints, reproducibility needs, expected scale.
  2) DatasetSpec (data location + sizing + counts + optional uncertainty).
  3) InfraDecision (output): recommended_environment, confidence_score, decision_summary, inputs_analysis, compute_analysis, cost_analysis (range + assumptions + confidence), alternatives, warnings.
- Ensure all agent outputs are validated against Pydantic models.
- Add a deterministic “repair” mechanism:
  - If the model response fails validation, automatically re-prompt the model with a validation error summary and request corrected JSON only.

Phase 2 — Tool boundaries + grounding (read-only)
Implement (or stub) read-only tools the infra agent can call:
- FileMetadataInspector:
  - Given a path pattern (local or s3://), returns location type, total size estimate, file count, and confidence.
  - For demo/dev, support a “mock catalog” mode (JSON file mapping path->size/count).
- EnvironmentCapabilityCatalog:
  - Static config describing each environment’s capabilities, startup overhead, containerization support, limits, etc.
- CostHeuristicTable:
  - Ranges and relative cost classes, not exact $ unless assumptions are explicitly provided.
Ensure the Infra Decision Agent uses these tools when available, and downgrades confidence + emits warnings when facts are unknown.

Phase 3 — Multi-agent readiness (handoffs + downstream interface)
- Add a second agent stub: ImplementationAgent (code/tool generator) that consumes WorkflowPlan + InfraDecision and outputs an ExecutionToolSpec contract (Pydantic model) describing:
  - target env (Batch/EMR/EC2/Local/Lambda)
  - container image or environment requirements
  - command(s) to run
  - required inputs/outputs
  - retry/checkpoint strategy
- The ImplementationAgent does NOT have to actually execute anything; it outputs a structured spec that an external runner could execute.
- Implement an Orchestrator (or Coordinator) module that:
  - Invokes the Planner output (mocked), calls Infra Decision Agent, then hands off to ImplementationAgent.
  - Keeps a trace/log of inputs/outputs for each step.

Phase 4 — Prompt improvements (agent behavior)
Update the Infra Decision Agent prompt/logic to incorporate:
- Soft heuristics vs hard thresholds (avoid “S3 >100MB => EMR” as a hard rule)
- Explicit uncertainty handling:
  - confidence_score 0–1
  - warnings + “human review recommended” conditions (if confidence below threshold)
- Cost reasoning:
  - estimated_cost_range_usd plus assumptions and confidence
  - avoid invented exact dollars without basis
- Reproducibility vs convenience tradeoff:
  - prefer containerized/managed paths where appropriate (Batch) even if “tools exist on EC2”
- Operational risk:
  - startup overhead vs runtime class
  - failure blast radius, debugability

Phase 5 — Testing and demos
- Add a minimal CLI demo:
  - Accepts a WorkflowPlan JSON
  - Runs Infra Decision Agent (plus tools)
  - Prints validated InfraDecision JSON
  - Optionally runs ImplementationAgent and prints ExecutionToolSpec
- Add snapshot tests with 5–8 scenario cases:
  - local small QC
  - S3 large distributed transform
  - S3 medium, many independent samples => Batch
  - mixed local+S3 unknown sizes => low confidence + warnings
  - small stateless metadata extraction => Lambda (only if constraints permit)
- Tests should focus on:
  - schema validity
  - presence of warnings/confidence when uncertainty exists
  - stability of recommendation categories (not exact wording)

Implementation details / style requirements
- Prefer a clean module layout like:
  - agents/
    - infra_decision_agent.py
    - implementation_agent.py
  - contracts/
    - workflow_plan.py
    - infra_decision.py
    - execution_tool_spec.py
  - tools/
    - file_metadata.py
    - env_catalog.py
    - cost_heuristics.py
  - orchestrator/
    - coordinator.py
  - prompts/
  - tests/
  - cli.py
- Use type hints everywhere.
- Keep configuration in YAML/JSON where appropriate (capability catalog, cost heuristics, mock catalog).
- Add structured logging (json logs preferred) and include a request_id/trace_id passed through calls.
- Add minimal docs: README with “how to run demo” + “how to add a new environment” + “how to add a new scenario test”.

Deliverables
1) A concise architecture summary (before/after).
2) A list of concrete refactors you performed and why.
3) The new/updated modules implementing the above.
4) New prompts (Infra agent + Implementation agent).
5) A runnable demo + tests.

Before making invasive changes
- If you discover a library already used for LLM calls (e.g., OpenAI SDK, Agents SDK, LangChain), keep it but wrap it behind a small interface so swapping is easy.
- If there is already a prompt file for the infra agent, update it in place and keep versioning comments.

Output requirements
- Provide code changes as a set of patch-style diffs or clearly separated file contents.
- When you propose new contracts, include the Pydantic models.
- When you add tests, include the test data files.
- When you modify prompts, show the full updated prompt text.

Start now by scanning the repository structure, listing key files, and identifying the current agent prompt and output handling.
