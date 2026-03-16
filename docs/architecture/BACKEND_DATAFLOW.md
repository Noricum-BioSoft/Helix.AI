# Backend Dataflow: Complete End-to-End Flow

This document describes the complete flow of a user command through the Helix.AI backend, from initial request to final response (synchronous or asynchronous).

## Overview

When a user sends a command to Helix.AI, it goes through multiple decision points and processing stages:

1. **Request Reception** - FastAPI endpoint receives the command
2. **Preflight Checks** - Rate limiting, prompt validation, session management
2.5. **Approval Gate** - Detects approval commands or stages new plans as `WorkflowCheckpoint`
3. **Intent Detection** - Classifies command as Q&A or Execute (session-context-aware)
4. **Tool Identification** - Determines which tool to use (existing or generate new)
5. **Infrastructure Decision** - Chooses execution environment (Local/EC2/EMR/Batch/Lambda)
6. **Execution Routing** - Routes to sync or async execution
7. **Tool Execution** - Runs on selected infrastructure
8. **Result Return** - Returns synchronous result or async job status

---

## Detailed Flow

### 1. Request Reception (`/execute` endpoint)

**Location**: `backend/main_with_mcp.py::execute()`

The FastAPI endpoint receives a POST request with:
- `command`: User's natural language command
- `session_id`: Optional session identifier

```python
@app.post("/execute")
async def execute(req: CommandRequest, request: Request):
    # Process command...
```

### 2. Preflight Checks

**Location**: `backend/main_with_mcp.py::execute()`

Before processing, the system performs several checks:

#### 2.1 Prompt Length Validation
- Validates command length to prevent abuse
- Rejects overly long commands

#### 2.2 Session Management
- Creates new session if `session_id` not provided
- Ensures session exists (thread-safe)
- Retrieves session context (uploaded files, previous results)

#### 2.3 Rate Limiting
- Tracks requests per identity (session ID or IP address)
- Enforces daily rate limits
- Prevents abuse

**Code Flow**:
```python
# Validate prompt length
_validate_prompt_length(req.command)

# Create/ensure session exists
if not req.session_id:
    req.session_id = history_manager.create_session()
else:
    history_manager.ensure_session_exists(req.session_id)

# Rate limit check
identity = _get_request_identity(request, req.session_id)
_check_and_increment_daily_counter(identity)

# Get session context
session_context = history_manager.sessions[req.session_id]
```

### 2.5. Approval Gate and Session-Aware State Machine

**Location**: `backend/main_with_mcp.py` (approval gate block) + `backend/workflow_checkpoint.py`

Before intent classification, the system checks whether this turn is part of an ongoing approval workflow.

#### 2.5.1 WorkflowCheckpoint

Every session maintains a `WorkflowCheckpoint` persisted by `history_manager.save_checkpoint()`.  It holds:
- `workflow_state`: Current lifecycle state (`IDLE`, `WAITING_FOR_APPROVAL`, `EXECUTING`, `COMPLETED`, etc.)
- `pending_plan`: The staged plan awaiting user approval
- `pending_command`: The original command that produced the plan

The checkpoint is loaded at the start of every `/execute` request so interrupted workflows resume correctly across turns.

#### 2.5.2 Approval command detection

If `_is_approval_command(command)` returns `True` (e.g. "Approve", "I approve", "Yes, proceed"):
- The pending `WorkflowCheckpoint` is loaded
- The staged plan is executed immediately
- No re-staging occurs — the response reflects execution, not another plan

#### 2.5.3 Staging new plans

For all other commands, `CommandRouter.route_command()` is called (with a `HELIX_GATE_ROUTE_TIMEOUT_S` async timeout to prevent LLM routing hangs).  If `_should_stage_for_approval(tool, command, params)` returns `True`:
- A `WorkflowCheckpoint(WAITING_FOR_APPROVAL)` is saved
- The response always has `status=workflow_planned` and `workflow_state=WAITING_FOR_APPROVAL`
- Missing file bindings are noted in the plan text but do **not** downgrade the status to `needs_inputs`
- The frontend renders the "I approve" button based on `workflow_state`

**Code Flow**:
```python
# Load checkpoint from previous turn
_cp = history_manager.load_checkpoint(req.session_id)

# Approval path
if _is_approval_command(req.command) and _cp:
    pending_plan = _cp.pending_plan
    # ... execute the plan ...

# Staging path
_approval_tool, _approval_params = await asyncio.wait_for(
    run_in_executor(lambda: _approval_router.route_command(req.command, session_context)),
    timeout=HELIX_GATE_ROUTE_TIMEOUT_S,
)
if _should_stage_for_approval(_approval_tool, req.command, _approval_params):
    history_manager.save_checkpoint(req.session_id,
        WorkflowCheckpoint.waiting_for_approval(pending_plan={...}))
    return build_standard_response(..., result={"status": "workflow_planned", ...})
```

---

### 3. Intent Detection

**Location**: `backend/intent_classifier.py::classify_intent()`

The system classifies the user's intent into one of two categories:

- **`"qa"`** - Question/Answer: User is asking for information, explanation, or guidance
- **`"execute"`** - Execution: User wants to perform an action or run a tool

`classify_intent` receives `session_context` and the current `workflow_state` to enable session-aware short-circuit rules (e.g. an approval command in `WAITING_FOR_APPROVAL` state always resolves to `"execute"` without an LLM call).

#### 3.1 Classification Methods

**Primary Method: LLM-based Classification**
- Uses `agents/intent-detector-agent.md` prompt
- LLM analyzes the command and returns multi-label classification:
  - `"question"` - Informational question
  - `"action"` - Instruction to perform something
  - `"data"` - Involves datasets/dataframes
  - `"workflow"` - Multi-step procedure

**Mapping to Binary Decision**:
- If `"action"` present → `"execute"`
- If only `"question"` (no `"action"`) → `"qa"`
- Otherwise → `"execute"` (safer default)

**Fallback Method: Heuristic Classification**
- Used when LLM is unavailable (mock mode, API errors)
- Pattern matching on:
  - Question starters: "what", "why", "how", "explain"
  - Execution verbs: "run", "execute", "analyze", "align"
  - File references: S3 URIs, file paths, file extensions
  - FASTA headers: `>sequence_id`

**Code Flow**:
```python
from backend.intent_classifier import classify_intent

intent = classify_intent(req.command)
# Returns: IntentDecision(intent="execute" | "qa", reason="...")
```

#### 3.2 Intent-Based Routing

**If `intent == "qa"`**:
- Routes to Q&A handler in `CommandProcessor`
- Uses agent to generate informational response
- No tool execution occurs
- Returns structured JSON response

**If `intent == "execute"`**:
- Continues to tool identification phase
- May generate/execute tools

### 4. Tool Identification

**Location**: `backend/agent.py::CommandProcessor.process()`

For execute intents, the system attempts to identify which tool to use through multiple strategies:

#### 4.1 Strategy 1: Agent Tool Mapping (Primary)

**Location**: `backend/agent.py::CommandProcessor._extract_tool_mapping_from_stream()`

The BioAgent (LangGraph ReAct agent) analyzes the command and attempts to map it to an existing tool:

1. **Streaming Extraction** (Fast):
   - Uses `agent.stream()` to watch execution in real-time
   - Looks for tool calls in streaming events
   - Stops early if tool call found
   - May miss tool calls if not in expected event format

2. **Invoke Extraction** (More Reliable):
   - Uses `agent.invoke()` to run agent to completion
   - Checks all result messages for tool calls
   - More reliable but slower

**If tool mapping found**:
- Returns `{"status": "tool_mapped", "tool_name": "...", "parameters": {...}}`
- Execution continues to infrastructure decision phase

**If no tool mapping found**:
- Falls back to Strategy 2 (Router)

#### 4.2 Strategy 2: Command Router (Fallback)

**Location**: `backend/command_router.py::CommandRouter.route_command()`

Rule-based keyword matching for known command patterns:

- Matches keywords to tool names
- Extracts parameters from command
- Fast but only works for predefined patterns

**Example Mappings**:
- `"align sequences"` → `sequence_alignment`
- `"mutate sequence"` → `mutate_sequence`
- `"phylogenetic tree"` → `phylogenetic_tree`

**If router matches**:
- Returns `(tool_name, parameters)`
- Execution continues to infrastructure decision phase

**If router doesn't match**:
- Falls back to Strategy 3 (Tool Generator)

#### 4.3 Strategy 3: Tool Generator (Dynamic Creation)

**Location**: `backend/tool_generator_agent.py::generate_and_execute_tool()`

If no existing tool matches, the system generates a new tool dynamically:

1. **Intent Check** (Safety):
   - Re-checks intent (prevents tool generation for Q&A)
   - Only generates tools for execute intent

2. **Input Discovery**:
   - Discovers input files from arguments and session context
   - Determines file sizes and locations (S3 vs local)

3. **Infrastructure Decision**:
   - Calls `infrastructure_decision_agent` to determine execution environment
   - Considers file sizes, locations, computational requirements

4. **Tool Generation**:
   - Uses `agents/tool-generator-agent.md` prompt
   - LLM researches appropriate bioinformatics tools
   - Generates Python code for the operation
   - Code is infrastructure-aware (EC2/EMR/Local patterns)

5. **Code Execution**:
   - Executes generated code on selected infrastructure
   - Returns results

**Code Flow**:
```python
# In main_with_mcp.py::call_mcp_tool()
if tool_name == "tool_generator":
    # Discover inputs
    discovered_inputs = _discover_inputs_from_args(arguments, session_context)
    
    # Generate and execute tool
    result = await generate_and_execute_tool(
        command=command,
        inputs=discovered_inputs,
        session_context=session_context
    )
```

### 5. Infrastructure Decision

**Location**: `backend/infrastructure_decision_agent.py::decide_infrastructure()`

For tool execution (especially tool generator), the system determines the optimal execution environment.

#### 5.1 Decision Factors

The infrastructure decision agent considers:

1. **File Location**:
   - S3: Data in AWS S3 buckets
   - Local: Data on local filesystem
   - Other: GCS, Azure Blob, etc.

2. **File Size**:
   - Small (<100MB): Suitable for local/EC2
   - Medium (100MB-10GB): May require cloud resources
   - Large (>10GB): Typically requires distributed processing (EMR)

3. **Computational Requirements**:
   - CPU: Single-threaded vs multi-threaded
   - Memory: RAM requirements
   - Runtime: Expected execution time
   - Parallelization: Whether operation can be distributed

4. **Tool Dependencies**:
   - Pre-installed tools (EC2 has bioinformatics tools)
   - Containerization needs
   - Installation complexity

5. **Cost Considerations**:
   - Data transfer costs (S3 egress ~$0.09/GB)
   - Compute costs (EC2, EMR, Batch pricing)
   - Time efficiency trade-offs

#### 5.2 Decision Matrix

| File Location | File Size | Computational Needs | Recommended Infrastructure |
|---------------|-----------|---------------------|---------------------------|
| S3 | >100MB | Any | **AWS EMR** |
| S3 | <100MB | Low/Medium | **Local or EC2** |
| S3 | <100MB | High/Parallel | **AWS Batch** |
| S3 | Any | Distributed/Spark | **AWS EMR** |
| Local | <100MB | Any | **Local execution** |
| Local | 100MB-10GB | Low/Medium | **Local or EC2** |
| Local | >10GB | Any | **Upload to S3 + EMR** |
| Any | Any | <15min, <10GB RAM | **Lambda** |

#### 5.3 Decision Methods

**Primary Method: LLM-based Decision**
- Uses `agents/infrastructure-decision-agent.md` prompt
- LLM analyzes operation and returns structured decision with reasoning

**Fallback Method: Heuristic Decision**
- Simple rules based on file sizes and locations
- Uses `HELIX_ASYNC_BYTES_THRESHOLD` (default: 100MB)
- Checks `HELIX_USE_EC2` environment variable

**Code Flow**:
```python
from backend.infrastructure_decision_agent import decide_infrastructure, InputAsset

inputs = [
    InputAsset(uri="s3://bucket/file.fastq", size_bytes=250*1024*1024, source="args")
]

decision = await decide_infrastructure(
    command="merge reads",
    inputs=inputs,
    session_context=session_context
)
# Returns: InfrastructureDecision(
#     infrastructure="EMR",
#     reasoning="Files on S3 are large (>100MB). EMR recommended...",
#     ...
# )
```

### 6. Execution Routing (Sync vs Async)

**Location**: `backend/execution_broker.py::ExecutionBroker.execute_tool()`

The execution broker determines whether to execute synchronously or asynchronously.

#### 6.1 Input Discovery

Before routing, the broker discovers all input files:
- From tool arguments
- From session context (uploaded files, previous results)
- Determines file sizes and locations

**Code Flow**:
```python
inputs = self._discover_inputs(req.arguments or {}, req.session_context or {})
estimated_bytes, unknown = self._estimate_total_bytes(inputs)
```

#### 6.2 Routing Policy

The broker evaluates routing policy based on:

1. **File Size Threshold**:
   - Default: `HELIX_ASYNC_BYTES_THRESHOLD` (100MB)
   - If total input size > threshold → async
   - If total input size < threshold → sync

2. **Tool Overrides**:
   - `FORCE_ASYNC_TOOLS = {"fastqc_quality_analysis"}`
   - `FORCE_SYNC_TOOLS = set()`
   - Overrides size-based routing

3. **Timeout Promotion** (Future):
   - If operation expected to take >X minutes → async

**Code Flow**:
```python
decision = await self._evaluate_routing_policy(
    tool_name=req.tool_name,
    estimated_bytes=estimated_bytes,
    unknown_inputs=unknown,
    inputs=inputs,
    command=req.original_command,
    session_context=req.session_context,
)
# Returns: RoutingDecision(mode="sync" | "async", reason="...", ...)
```

#### 6.3 Routing Decision

**If `decision.mode == "async"`**:
- Submits job to EMR (or other async infrastructure)
- Returns job ID and status endpoint
- User polls for results

**If `decision.mode == "sync"`**:
- Executes tool immediately
- Waits for completion
- Returns results directly

### 7. Tool Execution

#### 7.1 Synchronous Execution

**Location**: `backend/main_with_mcp.py::call_mcp_tool()`

For sync execution, tools run in one of these environments:

**Local Execution**:
- Runs on the backend server
- Direct file I/O, subprocess calls
- May download S3 files first

**EC2 Execution**:
- Requires `HELIX_USE_EC2=true`
- Uses `backend/ec2_executor.py`
- SSH to EC2 instance with pre-installed tools
- Tools available: `bbmerge.sh`, `samtools`, `bcftools`, etc.

**Tool Generator (Local/EC2)**:
- Generated code executes on selected infrastructure
- Infrastructure-aware code patterns

**Code Flow**:
```python
# In execution_broker.py
if decision.mode == "sync":
    output = await self._tool_executor(req.tool_name, tool_args)
    # _tool_executor routes to:
    # - Local execution
    # - EC2 execution (if HELIX_USE_EC2=true)
    # - Tool generator (if tool doesn't exist)
```

#### 7.2 Asynchronous Execution

**Location**: `backend/execution_broker.py::ExecutionBroker._submit_*_job()`

For async execution, jobs are submitted to AWS EMR:

**FastQC Jobs**:
- Special handling for `fastqc_quality_analysis`
- Uses `backend/job_manager.py` for EMR step submission

**Universal EMR Jobs**:
- Generic EMR job submission for any tool
- Creates EMR cluster if needed
- Submits Spark/PySpark job
- Stores results in S3

**Plan IR Jobs** (Multi-step workflows):
- Executes multi-step workflows on EMR
- Each step runs as EMR step
- Results flow between steps via S3

**Code Flow**:
```python
if decision.mode == "async":
    if req.tool_name == "fastqc_quality_analysis":
        result = await self._submit_fastqc_job(req)
    else:
        result = await self._submit_universal_emr_job(req)
    
    # Returns: {
    #     "status": "submitted",
    #     "job_id": "...",
    #     "status_endpoint": "/jobs/{job_id}/status",
    #     ...
    # }
```

### 8. Result Return

#### 8.1 Synchronous Results

**Location**: `backend/main_with_mcp.py::execute()`

For sync execution, results are returned immediately:

1. **History Storage**:
   - Stores command, tool, and result in session history
   - Updates session context with results

2. **Response Building**:
   - Builds standardized JSON response
   - Includes tool name, result, status, artifacts

3. **Response Return**:
   - Returns JSON to frontend
   - Frontend displays results immediately

**Code Flow**:
```python
# Store in history
history_manager.add_history_entry(
    req.session_id,
    req.command,
    tool_name,
    result
)

# Build response
standard_response = build_standard_response(
    prompt=req.command,
    tool=tool_name,
    result=result,
    session_id=req.session_id,
    mcp_route="/execute",
    success=True
)

# Return to frontend
return CustomJSONResponse(standard_response)
```

#### 8.2 Asynchronous Results

**Location**: `backend/main_with_mcp.py::execute()` + `/jobs/{job_id}/status` endpoint

For async execution, results are returned via job status polling:

1. **Initial Response**:
   - Returns job submission confirmation
   - Includes `job_id` and `status_endpoint`
   - Status: `"submitted"` or `"running"`

2. **Job Status Polling**:
   - Frontend polls `/jobs/{job_id}/status`
   - Returns current job status:
     - `"submitted"` - Job queued
     - `"running"` - Job executing
     - `"completed"` - Job finished successfully
     - `"failed"` - Job failed

3. **Result Retrieval**:
   - When status is `"completed"`, frontend retrieves results
   - Results stored in S3, referenced in job status
   - Frontend downloads/displays results

**Code Flow**:
```python
# Initial async response
{
    "status": "submitted",
    "type": "execution_result",
    "mode": "async",
    "tool_name": "fastqc_quality_analysis",
    "result": {
        "job_id": "job_123",
        "status": "submitted",
        "status_endpoint": "/jobs/job_123/status",
        "emr_cluster_id": "j-xxx",
        "s3_output_path": "s3://bucket/results/job_123/"
    },
    "routing": {
        "mode": "async",
        "reason": "Input size (250MB) exceeds threshold (100MB)",
        ...
    }
}

# Status polling endpoint
@app.get("/jobs/{job_id}/status")
async def get_job_status(job_id: str):
    # Query EMR for job status
    # Return current status and results if completed
```

---

## Flow Diagrams

### Complete Flow Diagram

See `docs/architecture/complete-dataflow.mmd` for a visual representation of the entire flow.

### Key Decision Points

```
User Command
    ↓
[Preflight Checks]  (rate limit, prompt length, session init)
    ↓
[Approval Gate]
    ├─→ _is_approval_command? → Load WorkflowCheckpoint → Execute Staged Plan
    └─→ CommandRouter + _should_stage_for_approval?
        ├─→ Yes → Save WorkflowCheckpoint(WAITING_FOR_APPROVAL)
        │         → Return workflow_planned response → Frontend shows "I approve" button
        └─→ No  → Continue ↓
    ↓
[Intent Detection]  (session-context-aware: qa vs execute)
    → "qa" → [Q&A Handler] → Response
    ↓ "execute"
[Tool Identification]
    ├─→ Agent Mapping (LangGraph ReAct) → Tool Found
    ├─→ Router Fast-Path (_use_historical, bio_diff_runs, patch_and_rerun, …) → Tool Found
    └─→ Tool Generator → Generate Tool
    ↓
[Infrastructure Decision] → Local/EC2/EMR/Batch/Lambda
    ↓
[Execution Routing] → Sync/Async
    ↓
[Tool Execution]
    ├─→ Sync → Local/EC2 → Immediate Results
    └─→ Async → EMR → Job Submission → Polling
    ↓
[Result Return] → build_standard_response → Frontend
    (workflow_state field drives button rendering)
```

---

## Example Scenarios

### Scenario 1: Simple Q&A

**User Command**: "What is sequence alignment?"

**Flow**:
1. Preflight checks pass
2. Intent detection: `"qa"` (question starter detected)
3. Routes to Q&A handler
4. Agent generates informational response
5. Returns response immediately

**Result**: Synchronous JSON response with explanation

---

### Scenario 2: Simple Tool Execution (Sync)

**User Command**: "Align these sequences: >seq1\nATCG..."

**Flow**:
1. Preflight checks pass
2. Intent detection: `"execute"` (FASTA header detected)
3. Tool identification: Agent maps to `sequence_alignment`
4. Infrastructure decision: Local (small input, no S3 files)
5. Execution routing: Sync (small input <100MB)
6. Tool execution: Runs locally
7. Result return: Immediate results

**Result**: Synchronous JSON response with alignment results

---

### Scenario 3: Large File Tool Execution (Async)

**User Command**: "Run quality analysis on s3://bucket/large_R1.fastq and s3://bucket/large_R2.fastq"

**Flow**:
1. Preflight checks pass
2. Intent detection: `"execute"` (execution verb + S3 URIs)
3. Tool identification: Agent maps to `fastqc_quality_analysis`
4. Input discovery: Discovers S3 files, checks sizes (250MB each)
5. Infrastructure decision: EMR (S3 files >100MB)
6. Execution routing: Async (total input 500MB > 100MB threshold)
7. Tool execution: Submits EMR job
8. Result return: Job ID and status endpoint

**Result**: Asynchronous response with job ID, frontend polls for status

---

### Scenario 4: Tool Generation

**User Command**: "Merge forward reads from s3://bucket/R1.fastq and reverse reads from s3://bucket/R2.fastq"

**Flow**:
1. Preflight checks pass
2. Intent detection: `"execute"` (execution verb + S3 URIs)
3. Tool identification:
   - Agent mapping: No existing tool matches
   - Router fallback: No keyword match
   - Tool generator: Generates new tool
4. Input discovery: Discovers S3 files, checks sizes (250MB each)
5. Infrastructure decision: EMR (S3 files >100MB, read merging is distributed operation)
6. Tool generation: LLM generates PySpark code for read merging using BBMerge
7. Execution routing: Async (large files)
8. Tool execution: Submits EMR job with generated code
9. Result return: Job ID and status endpoint

**Result**: Asynchronous response with job ID, generated tool executes on EMR

---

### Scenario 5: Multi-Step Workflow

**User Command**: "Trim reads from s3://bucket/reads.fastq, then align them, then create a phylogenetic tree"

**Flow**:
1. Preflight checks pass
2. Intent detection: `"execute"` (workflow detected: "then")
3. Workflow detection: `_looks_like_workflow()` returns True
4. Plan IR creation: Command router builds Plan IR with steps:
   - Step 1: `read_trimming` (trim reads)
   - Step 2: `sequence_alignment` (align trimmed reads)
   - Step 3: `phylogenetic_tree` (create tree from alignment)
5. Infrastructure decision: EMR (large S3 files, multi-step workflow)
6. Execution routing: Async (workflow + large files)
7. Plan execution: Submits EMR job with Plan IR
8. Workflow execution: EMR executes steps sequentially, passing results via S3
9. Result return: Job ID and status endpoint

**Result**: Asynchronous response with job ID, workflow executes on EMR

---

## Key Components Reference

### Core Files

- **`backend/main_with_mcp.py`**: FastAPI endpoints, approval gate, tool execution orchestration, `build_standard_response`
- **`backend/workflow_checkpoint.py`**: `WorkflowState` enum and `WorkflowCheckpoint` dataclass — the session state machine definition
- **`backend/history_manager.py`**: Session management, history storage, checkpoint persistence (`save_checkpoint`, `load_checkpoint`, `clear_checkpoint`)
- **`backend/agent.py`**: BioAgent (LangGraph ReAct), tool mapping, Q&A handling
- **`backend/intent_classifier.py`**: Session-aware intent detection (Q&A vs Execute)
- **`backend/command_router.py`**: Deterministic fast-path routing + LLM-based fallback routing
- **`backend/artifact_resolver.py`**: Semantic artifact/run resolution (`resolve_semantic_reference`) — resolves "first DEG results", "before batch exclusion", etc.
- **`backend/action_plan.py`**: Maps analytical intents to tool classes (`infer_action_type`, `map_action_to_tool`)
- **`backend/plan_binding.py`**: Validates tool input bindings (`validate_tool_bindings`)
- **`backend/orchestration/approval_policy.py`**: Determines which commands require user approval before execution
- **`backend/infrastructure_decision_agent.py`**: Infrastructure decision (Local/EC2/EMR/Batch/Lambda)
- **`backend/tool_generator_agent.py`**: Dynamic tool generation
- **`backend/execution_broker.py`**: Execution routing (sync/async), job submission
- **`backend/job_manager.py`**: EMR job management, cluster creation, step submission

### Iterative Bioinformatics Tools

These tools operate on the session's run history and are central to the iterative workflow:

- **`patch_and_rerun`** (`backend/agent_tools.py`): Applies a change request to the most recent analysis script and re-executes it
- **`bio_rerun`**: Re-runs an analysis with updated parameters (e.g. after sample exclusion)
- **`bio_diff_runs`**: Compares two historical runs — parameters, outputs, DEG tables
- **`go_enrichment_analysis`**: Runs Gene Ontology enrichment on a gene list artifact
- **`bulk_rnaseq_analysis`**: Full bulk RNA-seq pipeline (DESeq2-style differential expression)
- **`single_cell_analysis`**: scRNA-seq analysis pipeline

### Agent Prompts

- **`agents/intent-detector-agent.md`**: Intent classification prompt
- **`agents/infrastructure-decision-agent.md`**: Infrastructure decision prompt
- **`agents/tool-generator-agent.md`**: Tool generation prompt
- **`agent.md`**: Main BioAgent prompt (tool mapping, Q&A)

### Configuration

- **`HELIX_MOCK_MODE`**: Disables LLM (for testing)
- **`HELIX_USE_EC2`**: Enables EC2 execution
- **`HELIX_ASYNC_BYTES_THRESHOLD`**: Size threshold for async routing (default: 100MB)
- **`OPENAI_API_KEY` / `DEEPSEEK_API_KEY`**: LLM API keys

---

## Error Handling

### Intent Classification Errors

- **LLM unavailable**: Falls back to heuristic classification
- **Invalid response**: Falls back to heuristic classification
- **Default**: Treats as `"execute"` (safer for command-style UI)

### Tool Identification Errors

- **Agent timeout**: Falls back to router
- **Agent error**: Falls back to router
- **Router no match**: Falls back to tool generator
- **Tool generator error**: Returns error response

### Infrastructure Decision Errors

- **LLM unavailable**: Falls back to heuristic decision
- **File size unknown**: Uses conservative defaults (assumes large for S3, small for local)
- **Invalid decision**: Defaults to Local

### Execution Errors

- **Sync execution failure**: Returns error in response
- **Async job submission failure**: Returns error in response
- **EMR cluster creation failure**: Returns error with retry suggestion
- **Job execution failure**: Status endpoint returns error details

---

## Performance Considerations

### Latency

- **Intent detection**: ~100-500ms (LLM call)
- **Tool mapping**: ~1-5s (agent execution)
- **Infrastructure decision**: ~100-500ms (LLM call)
- **Sync execution**: Varies by tool (seconds to minutes)
- **Async job submission**: ~1-2s (EMR API calls)

### Optimization Strategies

1. **Caching**: Session context caches file metadata
2. **Early Returns**: Streaming tool extraction stops early when tool found
3. **Parallel Operations**: Input discovery and infrastructure decision can be parallelized
4. **Async Routing**: Large operations routed to async to avoid timeouts

---

## Future Enhancements

### Planned Improvements

1. **Enhanced Routing Policy**:
   - Timeout-based promotion to async
   - Cost-based routing decisions
   - User preference overrides

2. **Better Tool Discovery**:
   - Semantic tool matching
   - Tool capability descriptions
   - Tool recommendation system

3. **Workflow Optimization**:
   - Parallel step execution where possible
   - Intermediate result caching
   - Workflow visualization

4. **Infrastructure Optimization**:
   - Spot instance support for EMR
   - Auto-scaling based on workload
   - Cost tracking and optimization

---

## Related Documentation

- **`docs/architecture/current-dataflow.mmd`**: Current dataflow diagram (simplified)
- **`docs/architecture/target-dataflow.mmd`**: Target architecture diagram
- **`docs/architecture/ARCHITECTURE_EXPLANATION.md`**: System architecture overview
- **`agents/intent-detector-agent.md`**: Intent detection agent details
- **`agents/infrastructure-decision-agent.md`**: Infrastructure decision agent details
- **`agents/tool-generator-agent.md`**: Tool generator agent details



