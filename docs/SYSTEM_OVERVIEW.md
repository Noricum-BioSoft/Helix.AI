# Helix.AI System Overview

**Version:** 2.0  
**Last Updated:** 2026-01-18

---

## ⚠️ IMPORTANT: Two Orchestration Systems

Helix.AI currently has **two orchestration systems**:

### 1. PRIMARY System (Production) 🟢
- **File:** `backend/agent.py`
- **Entry Point:** `handle_command()`
- **Used By:** HTTP API (`/execute`, `/chat`, `/mcp/call_tool`)
- **Architecture:** Multi-agent graph (6+ agents)
- **Agents:** IntentDetector, BioinformaticsGuru, BioinformaticsExecutor, InfrastructureExpert, ExecutionBroker, DataVisualizer
- **Status:** Stable, serving all production traffic

### 2. EXPERIMENTAL System (Prototype) 🟡
- **File:** `backend/orchestrator.py`
- **Entry Point:** `Orchestrator.run_pipeline()`
- **Used By:** CLI demos (`cli_demo.py`) and unit tests only
- **Architecture:** Clean 2-agent pipeline
- **Agents:** InfrastructureDecisionAgent, ImplementationAgent
- **Status:** Prototype, **NOT used in production**

**The diagrams and descriptions below refer to the PRIMARY system unless otherwise noted.**

For detailed comparison and migration plan, see **[docs/ORCHESTRATION_DUALITY.md](./ORCHESTRATION_DUALITY.md)**

---

## What is Helix.AI?

Helix.AI is an autonomous bioinformatics assistant that interprets user requests, plans workflows, selects infrastructure, and executes analyses using real computational tools. The system operates on a **multi-agent architecture** where specialized agents handle distinct responsibilities with clear handoff contracts.

---

## Core Architecture (PRIMARY System)

### Multi-Agent System

Helix.AI uses a **single-owner rule**: every decision type has exactly one agent responsible for it. Agents exchange typed contracts and cannot perform actions outside their scope.

```
User Request
    ↓
Intent Detector → "execute" or "qa"?
    ↓
    ├─→ [qa] → Bioinformatics Guru (answers questions)
    │
    └─→ [execute] → Workflow Planner → "WHAT workflow to run?"
                          ↓
                    Infrastructure Agent → "WHERE to run?"
                          ↓
                    Implementation Agent → "HOW to package?"
                          ↓
                    Execution Broker → "Run it"
                          ↓
                    Results + Visualization
```

### Key Principle: Separation of Concerns

- **Intent Detector**: Classifies user intent (ask vs. execute)
- **Workflow Planner**: Selects workflow, decomposes into tasks, asks clarifying questions
- **Infrastructure Agent**: Decides execution environment (Local/EC2/EMR/Batch/Lambda)
- **Implementation Agent**: Creates execution specification (packaging, containers, parameters)
- **Tool Generator**: Dynamically creates tools when capabilities are missing
- **Execution Broker**: Actually executes tools and manages jobs
- **Data Visualizer**: Produces plots, tables, and reports

**No agent can skip another or perform actions outside its scope.**

---

## Foundational Concepts

### 1. Micro-/Macroflow Pattern

**Core architectural pattern defining task structure:**

- **Atomic Tasks (Microflows)**: Single operations (align, trim, QC)
- **Composite Tasks (Macroflows)**: Multi-step workflows (full RNA-seq pipeline)

Tasks can be user-specified or agent-proposed. The system maintains provenance and supports session continuity.

📖 **See:** [docs/TASK_ARCHITECTURE.md](./TASK_ARCHITECTURE.md)

---

### 2. Session Awareness

Helix.AI maintains persistent session context across interactions:
- Previous results and intermediate artifacts
- Uploaded files and dataset references
- Workflow history with provenance

**Golden Rule:** Always check session context before requesting data from the user.

📖 **See:** [docs/SESSION_MANAGEMENT.md](./SESSION_MANAGEMENT.md)

---

### 3. Structured Output

All responses follow a strict JSON schema with:
- Status codes (success, partial_success, failed, declined)
- Artifacts (sequences, tables, charts, trees, images)
- Provenance (tools used, versions, commands, parameters)

📖 **See:** [docs/OUTPUT_SCHEMA.md](./OUTPUT_SCHEMA.md)

---

## Agent Specifications

### Intent Detector
**Responsibility:** Classify user intent ("ask" or "execute")  
**Spec:** [agents/intent-detector-agent.md](../agents/intent-detector-agent.md)  
**Implementation:** `backend/intent_classifier.py`

---

### Workflow Planner
**Responsibility:** Decide WHAT workflow to execute  
**Spec:** [agents/workflow-planner-agent.md](../agents/workflow-planner-agent.md)  
**Implementation:** `backend/workflow_planner_agent.py`

**Owns:**
- Workflow playbooks (RNA-seq, scRNA-seq, WGS/WES, metagenomics, phylogenetics)
- Task decomposition (atomic vs. composite)
- Clarification protocol (ask for reference genome, organism, etc.)
- Feasibility assessment

---

### Infrastructure Decision Agent
**Responsibility:** Decide WHERE to execute (environment selection)  
**Spec:** [agents/infrastructure-decision-agent.md](../agents/infrastructure-decision-agent.md)  
**Implementation:** `backend/infrastructure_decision_agent.py`

**Decides:**
- Local (< 100MB, lightweight tools)
- EC2 (100MB - 10GB, moderate compute)
- EMR (> 10GB, distributed processing)
- AWS Batch (containerized jobs)
- Lambda (serverless, small tasks)

---

### Implementation Agent
**Responsibility:** Decide HOW to package and execute  
**Spec:** [agents/implementation-agent.md](../agents/implementation-agent.md)  
**Implementation:** `backend/implementation_agent.py`

**Creates:**
- Execution specifications (tool parameters, container images)
- Retry and checkpoint strategies
- Resource requirements

---

### Tool Generator Agent
**Responsibility:** Dynamically create tools when capabilities are missing  
**Spec:** [agents/tool-generator-agent.md](../agents/tool-generator-agent.md)  
**Implementation:** `backend/tool_generator_agent.py`

**Process:**
1. Detect tool gaps
2. Research appropriate tools and methods
3. Generate tool implementation
4. Validate and test

---

### Execution Broker
**Responsibility:** Execute tools and manage job lifecycle  
**Implementation:** `backend/execution_broker.py`

**Handles:**
- Routing (sync vs. async execution)
- Job submission and monitoring
- Input/output management
- Session state updates

---

## Workflow Playbooks

Pre-defined workflows for common bioinformatics analyses:

| Analysis Type | Workflow Steps | Clarifications Needed |
|---|---|---|
| **RNA-seq** | QC → trim → align → quantify → DE → viz | Organism, reference genome, strandedness |
| **scRNA-seq** | QC → normalize → HVGs → PCA → cluster → markers | Input format, metadata columns |
| **WGS/WES** | QC → align → mark dups → call variants → annotate | Reference genome, sample type |
| **Metagenomics** | QC → classify → abundance → diversity | Shotgun vs. amplicon, database |
| **Phylogenetics** | MSA → model selection → tree inference → viz | (Optional trimming/model selection) |

📖 **See:** [agents/workflow-planner-agent.md](../agents/workflow-planner-agent.md) for complete playbook details

---

## Guardrails and Policies

### Statistical Validity
- Report effect sizes with confidence intervals
- Apply multiple-testing correction (FDR, Bonferroni)
- Surface assumptions and confounders
- Include diagnostic plots

📖 **See:** [docs/STATISTICAL_GUIDELINES.md](./STATISTICAL_GUIDELINES.md)

---

### Safety and Privacy
- Treat human genomic data as sensitive PII
- Redact identifiable metadata
- Refuse harmful dual-use requests (pathogen enhancement, bioweapons)
- Comply with HIPAA, GDPR regulations

📖 **See:** [docs/SAFETY_POLICY.md](./SAFETY_POLICY.md)

---

## Agent Handoff Policy

**Enforced routing rules:**

1. **Intent Detector** must always be first
2. **Ask intent** → Bioinformatics Guru only
3. **Execute intent** → Workflow Planner → Infrastructure → Implementation → Execution Broker
4. **Tool Generator** is optional (invoked only when tool gaps exist)
5. **Data Visualizer** is terminal (no further handoffs)

**Illegal transitions:**
- Guru cannot call Execution Broker (read-only agent)
- Planner cannot skip Infrastructure Agent
- No agent can skip steps or call agents out of order

📖 **See:** [agents/handoff-policy.md](../agents/handoff-policy.md) and [agents/agent-responsibilities.md](../agents/agent-responsibilities.md)

---

## Non-Negotiable Principles

1. **Never fabricate computational results.** If a tool was not run, do not present derived outputs as real.
2. **Prefer real computation over reasoning-only** for algorithmic results.
3. **Reproducibility is mandatory:** Capture commands, parameters, versions, file hashes.
4. **Fail fast on invalid inputs** with actionable errors.
5. **Be transparent** about assumptions, defaults, limitations, and uncertainty.
6. **Protect privacy** and prevent harmful misuse.
7. **Leverage session context:** Check for existing data before requesting new inputs.

---

## Key Contracts

### IntentResult
```python
IntentResult(
    intent: "execute" | "qa",
    confidence: float,
    reason: str
)
```

### WorkflowPlan
```python
WorkflowPlan(
    workflow_id: str,
    description: str,
    data_inputs: List[DataInput],
    operations: List[OperationSpec],
    constraints: ConstraintSpec,
    expected_compute_intensity: "Low" | "Medium" | "High"
)
```

### InfraDecision
```python
InfraDecision(
    infrastructure: "Local" | "EC2" | "EMR" | "Batch" | "Lambda",
    confidence_score: float,  # 0-1, <0.5 suggests human review
    decision_summary: str,  # 1-2 sentence justification
    reasoning: str,  # Detailed explanation
    file_analysis: FileAnalysis,  # Input file sizes, locations, counts
    computational_requirements: ComputationalRequirements,  # CPU, memory, runtime
    cost_analysis: CostAnalysis(
        estimated_cost_range_usd: Tuple[float, float],  # (min, max) range
        cost_assumptions: str,  # Explicit assumptions (region, instance type, etc.)
        cost_confidence: float,  # Confidence in cost estimate (0-1)
        data_transfer_cost_usd: Optional[float],
        breakdown: Optional[dict]  # Cost breakdown by component
    ),
    alternatives: List[InfraAlternative],  # Alternative options with tradeoffs
    warnings: List[str],  # Warnings about decision (unknown sizes, missing data, etc.)
    inputs_analyzed: int
)
```

### ExecutionToolSpec
```python
ExecutionToolSpec(
    tool_name: str,
    infrastructure: "Local" | "EC2" | "EMR" | "Batch" | "Lambda",
    container_spec: Optional[ContainerSpec](  # None for native execution
        image: str,  # e.g., "docker.io/biocontainers/fastqc:0.11.9"
        image_type: "docker" | "singularity" | "conda" | "native",
        pull_policy: "Always" | "IfNotPresent" | "Never",
        env_vars: Dict[str, str],
        mount_paths: List[str],
        working_dir: Optional[str]
    ),
    commands: List[CommandSpec](  # Commands to execute (in order)
        name: str,  # Human-readable command name
        command: str,  # Shell command to execute
        inputs: List[str],  # Input file URIs
        outputs: List[str],  # Output file URIs
        success_criteria: Optional[str],  # e.g., "exit_code==0"
        timeout_minutes: Optional[float]
    ),
    retry_policy: RetryPolicy(
        max_retries: int,  # 0 = no retries
        retry_on: List["exit_code" | "timeout" | "oom" | "network" | "all"],
        backoff_multiplier: float,  # 2.0 = double wait time
        initial_delay_seconds: float
    ),
    resource_requirements: ResourceRequirements(
        min_cpu_cores: Optional[float],
        min_memory_gb: Optional[float],
        min_disk_gb: Optional[float],
        gpu_required: bool,
        gpu_count: Optional[int]
    ),
    expected_outputs: List[OutputSpec](
        uri: str,  # Output file URI (s3:// or local)
        format: Optional[str],  # e.g., "fastq", "bam", "vcf"
        required: bool,  # Whether required for success
        size_estimate_mb: Optional[float]
    ),
    confidence_score: float,  # 0-1, confidence in execution plan
    reasoning: str,  # Why this execution plan was chosen
    warnings: List[str],  # Warnings or caveats
    estimated_runtime_minutes: Optional[float],
    request_id: Optional[str],  # For tracing
    workflow_plan_hash: Optional[str],  # For reproducibility
    infra_decision_hash: Optional[str]  # For reproducibility
)
```

📖 **See:** `backend/contracts/` for complete contract definitions

---

## Tool Categories

### Pre-built Tools
- Sequence alignment (MAFFT, Clustal, MUSCLE)
- Phylogenetic tree inference (IQ-TREE, RAxML, FastTree)
- Quality control (FastQC)
- Differential expression (DESeq2, edgeR)
- Single-cell analysis (Scanpy)
- NCBI/UniProt/GO database queries

### Dynamically Generated Tools
- Read merging (BBMerge)
- Custom bioinformatics operations
- Dataset-specific analyses

📖 **See:** [agents/tool-generator-agent.md](../agents/tool-generator-agent.md)

---

## Testing Infrastructure

Helix.AI uses a **three-tier testing strategy** to ensure reliability at component, integration, and system levels.

### 1. Unit Tests (`tests/unit/backend/`)

**Purpose:** Validate individual agent logic and contracts in isolation

**Coverage:**
- Agent decision-making (Intent Detector, Workflow Planner, Infrastructure Agent)
- Contract validation (WorkflowPlan, InfraDecision, ExecutionToolSpec)
- Tool routing and execution policies
- Session management and history tracking
- Handoff policy enforcement

**Approach:** Fully mocked—no external dependencies (LLMs, AWS, S3)

**Example:**
```bash
# Run all unit tests
pytest tests/unit/backend/ -v

# Test specific agent
pytest tests/unit/backend/test_workflow_planner_agent.py -v

# Current status: 290/305 passing (94.8%)
```

**Key Test Files:**
- `test_workflow_planner_agent.py` - Workflow playbook matching and planning (27 tests)
- `test_intent_classifier.py` - Intent detection and confidence scoring
- `test_execution_broker_policy.py` - Routing decisions (sync vs. async)
- `test_infrastructure_based_async_routing.py` - Infrastructure-aware routing
- `test_agent_responsibilities_policy.py` - Agent ownership validation
- `test_handoff_policy.py` - Legal agent transitions

---

### 2. End-to-End Tests with Mocks (`tests/demo_scenarios/`)

**Purpose:** Validate multi-agent orchestration without real computation

**Coverage:**
- Complete request-to-response pipeline
- Agent handoffs and contract exchanges
- Tool mapping from natural language
- Error handling and edge cases
- Session context propagation

**Approach:** Mocked execution (LLMs may be real or mocked, but no actual tool execution)

**Scenarios:**
| Scenario | Test Focus |
|----------|-----------|
| `ask_basic_question.yaml` | QA intent routing to Bioinformatics Guru |
| `execute_fastqc_basic.yaml` | FastQC workflow planning + infra decision |
| `execute_single_cell_analysis.yaml` | Complex scRNA-seq workflow decomposition |
| `execute_with_codegen.yaml` | Tool Generator invocation for missing tools |
| `edge_low_confidence_intent.yaml` | Uncertainty handling and clarification |

**Running:**
```bash
# Run all demo scenarios
pytest tests/demo_scenarios/test_real_execution.py -v

# Run specific scenario
python tests/demo_scenarios/demo_cli.py scenarios/execute_fastqc_basic.yaml

# Compare against baselines
pytest tests/demo_scenarios/ --compare-baselines
```

**Output Validation:**
- ✅ Correct intent classification
- ✅ Appropriate workflow plan generated
- ✅ Valid infrastructure recommendation
- ✅ Proper agent sequence (no illegal handoffs)
- ✅ Tool parameters extracted correctly

---

### 3. End-to-End Tests with Real Execution (`tests/workflows/`)

**Purpose:** Validate actual workflow execution on real infrastructure

**Coverage:**
- Complete RNA-seq, scRNA-seq, phylogenetics workflows
- Real data processing on AWS (Local, EC2, EMR)
- S3 input/output management
- Job submission, monitoring, and result retrieval
- Compute cost and performance validation

**Approach:** Real execution—uses actual AWS infrastructure and computational tools

**Test Categories:**

**🚀 Fast Tests (`-m fast`):**
- Small datasets (< 100MB)
- Local execution only
- < 5 minutes per test
- Safe for CI/CD

**🐌 Slow Tests (`-m slow`):**
- Large datasets (> 10GB)
- EMR cluster execution
- 10-30 minutes per test
- Manual or nightly builds

**Running:**
```bash
# Run fast E2E tests
pytest tests/workflows/test_e2e_workflow_execution.py -v -m fast

# Run specific workflow
pytest tests/workflows/test_rnaseq_workflow.py::test_rnaseq_complete_pipeline -v -s

# Skip slow tests
pytest tests/workflows/ -v -m "not slow"

# Run all (WARNING: 30+ minutes, incurs AWS costs)
pytest tests/workflows/ -v -s
```

**Prerequisites:**
- AWS credentials configured (`aws sts get-caller-identity`)
- S3 bucket access (read: input data, write: results)
- EMR cluster permissions (for slow tests)

**Validation:**
- ✅ Job completes successfully (exit code 0)
- ✅ Output files exist in S3 at expected paths
- ✅ File types match expectations (FASTQ, BAM, HTML, etc.)
- ✅ File sizes within expected ranges
- ✅ Provenance captured (commands, tool versions, timestamps)

**Example E2E Test:**
```python
def test_fastqc_small_dataset():
    """Run FastQC on small test dataset locally."""
    response = submit_job(
        command="Run FastQC quality control",
        files=["s3://bucket/test_R1.fq", "s3://bucket/test_R2.fq"]
    )
    
    # Wait for completion (max 5 minutes)
    result = wait_for_job(response["job_id"], timeout=300)
    
    # Verify outputs
    assert result["status"] == "completed"
    assert len(result["output_files"]) >= 2  # HTML + ZIP per file
    assert all(f.endswith((".html", ".zip")) for f in result["output_files"])
```

---

### Test Execution Summary

| Test Type | Count | Duration | Mocked | Real Compute | Use Case |
|-----------|-------|----------|--------|--------------|----------|
| **Unit** | 305 | < 1 min | ✅ Full | ❌ None | Agent logic validation |
| **E2E (Mock)** | 14 scenarios | 2-5 min | ⚠️ Partial | ❌ None | Pipeline integration |
| **E2E (Real-Fast)** | 8 tests | 5-10 min | ❌ None | ✅ Local | Quick validation |
| **E2E (Real-Slow)** | 5 tests | 30-60 min | ❌ None | ✅ EMR | Full system test |

---

### Evaluation Tests (`tests/evals/`)

**Purpose:** Track regression in agent decision quality over time

**Coverage:**
- Intent classification accuracy (recall, precision, F1)
- Tool mapping correctness for known scenarios
- Infrastructure routing decisions vs. ground truth

**Running:**
```bash
# Run evaluation tests
pytest tests/evals/ -v

# Generate evaluation report
pytest tests/evals/ --generate-report
```

📖 **See:** `tests/README.md`, `tests/workflows/README_E2E_TESTS.md`, and `tests/demo_scenarios/GETTING_REAL_RESULTS.md`

---

## Development Workflow

### Adding a New Workflow Playbook

1. Define playbook in `backend/workflow_planner_agent.py`
2. Add playbook class with `matches()` and `create_workflow_plan()` methods
3. Register in `PLAYBOOKS` list
4. Write unit tests in `tests/unit/backend/test_workflow_planner_agent.py`
5. Update [agents/workflow-planner-agent.md](../agents/workflow-planner-agent.md) documentation

### Adding a New Agent

1. Define responsibility in [agents/agent-responsibilities.md](../agents/agent-responsibilities.md)
2. Create agent spec in `agents/<agent-name>-agent.md`
3. Implement in `backend/<agent_name>.py`
4. Define contracts in `backend/contracts/`
5. Write unit tests
6. Update [agents/handoff-policy.md](../agents/handoff-policy.md) with new transitions

---

## Related Documentation

### System Architecture
- [docs/TASK_ARCHITECTURE.md](./TASK_ARCHITECTURE.md) - Micro-/Macroflow pattern
- [docs/SESSION_MANAGEMENT.md](./SESSION_MANAGEMENT.md) - Session context and continuity
- [docs/OUTPUT_SCHEMA.md](./OUTPUT_SCHEMA.md) - JSON response format

### Agent Specifications
- [agents/workflow-planner-agent.md](../agents/workflow-planner-agent.md) - Workflow planning
- [agents/infrastructure-decision-agent.md](../agents/infrastructure-decision-agent.md) - Infrastructure selection
- [agents/implementation-agent.md](../agents/implementation-agent.md) - Execution packaging
- [agents/tool-generator-agent.md](../agents/tool-generator-agent.md) - Dynamic tool generation
- [agents/intent-detector-agent.md](../agents/intent-detector-agent.md) - Intent classification

### Policies and Guidelines
- [agents/agent-responsibilities.md](../agents/agent-responsibilities.md) - Agent ownership matrix
- [agents/handoff-policy.md](../agents/handoff-policy.md) - Routing rules and transitions
- [docs/STATISTICAL_GUIDELINES.md](./STATISTICAL_GUIDELINES.md) - Statistical validity
- [docs/SAFETY_POLICY.md](./SAFETY_POLICY.md) - Privacy and dual-use prevention

### Implementation Details
- `backend/contracts/` - Contract definitions (WorkflowPlan, InfraDecision, etc.)
- `backend/orchestrator.py` - Multi-agent coordination
- `backend/execution_broker.py` - Execution routing and job management
- `backend/plan_ir.py` - Plan intermediate representation

---

## Quick Start

**For Users:**
1. Submit bioinformatics request (e.g., "Analyze these RNA-seq files")
2. System classifies intent and plans workflow
3. Infrastructure is selected automatically
4. Jobs execute on appropriate platform
5. Results returned with visualizations and provenance

**For Developers:**
1. Read [agents/agent-responsibilities.md](../agents/agent-responsibilities.md)
2. Review contract definitions in `backend/contracts/`
3. Study [agents/handoff-policy.md](../agents/handoff-policy.md)
4. Run tests: `pytest tests/unit/backend/ -v`
5. See `tests/demo_scenarios/` for example workflows

---

## Version History

- **v2.0** (2026-01-18): Multi-agent architecture with Workflow Planner
- **v1.5** (2026-01): Infrastructure decision agent integration
- **v1.0** (2025): Initial release with tool-based execution

---

## Contact and Support

- **Issues:** GitHub Issues
- **Documentation:** `docs/` and `agents/`
- **Tests:** `tests/`
- **Examples:** `tests/demo_scenarios/`

---

**For detailed information on any component, see the linked documentation above.**
