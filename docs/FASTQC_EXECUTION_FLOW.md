# FastQC Execution Flow: A Multi-Agent Story

This document explains the multi-agent execution flow for a FastQC quality analysis request in the Helix.AI bioinformatics platform.

## Example Request

```
Run FastQC quality analysis on paired-end reads:
- s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq
- s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq
```

## The Journey of a FastQC Request

### Act 1: Understanding Intent (IntentDetector Agent)

**When**: 08:26:57  
**Duration**: ~5 seconds  
**Agent**: `IntentDetector`

The journey begins when the system receives the FastQC command with two S3 file paths.

The **IntentDetector Agent** is called first. Its job is to understand *what kind of request* this is. It examines the command and classifies it as an **"execute"** intent with labels `['action', 'data']`.

**Classification Result**:
- Intent: `execute`
- Labels: `['action', 'data']`
- Reason: `llm_classified_action_action,data`

**Why this agent?** Before doing anything, the system needs to understand if you're asking a question, requesting an action, or doing something else entirely. This prevents the system from executing commands when you're just asking about them.

---

### Act 2: Handoff to the Specialist (BioinformaticsExecutor)

**When**: 08:27:02  
**Agent**: `BioinformaticsExecutor` (Planner)

Once the intent is clear, the **HandoffPolicy** validates and approves a handoff:

```
IntentDetector → BioinformaticsExecutor
```

This creates an agent sequence:
```python
['IntentDetector', 'BioinformaticsExecutor']
```

The **BioinformaticsExecutor** takes over with a clear mission: *identify which specific tool should handle this FastQC request*.

**Strategy**: `stream → invoke → router → tool-generator`

---

### Act 3: Tool Identification (Streaming Strategy)

**When**: 08:27:03  
**Duration**: 11.16 seconds  
**LLM**: DeepSeek API

The BioinformaticsExecutor uses a clever streaming approach to quickly identify the right tool without waiting for full execution.

**Process**:
1. Starts agent.stream() with stream_mode='updates'
2. Streams events from the LLM (DeepSeek API)
3. Captures tool call in stream event #1 (0.6 seconds)
4. Extracts tool: `fastqc_quality_analysis`
5. Extracts parameters:
   ```python
   {
       'input_r1': 's3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq',
       'input_r2': 's3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq'
   }
   ```
6. Stops early—no need to wait for full execution

**Result**: Tool mapped in 11.16 seconds total.

**Why streaming?** This approach captures the tool decision immediately without waiting for the full graph execution, significantly improving response time.

---

### Act 4: Infrastructure Decision (Infrastructure Decision Agent)

**When**: 08:27:09  
**Duration**: ~22 seconds  
**Agent**: `InfrastructureDecisionAgent`

Now comes a critical decision point. Before executing, the system needs to decide: *Should this run locally or in the cloud (EMR)?*

The **Infrastructure Decision Agent** analyzes:
- File sizes (small test files)
- Computational requirements
- Cost-effectiveness
- Available resources
- Tool characteristics

**Decision**: Run **LOCAL** with 95% confidence

**Reasoning**:
- Files are small test files (test_mate_R1.fq and test_mate_R2.fq)
- Local execution in Docker is faster than spinning up cloud infrastructure
- More cost-effective for small-scale analysis
- EMR would add unnecessary overhead for this workload

---

### Act 5: Execution (Sandbox Executor)

**When**: 08:27:31  
**Duration**: 29.50 seconds  
**Executor**: Docker Sandbox

With the infrastructure decision made, execution begins:

#### Step 1: Download (2.57 seconds)
```
08:27:31 - 📥 Downloading test_mate_R1.fq from S3...
08:27:32 - 📥 Downloading test_mate_R2.fq from S3...
```

#### Step 2: Docker Sandbox Execution (1.39 seconds)
```
08:27:34 - 🐳 Running fastqc in Docker sandbox...
08:27:36 - ✅ fastqc completed in 1.39s
```

**Why Docker?** Provides:
- Isolated execution environment
- Consistent tool versions
- Security (sandboxing)
- No impact on host system

#### Step 3: Upload Results (1.63 seconds)
```
08:27:36 - 📤 Uploading results to S3...
08:27:38 - ✅ Results uploaded to s3://noricum-ngs-data/.../fastqc-results/
```

**Total execution time**: 29.50 seconds

---

### Epilogue: Response Building

**When**: 08:27:38  
**Duration**: < 1 millisecond

Finally, the system:
1. Saves the session history (0.73ms)
2. Adds history entry for `fastqc_quality_analysis`
3. Builds a standardized response with all metadata
4. Returns success to the user

**Total end-to-end time**: 40.66 seconds

---

## The Agent Chain

```
User Request
    ↓
IntentDetector (classify intent)
    ↓ "execute" with data
BioinformaticsExecutor (identify tool)
    ↓ "fastqc_quality_analysis"
Infrastructure Decision Agent (local vs cloud)
    ↓ "LOCAL" (95% confidence)
Sandbox Executor (Docker execution)
    ↓ Results to S3
Response Builder
    ↓
User receives results
```

---

## Why This Multi-Agent Approach?

### 1. IntentDetector
**Purpose**: Ensures the system understands *what you want* before doing anything

**Benefits**:
- Prevents accidental execution of commands
- Routes questions to knowledge base instead of execution
- Enables appropriate response for different intent types

### 2. BioinformaticsExecutor
**Purpose**: Maps natural language to specific bioinformatics tools

**Benefits**:
- Users don't need to know exact tool names
- Handles parameter extraction automatically
- Supports multiple ways of describing the same operation

### 3. Infrastructure Decision Agent
**Purpose**: Optimizes cost and performance by choosing the right execution environment

**Benefits**:
- Automatic cost optimization
- Performance optimization (local for small, EMR for large)
- Resource-aware decision making
- Confidence scoring for transparency

### 4. Sandbox Executor
**Purpose**: Provides safe, isolated execution

**Benefits**:
- Security through isolation
- Consistent execution environment
- No impact on host system
- Easy to scale and replicate

---

## Performance Breakdown

| Phase | Agent/Component | Duration | % of Total |
|-------|----------------|----------|------------|
| Intent Detection | IntentDetector | 5.3s | 13% |
| Tool Mapping | BioinformaticsExecutor | 11.2s | 28% |
| Infrastructure Decision | Infrastructure Decision Agent | 22.0s | 54% |
| Execution | Sandbox Executor | 6.6s | 16% |
| File Download | S3 Download | 2.6s | 6% |
| FastQC Processing | Docker/FastQC | 1.4s | 3% |
| File Upload | S3 Upload | 1.6s | 4% |
| Response Building | Response Builder | <0.1s | <1% |
| **Total** | **End-to-End** | **40.7s** | **100%** |

---

## HandoffPolicy Enforcement

The system uses a **HandoffPolicy** to ensure agents collaborate in the correct order with proper validation:

```python
Agent Sequence: ['IntentDetector', 'BioinformaticsExecutor']
```

**Validation Points**:
1. `IntentDetector` is registered first
2. Handoff validated: `IntentDetector → BioinformaticsExecutor`
3. `BioinformaticsExecutor` is registered
4. Sequence is tracked throughout execution

This modular design makes the system:
- **Flexible**: Easy to add new agents or modify existing ones
- **Maintainable**: Each agent has a single, clear responsibility
- **Intelligent**: Resource allocation and routing based on actual analysis
- **Auditable**: Complete trace of agent decisions and handoffs

---

## Example Log Trace

For reference, here's the actual log trace showing the agent sequence:

```
2026-01-23 08:26:57 - [HandoffPolicy] Registered agent call: IntentDetector
2026-01-23 08:26:57 - [HandoffPolicy] Current sequence: ['IntentDetector']
2026-01-23 08:27:02 - Intent classified: execute (labels: ['action', 'data'])
2026-01-23 08:27:02 - [HandoffPolicy] Validated handoff: IntentDetector → BioinformaticsExecutor
2026-01-23 08:27:02 - [HandoffPolicy] Registered agent call: BioinformaticsExecutor
2026-01-23 08:27:02 - [HandoffPolicy] Current sequence: ['IntentDetector', 'BioinformaticsExecutor']
2026-01-23 08:27:03 - ✅ Tool call found: fastqc_quality_analysis
2026-01-23 08:27:31 - ✅ Infrastructure Decision: Local (confidence: 0.95)
2026-01-23 08:27:38 - ✅ FastQC execution completed
```

---

## Session Information

- **Session ID**: `c3ee6a71-1474-4bb6-9ac7-54882cb91580`
- **Agent Mode**: `TOOL MAPPING ONLY` (no execution in agent graph)
- **Execution Mode**: `LOCAL` (Docker sandbox)
- **Total Time**: 40.66 seconds
- **Status**: ✅ Success
