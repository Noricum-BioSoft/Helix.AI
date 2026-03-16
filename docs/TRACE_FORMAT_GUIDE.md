# Orchestrator Trace Format Guide

**Date:** 2026-01-18  
**Related:** P1-1 Provenance Hashing Fix  
**Sample:** `docs/SAMPLE_TRACE.json`

---

## Overview

The `OrchestratorTrace` captures the complete execution history of a multi-agent pipeline, including:
- Request metadata (ID, command, timing)
- Agent invocations (inputs, outputs, hashes, performance)
- Final outputs (InfraDecision, ExecutionToolSpec)

This trace enables:
- **Provenance tracking** (what changed between runs)
- **Performance monitoring** (agent latency, bottlenecks)
- **Debugging** (error tracking, confidence scores)
- **Caching/deduplication** (detect repeated work via hashes)
- **Baseline regression testing** (compare traces across versions)

---

## Trace Structure

### Top-Level Fields

```json
{
  "request_id": "demo_small_local_20260118",
  "user_command": "Run FastQC quality control on local FASTQ files",
  "start_time": 1737244800.123456,
  "end_time": 1737244805.789012,
  "duration_ms": 5665.556,
  "success": true,
  "error": null,
  "timestamp": "2026-01-18T12:00:00.123456",
  "invocations": [...],
  "infra_decision": {...},
  "execution_spec": {...}
}
```

| Field | Type | Description |
|-------|------|-------------|
| `request_id` | string | Unique identifier for this request |
| `user_command` | string | Original user command/prompt |
| `start_time` | float | Unix timestamp when pipeline started |
| `end_time` | float | Unix timestamp when pipeline completed |
| `duration_ms` | float | Total pipeline duration in milliseconds |
| `success` | boolean | Whether the pipeline completed successfully |
| `error` | string\|null | Error message if pipeline failed |
| `timestamp` | string | ISO 8601 timestamp for start_time |
| `invocations` | array | List of agent invocations (see below) |
| `infra_decision` | object | Final infrastructure decision contract |
| `execution_spec` | object | Final execution specification contract |

---

## Agent Invocation Structure

Each agent invocation in the `invocations` array has this structure:

```json
{
  "agent_name": "InfrastructureDecisionAgent",
  "request_id": "demo_small_local_20260118",
  "start_time": 1737244800.123456,
  "end_time": 1737244802.456789,
  "duration_ms": 2333.333,
  "success": true,
  "error": null,
  "input_hash": "a1b2c3d4e5f67890...",
  "output_hash": "9876543210fedcba...",
  "confidence_score": 0.92,
  "warnings": [],
  "timestamp": "2026-01-18T12:00:00.123456"
}
```

### Fields

| Field | Type | Description |
|-------|------|-------------|
| `agent_name` | string | Canonical agent name (from `AgentName` enum) |
| `request_id` | string | Request ID (same as parent trace) |
| `start_time` | float | Unix timestamp when agent started |
| `end_time` | float | Unix timestamp when agent completed |
| `duration_ms` | float | Agent execution time in milliseconds |
| `success` | boolean | Whether the agent succeeded |
| `error` | string\|null | Error message if agent failed |
| `input_hash` | string | SHA256 hash of input contract (64 chars) |
| `output_hash` | string | SHA256 hash of output contract (64 chars) |
| `confidence_score` | float | Agent's confidence in its output (0.0-1.0) |
| `warnings` | array | List of warning messages from agent |
| `timestamp` | string | ISO 8601 timestamp for start_time |

### Key Features

**1. Provenance Hashing (P1-1 Fix)**
- `input_hash` and `output_hash` are SHA256 hashes of contract JSON
- Enables fast comparison: "Did the contract change?"
- Use case: Baseline regression testing

**2. Performance Metrics**
- `duration_ms` tracks agent latency
- Helps identify bottlenecks in multi-agent pipelines

**3. Confidence Tracking**
- `confidence_score` (0.0-1.0) indicates agent certainty
- Low confidence (<0.5) should trigger user interaction
- Very low confidence (<0.3) should block execution

**4. Warning Collection**
- `warnings` array captures non-fatal issues
- Example: "Container image version not specified"

---

## Hash Computation

Hashes are computed using `compute_contract_hash()`:

```python
import hashlib

def compute_contract_hash(contract_json: str) -> str:
    """Compute SHA256 hash of contract JSON."""
    hash_obj = hashlib.sha256(contract_json.encode('utf-8'))
    return hash_obj.hexdigest()  # 64-character hex string
```

### Example

**Input Contract:**
```json
{
  "description": "Quality control analysis...",
  "operations": [...],
  "data_inputs": [...]
}
```

**Hash:**
```
a1b2c3d4e5f67890abcdef1234567890abcdef1234567890abcdef1234567890
```

### Properties

- **Deterministic:** Same input always produces same hash
- **Stable:** Pydantic `model_dump_json()` uses sorted keys
- **Fast:** SHA256 computation is O(n) in JSON length
- **Collision-resistant:** SHA256 is cryptographically secure

---

## Use Cases

### 1. Deduplication

**Problem:** Avoid re-running identical requests

**Solution:** Check `input_hash` before invoking agent

```python
# Check if we've seen this input before
if input_hash in cache:
    return cache[input_hash]  # Skip agent call!

# Otherwise, invoke agent and cache result
result = await agent.run(input)
cache[input_hash] = result
```

**Benefit:** Save API costs, reduce latency

---

### 2. Provenance Tracking

**Problem:** Detect unexpected behavior changes

**Solution:** Compare `output_hash` against baseline

```python
# Compare current trace to baseline
baseline_trace = load_baseline("execute_fastqc_basic.json")
current_trace = await orchestrator.run_pipeline(...)

if not current_trace.has_same_contracts_as(baseline_trace):
    diff = current_trace.get_contract_diff_summary(baseline_trace)
    print(f"Contract changed: {diff}")
    # Example: {"InfrastructureDecisionAgent": "output_changed"}
```

**Benefit:** Catch regressions, audit changes

---

### 3. Baseline Regression Testing

**Problem:** Ensure agents produce consistent outputs

**Solution:** Store trace as baseline, compare on each test run

```python
# Baseline test
def test_fastqc_baseline():
    trace = await orchestrator.run_pipeline(command, workflow_plan)
    baseline = load_baseline("fastqc_baseline.json")
    
    # Fast comparison using hashes
    assert trace.has_same_contracts_as(baseline), \
        f"Contract changed: {trace.get_contract_diff_summary(baseline)}"
```

**Benefit:** Detect unintended behavior changes

---

### 4. Performance Monitoring

**Problem:** Identify slow agents

**Solution:** Analyze `duration_ms` per agent

```python
# Analyze agent performance
for inv in trace.invocations:
    if inv.duration_ms > 5000:  # >5 seconds
        print(f"Slow agent: {inv.agent_name} took {inv.duration_ms}ms")
```

**Benefit:** Optimize pipeline latency

---

### 5. Confidence Monitoring

**Problem:** Detect low-quality outputs

**Solution:** Check `confidence_score` and `warnings`

```python
# Check for low confidence
for inv in trace.invocations:
    if inv.confidence_score < 0.5:
        print(f"Low confidence: {inv.agent_name} = {inv.confidence_score}")
        # Trigger user interaction
    
    if len(inv.warnings) >= 3:
        print(f"Many warnings: {inv.agent_name} has {len(inv.warnings)} warnings")
        # Surface to user
```

**Benefit:** Improve user experience, prevent errors

---

## Complete Example

See `docs/SAMPLE_TRACE.json` for a full example trace including:

1. **Request metadata:**
   - ID: `demo_small_local_20260118`
   - Command: "Run FastQC quality control on local FASTQ files"
   - Duration: 5.67 seconds

2. **Two agent invocations:**
   - **InfrastructureDecisionAgent** (2.33s, confidence=0.92)
     - Input hash: `a1b2c3d4...`
     - Output hash: `98765432...`
     - Decision: Local execution (small files)
   
   - **ImplementationAgent** (3.22s, confidence=0.88)
     - Input hash: `f0e1d2c3...`
     - Output hash: `12345678...`
     - Spec: FastQC with biocontainers image
     - Warning: "Container image version not specified"

3. **Final outputs:**
   - **InfraDecision:** Local, cost $0.00-$0.01, confidence 0.92
   - **ExecutionToolSpec:** FastQC command, Docker container, 2 cores, 2GB RAM

---

## Trace File Organization

### Recommended Structure

```
traces/
├── baselines/              # Baseline traces for regression testing
│   ├── execute_fastqc_basic.json
│   ├── execute_read_merging_emr.json
│   └── ...
├── production/             # Production request traces
│   ├── 2026-01-18/
│   │   ├── req_abc123.json
│   │   ├── req_def456.json
│   │   └── ...
│   └── ...
└── experiments/            # Experimental/debug traces
    ├── agent_comparison_v1_v2.json
    └── ...
```

---

## Trace Utilities

The `OrchestratorTrace` class provides utility methods:

### `get_invocation_hashes()`

**Returns:** Dict mapping agent_name → (input_hash, output_hash)

```python
hashes = trace.get_invocation_hashes()
# {
#   "InfrastructureDecisionAgent": ("a1b2c3d4...", "98765432..."),
#   "ImplementationAgent": ("f0e1d2c3...", "12345678...")
# }
```

---

### `has_same_contracts_as(other_trace)`

**Returns:** True if all agents have matching hashes

```python
baseline = load_baseline("fastqc_baseline.json")
current = await orchestrator.run_pipeline(...)

if current.has_same_contracts_as(baseline):
    print("✅ Behavior unchanged")
else:
    print("⚠️ Behavior changed!")
```

**Complexity:** O(n) where n = number of agents (fast!)

---

### `get_contract_diff_summary(other_trace)`

**Returns:** Dict mapping agent_name → diff status

**Possible statuses:**
- `"added"` - Agent in current but not baseline
- `"removed"` - Agent in baseline but not current
- `"unchanged"` - Hashes match
- `"input_changed"` - Input hash differs
- `"output_changed"` - Output hash differs
- `"both_changed"` - Both hashes differ

```python
diff = current.get_contract_diff_summary(baseline)
# {
#   "InfrastructureDecisionAgent": "unchanged",
#   "ImplementationAgent": "output_changed"
# }
```

---

## Integration with Testing

### Scenario Baseline

Each scenario can have a baseline trace:

```yaml
# scenarios/execute_fastqc_basic.yaml
name: execute_fastqc_basic
command: "Run FastQC on small local files"
baseline_trace: baselines/execute_fastqc_basic.json
```

### Test

```python
def test_scenario_matches_baseline():
    scenario = load_scenario("execute_fastqc_basic.yaml")
    baseline = load_trace(scenario["baseline_trace"])
    
    trace = await orchestrator.run_pipeline(
        command=scenario["command"],
        workflow_plan=scenario["workflow_plan"]
    )
    
    # Fast hash-based comparison
    assert trace.has_same_contracts_as(baseline), \
        f"Behavior changed: {trace.get_contract_diff_summary(baseline)}"
```

---

## Logging Format

When agents complete, logs include truncated hashes:

```
INFO [demo_small_local_20260118] InfrastructureDecisionAgent completed: 
  Local (confidence=0.92) [input_hash=a1b2c3d4..., output_hash=98765432...]

INFO [demo_small_local_20260118] ImplementationAgent completed: 
  fastqc on Local (confidence=0.88) [input_hash=f0e1d2c3..., output_hash=12345678...]
```

**Format:** First 8 characters of hash (for brevity)

---

## Future Enhancements

### 1. Trace Aggregation

**Problem:** Analyze patterns across many traces

**Solution:** Aggregate traces into analytics database

```sql
SELECT agent_name, AVG(duration_ms), AVG(confidence_score)
FROM agent_invocations
WHERE DATE(timestamp) = '2026-01-18'
GROUP BY agent_name;
```

---

### 2. Hash-Based Caching Layer

**Problem:** Repeated identical requests waste resources

**Solution:** Build cache keyed by input_hash

```python
class HashCache:
    def get(self, input_hash: str) -> Optional[Contract]:
        return redis.get(f"cache:{input_hash}")
    
    def set(self, input_hash: str, output: Contract):
        redis.setex(f"cache:{input_hash}", ttl=3600, value=output.json())
```

---

### 3. Automatic Baseline Updates

**Problem:** Manually updating baselines is tedious

**Solution:** Auto-update if all tests pass

```python
if all_tests_passed and no_confidence_warnings:
    # Update baseline
    trace.to_json("baselines/execute_fastqc_basic.json")
    git_commit("Update baseline: all tests passing")
```

---

## Summary

The `OrchestratorTrace` format provides:

✅ **Provenance hashing** (P1-1 fix) - Track what changed  
✅ **Performance monitoring** - Identify bottlenecks  
✅ **Confidence tracking** - Detect low-quality outputs  
✅ **Baseline testing** - Fast hash-based comparison  
✅ **Debugging** - Complete audit trail of execution  

**Key Files:**
- `backend/orchestrator.py` - Trace generation
- `docs/SAMPLE_TRACE.json` - Example trace
- `docs/PROVENANCE_HASHING_FIX.md` - P1-1 implementation details

🎉 **Complete observability for multi-agent pipelines!**
