# Orchestration Duality - P1 Architectural Issue

**Date:** 2026-01-18  
**Issue:** Two parallel orchestration systems exist, unclear which is authoritative  
**Priority:** P1 (High)  
**Status:** 🔴 Documented, Needs Resolution

---

## Problem Statement

Helix.AI has **TWO PARALLEL ORCHESTRATION SYSTEMS** operating simultaneously, creating architectural confusion and maintenance burden:

1. **OLD System** (`backend/agent.py`) - 6+ agent multi-agent graph
2. **NEW System** (`backend/orchestrator.py`) - 2-agent clean pipeline

**Impact:**
- Unclear which system is authoritative
- Different entry points use different orchestrators
- Tests reference agent names from OLD system
- Documentation describes both without clear distinction
- Difficult to reason about system behavior
- Maintenance burden (bug fixes must target both systems)

---

## System 1: OLD Multi-Agent Orchestrator (PRODUCTION)

### Location
**`backend/agent.py`** (1,798 lines)

### Architecture
Multi-agent graph with 6+ specialized agents:

```
User Command
    ↓
IntentDetector ("ask" or "execute"?)
    ↓
    ├─→ [ask] → BioinformaticsGuru (answers questions)
    │
    └─→ [execute] → BioinformaticsExecutor (Planner)
                        ↓
                    InfrastructureExpert (selects environment)
                        ↓
                    (optional) CodeGenerator (tool gaps)
                        ↓
                    ExecutionBroker (executes jobs)
                        ↓
                    DataVisualizer (plots/reports)
```

### Key Components
```python
# Entry point
from backend.agent import handle_command

result = await handle_command(
    command="Run FastQC",
    session_id="user123",
    session_context={...}
)
```

### Agent Names
- `IntentDetector`
- `BioinformaticsGuru` (QA agent)
- `BioinformaticsExecutor` (Planner)
- `InfrastructureExpert` (Infrastructure selector)
- `CodeGenerator` (Tool generator)
- `ExecutionBroker` (Job executor - **not actually an agent!**)
- `DataVisualizer` (Visualization)

### Policy Enforcement
- `HandoffPolicy` class validates agent transitions
- `check_plan_not_mutated()` ensures InfrastructureExpert doesn't change plan
- Staged functions with integrated policy checks

### Used By (PRODUCTION)
✅ **HTTP API** (`backend/main_with_mcp.py`):
- `/execute` endpoint → `handle_command()`
- `/chat` endpoint → `handle_command()`
- `/mcp/call_tool` with `bioinformatics_agent` → `handle_command()`

✅ **Scripts:**
- `tests/demo_scenarios/run_fastqc_full_execution.py` → HTTP API → `handle_command()`

✅ **Tests:**
- `tests/demo_scenarios/test_scenarios.py` (asserts agent names like "BioinformaticsExecutor")
- Most existing tests

### Status
🟢 **PRODUCTION** - Currently serving all HTTP requests

---

## System 2: NEW Clean Orchestrator (PROTOTYPE)

### Location
**`backend/orchestrator.py`** (443 lines)

### Architecture
Clean 2-agent pipeline with Pydantic contracts:

```
User Command + WorkflowPlan
    ↓
InfrastructureDecisionAgent (WHERE to execute)
    │
    │ Input: WorkflowPlan (with size_bytes, location_type)
    │ Output: InfraDecision (with confidence, cost ranges)
    │
    ↓
ImplementationAgent (HOW to execute)
    │
    │ Input: WorkflowPlan + InfraDecision
    │ Output: ExecutionToolSpec (with commands, retry policy)
    │
    ↓
(External runner executes - NOT part of orchestrator)
```

### Key Components
```python
# Entry point
from backend.orchestrator import Orchestrator

orchestrator = Orchestrator()
trace = await orchestrator.run_pipeline(
    command="Run FastQC",
    workflow_plan=workflow_plan  # Pydantic model
)
```

### Agent Names (Canonical)
- `InfrastructureDecisionAgent`
- `ImplementationAgent`

### Features
✅ **Strict Pydantic contracts** (InfraDecision, ExecutionToolSpec)  
✅ **Provenance hashing** (SHA256 for all inputs/outputs)  
✅ **Rich metadata flow** (size_bytes, location_type)  
✅ **Trace utilities** (hash-based comparison, diff summaries)  
✅ **Clean separation** (planning vs execution)  
✅ **Observable** (structured logging with hashes)

### Used By (PROTOTYPE ONLY)
⚠️ **CLI Demo** (`backend/cli_demo.py`):
- Development/testing only, not production

⚠️ **Unit Tests:**
- `tests/unit/backend/test_phase3_agents.py`
- `tests/unit/backend/test_phase5_snapshots.py`

### Status
🟡 **PROTOTYPE** - Not used in production HTTP API

---

## Comparison Matrix

| Aspect | OLD System (agent.py) | NEW System (orchestrator.py) |
|--------|----------------------|------------------------------|
| **Status** | 🟢 Production | 🟡 Prototype |
| **Entry Point** | `handle_command()` | `Orchestrator.run_pipeline()` |
| **HTTP API** | ✅ Used | ❌ Not used |
| **Agent Count** | 6+ agents | 2 agents |
| **Architecture** | Multi-agent graph | Linear pipeline |
| **Contracts** | Mixed (dict + dataclass) | Strict Pydantic |
| **Provenance** | ❌ No hashing | ✅ SHA256 hashing |
| **Metadata Flow** | ⚠️ Partial | ✅ Complete |
| **Policy Checks** | ✅ Extensive | ❌ Minimal |
| **Agent Names** | Legacy (BioinformaticsExecutor) | Canonical (ImplementationAgent) |
| **Lines of Code** | 1,798 | 443 |
| **Complexity** | High | Low |
| **Tests** | Extensive | Limited |

---

## Root Cause of Confusion

### 1. Entry Point Ambiguity

**Question:** "Where does a request go?"

**Answer depends on entry point:**

```python
# Via HTTP API (/execute) → OLD system
POST /execute
└─> handle_command() from backend.agent
    └─> IntentDetector → BioinformaticsExecutor → InfrastructureExpert → ...

# Via CLI Demo → NEW system
python backend/cli_demo.py
└─> Orchestrator.run_pipeline()
    └─> InfrastructureDecisionAgent → ImplementationAgent

# Via Tests → Both!
test_scenarios.py → OLD system (asserts "BioinformaticsExecutor")
test_phase3_agents.py → NEW system (uses Orchestrator)
```

### 2. Documentation Ambiguity

**`docs/SYSTEM_OVERVIEW.md` describes both systems without clear labels:**

```markdown
## Core Architecture

[Shows 6+ agent graph: IntentDetector → Guru → Planner → Infra → Broker → Visualizer]

## Key Contracts

### InfraDecision
[Shows NEW system Pydantic contract]

### ExecutionToolSpec
[Shows NEW system Pydantic contract]
```

**Problem:** Diagram shows OLD system, contracts show NEW system!

### 3. Agent Name Confusion

**Tests assert OLD names:**
```python
assert trace.agent_sequence[1] == "BioinformaticsExecutor"
```

**Orchestrator uses NEW names:**
```python
agent_name=AgentName.IMPLEMENTATION.value  # "ImplementationAgent"
```

**Result:** Tests can pass while validating wrong system!

### 4. Missing Workflow Planner

**NEW orchestrator** expects a `WorkflowPlan` to already exist:
```python
trace = await orchestrator.run_pipeline(
    command="Run FastQC",
    workflow_plan=workflow_plan  # Where does this come from?
)
```

**OLD system** generates the plan internally:
```python
result = await handle_command("Run FastQC")  # Plan created inside
```

**Problem:** NEW system is incomplete - no equivalent of "BioinformaticsExecutor" to generate WorkflowPlan!

---

## Impact Analysis

### For Development
- 🔴 **High confusion** about which system to modify
- 🔴 **Duplicate work** (bug fixes must target both systems)
- 🔴 **Inconsistent behavior** between CLI and HTTP
- 🟡 **Test maintenance** (must update both OLD and NEW tests)

### For Production
- 🟢 **No immediate impact** (OLD system is stable and working)
- 🟡 **Future migration risk** (unclear path to adopt NEW system)
- 🟡 **Feature parity** (NEW system missing IntentDetector, Planner, Visualizer)

### For Observability
- 🔴 **Inconsistent traces** (OLD system: no hashes, NEW system: SHA256 hashes)
- 🔴 **Mixed agent names** (OLD: BioinformaticsExecutor, NEW: ImplementationAgent)
- 🟡 **Hard to compare** behavior between systems

---

## Recommendations

### Option A: Deprecate NEW System (Conservative) ⚠️

**Approach:**
1. Label `orchestrator.py` as **EXPERIMENTAL - DO NOT USE IN PRODUCTION**
2. Continue using OLD system (`agent.py`) for all production traffic
3. Gradually incorporate ideas from NEW system into OLD system:
   - Add Pydantic contracts to OLD system
   - Add provenance hashing to OLD system
   - Improve metadata flow in OLD system
4. Eventually remove `orchestrator.py` when features absorbed

**Pros:**
✅ Minimal immediate disruption  
✅ Proven system stays in production  
✅ Can cherry-pick best ideas from NEW system

**Cons:**
❌ Continues complexity of OLD system  
❌ Loses clean architecture of NEW system  
❌ Work invested in NEW system wasted

**Timeline:** 1-2 months to absorb features, 3-6 months to stabilize

---

### Option B: Complete NEW System & Migrate (Progressive) 🚀

**Approach:**
1. **Phase 1:** Complete NEW orchestrator (2-4 weeks)
   - Add WorkflowPlannerAgent (generates WorkflowPlan from user command)
   - Add IntentDetector (routes ask vs execute)
   - Add DataVisualizer (produces plots/reports)
   - Port policy checks from OLD system

2. **Phase 2:** HTTP API dual-mode (1-2 weeks)
   - Add `/execute_v2` endpoint using NEW orchestrator
   - Keep `/execute` using OLD system
   - Run both in parallel (shadow mode)

3. **Phase 3:** Validation (2-3 weeks)
   - Compare traces from OLD vs NEW system
   - Verify contract hashes match
   - Performance testing
   - Fix discrepancies

4. **Phase 4:** Migration (1-2 weeks)
   - Switch `/execute` to NEW orchestrator
   - Deprecate `/execute_legacy`
   - Update all tests to NEW agent names

5. **Phase 5:** Cleanup (1 week)
   - Remove OLD system (`agent.py`)
   - Update documentation
   - Archive OLD tests

**Pros:**
✅ Clean architecture with Pydantic contracts  
✅ Provenance hashing built-in  
✅ Better metadata flow  
✅ Simpler codebase (443 vs 1798 lines)  
✅ Easier to reason about and maintain

**Cons:**
❌ Significant engineering effort (8-12 weeks)  
❌ Migration risk (potential bugs)  
❌ Tests must be rewritten  
❌ Documentation must be updated

**Timeline:** 2-3 months for complete migration

---

### Option C: Hybrid Approach (Pragmatic) ✅ RECOMMENDED

**Approach:**
1. **Immediate (1 week):**
   - Label systems clearly in documentation:
     - `agent.py` → **"PRIMARY ORCHESTRATOR (Production)"**
     - `orchestrator.py` → **"EXPERIMENTAL (Prototype only)"**
   - Add comments in code marking each system
   - Update `SYSTEM_OVERVIEW.md` with clear "Production vs Prototype" sections

2. **Short-term (4-6 weeks):**
   - Add **best features** from NEW system to OLD system:
     - Provenance hashing (SHA256)
     - Rich metadata flow (pass size_bytes, location_type)
     - Structured logging with hashes
   - Keep OLD system as primary
   - Use NEW system for experimentation only

3. **Medium-term (3-6 months):**
   - Gradually refactor OLD system toward NEW architecture:
     - Extract agents into separate modules
     - Add Pydantic contracts alongside existing dataclasses
     - Simplify agent graph (reduce from 6+ to 4-5 agents)
   - Maintain backward compatibility throughout

4. **Long-term (6-12 months):**
   - Decide on full migration based on:
     - Production stability of OLD system with NEW features
     - Business value of clean architecture
     - Engineering bandwidth
   - If migration justified, follow Option B phases

**Pros:**
✅ Immediate clarity (labels systems clearly)  
✅ Low risk (incremental improvements)  
✅ Best of both worlds (OLD stability + NEW features)  
✅ Flexible timeline (adapt based on feedback)

**Cons:**
⚠️ Continues dual-system complexity short-term  
⚠️ Requires discipline (resist temptation to diverge)

**Timeline:** 1 week for clarity, 4-6 weeks for feature absorption

---

## Immediate Action Items (This Week)

### 1. Label Systems Clearly ✅

**In `backend/agent.py`:**
```python
"""
PRIMARY ORCHESTRATOR - PRODUCTION

This is the primary multi-agent orchestrator used by the HTTP API.
All production traffic flows through handle_command().

Status: 🟢 PRODUCTION
Entry Point: handle_command()
Used By: /execute, /chat, /mcp/call_tool

For experimental clean orchestrator, see backend/orchestrator.py
"""
```

**In `backend/orchestrator.py`:**
```python
"""
EXPERIMENTAL ORCHESTRATOR - PROTOTYPE ONLY

This is a clean 2-agent pipeline with Pydantic contracts and provenance hashing.
Currently used only for CLI demos and unit tests.

Status: 🟡 PROTOTYPE
Entry Point: Orchestrator.run_pipeline()
Used By: cli_demo.py, test_phase3_agents.py

For production orchestrator, see backend/agent.py
"""
```

### 2. Update Documentation ✅

**In `docs/SYSTEM_OVERVIEW.md`:**

Add section at top:
```markdown
## ⚠️ IMPORTANT: Two Orchestration Systems

Helix.AI currently has two orchestration systems:

1. **PRIMARY (Production):** Multi-agent graph in `backend/agent.py`
   - Entry point: `handle_command()`
   - Used by: HTTP API (/execute, /chat)
   - Agents: IntentDetector, BioinformaticsGuru, BioinformaticsExecutor, etc.
   - Status: 🟢 Stable, in production

2. **EXPERIMENTAL (Prototype):** Clean pipeline in `backend/orchestrator.py`
   - Entry point: `Orchestrator.run_pipeline()`
   - Used by: CLI demos, unit tests only
   - Agents: InfrastructureDecisionAgent, ImplementationAgent
   - Status: 🟡 Prototype, not used in production

The diagram below shows the PRIMARY system. For EXPERIMENTAL system details,
see `docs/ORCHESTRATION_DUALITY.md`.
```

### 3. Add System Markers in Tests ✅

**In `tests/demo_scenarios/test_scenarios.py`:**
```python
"""
Tests for PRIMARY orchestrator (backend/agent.py).

These tests validate the production multi-agent system used by the HTTP API.
For tests of the EXPERIMENTAL orchestrator (backend/orchestrator.py),
see tests/unit/backend/test_phase3_agents.py.
"""
```

**In `tests/unit/backend/test_phase3_agents.py`:**
```python
"""
Tests for EXPERIMENTAL orchestrator (backend/orchestrator.py).

These tests validate the prototype 2-agent pipeline with Pydantic contracts.
For tests of the PRIMARY orchestrator (backend/agent.py),
see tests/demo_scenarios/test_scenarios.py.
"""
```

---

## Decision Matrix

| Criterion | Option A: Deprecate NEW | Option B: Migrate | Option C: Hybrid ✅ |
|-----------|------------------------|-------------------|---------------------|
| **Immediate Clarity** | ⚠️ Medium | ❌ Low | ✅ High |
| **Engineering Effort** | 🟢 Low (1-2 months) | 🔴 High (3 months) | 🟡 Medium (4-6 weeks) |
| **Migration Risk** | 🟢 Low | 🔴 High | 🟡 Low-Medium |
| **Clean Architecture** | ❌ No | ✅ Yes | ⚠️ Eventually |
| **Backward Compat** | ✅ Full | ❌ Breaking | ✅ Full |
| **Provenance Hashing** | 🟡 Must add | ✅ Built-in | 🟡 Must add |
| **Timeline** | 1-2 months | 3 months | 1 week → 4-6 weeks → 6-12 months |
| **Recommended** | ❌ | ❌ | ✅ |

---

## Conclusion

**Recommendation:** **Option C: Hybrid Approach**

**Rationale:**
1. **Immediate clarity** (1 week) - Label systems, update docs, no confusion
2. **Low risk** - Incremental improvements, no big bang migration
3. **Best features** - Absorb provenance hashing, metadata flow from NEW system
4. **Flexible** - Can adapt timeline based on production needs
5. **Pragmatic** - Doesn't waste work on either system

**Next Steps:**
1. ✅ **This week:** Label systems clearly (code comments + docs)
2. ✅ **Next 4-6 weeks:** Add provenance hashing + metadata flow to PRIMARY system
3. ⏳ **Next 3-6 months:** Gradually refactor PRIMARY toward cleaner architecture
4. ⏳ **Next 6-12 months:** Decide on full migration vs convergence

---

## References

- **PRIMARY Orchestrator:** `backend/agent.py`
- **EXPERIMENTAL Orchestrator:** `backend/orchestrator.py`
- **Agent Registry:** `backend/config/agent_registry.py`
- **HTTP API:** `backend/main_with_mcp.py`
- **CLI Demo:** `backend/cli_demo.py`
- **System Overview:** `docs/SYSTEM_OVERVIEW.md`

---

**Document Status:** 🔴 **Action Required** - System labeling and documentation updates needed this week
