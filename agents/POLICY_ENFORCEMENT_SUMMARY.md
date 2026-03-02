# Agent Responsibilities Policy Enforcement - Implementation Summary

This document summarizes the policy enforcement system implemented to make `agents/agent-responsibilities.md` a binding protocol.

## What Was Implemented

### Phase 1: Contract-Level Enforcement ("Must Not" Checks)

Created `backend/policy_checks.py` with validation functions that enforce the "must not" rules from the responsibilities spec:

#### 1. Planner Output Validation (`check_planner_output`)
**Rule**: Bioinformatics Executor (Planner) must not choose infrastructure.

**Checks**:
- Scans for infrastructure keywords (EC2, EMR, instance types, cluster configs)
- Rejects plans with `environment`, `recommended_environment`, or `target` fields
- Prevents planner from making infrastructure decisions

**Example violation**:
```python
invalid_plan = {
    "steps": [...],
    "environment": "EMR"  # ❌ NOT ALLOWED
}
check_planner_output(invalid_plan)
# Raises: PolicyViolationError with fix instructions
```

#### 2. Plan Mutation Detection (`check_plan_not_mutated`)
**Rule**: Infrastructure Expert must not mutate workflow plan steps.

**Checks**:
- Compares plan steps before and after Infrastructure Expert runs
- Detects added, removed, or modified steps
- Enforces the No-Mutation Rule (agents may only add annotations, not change upstream contracts)

**Example violation**:
```python
before = [{"id": "step1", "tool_name": "fastqc", ...}]
after = [{"id": "step1", "tool_name": "DIFFERENT_TOOL", ...}]  # ❌ Changed!
check_plan_not_mutated(before, after, "InfrastructureExpert")
# Raises: PolicyViolationError with diff
```

#### 3. Code Generator Output Validation (`check_codegen_output`)
**Rule**: Code Generator must not change scientific intent.

**Checks**:
- Ensures ExecutionSpec doesn't rewrite workflow steps
- Prevents CodeGen from adding infrastructure decision fields
- Scientific intent remains owned by Planner

**Example violation**:
```python
execution_spec = {
    "steps": [...],  # ❌ ExecutionSpec shouldn't have steps
    "instance_type": "m5.xlarge"  # ❌ Infrastructure decision
}
check_codegen_output(original_plan, execution_spec)
# Raises: PolicyViolationError
```

#### 4. Visualizer Output Validation (`check_visualizer_output`)
**Rule**: Visualizer must not introduce workflow steps or infrastructure changes.

**Checks**:
- Ensures output doesn't contain `steps`, `workflow`, or `environment` fields
- Prevents Visualizer from generating execution code
- Output should be VisualizationArtifacts only

**Example violation**:
```python
invalid_output = {
    "summary": "Results",
    "steps": [...]  # ❌ NOT ALLOWED
}
check_visualizer_output(invalid_output)
# Raises: PolicyViolationError
```

### Phase 2: Contract Tightening

**Outcome**: Contracts already exist in `shared/contracts.py` with proper types:
- `IntentResult`: Intent classification output
- `WorkflowPlan`: Wrapper around Plan IR with metadata
- `InfraDecisionContract`: Infrastructure selection output (with constraints/warnings)
- `ExecutionSpec`: Runnable specification
- `ExecutionResult`: Execution outcomes
- `VisualizationArtifacts`: Visualization outputs

The `InfraDecision` model in `backend/contracts/infra_decision.py` already includes:
- `warnings: List[str]` for uncertainty tracking
- `alternatives: List[InfraAlternative]` for tradeoff analysis
- `confidence_score: float` for decision confidence

### Phase 3: Policy Unit Tests

Created `tests/unit/backend/test_agent_responsibilities_policy.py` with **39 comprehensive tests**:

#### Test Categories

**1. Handoff Policy Legality (9 tests)**
- Intent Detector must be first
- Ask intent routes to Guru
- Execute intent routes to Planner
- Guru cannot handoff to Broker
- Infra must not be invoked before Planner
- Broker must be after Infra or CodeGen
- Visualizer must be after Broker
- CodeGen is optional
- Guru escalation requires user consent

**2. Contract "Must Not" Checks (19 tests)**
- Planner cannot choose environment
- Planner cannot include instance types
- Planner cannot include EMR config
- Planner cannot set execution target
- Infra cannot mutate steps
- Infra cannot add/remove steps
- CodeGen cannot rewrite steps
- CodeGen cannot choose infrastructure
- Visualizer cannot include steps
- Visualizer cannot choose environment
- Visualizer cannot include commands

**3. Happy-Path Sequences (4 tests)**
- Ask workflow: Intent → Guru
- Execute workflow (minimal): Intent → Planner → Infra → Broker → Visualizer
- Execute workflow (with CodeGen): Intent → Planner → Infra → CodeGen → Broker → Visualizer
- Full end-to-end validation

**4. Edge Cases (4 tests)**
- Error messages include fix instructions
- Empty agent sequence fails
- Unknown intent fails
- Plan mutation errors include diffs

**5. Integration Tests (3 tests)**
- HandoffPolicy accessible from agent module
- Policy checks use consistent error types
- Contracts importable from shared module

**Run tests**:
```bash
pytest tests/unit/backend/test_agent_responsibilities_policy.py -v
# Result: 39 passed in 0.24s ✓
```

### Phase 4: Orchestrator Refactor

Added staged functions to `backend/agent.py` (`CommandProcessor` class):

```python
# Stage 1: Intent Detection (always first)
async def _run_intent_detector(command: str) -> IntentResult

# Stage 2a: Ask path
async def _run_guru(command: str, session_context: Dict) -> Dict

# Stage 2b: Execute path
async def _run_planner(command: str, session_id: str, session_context: Dict) -> Dict

# Stage 3: Infrastructure selection
async def _run_infra(plan_or_mapping: Dict) -> Optional[Dict]

# Stage 4: Optional code generation
async def _run_codegen_if_needed(plan: Dict, infra: Dict) -> Optional[Dict]

# Stage 5: Execution handoff
def _handoff_to_broker(execution_spec: Dict) -> Dict

# Stage 6: Visualization
async def _run_visualizer(execution_result: Dict) -> Optional[Dict]
```

Each stage:
1. Validates handoff with `_validate_next_agent()`
2. Registers agent call with `_register_agent_call()`
3. Runs the agent logic
4. Validates output schema (if applicable)
5. Runs relevant policy checks

**Integration points**:
- `CommandProcessor.agent_sequence` tracks execution order
- `HandoffPolicy` validates transitions
- Policy checks validate outputs
- Existing `process()` method coordinates stages

## What Failures Look Like

### Handoff Policy Violation

```
PolicyViolationError: Illegal handoff: BioinformaticsGuru → ExecutionBroker.
Allowed transitions: ['BioinformaticsExecutor']

Valid workflows:
- Ask questions: IntentDetector → BioinformaticsGuru
- Execute workflows: IntentDetector → BioinformaticsExecutor → 
  InfrastructureExpert → ExecutionBroker → DataVisualizer
```

### Contract Violation (Planner choosing infrastructure)

```
PolicyViolation: Bioinformatics Executor (Planner) must not choose infrastructure.

The Planner's role is to create a workflow plan (steps, tools, inputs/outputs) 
without making infrastructure decisions. Infrastructure selection is owned by the 
Infrastructure Expert agent.

Violations found:
  - Infrastructure keyword 'emr' found in field name: plan.environment
  - Field 'environment' found in plan. Planner must not choose infrastructure; 
    that's the Infrastructure Expert's role.

How to fix:
  - Remove infrastructure-related fields from the plan
  - Keep the plan focused on WHAT to do (steps/tools), not WHERE to do it (infra)
  - Let the Infrastructure Expert decide environment, instance types, etc.
```

### Contract Violation (Infra mutating plan)

```
PolicyViolation: InfrastructureExpert must not mutate workflow plan steps.

Per the No-Mutation Rule in agent-responsibilities.md, agents may only:
  - Add annotations (warnings, assumptions, constraints)
  - Produce a new downstream contract

But InfrastructureExpert modified the plan steps:
  Line 3: '  "tool_name": "fastqc"' → '  "tool_name": "fastqc_distributed"'
  Line 8: '  "arguments": {}' → '  "arguments": {"cluster_mode": true}'

How to fix:
  - InfrastructureExpert should return its output separately (e.g., InfraDecision)
  - InfrastructureExpert should not modify the input WorkflowPlan.steps
  - Use constraints/advice fields to suggest non-invasive changes
```

## How to Extend the Policy

### Adding a New Agent

1. **Define agent in `agents/agent-responsibilities.md`**:
   ```markdown
   ## New Agent Name
   **Owner of:** [specific decision type]
   
   **Input:** [contract types]
   **Output:** [contract type]
   
   **May:** [allowed actions]
   **Must not:** [forbidden actions]
   ```

2. **Add to `AgentRole` enum** in `backend/agent.py`:
   ```python
   class AgentRole(str, Enum):
       # ...existing...
       NEW_AGENT = "NewAgentName"
   ```

3. **Update `HandoffPolicy.ALLOWED_HANDOFFS`**:
   ```python
   ALLOWED_HANDOFFS = {
       # ...existing...
       AgentRole.SOME_AGENT: [AgentRole.NEW_AGENT],
       AgentRole.NEW_AGENT: [AgentRole.NEXT_AGENT],
   }
   ```

4. **Create contract in `shared/contracts.py`**:
   ```python
   class NewAgentOutput(BaseModel):
       """OWNER: New Agent Name"""
       # fields...
   ```

5. **Add policy check in `backend/policy_checks.py`**:
   ```python
   def check_new_agent_output(output: Union[Dict, BaseModel]) -> None:
       """
       Validate that New Agent output does NOT [forbidden action].
       
       Per agent-responsibilities.md:
           "New Agent - Must not: [rule]"
       """
       # validation logic
       if violation:
           raise PolicyViolationError("...")
   ```

6. **Add staged function in `backend/agent.py`**:
   ```python
   async def _run_new_agent(self, input_data: Dict) -> Dict:
       """Stage N: New Agent."""
       from backend.policy_checks import check_new_agent_output
       
       self._validate_next_agent(AgentRole.NEW_AGENT)
       self._register_agent_call(AgentRole.NEW_AGENT)
       
       # Run agent logic
       result = await some_agent_function(input_data)
       
       # Validate output
       check_new_agent_output(result)
       
       return result
   ```

7. **Add tests in `tests/unit/backend/test_agent_responsibilities_policy.py`**:
   ```python
   class TestNewAgentContractChecks:
       def test_new_agent_valid_output(self):
           # Test happy path
       
       def test_new_agent_cannot_do_forbidden_thing(self):
           # Test policy violation
   ```

### Adding a New "Must Not" Rule

1. **Document rule in `agents/agent-responsibilities.md`** under agent's "Must not:" section

2. **Implement check in `backend/policy_checks.py`**:
   ```python
   def check_agent_new_rule(output: Union[Dict, BaseModel]) -> None:
       """Validate agent does not [forbidden action]."""
       # Check logic
       if violation_detected:
           message = (
               f"PolicyViolation: {AgentName} must not [rule].\n\n"
               f"Violations found:\n  - {violation_details}\n\n"
               f"How to fix:\n  - {fix_instruction}"
           )
           raise PolicyViolationError(message)
   ```

3. **Integrate check into staged function**:
   ```python
   async def _run_agent(self, ...):
       # ... existing code ...
       check_agent_new_rule(output)
       return output
   ```

4. **Add test case**:
   ```python
   def test_agent_cannot_do_new_forbidden_thing(self):
       invalid_output = {
           # ... construct violation ...
       }
       with pytest.raises(PolicyCheckError, match="must not"):
           check_agent_new_rule(invalid_output)
   ```

### Adding a New Handoff Route

1. **Update `HandoffPolicy.ALLOWED_HANDOFFS`**:
   ```python
   ALLOWED_HANDOFFS = {
       AgentRole.FROM_AGENT: [
           # ...existing...
           AgentRole.NEW_TARGET,  # Add new allowed target
       ],
   }
   ```

2. **Document rationale in `agents/agent-responsibilities.md`**

3. **Add test case**:
   ```python
   def test_from_agent_can_handoff_to_new_target(self):
       policy = HandoffPolicy()
       # Should not raise
       policy.validate_handoff(
           AgentRole.FROM_AGENT,
           AgentRole.NEW_TARGET
       )
   ```

## Integration with Existing Code

### Current State
- **Handoff policy**: Fully integrated in `CommandProcessor.process()`
- **Policy checks**: Integrated in staged functions (partial - TODOs noted)
- **Contracts**: Defined in `shared/contracts.py` and `backend/contracts/`
- **Tests**: Comprehensive coverage (39 tests pass)

### Next Steps for Full Integration

1. **Complete Planner integration**:
   - Parse tool_mapping into WorkflowPlan
   - Run full `check_planner_output()` validation

2. **Complete Infrastructure Expert integration**:
   - Call `infrastructure_decision_agent_v2.py` from `_run_infra()`
   - Capture plan steps before/after and validate no mutation

3. **Complete Code Generator integration**:
   - Add conditional logic to invoke CodeGen only when needed
   - Parse output into ExecutionSpec and validate

4. **Complete Visualizer integration**:
   - Add visualization agent call after execution
   - Validate output with `check_visualizer_output()`

5. **End-to-end workflow**:
   - Wire up all stages in `CommandProcessor.process()`
   - Replace current tool-mapping flow with staged workflow
   - Maintain backward compatibility

## Files Modified/Created

### Created
- `backend/policy_checks.py` (520 lines) - Contract-level enforcement
- `tests/unit/backend/test_agent_responsibilities_policy.py` (717 lines) - Comprehensive tests
- `agents/POLICY_ENFORCEMENT_SUMMARY.md` (this file) - Documentation

### Modified
- `backend/agent.py` - Added:
  - Documentation header with policy overview
  - Staged orchestration functions (_run_intent_detector, _run_planner, etc.)
  - Integration with policy checks (TODOs for full integration)

## Testing

**Run all policy tests**:
```bash
pytest tests/unit/backend/test_agent_responsibilities_policy.py -v
```

**Expected output**:
```
39 passed in 0.24s ✓
```

**Run specific test category**:
```bash
# Handoff policy tests
pytest tests/unit/backend/test_agent_responsibilities_policy.py::TestHandoffPolicyLegality -v

# Contract checks
pytest tests/unit/backend/test_agent_responsibilities_policy.py::TestPlannerContractChecks -v

# Happy paths
pytest tests/unit/backend/test_agent_responsibilities_policy.py::TestHappyPathSequences -v
```

## Benefits

1. **Enforceable Boundaries**: Agents cannot accidentally violate their responsibilities
2. **Early Detection**: Policy violations caught at runtime with actionable errors
3. **Test Coverage**: 39 tests ensure policy remains enforced
4. **Maintainability**: Clear separation of concerns per agent
5. **Extensibility**: Well-defined patterns for adding agents/rules
6. **Documentation**: Responsibilities spec is now the source of truth

## Summary

The implementation makes `agents/agent-responsibilities.md` a **binding protocol** through:
- **Handoff routing** enforcement (HandoffPolicy)
- **Contract-level** enforcement (policy_checks.py)
- **Comprehensive testing** (39 unit tests)
- **Staged orchestration** (refactored CommandProcessor)

All new checks produce **actionable error messages** with fix instructions, making it easy to maintain correctness as the system evolves.
