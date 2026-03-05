# Agent Handoff Policy

This document describes the enforced agent handoff policy implemented in `backend/agent.py`.

## Overview

The `HandoffPolicy` class enforces the routing rules defined in [agent-responsibilities.md](./agent-responsibilities.md) to ensure:

1. **Predictable workflows**: Agents are invoked in a well-defined order
2. **Clear ownership**: Each decision type has exactly one owner
3. **No illegal calls**: Agents cannot skip steps or call agents they shouldn't

## Visual Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                     User Request                                 │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
                  ┌────────────────────┐
                  │ Intent Detector    │ (Always first)
                  └────────┬───────────┘
                           │
                ┌──────────┴──────────┐
                │                     │
                ▼                     ▼
      ┌─────────────────┐   ┌─────────────────┐
      │ Ask Intent      │   │ Execute Intent  │
      │ (Q&A)           │   │ (Workflow)      │
      └────────┬────────┘   └────────┬────────┘
               │                     │
               ▼                     ▼
      ┌─────────────────┐   ┌─────────────────┐
      │ Bioinformatics  │   │ Bioinformatics  │
      │ Guru            │   │ Executor        │
      │ (Answer only)   │   │ (Planner)       │
      └────────┬────────┘   └────────┬────────┘
               │                     │
               │ (escalate           ▼
               │  with user   ┌─────────────────┐
               │  consent)    │ Infrastructure  │
               └──────────────▶ Expert          │
                              └────────┬────────┘
                                       │
                            ┌──────────┴──────────┐
                            │                     │
                            ▼                     ▼
                   ┌─────────────────┐   ┌─────────────────┐
                   │ Code Generator  │   │ Execution       │
                   │ (Optional)      │   │ Broker          │
                   └────────┬────────┘   └────────┬────────┘
                            │                     │
                            └──────────┬──────────┘
                                       │
                                       ▼
                              ┌─────────────────┐
                              │ Execution       │
                              │ Broker          │
                              └────────┬────────┘
                                       │
                                       ▼
                              ┌─────────────────┐
                              │ Data            │
                              │ Visualizer      │
                              │ (Terminal)      │
                              └─────────────────┘
```

## Policy Rules

### 1. Intent Detector is Always First

Every user request must start with the Intent Detector agent, which classifies the intent as either:
- `ask` (Q&A request) → Routes to Bioinformatics Guru
- `execute` (Workflow execution) → Routes to Bioinformatics Executor (Planner)

```python
# Valid
sequence = [IntentDetector, Guru]

# Invalid - raises PolicyViolationError
sequence = [Guru, IntentDetector]
```

### 2. Ask Path: IntentDetector → Guru

For Q&A requests (`intent=ask`):
- Only the Guru agent can answer
- Guru may escalate to Planner, but **only with explicit user consent**

```python
# Valid ask workflow
sequence = [IntentDetector, Guru]

# Valid escalation (with user consent)
sequence = [IntentDetector, Guru, Planner, Infra, Broker, Visualizer]

# Invalid - escalation without user consent
policy.validate_handoff(from_agent=Guru, to_agent=Planner, user_consent=False)
# Raises: PolicyViolationError("Guru → Planner escalation requires explicit user consent")
```

### 3. Execute Path: IntentDetector → Planner → Infra → Broker → Visualizer

For execution requests (`intent=execute`):

**Mandatory sequence:**
1. **Intent Detector** (always first)
2. **Planner** (creates workflow plan)
3. **Infrastructure Expert** (decides where to execute)
4. **Code Generator** (optional - fills tool gaps)
5. **Execution Broker** (performs side effects)
6. **Data Visualizer** (creates visualizations)

```python
# Valid execute workflow (without CodeGen)
sequence = [
    IntentDetector,
    Planner,
    Infra,
    Broker,
    Visualizer
]

# Valid execute workflow (with CodeGen)
sequence = [
    IntentDetector,
    Planner,
    Infra,
    CodeGen,
    Broker,
    Visualizer
]

# Invalid - skipping Infra
sequence = [IntentDetector, Planner, Broker]
# Raises: PolicyViolationError("Illegal handoff: Planner → Broker")
```

### 4. Disallowed Transitions (Illegal Calls)

The following transitions are explicitly forbidden:

| From Agent | To Agent | Why Forbidden |
|------------|----------|---------------|
| Guru | Broker | Guru is read-only, cannot trigger execution |
| Guru | Infra | Guru doesn't make infrastructure decisions |
| Planner | Broker | Must go through Infra first |
| Planner | CodeGen | Must go through Infra first |
| Infra | Visualizer | Must go through Broker first |
| CodeGen | Visualizer | Must go through Broker first |
| Visualizer | Any | Visualizer is terminal (no further handoffs) |

```python
# These all raise PolicyViolationError
policy.validate_handoff(from_agent=Guru, to_agent=Broker)
policy.validate_handoff(from_agent=Planner, to_agent=Broker)
policy.validate_handoff(from_agent=Infra, to_agent=Visualizer)
```

## Usage

### In Code

The policy is automatically enforced by `CommandProcessor`:

```python
from backend.agent import CommandProcessor

processor = CommandProcessor()

# Process a command - policy is automatically enforced
result = await processor.process("What is DNA?", session_id="123")

# View the agent execution sequence
print(processor.agent_sequence)
# Output: [AgentRole.INTENT_DETECTOR, AgentRole.GURU]
```

### Validating Custom Workflows

You can validate custom agent sequences:

```python
from backend.agent import HandoffPolicy, AgentRole, PolicyViolationError

policy = HandoffPolicy()

# Validate a complete workflow
try:
    sequence = [
        AgentRole.INTENT_DETECTOR,
        AgentRole.PLANNER,
        AgentRole.INFRA,
        AgentRole.BROKER,
        AgentRole.VISUALIZER,
    ]
    policy.validate_workflow_sequence(sequence, intent="execute")
    print("✅ Workflow is valid")
except PolicyViolationError as e:
    print(f"❌ Invalid workflow: {e}")
```

### Getting Allowed Next Agents

Query which agents can be called next:

```python
from backend.agent import get_handoff_policy, AgentRole

policy = get_handoff_policy()

# From Planner, where can we go?
allowed = policy.get_allowed_next_agents(AgentRole.PLANNER)
print(allowed)
# Output: [AgentRole.INFRA]

# From Infra, where can we go?
allowed = policy.get_allowed_next_agents(AgentRole.INFRA)
print(allowed)
# Output: [AgentRole.CODEGEN, AgentRole.BROKER]
```

## Error Handling

When a policy violation occurs, the system returns a standardized error response:

```json
{
  "status": "error",
  "success": false,
  "error": "POLICY_VIOLATION",
  "message": "Illegal handoff: Planner → Broker. Allowed transitions: [InfrastructureExpert]",
  "agent_sequence": ["IntentDetector", "BioinformaticsExecutor"],
  "errors": [
    {
      "code": "POLICY_VIOLATION",
      "message": "Illegal handoff: Planner → Broker. Allowed transitions: [InfrastructureExpert]",
      "severity": "error"
    }
  ]
}
```

## Design Rationale

### Why Enforce These Rules?

1. **Separation of Concerns**: Each agent has a single, well-defined responsibility
2. **Testability**: Predictable workflows are easier to test and validate
3. **Observability**: Clear agent sequence makes debugging and tracing easier
4. **Safety**: Prevents agents from performing actions outside their scope
5. **Maintainability**: Changes to one agent don't affect others

### Why CodeGen is Optional

The Code Generator is only invoked when:
- The Planner requires a capability not in the existing toolbox, OR
- Infrastructure constraints require packaging not available

This keeps the common path fast (skip CodeGen when tools exist) while allowing dynamic tool generation when needed.

### Why Guru Can't Execute

The Bioinformatics Guru is a read-only agent that answers questions. Allowing it to trigger execution would:
- Blur responsibilities (answering vs. executing)
- Create security concerns (education agent performing side effects)
- Make it harder to reason about system behavior

## Testing

The policy is thoroughly tested in `tests/unit/backend/test_handoff_policy.py`:

```bash
# Run policy tests
pytest tests/unit/backend/test_handoff_policy.py -v

# Run with coverage
pytest tests/unit/backend/test_handoff_policy.py --cov=backend.agent --cov-report=html
```

## Future Enhancements

Potential improvements to the policy:

1. **Dynamic policy configuration**: Load rules from a config file
2. **Policy versioning**: Support multiple policy versions for A/B testing
3. **Agent capability checks**: Verify agents have required permissions before handoff
4. **Audit logging**: Record all handoffs for compliance and debugging
5. **Circuit breakers**: Prevent agent loops or infinite recursion

## Related Documents

- [agent-responsibilities.md](./agent-responsibilities.md) - Full agent specification
- [intent-detector-agent.md](./intent-detector-agent.md) - Intent classification
- [infrastructure-decision-agent.md](./infrastructure-decision-agent.md) - Infrastructure selection
- [tool-generator-agent.md](./tool-generator-agent.md) - Dynamic tool generation
