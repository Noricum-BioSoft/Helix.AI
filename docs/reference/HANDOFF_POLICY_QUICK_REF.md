# Handoff Policy Quick Reference

## Valid Agent Sequences

### Ask Intent (Q&A)
```
IntentDetector → Guru
```

### Ask Intent with Escalation
```
IntentDetector → Guru → Planner → Infra → Broker → Visualizer
                        ↑
                 (requires user consent)
```

### Execute Intent (Basic)
```
IntentDetector → Planner → Infra → Broker → Visualizer
```

### Execute Intent (with CodeGen)
```
IntentDetector → Planner → Infra → CodeGen → Broker → Visualizer
```

## Allowed Transitions Table

| Current Agent | Can Call Next |
|--------------|---------------|
| IntentDetector | Guru, Planner |
| Guru | Planner (with user consent) |
| Planner | Infra |
| Infra | CodeGen, Broker |
| CodeGen | Broker |
| Broker | Visualizer |
| Visualizer | (none - terminal) |

## Common Errors

### ❌ Skipping Infra
```python
# WRONG: Planner → Broker (skips Infra)
IntentDetector → Planner → Broker
```
**Error**: `Illegal handoff: BioinformaticsExecutor → ExecutionBroker`

### ❌ Guru Calling Broker
```python
# WRONG: Guru cannot trigger execution
IntentDetector → Guru → Broker
```
**Error**: `Illegal handoff: BioinformaticsGuru → ExecutionBroker`

### ❌ Starting Without IntentDetector
```python
# WRONG: Must start with IntentDetector
Planner → Infra → Broker
```
**Error**: `First agent must be IntentDetector`

### ❌ Wrong Agent for Intent
```python
# WRONG: Ask intent must go to Guru, not Planner
IntentDetector → Planner  # when intent="ask"
```
**Error**: `For intent 'ask', second agent must be BioinformaticsGuru`

## Code Examples

### Query Allowed Next Agents
```python
from backend.agent import get_handoff_policy, AgentRole

policy = get_handoff_policy()

# What can I call from Planner?
allowed = policy.get_allowed_next_agents(AgentRole.PLANNER)
print(allowed)  # [AgentRole.INFRA]
```

### Validate a Handoff
```python
from backend.agent import HandoffPolicy, AgentRole, PolicyViolationError

policy = HandoffPolicy()

try:
    policy.validate_handoff(
        from_agent=AgentRole.PLANNER,
        to_agent=AgentRole.BROKER
    )
except PolicyViolationError as e:
    print(f"Invalid handoff: {e}")
```

### Check Intent Routing
```python
from backend.agent import get_handoff_policy

policy = get_handoff_policy()

# What agent handles "ask" intent?
next_agent = policy.get_next_agent_for_intent("ask")
print(next_agent)  # AgentRole.GURU

# What agent handles "execute" intent?
next_agent = policy.get_next_agent_for_intent("execute")
print(next_agent)  # AgentRole.PLANNER
```

## Decision Tree

```
User Request
    │
    ▼
Is IntentDetector first?
    │
    ├─ NO → PolicyViolationError ❌
    │
    └─ YES → Classify Intent
              │
              ├─ "ask" → Route to Guru ✅
              │           │
              │           └─ Guru wants to escalate?
              │                 │
              │                 ├─ User consent? YES → Continue to Planner ✅
              │                 └─ User consent? NO → PolicyViolationError ❌
              │
              └─ "execute" → Route to Planner ✅
                             │
                             └─ Planner → Infra → (CodeGen?) → Broker → Visualizer ✅
```

## Agent Responsibilities (One-line Summary)

| Agent | Responsibility | Can Execute? |
|-------|---------------|-------------|
| IntentDetector | Classify user intent | No |
| Guru | Answer questions, provide guidance | No |
| Planner | Create workflow plans | No |
| Infra | Select execution environment | No |
| CodeGen | Generate missing tools (optional) | No |
| Broker | Execute workflows, perform side effects | **Yes** |
| Visualizer | Create visualizations and reports | No |

**Key Rule**: Only the Execution Broker can perform side effects (job submission, I/O, etc.)
