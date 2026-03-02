# Helix.AI Demo & Eval Scenarios

This directory contains the **scenario-driven demo and evaluation framework** for the Helix.AI multi-agent system.

## Purpose

This framework serves **three purposes simultaneously**:
1. **Developer regression tests** - automated testing via pytest
2. **Automated evals** - agent behavior validation and tracking over time
3. **Human-readable demos** - CLI-based demonstrations of system capabilities

## Architecture

```
tests/demo_scenarios/
├── scenarios/              # YAML scenario definitions (source of truth)
│   ├── ask_*.yaml         # Question-answering scenarios
│   ├── execute_*.yaml     # Workflow execution scenarios
│   └── edge_*.yaml        # Edge cases and error conditions
├── framework/             # Core evaluation framework
│   ├── scenario.py        # Scenario models and loader
│   ├── executor.py        # Scenario execution engine
│   ├── tracer.py          # Agent call trace recorder
│   ├── validator.py       # Contract and policy validators
│   └── reporter.py        # Results formatter and diffing
├── test_scenarios.py      # Pytest integration
├── demo_cli.py            # Interactive CLI demo mode
└── README.md              # This file
```

## Quick Start

### Running Pytest Tests
```bash
# Run all scenario tests
pytest tests/demo_scenarios/test_scenarios.py -v

# Run specific category
pytest tests/demo_scenarios/test_scenarios.py -k "ask"

# Run with detailed trace output
pytest tests/demo_scenarios/test_scenarios.py -v --show-traces
```

### Running CLI Demo Mode
```bash
# Interactive demo with all scenarios
python tests/demo_scenarios/demo_cli.py

# Run specific scenario
python tests/demo_scenarios/demo_cli.py --scenario ask_basic_question

# Run and show diff against baseline
python tests/demo_scenarios/demo_cli.py --scenario execute_fastqc --show-diff
```

## Scenario Definition Format

Scenarios are defined in YAML with the following structure:

```yaml
# scenarios/ask_basic_question.yaml
metadata:
  id: ask_basic_question
  category: ask
  description: "Simple bioinformatics question requiring Guru response"
  tags: [baseline, guru-only]

input:
  user_prompt: "What is the difference between RNA-seq and DNA-seq?"
  session_context: {}

expected_behavior:
  agent_sequence:
    - IntentDetector
    - BioinformaticsGuru
  
  intent:
    type: "ask"
    confidence_min: 0.7
  
  contracts:
    - agent: IntentDetector
      output_type: IntentResult
      validations:
        - field: intent
          equals: "ask"
        - field: confidence
          greater_than: 0.7
    
    - agent: BioinformaticsGuru
      output_type: Answer
      validations:
        - field: text
          contains: ["RNA", "transcriptome", "DNA"]
        - field: text
          min_length: 100

  policy_checks:
    - no_execution_side_effects
    - no_infrastructure_decisions
    - strict_role_separation
```

## Validation Layers

The framework validates multiple layers:

### 1. Agent Sequence Validation
- Correct agent invocation order
- Adherence to handoff policy
- Intent-based routing

### 2. Contract Validation
- Each agent outputs the correct contract type
- Contract fields meet expectations
- No unauthorized contract mutations

### 3. Policy Enforcement
- Role separation (e.g., Guru doesn't execute)
- Infrastructure decisions only from Infrastructure Expert
- No side effects except from Execution Broker

### 4. Confidence & Uncertainty
- Agents express appropriate confidence levels
- Low confidence triggers clarifying questions
- Uncertainty is surfaced, not hidden

## Scenario Categories

### Ask Scenarios (`ask_*.yaml`)
Test the question-answering pathway:
- Intent Detector → Guru
- No execution, no infrastructure selection
- Knowledge retrieval and explanation

### Execute Scenarios (`execute_*.yaml`)
Test the workflow execution pathway:
- Intent Detector → Planner → Infra → Broker → Visualizer
- Optional CodeGen insertion
- Full workflow lifecycle

### Edge Cases (`edge_*.yaml`)
Test boundary conditions:
- Low confidence handling
- Missing information
- Ambiguous requests
- Policy violations

### Multi-Turn Scenarios (`multi_*.yaml`)
Test conversation context:
- Follow-up questions
- Guru escalation to execution (with consent)
- Context preservation

## Output Artifacts

Each scenario run produces:

### 1. Execution Trace
```json
{
  "scenario_id": "ask_basic_question",
  "timestamp": "2026-01-17T10:30:00Z",
  "agent_sequence": ["IntentDetector", "BioinformaticsGuru"],
  "contracts": [
    {
      "agent": "IntentDetector",
      "contract_type": "IntentResult",
      "data": {...}
    },
    {
      "agent": "BioinformaticsGuru",
      "contract_type": "Answer",
      "data": {...}
    }
  ],
  "validation_results": {...},
  "duration_ms": 2341
}
```

### 2. Validation Report
- Pass/fail for each validation
- Contract diff (if baseline exists)
- Policy violation details (if any)

### 3. Baseline Snapshot
- Stored in `baselines/` for regression detection
- Can be updated with `--update-baseline` flag

## Best Practices

### Writing Scenarios

1. **Be specific**: Each scenario should test one clear behavior
2. **Use real prompts**: Based on actual user requests when possible
3. **Validate contracts, not prose**: Focus on structured outputs
4. **Set meaningful thresholds**: Confidence levels, response length, etc.

### Maintaining Baselines

1. **Review diffs carefully** before updating baselines
2. **Document intentional changes** in commit messages
3. **Run full suite** before releasing changes

### Debugging Failures

```bash
# Show detailed trace
pytest tests/demo_scenarios/test_scenarios.py::test_ask_basic_question -vv

# Run in demo mode for interactive inspection
python tests/demo_scenarios/demo_cli.py --scenario ask_basic_question --verbose

# Compare against baseline
python tests/demo_scenarios/demo_cli.py --scenario ask_basic_question --show-diff
```

## Integration with CI/CD

Scenarios can be integrated into CI:

```yaml
# .github/workflows/test-agents.yml
- name: Run Agent Demo Scenarios
  run: |
    pytest tests/demo_scenarios/test_scenarios.py \
      --junitxml=test-results/demo-scenarios.xml \
      --scenario-report=test-results/scenario-report.json
```

## Design Principles

1. **Observable**: Agent behavior captured in structured traces
2. **Testable**: Automated assertions on contracts and sequences
3. **Diffable**: Changes over time are visible and trackable
4. **Reproducible**: Same input produces consistent agent behavior
5. **Human-readable**: Can be run interactively for demos

## Related Documentation

- `agents/agent-responsibilities.md` - Agent role specifications
- `agents/HANDOFF_POLICY_QUICK_REF.md` - Valid agent transitions
- `shared/contracts.py` - Contract type definitions
- `backend/agent.py` - Handoff policy implementation
