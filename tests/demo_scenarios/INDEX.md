# Demo/Eval System Documentation Index

Quick navigation to all documentation files.

## 📚 Documentation Files

### Getting Started
- **[QUICK_START.md](QUICK_START.md)** - Start here! 5-minute introduction
- **[README.md](README.md)** - Complete system documentation

### For Developers
- **[SCENARIO_AUTHORING.md](SCENARIO_AUTHORING.md)** - How to write scenarios
- **[INTEGRATION.md](INTEGRATION.md)** - CI/CD integration guide
- **[IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md)** - Technical architecture

### For Project Root
- **[../DEMO_EVAL_SYSTEM_COMPLETE.md](../../DEMO_EVAL_SYSTEM_COMPLETE.md)** - Implementation summary at project root

## 🚀 Quick Commands

### Verification
```bash
python tests/demo_scenarios/verify_installation.py
```

### Interactive Demo
```bash
# List scenarios
python tests/demo_scenarios/demo_cli.py --list

# Run specific scenario
python tests/demo_scenarios/demo_cli.py --scenario ask_basic_question

# Interactive mode
python tests/demo_scenarios/demo_cli.py
```

### Pytest Tests
```bash
# All tests
pytest tests/demo_scenarios/test_scenarios.py -v

# Specific category
pytest tests/demo_scenarios/test_scenarios.py -k "ask" -v

# With traces
pytest tests/demo_scenarios/test_scenarios.py -v --show-traces

# Create baselines
pytest tests/demo_scenarios/test_scenarios.py --update-baseline

# Compare with baselines
pytest tests/demo_scenarios/test_scenarios.py --compare-baseline
```

## 📁 Directory Structure

```
tests/demo_scenarios/
├── framework/           # Core evaluation engine
├── scenarios/           # Scenario definitions (YAML)
├── baselines/           # Baseline traces (generated)
├── test_scenarios.py    # Pytest integration
├── demo_cli.py          # Interactive CLI
└── *.md                 # Documentation
```

## 🎯 Where to Start

### As a New User
1. Read **QUICK_START.md** (5 min)
2. Run `verify_installation.py`
3. Try `demo_cli.py --list`
4. Run your first scenario

### As a Developer
1. Read **QUICK_START.md**
2. Review example scenarios in `scenarios/`
3. Read **SCENARIO_AUTHORING.md**
4. Create your first custom scenario

### As a DevOps/SRE
1. Read **QUICK_START.md**
2. Read **INTEGRATION.md**
3. Set up CI/CD integration
4. Configure baseline management

### As an Architect
1. Read **IMPLEMENTATION_SUMMARY.md**
2. Review `framework/` source code
3. Understand validation layers
4. Plan extensions

## 📊 System Statistics

- **Total Lines of Code**: ~3,000
- **Framework Modules**: 5
- **Example Scenarios**: 10
- **Documentation Files**: 6
- **Validation Operators**: 8
- **Test Coverage**: All agent paths

## 🔗 Related Documentation

- `agents/agent-responsibilities.md` - Agent role specifications
- `agents/HANDOFF_POLICY_QUICK_REF.md` - Valid agent transitions
- `shared/contracts.py` - Contract type definitions
- `backend/agent.py` - Handoff policy implementation

## 💡 Tips

- Start with existing scenarios as templates
- Use `--show-traces` for debugging
- Create baselines early to detect regressions
- Run scenarios in CI/CD on every PR
- Keep scenarios focused on one behavior each
