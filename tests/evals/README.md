## Helix.AI agent evaluation suite (`tests/evals/`)

This folder complements the existing unit/integration tests with **repeatable, metrics-driven agent evaluations**.

### Goals
- **Functionality**: tool routing + parameter extraction match expectations on diverse inputs.
- **Consistency**: the same inputs produce the same outputs (in deterministic/mock tiers).
- **Safety**: Q&A prompts must not trigger tool generation; disallowed patterns should be blocked/regressed quickly.
- **Quality**: capture accuracy/relevance/safety judgments for optional live-LLM runs, with a human feedback loop.

### Evaluation tiers
- **Tier 0 (CI / deterministic)**: runs with `HELIX_MOCK_MODE=1` (default in `tests/conftest.py`).
  - Tests `backend.intent_classifier.classify_intent`
  - Tests `backend.command_router.CommandRouter` tool + parameter mapping
  - Tests safety gates that are independent of LLM behavior (e.g., toolgen gating)
- **Tier 1 (CI / deterministic regression set)**: curated “golden” cases stored as JSONL in `tests/evals/cases/`.
- **Tier 2 (optional / live-LLM)**: requires API keys. Produces artifacts (JSON) for human scoring and regression promotion.

### Metrics (Tier 0/1)
- **Tool accuracy**: expected tool name matches the router’s tool selection.
- **Parameter accuracy**: expected parameter **subset** matches the extracted parameters.
- **Intent accuracy**: expected `execute|qa` matches classifier output.
- **Safety checks**:
  - Q&A intent does **not** trigger tool generation (see also `tests/backend/test_toolgen_gating.py`)
  - Unknown/ambiguous prompts route to safe fallbacks (no “magic” privileged tools)

### Running
- **All evals**:

```bash
python -m pytest tests/evals -v
```

- **Only router evals**:

```bash
python -m pytest tests/evals/test_eval_router_tool_mapping.py -v
```

### Adding a new case
Add a line to one of the JSONL files in `tests/evals/cases/` and include:
- `id`: stable identifier
- `input`: the prompt/command
- `expect`: expected tool/params/intent

See the existing case files for examples.





