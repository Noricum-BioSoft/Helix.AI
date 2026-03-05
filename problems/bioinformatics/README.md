# Bioinformatics problem library

This folder is a **curated repository of bioinformatics problems** used to test Helix.AI behavior (intent detection, tool routing, parameter extraction, and multi-step workflow planning).

## What goes in here

- **Problems**: realistic user prompts that represent common bioinformatics tasks (QC, trimming, read merging, alignment, phylogenetics, bulk RNA-seq, single-cell, metadata lookups).
- **Expectations**: the *minimum* expected behavior for deterministic test tiers (e.g. which tool should be selected; which parameters must be extracted).

## Where the tests live

Deterministic evaluation cases are stored as JSONL in:

- `tests/evals/cases/bioinformatics_router_tool_mapping.jsonl`
- `tests/evals/cases/bioinformatics_intent_classification.jsonl`

They are executed by:

- `tests/evals/test_eval_router_tool_mapping_bio.py`
- `tests/evals/test_eval_intent_classifier_bio.py`

## Running locally

```bash
python -m pytest tests/evals -v
```

## Adding a new problem

1. Pick the right JSONL file:
   - **Intent-only** cases → `bioinformatics_intent_classification.jsonl`
   - **Router / parameters / plan** cases → `bioinformatics_router_tool_mapping.jsonl`
2. Add a single JSON object per line with:
   - `id`: stable identifier (use `bio.<area>.<short_name>`)
   - `tags`: list of tags (area + goal like `safety`, `plan`, `s3`, `paths`)
   - `input`: `{ "text": ... }` (intent) or `{ "command": ..., "session_context": ... }` (router)
   - `expect`: expected `intent`, or expected `tool_name` + optional `parameters_subset`, or `plan_tool_names`
3. Keep `parameters_subset` strict and minimal: assert only the fields that must be correct.

