# Generate use cases for Helix.AI

## Purpose
Generate realistic, testable use cases from Helix.AI capabilities and existing code paths.

## When to use
Use this skill when:
- creating workflow tests,
- creating eval cases,
- expanding demo scenarios,
- increasing benchmark coverage.

## Use-case families
Always cover:
1. user provides explicit plan,
2. system generates plan,
3. iterative refinement after initial output,
4. plot update request,
5. additional file integration,
6. bug-fix / recovery request,
7. malformed input or inconsistent files,
8. workflow resume from checkpoint,
9. frontend retry / resume flow,
10. partial agent/tool failure.

## Procedure
1. Read current workflow and pipeline implementations.
2. Read current tests and demo scenarios.
3. Identify under-tested workflow branches.
4. Generate cases that are:
   - realistic,
   - reproducible,
   - automatable,
   - aligned with product behavior.
5. Prefer coverage over novelty.
6. Add or propose cases under:
   - `tests/demo_scenarios/`
   - `tests/workflows/`
   - `tests/evals/`
   - `benchmarks/cases/`

## Output format
For each case include:
- title,
- user prompt,
- input files or fixtures,
- expected workflow behavior,
- expected output properties,
- failure modes to check,
- pass criteria.

## Success criteria
- each major workflow family has coverage,
- cases are executable or near-executable,
- cases can be used in benchmarks or regression suites.
