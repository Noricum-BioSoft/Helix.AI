# Release gate for Helix.AI

## Purpose
Determine whether Helix.AI is ready for release based on fresh tests and benchmarks.

## When to use
Use this skill when:
- preparing a release candidate,
- stabilization work is nearly complete,
- benchmark results have been refreshed.

## Procedure
1. Read `benchmarks/release_thresholds.yaml`.
2. Read latest results from:
   - `artifacts/test_results/`
   - `artifacts/benchmark_results/`
   - `artifacts/failure_summaries/`
3. Determine:
   - which suites are fresh,
   - which thresholds passed,
   - which blockers remain.
4. Write `artifacts/release_readiness.json`.

## Required output fields
- `timestamp`
- `release_ready`
- `suites`
- `metrics`
- `threshold_results`
- `blockers`
- `notes`

## Decision rules
- any failed release threshold => `release_ready: false`
- missing required suite => `release_ready: false`
- unresolved critical regression => `release_ready: false`

## Success criteria
- readiness status is explicit,
- blockers are concrete,
- result is machine-readable and traceable.
