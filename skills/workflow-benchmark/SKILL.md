# Workflow benchmark for Helix.AI

## Purpose
Run benchmark and eval cases for Helix.AI workflows and compare results to release thresholds.

## When to use
Use this skill when:
- planner behavior changes,
- execution/orchestration changes,
- frontend workflow UX changes materially,
- preparing release readiness evaluation.

## Benchmark families
- planning quality
- execution success
- iterative refinement correctness
- bug recovery
- file integration
- frontend workflow stability
- latency and reliability

## Procedure
1. Read `benchmarks/release_thresholds.yaml`.
2. Select relevant benchmark cases.
3. Execute the benchmark suite or targeted subset.
4. Record:
   - suite name,
   - command used,
   - pass/fail counts,
   - timing if available,
   - notable regressions.
5. Update:
   - `artifacts/benchmark_results/`
   - `artifacts/release_readiness.json`

## Rules
- use fixed cases where possible,
- compare to baseline when available,
- do not mark success without writing machine-readable results.

## Success criteria
- benchmark results are fresh,
- threshold comparison is explicit,
- release readiness status is updated.
