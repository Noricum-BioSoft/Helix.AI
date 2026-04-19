# Holdout Access Policy

This policy prevents benchmark overfitting by limiting `release_holdout` case visibility.

## Rules

- `release_holdout` prompts, expected routes, and scorer targets are not visible to implementation contributors.
- Only QA/Benchmark Lead and delegated release reviewers can execute holdout suites directly.
- Holdout execution artifacts are retained in `artifacts/benchmark_results/holdout/`.
- Implementation contributors can only see aggregate holdout status (`pass/fail`, counts, blocker classes), not per-case details.
- Any exception requires written approval from QA/Benchmark Lead and Engineering Manager with expiry date.

## Change Control

- Changes to holdout case definitions require:
  - issue ticket,
  - rationale,
  - risk note,
  - approval from QA/Benchmark Lead.
- Split changes in `benchmarks/splits.yaml` must be reviewed by Security Lead for leakage risk.

