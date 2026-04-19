# Fix until green for Helix.AI

## Purpose
Take a failing behavior, reproduce it, fix it at the root cause, add regression coverage, and rerun the smallest relevant suite until green.

## When to use
Use this skill when:
- tests fail,
- workflow cases fail,
- frontend regressions are found,
- iterative workflow behavior is broken.

## Procedure
1. Reproduce the failure.
2. Identify whether the fault is in:
   - routing,
   - orchestration,
   - contract/schema handling,
   - checkpoint/history,
   - pipeline logic,
   - frontend state or API integration,
   - test fixture drift.
3. Fix the root cause.
4. Add or update the most specific regression test.
5. Rerun the smallest relevant verification command.
6. If green, expand to the next larger suite.
7. Update:
   - `artifacts/worklog.md`
   - `artifacts/failure_summaries/latest.md`
   - `artifacts/test_results/`

## Rules
- do not remove assertions to make tests pass,
- do not weaken release gates,
- do not hide failures behind fallback logic unless product behavior requires it,
- prefer structural fixes over prompt-only patches when the issue is architectural.

## Success criteria
- failure is reproducible before fix,
- targeted regression passes after fix,
- broader impacted suite passes,
- root cause documented.
