# Fixture Versioning and Governance

This document defines how benchmark fixtures are versioned, reviewed, and updated.

## Versioning Rules

- Every fixture must have:
  - immutable path/version identity (for example `counts_matrix.v1.tsv`), or
  - content hash tracking in the approval registry.
- Any material fixture content change requires a new version identifier.
- Existing scenario files may only be repointed to a new fixture version via review.

## Required Metadata

Each fixture entry in `fixture_approval_registry.yaml` must include:

- source provider and source URL,
- license type and redistribution status,
- human-data/sensitivity classification,
- approval status and restrictions,
- retention profile,
- reviewer role and review notes.

## Change Workflow

1. Propose fixture addition/update in a PR.
2. Update:
   - `benchmarks/compliance/fixture_approval_registry.yaml`,
   - affected scenario files in `benchmarks/cases/`,
   - this governance policy if rules changed.
3. Run compliance validation:
   - all scenario fixture paths must exist in registry,
   - no duplicate fixture registry keys.
4. Obtain approvals:
   - Security Lead for privacy/sensitivity,
   - QA/Benchmark Lead for benchmark fitness.
5. Merge only after CI compliance checks pass.

## Release Guardrails

- `pending` or `rejected` fixtures are not eligible for release holdout execution.
- Placeholder fixtures (`*.placeholder`) are allowed only for development planning;
  they block release-gate readiness until replaced.
- Any fixture with `restricted_human_data` must be run in approved secure environments
  and must not be exported as raw data in CI artifacts.

## Audit Trail

- Keep fixture changes reviewable through git history.
- Include rationale and risk notes in each fixture-change PR description.
