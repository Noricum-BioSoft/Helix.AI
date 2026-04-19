# Benchmark Artifact Retention and Redaction Policy

This policy governs benchmark artifacts and logs produced by Helix benchmark runs.

## Data Classes

- `public`: synthetic/non-sensitive benchmark content.
- `internal`: internal fixtures, prompts, intermediate run traces.
- `restricted_human_data`: any fixture or artifact derived from human genomic/clinical data.

## Retention Windows

- `public`:
  - benchmark reports: 365 days
  - debug logs: 90 days
- `internal`:
  - benchmark reports: 180 days
  - debug logs: 60 days
- `restricted_human_data`:
  - benchmark reports (redacted only): 90 days
  - raw artifacts/logs: 30 days maximum

## Redaction Requirements

- Never store or publish direct identifiers in benchmark artifacts.
- For restricted human data:
  - store aggregate metrics only in shared reports,
  - redact sample-level identifiers before artifact persistence,
  - do not upload raw fixture files as CI artifacts.
- Error traces should redact file payload snippets and any unique sample identifiers.

## Allowed CI Artifacts

- Allowed:
  - `artifacts/benchmark_results/gate_report.json` (summary-level only)
  - aggregate score summaries and blocker classes
- Disallowed:
  - raw fixture content,
  - full sample-level outputs for restricted human data,
  - unredacted execution logs containing sensitive fields.

## Deletion and Enforcement

- A scheduled cleanup job must enforce retention windows by class.
- Security Lead owns redaction policy audits.
- QA/Benchmark Lead owns benchmark artifact scope audits.
- Any policy violation blocks release readiness until remediated and documented.
