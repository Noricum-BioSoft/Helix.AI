# P3A — Execution Fabric core (test results)

Date: 2026-10-08
Gate flag: `HELIX_EXECUTION_FABRIC_V1` (default off)
Profile: `HELIX_EXECUTION_PROFILE` (default `local-only`)
Python: `.venv/bin/python` (3.9.6), `HELIX_MOCK_MODE=1`

## What P3A adds

A vendor-neutral path from a scientific requirement to an executed result:

    requirement → capability → compatible providers → selected provider
                → approved execution → result

Three distinct layers are kept separate on purpose:

- **Capability descriptors** (`backend/config/capabilities/*.yaml`) — static
  facts about what a capability needs and produces, manually curated, merged
  with `_defaults.yaml`. Namespaced ids (`<provider>:<name>`).
- **Execution profile** (`backend/config/profiles/local-only.yaml`) — which
  providers *this deployment* turns on. `local-only` enables `local_compute`,
  `nextflow` (local executor) and the `mock_experimental_provider` (a local
  fake); no cloud. Experimental capabilities are still gated at runtime by the
  approval invariant, not by the profile.
- **Runtime health** — best-effort per-provider flag (defaults healthy).

`CapabilityRegistry.find()` returns only profile-enabled + healthy descriptors;
a disabled provider's capabilities are never returned (but remain visible in the
catalog with `enabled: false`). Alias resolution is deterministic and prefers an
enabled provider.

**Providers** implement a single `ExecutionProvider` Protocol
(`describe_capabilities`, `validate_request`, `estimate`, `prepare_execution`,
`execute`, `get_status`, `retrieve_results`, `cancel`, `provenance`):

- `LocalComputeProvider` — runs existing tools synchronously via
  `backend.main.dispatch_tool`; `client_dedup` idempotency.
- `NextflowProvider` — async; injectable submit/status/results fns; raises
  `ProviderNotConfigured` with no backend wired (real wiring is P3B).
- `MockExperimentalProvider` — deterministic seeded fake lab; `native`
  idempotency (first write wins); refuses to run without an approval.
- `LegacyBrokerProvider` — wraps `dispatch_tool` for data fetches.

**ExecutionFabric** enforces invariants (request↔intent match, approval is
approved + matches intent, experimental requires approval), persists the
`ExecutionRequest` to the ledger *before* the provider call, and on a duplicate
idempotency key re-attaches to the existing run instead of resubmitting. On a
`ProviderTimeout` it re-queries before any retry and records an `unknown` run
rather than blindly resubmitting. Read-only data/advisory capabilities skip the
approval requirement.

**Seam:** `GET /capabilities` exposes the catalog; `ExecutionBroker.execute_tool`
calls `try_fabric_execution(...)` first, which returns `None` (legacy path
unchanged) unless the flag is on *and* the session has an approved,
`READY_TO_EXECUTE` intent whose capability matches the tool.

## Targeted suites

| Suite | Result |
|---|---|
| `tests/unit/backend/execution/test_registry.py` | 5 passed |
| `tests/unit/backend/execution/test_provider_contract.py` | 8 passed |
| `tests/unit/backend/execution/test_fabric.py` | 7 passed |
| `tests/unit/backend/execution/test_broker_delegation.py` | 4 passed |
| `tests/unit/backend/execution/test_capabilities_endpoint.py` | 4 passed |
| `tests/acceptance/test_orch_001_golden_path.py` | 11 passed, 11 skipped (4 P3A steps un-skipped) |

## Flag on/off parity

Broker behaviour is unchanged with the flag on vs off:

- `pytest tests/unit/backend/test_execute_dispatch.py test_execution_broker_policy.py`
  → **39 passed** with `HELIX_EXECUTION_FABRIC_V1` unset.
- Same suite with `HELIX_EXECUTION_FABRIC_V1=1` → **39 passed**.

## Full unit baseline

`pytest tests/unit` → **1037 passed, 1 skipped, 5 failed** in ~115s.

The 5 failures are the same pre-existing, P3A-unrelated failures as the P1/P2
baseline:

- `test_analysis_executor_retries.py::test_execute_analysis_plan_exhausts_attempts_without_success`
- `test_benchmark_refactor_architecture.py::test_approval_command_accepts_prefixed_approve_phrase`
- `test_demo_data_integrity.py::TestS3DataExists::test_all_followup_s3_uris_exist` (needs S3)
- `test_orchestration_modules.py::test_approval_policy_is_action_based`
- `test_phase25_benchmark_turn_fixtures.py::test_turn_03_approve_command_is_recognized`

## Intentional test updates (documented)

- Three trace-propagation tests enumerate the full ledger record-kind set; the
  fabric adds `execution_requests` and `execution_runs`, so their expected count
  maps now include those two kinds (both `0` in a pre-execution loop):
  `test_approval_endpoint`, `test_scientific_plan_persistence`,
  `test_trace_id_propagation`.

## Contracts

`python -m backend.contracts.schema_export --check` → up to date.

## Deferred (not P3A)

- Real Nextflow submit/status/results wiring and the ORCH-001
  `one_real_computational_backend` step (tagged P3B/P3A-nextflow).
- Cloud/AWS adapters (P4). Data providers exist as descriptors + a legacy
  adapter but are not enabled in the `local-only` profile yet.
