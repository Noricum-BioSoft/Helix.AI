# CarolinaCloud — reference hosted external execution provider

**Role in the platform:** reference implementation of a *hosted third‑party scientific execution service* behind the generic `ExecutionProvider` contract (plan Phase 3B). It exists to answer one architectural question — *can such a provider be integrated entirely through the generic contract without leaking provider assumptions into the orchestration core?* — not to make Noricum depend on CarolinaCloud. Everything provider‑specific belongs in this file, `backend/execution/providers/carolina_cloud.py`, `backend/config/capabilities/carolina_cloud.yaml` and `backend/config/profiles/carolina-cloud-dev.yaml`.

Status: **unverified notes from public documentation (2026‑10‑07).** See "Verify before production use".

## What it offers (per public OpenAPI v1.2.0, `https://api.carolinacloud.io/api/schema/`)

Vendored snapshot: `tests/fixtures/carolina_cloud/openapi.yaml`.

| Surface | Endpoints | Relevance |
|---|---|---|
| Managed Nextflow | `POST /api/pipeline-launch/`; `GET/DELETE /api/pipeline-run/{id}/`; `…/console/`; `…/tasks/`; `…/debug-log/`; `GET /api/pipeline-run/` | The execution path Helix uses |
| Compute instances | `/api/instance/*` (VMs/containers; flavors incl. `genomics`, `nextflow-head`; GPUs) | Not used by the adapter |
| Storage | `/api/buckets/` (S3‑compatible object storage) | `outdir` + input staging |
| Billing / limits | `/api/projects/{id}/spend/`, `/api/limits/`, `/api/whoami/` | `validate_request`, cost reconciliation |

Auth: single Bearer API key. "Predefined workflows" are nf‑core pipelines (or any git URL / raw `main.nf`), i.e. the same identifiers Helix already maps in `nextflow_executor.HELIX_TO_PIPELINE`. Their executor plugin (`nf-ccloud`) is also on the Nextflow plugin registry, so the generic `nextflow` provider can target CarolinaCloud from a local head with `-plugins nf-ccloud` (alternative shape; no adapter code).

Published pricing at time of writing: `$0.005/vCPU/hr`, NVMe scratch `$0.0001/GiB/hr`, object storage `$0.009/GiB/mo`, `$0` egress. Treat as input to `estimate()`, re‑read from the rate card, never hard‑coded as truth.

## Mapping onto `ExecutionProvider`

| Method | Call(s) | Notes |
|---|---|---|
| `describe_capabilities()` | static YAML (∩ catalog when P3C exists) | one capability per supported pipeline |
| `validate_request(req)` | `GET /api/buckets/`, `GET /api/limits/` | pipeline known; inputs resolvable; outdir bucket owned; within limits |
| `estimate(req)` | rate card × summed process directives | add inbound transfer cost when inputs are outside CC storage |
| `prepare_execution(req)` | stage inputs to CC bucket; samplesheet via `staged_files` (≤ 5 MiB) | pin `-r <tag>`; build `config` text for resource overrides |
| `execute(prepared)` | `POST /api/pipeline-launch/` → `{id, state, head_uuid}` | **client‑side dedup** on `idempotency_key` (see below) |
| `get_status(handle)` | `GET /api/pipeline-run/{id}/` → `status: active\|closed\|aborted`, `outcome: succeeded\|failed\|null` | polling only; no webhooks in spec |
| `retrieve_results(handle)` | list `outdir` | egress free |
| `cancel(handle)` | `DELETE /api/pipeline-run/{id}/` | idempotent per spec |
| `provenance(handle)` | `…/console/` (verdict, task counts, cost, `nf_run_name`, `nf_session_uuid`, `duration_ms`), `…/tasks/` | `debug-log/` attached as failure Observation |

`PipelineLaunch` fields: `pipeline`, `pipeline_type: ref|script`, `engine_flags[]` (`-resume`, `-profile`, `-r`; `-c`/`-plugins` rejected — platform‑managed), `pipeline_params[]`, `staged_files[{name, content}]`, `outdir` (required), `config` (nextflow.config text), `external_aws_creds?`.

## Idempotency

The spec exposes no idempotency key. Adapter declares `idempotency_support: client_dedup`:

1. Fabric persists `ExecutionRequest{idempotency_key}` before calling `execute()`.
2. Adapter encodes `idempotency_key[:12]` into the Nextflow run name via `engine_flags: ["-name", "helix-<key12>"]`.
3. On timeout/retry, adapter lists `GET /api/pipeline-run/` and matches `nf_run_name` before submitting again; a match re‑attaches instead of relaunching.

Residual risk: a submission that was accepted but whose run has not yet appeared in the list. Keep retry back‑off ≥ the observed listing latency; record both attempts if duplication is detected.

## Policy decisions baked into the adapter

- **Never populate `external_aws_creds`.** It would send an AWS secret to a third party. Outputs always go to CC‑owned buckets; results are pulled back by Helix.
- `security.data_classes_allowed: [public, internal]` until the compliance posture below is confirmed in writing. Sensitive/PHI classes are denied by provider authorization.
- Only CC‑owned buckets as `outdir`; input staging uses the storage layer's `s3_compatible` endpoint (CC buckets are `s3://` URIs on a non‑AWS endpoint — this is why `location_type`/`storage_provider` were generalised).
- Pipelines are launched only by pinned tag (`-r`), never `dev`.

## Alternative shape: local head + `nf-ccloud`

Profile entry for the generic `nextflow` provider: `config: {plugins: [nf-ccloud], executor: ccloud, api_key_env: HELIX_CAROLINA_CLOUD_API_KEY}`. Zero adapter code; Helix host must stay up for the run; loses `console/`/`tasks/`/`debug-log/` provenance. Both shapes may coexist in one profile.

## Verify before production use

Independently confirm with the provider (none of these are asserted by the API docs):

- API version and behaviour against the vendored snapshot (re‑diff `/api/schema/`)
- pricing / rate card and billing granularity
- credential model (API key scoping, rotation, org vs. user keys)
- rate limits and quotas (`/api/limits/`)
- storage semantics: `outdir` retention, bucket ownership, S3‑compatible endpoint details
- compliance posture, BAA availability, data residency (presumably US)
- whether pool workers can read **inputs** from external S3 (spec shows external creds for `outdir` only)
- webhooks vs. polling roadmap
- run listing latency (affects the dedup window above)

Record the answers here with a date; update `capabilities/carolina_cloud.yaml` (`data_classes_allowed`, `regions`) accordingly.

## Smoke test (manual, not CI)

`HELIX_EXECUTION_PROFILE=carolina-cloud-dev HELIX_CAROLINA_CLOUD_API_KEY=… python scripts/smoke_provider.py carolina_cloud nf-core/rnaseq --profile test` → record `console/` cost + duration in `artifacts/test_results/p3b_carolina_cloud_smoke_<date>.md`.
