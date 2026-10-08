---
name: Noricum scientific orchestration platform
overview: Incremental, flag-gated evolution of Helix.AI (reasoning/orchestration) and DataWeaver.AI (scientific context layer) into the Objective → Plan → Security → Capability → Provider → Human Approval → Execution → Evidence → Learning → Decision loop, without rewriting either product and without coupling Noricum to one infrastructure provider.
status: proposed (rev 2 — incorporates NORICUM_IMPLEMENTATION_PLAN_REVIEW_FEEDBACK, 2026-10-07)
depends_on:
  - docs/architecture/NORICUM_PLATFORM_GAP_ANALYSIS.md
  - .cursor/plans/system-structural-fixes.md (Slice A persistence contract)
  - docs/providers/carolina_cloud.md (reference external provider notes)
derived_plans:
  - DataWeaver.AI/.cursor/plans/noricum-scientific-context-layer.md (P5, P6, DW side of P7/P9 — keep in sync with this file)
todos:
  - id: p0-contracts
    content: "P0: contracts (Objective, Plan w/ structured rationale, Assessment+PolicyEnvelope, Recommendation, ProviderAuthorization, ExecutionIntent, Approval, ExecutionRequest w/ idempotency, Evidence), trace_id, domain-ID generation, ExecutionProfile, auth mode guard, invariants module, ORCH-001 skeleton, flags, fixtures"
    status: completed
  - id: p1-objective-plan
    content: "P1 (Helix): ScientificObjective + persisted ScientificPlan + ExecutionIntent + approval endpoint bound to execution_intent_hash"
    status: completed
  - id: p2-security-assessment
    content: "P2 (Helix): Secure Science pre-routing assessment → PolicyEnvelope; DualUseTriageCheck; new checkpoint states"
    status: completed
  - id: p3a-fabric-core
    content: "P3A (Helix): CapabilityDescriptor + static registry + execution profiles + provider factory + Local/Nextflow/MockExperimental providers + ExecutionFabric + conformance suite + broker compat"
    status: completed
  - id: p4-recommender-authz
    content: "P4 (Helix): Execution Recommender over registry + provider-specific authorization against the PolicyEnvelope"
    status: pending
  - id: p3b-external-provider
    content: "P3B (Helix): one real hosted external provider behind the generic contract (reference: CarolinaCloud) — fake client, conformance, async lifecycle, cost/provenance, manual smoke"
    status: pending
  - id: p5-dw-knowledge-model
    content: "P5 (DataWeaver): canonical knowledge tables + relationships + provenance_events; compatibility adapters over legacy Design/Build/Test/File; no legacy migration"
    status: pending
  - id: p6-dw-context-api
    content: "P6 (DataWeaver): Scientific Context API with idempotent PUT-by-domain-id writes; analyses emit AnalysisResult/Artifact/ProvenanceEvent; service auth"
    status: pending
  - id: p7-integration
    content: "P7: Helix KnowledgeStore protocol + DataWeaver client + local ledger store; dual-mode parity on stable IDs"
    status: pending
  - id: p8-learning
    content: "P8 (Helix): EvidenceAssessment (Observation/Interpretation/Decision) + NextDecision → proposed objective (never auto-activated)"
    status: pending
  - id: p9-ui-evals-gate
    content: "P9: UI cards, ORCH-001 green end-to-end, release-gate thresholds, docs"
    status: pending
  - id: p3c-discovery
    content: "P3C (independent, LATER/SHOULD): nf-core catalog + schema ingestion → curated descriptors, deterministic search, CI freshness"
    status: pending
---

# Noricum Scientific Orchestration Platform — Implementation Plan (rev 2)

Gap analysis and architectural decisions D1–D8: `docs/architecture/NORICUM_PLATFORM_GAP_ANALYSIS.md`. Review feedback incorporated in this revision: `~/Downloads/NORICUM_IMPLEMENTATION_PLAN_REVIEW_FEEDBACK.md` (items 1–18; mapping in the last section).

> **Helix owns reasoning and orchestration. DataWeaver owns persistent scientific context and memory. The Execution Fabric connects existing scientific infrastructure. The Secure Science layer governs whether and how execution occurs.**

## Principles

1. **Strangler, not rewrite** — in Helix *and* in DataWeaver. Every phase wraps existing modules; legacy paths and legacy tables stay until parity is proven.
2. **Flag‑gated.** `HELIX_SCIENCE_GATE_V1`, `HELIX_EXECUTION_FABRIC_V1`, `HELIX_EXECUTION_RECOMMENDER_V1`, `HELIX_KNOWLEDGE_STORE=off|local|dataweaver|dual`. Default off until each phase's exit criteria hold.
3. **One write path.** New records are written through one `record_*` helper per type, routed through `KnowledgeStore` from P7 on.
4. **Structured rationale, not reasoning traces.** Persist decisions, evidence refs, assumptions, confidence, constraints and concise explanations as `RationaleItem[]`. Never persist or depend on opaque model reasoning. Applies to plans, recommendations, security outcomes, interpretations and next‑decisions.
5. **IDs originate once, at the domain layer.** Every shared domain object gets its UUID (uuid7) in Helix's orchestration layer and keeps it across local ledger and DataWeaver. Stores never mint IDs for shared objects.
6. **One `trace_id` per scientific loop.** Every artifact of a loop carries the same `trace_id` (`orch_…`).
7. **Fail closed.** Unknown provider → not recommendable; no screening adapter → `REQUIRE_REVIEW`; approval mismatch → no execution; `DENIED`/`WAITING_FOR_SECURITY_REVIEW` → no executable intent; decisions never auto‑execute. Invariants live in code validators and tests, not prompts.
8. **Vendor neutrality.** No provider name appears outside a YAML descriptor, an adapter module, or `docs/providers/`. External providers are reference adapters, never platform dependencies.
9. **Baseline must not regress.** 876 unit tests pass today; 5 pre‑existing failures are tracked, not masked (`artifacts/test_results/unit_baseline_2026-10-07_platform_gap_analysis.md`).

## Target flow (enforced in `HandoffPolicy` and `WorkflowCheckpoint`)

```
ScientificObjective
  → ScientificPlan (persisted, versioned, hashed)
  → Secure Science Assessment  → PolicyEnvelope          [pre‑routing: can this class of action proceed?]
  → Capability resolution → ExecutionRecommendation
  → Provider‑specific Authorization                      [post‑routing: may it proceed via THIS provider under THESE conditions?]
  → ExecutionIntent (immutable digest)
  → HumanApproval (bound to execution_intent_hash)
  → ExecutionProvider (idempotent submission) → ExecutionRun → Artifacts
  → persist to DataWeaver
  → Observation → Interpretation → Decision (terminal; proposed objective stays `proposed`)
```

Agent roles: `IntentDetector → Planner → SecurityAssessor → Recommender → ProviderAuthorizer → [CodeGen] → HumanApproval → Fabric → Visualizer → Learner → Decision`. Users see one "security" step; internally the two checks are distinct.

---

## Phase 0 — Contracts, IDs, trace, flags, invariants, fixtures (Helix; ~3 days)

**Goal:** types first, zero behaviour change. Everything below is additive.

**Cross‑cutting fields on every contract:** `trace_id: str` (`orch_` + uuid7), `schema_version: int`, `created_at`. `backend/contracts/ids.py` — `new_id(kind) -> UUID` (uuid7), `new_trace_id()`; the only place IDs are minted.

Files (new):
- `backend/contracts/rationale.py` — `RationaleItem{evidence_refs: list[str] = [], statement: str, conclusion: str, confidence: float | None}`. Used everywhere a free‑text "reasoning"/"explanation" would otherwise be persisted.
- `backend/contracts/scientific_objective.py` — `ScientificObjective{objective_id, trace_id, version, objective, biological_system?, question, constraints{…, execution?: {allowed_providers[], preferred_provider?, data_residency?}}, desired_evidence[], known_inputs[ArtifactRef], expected_outputs[], success_criteria[], user_context{}, status: proposed|active|answered|abandoned, parent_decision_id?}`.
- `backend/contracts/scientific_plan.py` — `ScientificPlan{plan_id, trace_id, version, objective_id, plan_hash, plan_rationale: list[RationaleItem], steps: [PlanStep + kind: reasoning|computational|experimental, required_data[], required_capabilities[], assumptions[], depends_on[]], security_requirements[], expected_outputs[], status: draft|assessed|recommended|authorized|approved|rejected|executing|completed|failed, ir: plan_ir.Plan}`. Wraps — does not replace — `plan_ir.Plan`. `plan_hash` covers `steps + ir + required_*`, not rationale.
- `backend/contracts/security_assessment.py` — `SecurityOutcome = ALLOW|ALLOW_WITH_APPROVAL|REQUIRE_REVIEW|DENY`; `PolicyEnvelope{allowed_data_classes[], external_execution_allowed: bool, screening_required: bool, human_approval_required: bool, allowed_regions[], allowed_provider_categories[], max_cost_usd?}`; `SecurityAssessment{assessment_id, trace_id, plan_id, plan_hash, assessment_hash, outcome, rationale: list[RationaleItem], applied_policies[{policy_id, version, result}], screening_results[{provider, status, summary, reference_id}], required_approvals[{role, reason}], policy_envelope, audit{actor, timestamp, gate_version, profile}}`.
- `backend/contracts/execution_recommendation.py` — `ExecutionRecommendation{recommendation_id, trace_id, plan_id, step_id?, provider_id, capability_id, candidate_set[provider_id], criteria_scores{suitability, capability_match, inputs_ready, cost, turnaround, throughput, security, privacy, availability, prior_performance, provenance}, weights{}, rationale: list[RationaleItem], alternatives[], warnings[], confidence, legacy_infrastructure?: Literal[...]}`.
- `backend/contracts/provider_authorization.py` — `ProviderAuthorization{authorization_id, trace_id, recommendation_id, assessment_id, provider_id, capability_id, result: AUTHORIZED|DENIED|REQUIRES_SCREENING, checks[{check_id, result, detail}], rationale: list[RationaleItem]}`.
- `backend/contracts/execution_intent.py` — **immutable** `ExecutionIntent{execution_intent_id, trace_id, plan_id, plan_hash, assessment_id?, assessment_hash?, recommendation_id?, authorization_id?, provider_id, capability_id, input_manifest_hash, execution_parameters_hash, execution_intent_hash}` (`frozen=True`; hash computed over all other fields; `assessment_*`/`recommendation_*`/`authorization_*` nullable only while the corresponding flags are off — P1 populates `provider_id` from the legacy `InfraDecision`).
- `backend/contracts/human_approval.py` — `HumanApproval{approval_id, trace_id, execution_intent_id, execution_intent_hash, plan_id, plan_hash, decision: approved|rejected|changes_requested, principal: Principal, timestamp, note?, scope: plan|step[], selected_recommendation_id?}`. `Principal{subject_id, display_name, identity_provider: dev_header|oidc|natural_language, auth_method, roles[], tenant_id?, project_id?, audit_signature?}`.
- `backend/contracts/execution_request.py` — `ExecutionRequest{execution_request_id, trace_id, execution_intent_id, execution_intent_hash, approval_id, idempotency_key, provider_id, capability_id, inputs[], parameters{}, submitted_at?, provider_handle?}`; `idempotency_key = sha256(execution_intent_hash + nonce)`, generated once and **persisted before** any provider call. `ExecutionRun{execution_run_id, trace_id, execution_request_id, status, provider_handle, started_at, finished_at?, cost_actual?, outputs[ArtifactRef]}`.
- `backend/contracts/evidence_assessment.py` — `Observation{observation_id, trace_id, run_ids[], metric, value, unit?, artifact_refs[]}` (no claims), `Interpretation{supports[], contradicts[], uncertainties[], rationale: list[RationaleItem]}`, `NextDecision{decision_id, kind: stop|rerun_analysis|acquire_data|repeat_experiment|new_condition|switch_provider|escalate, rationale: list[RationaleItem], evidence_ids[], proposed_objective?: ScientificObjective}`, `EvidenceAssessment{objective_id, trace_id, run_ids[], observations[], interpretation, decision}`.
- `backend/contracts/provenance.py` — `ProvenanceRecord{event_id, trace_id, event_type, actor{kind: user|agent|system, id}, tool, tool_version, model_version?, agent_id?, agent_version?, prompt_template_id?, prompt_template_version?, prompt_template_hash?, code_commit?, schema_version, parameters{}, execution_backend?, provider_id?, inputs[], outputs[], timestamp}`. Prompts are referenced by template id/version/hash, never stored verbatim.
- `shared/capability_registry.py` — `CapabilityDescriptor{provider, capability_id, category: computational|experimental|data|advisory, inputs[], outputs[], supported_systems[], constraints{}, cost{model, unit_usd_range?}, turnaround{estimated_days?, estimated_minutes?}, integration{api_available, adapter, execution_mode: sync|async, idempotency_support: native|client_dedup|none}, security{screening_required, human_approval_required, data_classes_allowed[], regions[]}, prior_performance{}, provenance{source: manual|generated|provider, catalog_sha?}}`. No `available` field (enablement is a profile concern).
- `backend/execution/providers/base.py` — `ExecutionProvider` Protocol: `describe_capabilities()`, `validate_request(req)`, `estimate(req)`, `prepare_execution(req) -> PreparedExecution`, `execute(prepared) -> ExecutionHandle`, `get_status(handle)`, `retrieve_results(handle)`, `cancel(handle)`, `provenance(handle) -> ProvenanceRecord`, `idempotency_support` property. Adapter docstring must state how `idempotency_key` is honoured.
- `backend/config/execution_profile.py` — `ProviderConfig{id, enabled, adapter?: "module:Class", config{}, credentials?: {source: env|profile|role|api_key|none, env_var?}}`, `ExecutionPolicy{allow_cloud, allowed_regions[], data_residency?, max_cost_usd_per_run?}`, `ExecutionProfile{profile, default_provider, policy, providers[], storage?: reserved}`; loader keyed by `HELIX_EXECUTION_PROFILE` (default `local-only`), `${ENV}` interpolation, adapter import check. Shipped: `profiles/local-only.yaml`, `profiles/aws-dev.yaml`. (Provider‑specific profiles arrive with their adapters.)
- `backend/config/auth_mode.py` — `HELIX_AUTH_MODE=dev_header|oidc` (default `dev_header`), `HELIX_ENV=development|staging|production`. **Startup fails** (`main.py` lifespan) when `HELIX_ENV=production` and `HELIX_AUTH_MODE=dev_header`. `get_principal(request) -> Principal` is the only way code obtains a principal.
- `backend/config/feature_flags.py` — typed accessors.
- `backend/orchestration/invariants.py` — executable invariants (see "Execution‑state invariants"); each is a function raising `InvariantViolation`, called by the services that own the transition.
- `tests/acceptance/test_orch_001_golden_path.py` — **ORCH‑001 skeleton**: the full flow as sequential steps, each step `pytest.skip("phase N")` until implemented; the skip list shrinks each phase and must be empty by P9.
- Tests: `test_platform_contracts.py` (round‑trip, hash stability incl. `execution_intent_hash` and `assessment_hash`, JSON‑schema export to `shared/schemas/contracts/*.json` — copied into DataWeaver as `backend/app/schemas/contracts/` — `RationaleItem` required where `reasoning` used to be), `test_ids_and_trace.py`, `test_execution_profile.py`, `test_auth_mode_guard.py` (production + dev_header → startup error), `test_invariants.py`.

Exit: new tests green; imports clean; no production path touched; ORCH‑001 collected with all steps skipped.

### Execution‑state invariants (code‑level, grown per phase)

| Invariant | Enforced in | Phase |
|---|---|---|
| `DENIED` → no transition to `READY_TO_EXECUTE`/`EXECUTING` | `WorkflowCheckpoint.transition` | P2 |
| `WAITING_FOR_SECURITY_REVIEW` → `ExecutionIntent` cannot be created | `intent_builder` | P2 |
| `approval.plan_hash == current plan_hash` | `approve_plan()` | P1 |
| `approval.execution_intent_hash == request.execution_intent_hash` | `ExecutionFabric.run()` | P1 (ledger) / P3A (fabric) |
| `ExecutionIntent.provider_id ∈ recommendation.candidate_set` | `intent_builder` | P4 |
| `envelope.screening_required` → no run without accepted screening result or authorized review override | `ProviderAuthorizer`, `ExecutionFabric.run()` | P2/P4 |
| `HumanApproval.decision != approved` → `ExecutionFabric.run()` raises | `ExecutionFabric.run()` | P3A |
| `category == experimental` → `HumanApproval` required regardless of policy | `ProviderAuthorizer` | P4 |
| Same `idempotency_key` → at most one provider submission | `ExecutionFabric` + `JobManager` ledger | P3A |
| `NextDecision.proposed_objective.status == proposed` until explicit user action | `store_objective`, `/execute` | P8 |
| All records of a loop share one `trace_id` | `record_*` helpers | P1→ |

### Configurable execution providers (cross‑cutting: P0, P3A, P4)

**Problem (today).** AWS service names are the type system: `Literal["Local","EC2","EMR","Batch","Lambda"]` in `contracts/infra_decision.py`, `execution_spec.py` and the infra‑agent prompt; `infra == "EMR" → async` string dispatch in `execution_broker.py`; `available: … # HELIX_USE_EC2` baked into `environment_capabilities.yaml`; `AWS_REGION` read in ~15 modules, `boto3` imported in 12; `location_type: Literal["S3","Local","URL","Unknown"]`; `cost_heuristics.yaml` priced for us‑east‑1.

**Design.** Three layers kept separate — **capability descriptors** (static, what a provider can do), **execution profile** (per deployment: enabled providers, how to reach them, policy limits), **runtime health** (adapter ping at startup). Candidates = `profile.enabled ∩ registry.find(capability) ∩ envelope.allowed`. Precedence (top wins): PolicyEnvelope → profile `policy` → `objective.constraints.execution` → recommender scoring → human picks an alternative at approval (`selected_recommendation_id`, must be in `candidate_set`). Adapters receive a `ProviderConfig` and never read env; credentials resolve via `credentials.source`; no secrets in YAML. Nextflow executors (awsbatch / google‑batch / azurebatch / k8s / slurm / vendor plugins) selected by the profile are the primary multi‑cloud lever; bespoke adapters are for what Nextflow doesn't cover (EMR/Spark, legacy EC2) and for hosted job‑service APIs (P3B). `integration.execution_mode` replaces string dispatch. `location_type` gains `ObjectStore` + `storage_provider: s3|s3_compatible|gcs|azure_blob` inferred from scheme **and** the profile's known endpoints — `s3://` is not assumed to mean AWS. Cost tables per provider under `backend/config/cost/{provider_id}.yaml`. `InfraDecision.infrastructure` stays as `legacy_infrastructure` until P4 removes the string dispatch.

Out of scope, flagged: artifact *storage* abstraction (`history_manager` is S3‑only); profile reserves a `storage:` block.

---

## Phase 1 — Objective, persisted plan, execution intent, explicit approval (Helix; ~4 days)

**1.1 Objective.** `backend/orchestration/objective_builder.py::build_objective(command, session_context, intent_result) -> ScientificObjective` (LLM; deterministic in mock mode). Mints `objective_id` + `trace_id`; stored on `WorkflowCheckpoint.objective` and as ledger artifact `type="scientific_objective"`.

**1.2 Plan artifact.** In the staging branch of `/execute`: `ScientificPlan.from_plan_ir(plan_ir, objective, rationale: list[RationaleItem])` — rationale items are produced by the planner as structured output (statement/conclusion/evidence_refs), *not* by persisting the planner's prose. `record_plan()` → ledger artifact `sessions/{sid}/plans/{plan_id}.v{n}.json`. `checkpoint.pending_plan` keeps the IR dict (back‑compat) and gains `pending_plan_id`/`pending_plan_hash`/`trace_id`. Revisions create new versions with `supersedes`.

**1.3 Execution intent (P1 shape).** `backend/orchestration/intent_builder.py::build_intent(plan, infra_decision, inputs) -> ExecutionIntent` with `provider_id = legacy:{InfraDecision.infrastructure}`, `capability_id = local_compute:{tool_name}`, `input_manifest_hash` over resolved input URIs + sizes + content hashes where known, `execution_parameters_hash` over step arguments; assessment/recommendation/authorization fields `None` while their flags are off. Persisted as ledger artifact `type="execution_intent"`; `checkpoint.pending_execution_intent_id/hash`.

**1.4 Approval.** Router `backend/api/approvals.py`: `POST /session/{sid}/intents/{intent_id}/approve|reject|request-changes`. Principal via `get_principal()` (dev_header mode: `X-Helix-User`; must be present). Verifies `execution_intent_hash` **and** `plan_hash` against the staged values (409 on either mismatch); writes `HumanApproval` via `record_approval`; transitions checkpoint. The NL approval path (`is_approval_command`) calls the same `approve_intent()` service with `identity_provider="natural_language"` and the same hash checks — a plan or intent that changed since staging cannot be approved by "yes".

**1.5 Handoff.** `AgentRole.HUMAN_APPROVAL` in `agent_registry.py`/`HandoffPolicy`; `Infra → HumanApproval → Broker` legal; `Infra → Broker` only when `approval_policy` says none needed (read‑only tools). Broker, when invoked, re‑verifies `approval.execution_intent_hash == intent.execution_intent_hash` (invariants module) even before the fabric exists.

Tests: `test_objective_builder.py`, `test_scientific_plan_persistence.py` (version/hash/supersedes; rationale is `RationaleItem[]`; bundle includes plan), `test_execution_intent.py` (hash changes when provider/inputs/params change; frozen), `test_approval_endpoint.py` (approve/reject; intent‑hash mismatch 409; plan‑hash mismatch 409; principal recorded; NL path yields identical record and respects mismatch), `test_handoff_policy.py` update, `test_trace_id_propagation.py`. ORCH‑001: steps Objective → Plan → Intent → Approval un‑skipped. Re‑evaluate the 3 approval‑classifier baseline failures (should become mockable via the endpoint).

Exit: `/execute` unchanged with flags off; with `HELIX_KNOWLEDGE_STORE=local` every staged plan has objective, plan, intent and (when approved) approval records sharing a `trace_id`; changing provider or inputs after staging invalidates the approval.

---

## Phase 2 — Secure Science pre‑routing assessment (Helix; ~5 days)

Answers **"can this class of action proceed?"** before any provider is chosen. Provider‑specific checks are deliberately *not* here (see P4).

Package `backend/security/`:
- `assessor.py` — `SecureScienceAssessor(checks).assess(plan, objective, context) -> SecurityAssessment` incl. `PolicyEnvelope`. Combine rule `DENY > REQUIRE_REVIEW > ALLOW_WITH_APPROVAL > ALLOW`; envelope is the intersection of per‑check envelopes; `assessment_hash` over outcome + envelope + applied policies. Persisted via `record_assessment`.
- `checks/base.py` — `SecurityCheck` Protocol: `check_id`, `version`, `applies_to(plan)`, `evaluate(plan, objective, context) -> CheckResult{outcome, rationale: RationaleItem[], policies[], screening_results[], required_approvals[], envelope_constraints}`.
- `checks/action_classification.py` — step `kind`/capability category → risk category (`read_only`, `compute`, `external_compute`, `experimental`, `synthesis`); experimental/synthesis → `human_approval_required`, `screening_required` in the envelope.
- `checks/data_sensitivity.py` — reuses `upload_intake_policy` results (`sensitivity_class`, `scan_flags`) → `allowed_data_classes`, `external_execution_allowed`, `allowed_regions` from `backend/config/security/policies/*.yaml` (declarative; `docs/SAFETY_POLICY.md` categories become rules).
- `checks/dual_use_triage.py` — **`DualUseTriageCheck`**: keyword/entity triage over objective + plan text. Outputs only `ALLOW` or `REQUIRE_REVIEW`. Module docstring and `docs/SAFETY_POLICY.md`: *"conservative routing to manual review; not a biological‑risk classifier and not a substitute for specialised screening."* Never emits safe/unsafe labels.
- `checks/sequence_screening.py` — `SequenceScreeningAdapter.screen(sequences) -> ScreeningResult`; `MockScreeningAdapter` (`provider="mock"`), `NotConfiguredAdapter` → `REQUIRE_REVIEW` whenever `screening_required`. Real adapters (IBBIS Common Mechanism, SecureDNA, vendor APIs) are documented extension points.
- `checks/identity.py` — principal present and role allowed for the risk category (via `get_principal`); missing → `REQUIRE_REVIEW`.
- `checks/manual_review.py` — explicit `requires_manual_review` on plan/objective.

Wiring: `/execute` staging, after plan, before infra: `DENY` → checkpoint `DENIED` (no intent may be built); `REQUIRE_REVIEW` → `WAITING_FOR_SECURITY_REVIEW` + `POST /session/{sid}/assessments/{id}/review` (reviewer principal, `backend/api/security.py`); `ALLOW_WITH_APPROVAL` → envelope `human_approval_required=true`; `ALLOW` → continue. `ExecutionIntent` now records `assessment_id/hash`. `HandoffPolicy`: `SECURITY_ASSESSOR` between `PLANNER` and `INFRA`. `_emit_policy_audit_event` becomes ledger write + log, with `trace_id`.

Tests: unit per check; combination/envelope‑intersection matrix; `test_execute_security_flow.py` (mock): DENY blocks intent creation; REQUIRE_REVIEW survives reload; ALLOW_WITH_APPROVAL forces approval; no‑adapter → REQUIRE_REVIEW for screening‑required; `DualUseTriageCheck` never returns DENY/ALLOW_WITH_APPROVAL. Invariants rows for P2 enabled. ORCH‑001: Assessment step un‑skipped.

Exit: with `HELIX_SCIENCE_GATE_V1=1`, 100% of staged plans have a `SecurityAssessment` with envelope; eval `security_gate_no_bypass` in `benchmarks/cases/routing_safety/`.

---

## Phase 3A — Execution Fabric core (Helix; ~6 days) — MUST

**Goal:** scientific requirement → capability → compatible providers → selected provider → approved execution → result, with Local, Nextflow and MockExperimental only. Manually curated descriptors.

- **Registry.** `backend/config/capabilities/`: `local_compute.yaml` (generated from `tool_schemas.py` + `environment_capabilities.yaml` Local), `nextflow.yaml` (capabilities = current `HELIX_TO_PIPELINE` entries + `nf-core/rnaseq`, hand‑written; executor from profile), `mock_experimental_provider.yaml` (`protein_expression`, `activity_assay`; `human_approval_required`, `screening_required`, `turnaround.estimated_days: 10`, `cost.model: per_sample`), `ncbi_entrez.yaml`/`uniprot.yaml` (`category: data`), `_defaults.yaml`. `backend/execution/registry.py::CapabilityRegistry(profile)` — validate, index, `find(capability_id, constraints)` → profile‑enabled + healthy only; `GET /capabilities` (`backend/api/capabilities.py`) with `enabled`/`healthy`. Planner resolves `required_capabilities` by alias, fallback `local_compute:{tool_name}`.
- **Providers** (`backend/execution/providers/`): `local_compute.py` (wraps `_tool_executor`/`sandbox_executor`; `idempotency_support: client_dedup` via ledger), `nextflow.py` (wraps `nextflow_executor` + `JobManager.create_nextflow_job`; `executor/queue/workdir_uri/plugins/-c` from `ProviderConfig`; `client_dedup` keyed on `-name {idempotency_key[:12]}` + ledger), `mock_experimental.py` (in‑memory lab; synthetic `Measurement`s flagged `provider="mock"`; refuses to run without an approved `HumanApproval` for the exact intent; **`native` idempotency**: same key → same fake run), `legacy_broker.py` (wraps `ExecutionBroker.execute_tool` for unmapped tools). `factory.py::build_providers(profile)`; no `boto3` import outside `aws_*.py`.
- **Fabric.** `backend/execution/fabric.py::ExecutionFabric.run(intent, approval, request) -> ExecutionRun`: invariants (approval approved; `approval.execution_intent_hash == intent.hash == request.execution_intent_hash`; experimental ⇒ approval present); **idempotency**: look up `request.idempotency_key` in the ledger — if a submission exists, return/attach to it; otherwise persist the request *then* `validate → prepare → execute`; on timeout after `execute`, re‑query by key/handle before any retry; `provenance()` recorded with `trace_id`. `JobManager` gains `provider_id`, `provider_handle`, `idempotency_key`, `trace_id`. `ExecutionBroker.execute_tool` delegates when `HELIX_EXECUTION_FABRIC_V1=1`.
- **AWS adapters deferred:** `aws_emr.py`/`aws_ec2.py` remain reachable through `legacy_broker.py` in 3A; they move behind the contract in P4 alongside the dispatch cleanup (`aws-dev` profile).

Tests: conformance suite `tests/unit/backend/execution/test_provider_contract.py` parametrised over providers (fake Nextflow); profile tests (`local-only` lists no cloud providers; disabled never returned); fabric: approval binding, fail‑closed on `decision != approved`, **duplicate `execute()` with same key → one submission** (simulated timeout), experimental without approval raises; parity: `test_execution_broker_policy.py` + `test_execute_dispatch.py` flag on/off; mock‑lab path (plan → assessment `ALLOW_WITH_APPROVAL` → intent → approve → execute → measurements). ORCH‑001: Capability → Provider → Execution → Result steps un‑skipped (Local + Mock).

Exit: flag‑on parity on unit + `tests/workflows`; `GET /capabilities` lists local, nextflow, mock_experimental, data providers; mock‑lab scenario green; idempotency test green.

---

## Phase 4 — Execution Recommender + provider‑specific authorization (Helix; ~4 days) — MUST

- **Recommender.** `infrastructure_decision_agent.py::recommend_execution(step, registry, objective, assessment) -> ExecutionRecommendation`: candidates = `registry.find()` ∩ `envelope` (data classes, `external_execution_allowed`, regions, categories) ∩ profile policy ∩ `objective.constraints.execution`; deterministic criteria scoring (`backend/config/execution_recommender_weights.yaml`); `rationale: RationaleItem[]` generated from the score table (LLM may phrase statements; items must reference criteria ids; mock mode templated). `candidate_set` recorded. `legacy_infrastructure` derived for compute providers.
- **Provider‑specific authorization.** `backend/security/provider_authorizer.py::authorize(recommendation, assessment, descriptor, context) -> ProviderAuthorization`: provider accepts the data class; supports the capability; satisfies residency/regions; `screening_required` satisfied (accepted `ScreeningResult` or authorized review override); descriptor `human_approval_required` honoured; recommendation still complies with the envelope (re‑check — the envelope is authoritative). `DENIED` → back to recommender for next candidate or `WAITING_FOR_SECURITY_REVIEW`. `ExecutionIntent` now carries `recommendation_id` + `authorization_id`, and `provider_id ∈ candidate_set` is enforced.
- **Dispatch cleanup.** `execution_broker.py` reads `descriptor.integration.execution_mode`; `aws_emr.py`/`aws_ec2.py` behind the contract (`aws-dev` profile; cost from `backend/config/cost/aws_*.yaml`); `environment_capabilities.yaml` becomes the generated source for `aws_*`/`local_compute` descriptors; `DatasetSpec.location_type` generalised (`ObjectStore` + `storage_provider`).
- UI payload: `recommendation` + `authorization` in the `/execute` plan card; ledger artifacts `execution_recommendation`, `provider_authorization`.

Tests: scoring determinism; envelope excludes providers; profile `allow_cloud: false` → local/nextflow‑local only; objective `allowed_providers` respected; authorization denies on data class / region / missing screening; intent with provider outside `candidate_set` rejected; rationale items reference every criterion; `check_plan_not_mutated` holds; `test_infra_decision_validation.py` unchanged. ORCH‑001: Recommendation → Authorization steps un‑skipped.

Exit: every staged step has recommendation + authorization; `HELIX_EXECUTION_RECOMMENDER_V1=1` default on in dev.

---

## Phase 3B — One real hosted external provider (Helix; ~3 days) — SHOULD

**Architectural question:** can a hosted third‑party scientific execution service be integrated *entirely* through the generic `ExecutionProvider` contract, without provider assumptions leaking into the orchestration core?

**Reference implementation:** CarolinaCloud managed‑Nextflow API. All provider‑specific detail (endpoints, payload mapping, pricing, storage semantics, policy decisions, open questions to verify) lives in **`docs/providers/carolina_cloud.md`**, not here.

Deliverables: `backend/execution/providers/carolina_cloud.py` + `profiles/carolina-cloud-dev.yaml` (`credentials.source: api_key`); `capabilities/carolina_cloud.yaml` (descriptor; `execution_mode: async`, `idempotency_support: client_dedup`, `data_classes_allowed: [public, internal]` until compliance is verified); `FakeCarolinaCloudClient` built from the vendored OpenAPI snapshot `tests/fixtures/carolina_cloud/openapi.yaml` + recorded responses; async lifecycle via `JobManager` poller; cost/provenance mapped into `ProvenanceRecord` (`cost_actual`, per‑task records); adapter never forwards Helix cloud credentials to the provider.

Tests: conformance suite passes with the fake; idempotency (re‑submit after simulated timeout → one remote run via client‑side dedup on `idempotency_key`); authorization denies `sensitive` data class; one **manual** smoke run (`nf-core/rnaseq -profile test`) recorded in `artifacts/test_results/` — not CI.

Success criterion: zero `carolina`/provider‑specific identifiers outside `providers/carolina_cloud.py`, its YAML, its profile and `docs/providers/`. ORCH‑001 "one real computational backend" may be satisfied by Nextflow‑local *or* this provider.

Before any production use, independently verify the list in `docs/providers/carolina_cloud.md` (API version/behaviour, pricing, credentials, limits, storage semantics, compliance/BAA, residency, retention, input‑access semantics, webhooks vs polling). Unverified assumptions never enter core contracts.

---

## Phase 3C — Capability Discovery (independent; ~2 days) — SHOULD / LATER

Not required to prove the loop; the existing pipelines + `nf-core/rnaseq` are enough for MVP validation. Develop on a branch and merge when useful.

Build‑time, curated: `scripts/discover_capabilities.py` + `backend/execution/discovery/{nf_core.py, provider_self_describe.py}` → candidate descriptors (`capability_id: nf_core:{name}`, inputs from `nextflow_schema.json`, pinned tagged release, `components`, `provenance{source: nf_core, catalog_sha}`) → `capabilities/allowlist.yaml` promotion → committed `nf_core_generated.yaml` (PR‑reviewed); unpromoted candidates git‑ignored. Rules: tagged releases only, skip archived/schema‑less; discovery never produces experimental/synthesis capabilities. Deterministic retrieval `CapabilityRegistry.search(text, k)` over `_index.json` (mock mode: hashed bag‑of‑words); Tool Generator Agent becomes explicit last resort emitting candidates, not registrations. `--diff` CI freshness check. Fixtures already vendored: `tests/fixtures/nf_core/pipelines.json` + `schemas/{rnaseq,sarek,ampliseq}.json`. Later sources (WorkflowHub, Dockstore, bio.tools) are extension points only.

Tests: mapping, exclusions, allowlist exactness, `search()` ranks `nf_core:rnaseq` first for a bulk RNA‑seq DE query, `HELIX_TO_PIPELINE` ≡ registry alias view, Tool Generator not invoked above threshold.

---

## Phase 5 — DataWeaver knowledge model v1 (~5 days) — MUST

> Executed in the DataWeaver.AI repo. Working plan there: `DataWeaver.AI/.cursor/plans/noricum-scientific-context-layer.md` (same phase numbers; this section is the normative summary, that file carries repo‑level detail). Until a shared contracts package exists (revisit at P7), DataWeaver validates against JSON‑schema snapshots exported from `backend/contracts/` by `test_platform_contracts.py`.

**Strangler inside DataWeaver:** introduce the canonical model; expose legacy data through compatibility adapters; **do not migrate or retire legacy tables in this phase.**

Package `backend/app/models/knowledge/` (UUID PKs **supplied by the caller**, `trace_id`, `schema_version`, `attrs JSONB`, timestamps, `created_by`): `scientific_objectives`, `hypotheses`, `experiment_plans` (Helix `ScientificPlan` JSON + hash), `experiments`, `experimental_conditions`, `samples`, `biological_entities`, `artifacts`, `measurements`, `observations`, `analysis_results`, `conclusions`, `decisions`, `capabilities` (+ `source`, `catalog_sha`, `pinned_release`, `promoted_by/at`), `execution_providers`, `execution_intents`, `execution_requests` (+ `idempotency_key` unique), `security_assessments` (+ `policy_envelope JSONB`), `provider_authorizations`, `human_approvals`, `execution_runs`, `provenance_events` (fields per `ProvenanceRecord`, incl. agent/prompt‑template/code_commit/schema_version), `relationships` (polymorphic edge: `subject_type, subject_id, predicate, object_type, object_id, confidence?, provenance_event_id?, attrs`; unique 5‑tuple; both‑direction indexes; predicates v1: `motivated, produced, analyzed_by, supports, contradicts, based_on, generated, derived_from, measured_on, part_of, executed_by, assessed_by, authorized_by, approved_by, supersedes`).

Compatibility adapters (`backend/app/services/compat/`): `Design/Build → BiologicalEntity`, `Test → Measurement`, `File → Artifact`, `FileRelationship → relationships(derived_from)` as read adapters (SQL views where the dialect allows, Python adapters otherwise) so the Context API can serve legacy records without copying them. Existing `/api/bio/*` endpoints untouched. `workflow_phases` tables: leave as‑is, document as legacy.

Migration: one additive Alembic revision after `f7c123bd3c74` (new tables only); `database.create_tables()` registers metadata. Legacy data migration and table retirement → **LATER** (after P7 parity is proven).

Tests (SQLite + optional Postgres compose): round‑trips with caller‑supplied IDs (duplicate ID → integrity error, not silent new row), edge uniqueness, recursive `get_provenance(id, depth)`, compat adapters return legacy rows as canonical shapes, `alembic upgrade head` on empty and seeded DBs, all legacy tests still green.

Exit: new tables present; legacy endpoints unchanged; compat views serve seeded legacy data.

---

## Phase 6 — DataWeaver Scientific Context API (~5 days) — MUST

Router package `backend/app/api/context/` (prefix `/api/context`; `X-Service-Key` + `project_id` scoping; `backend/app/auth.py`).

**Writes are idempotent upserts keyed by caller‑supplied domain IDs:** `PUT /objectives/{objective_id}`, `PUT /plans/{plan_id}`, `PUT /assessments/{assessment_id}`, `PUT /authorizations/{id}`, `PUT /intents/{execution_intent_id}`, `PUT /approvals/{approval_id}`, `PUT /executions/{execution_request_id}` (+ `PATCH` status/results), `PUT /runs/{execution_run_id}`, `PUT /observations/{id}`, `PUT /analysis-results/{id}`, `PUT /conclusions/{id}`, `PUT /decisions/{id}`, `PUT /provenance-events/{event_id}`, `PUT /relationships` (natural key). Same body twice → 200, no duplicate; conflicting body for same id → 409 (immutable types: intent, approval, assessment) or new version (versioned types: objective, plan). Every write creates its `ProvenanceEvent` and edges in one transaction and carries `trace_id`.

Reads: `GET /objectives/{id}/context` (single payload for Helix), `/hypotheses`, `/evidence`, `/experiments`, `/measurements`, `/provenance/{type}/{id}?depth=`, `/decisions`, `/relationships`, and `GET /traces/{trace_id}` (everything in one loop, ordered).

Analysis integration: `DataAnalyzer`, `simple_visualizer`, `IntelligentMerger`, `DataQAService` gain `persist(project_id, inputs[]) -> AnalysisResult + Artifact(s) + ProvenanceEvent(tool, version, parameters, code_commit)`; `bio_matcher` endpoints call it when `project_id` supplied. Ingest `POST /api/context/datasets/ingest` (CSV/Excel); other formats as extension points.

Also: declare `plotly`, `openpyxl`, `numpy`; drop unused `celery`/`redis`; `Dockerfile` + CI running `backend/tests`; OpenAPI snapshot `docs/openapi.yaml`.

Tests: per‑endpoint (TestClient + SQLite); idempotency (double PUT, conflicting PUT); `/traces/{trace_id}` returns the full loop; contract test vs `docs/openapi.yaml`; analysis persistence links result to input via `analyzed_by`.

Exit: OpenAPI snapshot committed; Helix fixture `tests/fixtures/dataweaver_openapi.json` generated from it.

---

## Phase 7 — Helix ↔ DataWeaver integration (~4 days) — MUST

- `backend/knowledge/store.py` — `KnowledgeStore` Protocol (`get_objective_context`, `get_related_hypotheses`, `get_available_evidence`, `get_prior_experiments`, `get_measurements`, `get_provenance`, `get_prior_decisions`, `get_trace`, `store_*` for every contract — **all `store_*` take the pre‑minted domain object and never return a new id**).
- `backend/knowledge/local_ledger_store.py` over `history_manager`; `backend/knowledge/dataweaver_client.py` (httpx, retries with idempotent PUTs, typed from the OpenAPI fixture); `get_knowledge_store()` by `HELIX_KNOWLEDGE_STORE`; `dual` writes both and **asserts id/hash parity** per object (mismatch → logged `ParityViolation` artifact, never silent).
- `context_builder.py` pulls `get_objective_context()` into the planner prompt (bounded) when an objective id is known.
- `record_*` helpers route through `KnowledgeStore`; `trace_id` propagated.

Tests: Protocol conformance for both stores; recorded‑response contract tests (`respx`); dual mode: identical object ids, hashes and counts in both stores for a full mock loop; retry of a failed `PUT` creates no duplicate. ORCH‑001: "Persist to DataWeaver" + "IDs stable across stores" un‑skipped.

Exit: closed‑loop mock‑lab demo persists the whole loop in DataWeaver; `GET /api/context/traces/{trace_id}` returns objective → … → measurements; `provenance/measurement/{id}` reaches the objective.

---

## Phase 8 — Learning and decision (Helix; ~4 days) — MUST

- `backend/learning/evidence_assessor.py::assess(objective, run_ids, store) -> EvidenceAssessment`: observations extracted deterministically where possible (counts, thresholds from `success_criteria`) and schema‑constrained to contain no claims; interpretation and decision are LLM with Pydantic‑enforced `RationaleItem[]` + explicit uncertainty lists; `NextDecision.evidence_ids` must reference persisted observation/artifact ids (validator); mock‑mode templated.
- `Learner`, `Decision` roles in `HandoffPolicy` (Visualizer → Learner → Decision; Decision terminal). `proposed_objective` stored with `status=proposed`; activation requires an explicit user turn (invariant).
- "Decision Packet" = `EvidenceAssessment` + plan + assessment + authorization + approval + provenance links, via `bundle_generator`.
- Provider failure artifacts (e.g. a provider's debug report) attach as failure Observations.

Tests: observation ≠ interpretation (schema, not prose); decision kinds coverage; `evidence_ids` must resolve; proposed objective never transitions without user action; packet in bundle. ORCH‑001: Observation → Interpretation → Decision steps un‑skipped.

Exit: eval family `closed_loop_mock_experiment` passes with rationale citing evidence ids.

---

## Phase 9 — UI, ORCH‑001, release gate, docs (~4 days) — MUST

- Frontend: `SecurityAssessmentCard` (outcome + envelope + rationale items), `ExecutionRecommendationCard` (criteria table, alternatives, authorization result), `ApprovalCard` (principal, intent hash, approve/reject via new endpoint, alternative selection), `EvidenceDecisionCard`. Extend `helixApi.ts` + tests. No redesign.
- **ORCH‑001 green end‑to‑end** with deterministic fixtures in CI (Local + MockExperimental; Nextflow‑local where the runner has Nextflow; external provider via fake).
- Benchmarks: suites `security_gate`, `execution_fabric_parity`, `closed_loop`, `orch_001`; thresholds in `benchmarks/release_thresholds.yaml` (`security_gate.min_pass_rate: 1.0`, `orch_001.min_pass_rate: 1.0`, `closed_loop.min_pass_rate: 0.9`); `release_readiness.json` includes `platform_phase_flags` and `auth_mode`.
- Docs: `docs/architecture/ROUTING_FLOW.md`, `BACKEND_DATAFLOW.md`, `agents/agent-responsibilities.md`, `handoff-policy.md` (new roles), `docs/SAFETY_POLICY.md` (executable policies; triage disclaimer), `docs/providers/README.md` (how to add a provider), DataWeaver `docs/API.md`.

Exit: full release gate with all flags on; readiness `true` or blockers documented.

---

## ORCH‑001 — Scientific Orchestration Golden Path (acceptance test)

`tests/acceptance/test_orch_001_golden_path.py`, created in P0 as a skipped skeleton, fully green by P9. Flow: Objective → persisted Plan → Security Assessment (envelope) → Capability resolution → Recommendation → Provider Authorization → ExecutionIntent → Human Approval → ExecutionProvider → Result → persist to DataWeaver → Observation → Interpretation → Decision.

Assertions: every major state persisted; all records share one `trace_id`; domain ids identical across local and DataWeaver stores; plan/assessment/intent hashes stable across re‑serialisation; approval binds to the exact intent (mutate provider or inputs → approval invalid); no provider executes without authorization; no experimental execution without approval; retrying `execute()` after a simulated timeout creates no duplicate work; provenance reaches from the final decision to the original inputs; DataWeaver `GET /traces/{trace_id}` returns the loop; Observation contains no claims; Decision references evidence ids; proposed objective remains `proposed`; MockExperimental path works; one real computational backend works; runs on deterministic fixtures in CI.

---

## Sequencing

```
P0 ─► P1 ─► P2 ─► P3A ─► P4 ─► P3B ─► P7 ─► P8 ─► P9
             │                          ▲
             └── (DataWeaver) P5 ─► P6 ─┘
P3C (discovery) — independent branch, merge when useful
```

P5/P6 run in parallel with P2–P4 by a second developer/agent; P7 needs both tracks.

**Schedule guidance (not calendar guarantees).** Engineering estimate ≈ 45 dev‑days single‑track. Realistic MVP golden path (ORCH‑001 green): **8–12 weeks**, driven by repository condition, integration friction, provider behaviour, migrations, security semantics and state‑management complexity. Hardening / pilot readiness: **+4–8 weeks** (real user feedback, provider failures, retry/idempotency edge cases, data compatibility, OIDC, auditability, operational tooling, benchmark stabilisation). The dominant risk is not the Pydantic models; it is proving that distributed scientific workflows stay reliable, auditable, secure and scientifically coherent across real external systems.

## Priorities

**MUST (prove the thesis):** structured ScientificObjective; persisted ScientificPlan; structured rationale; pre‑routing Secure Science assessment + PolicyEnvelope; capability registry; execution profile; Local, Nextflow, MockExperimental providers; execution recommender; provider‑specific authorization; ExecutionIntent + immutable approval binding; idempotent submission; DataWeaver knowledge model + Context API; Helix `KnowledgeStore`; provenance (incl. agent/prompt/code/schema versions); trace ids; Observation/Interpretation/Decision separation; ORCH‑001; benchmark/release gate; dev‑header auth cannot run in production.

**SHOULD (pilot readiness):** one real external provider adapter (P3B); OIDC principal; capability search (P3C); provider cost estimates; richer security plugin integrations; connector conformance automation; provider health monitoring; real dual‑write verification against a deployed DataWeaver.

**LATER:** generalised nf‑core ingestion beyond the allowlist; WorkflowHub/Dockstore discovery; native GCP/Azure Batch adapters; capability marketplace; multi‑tenant execution profiles; artifact‑storage abstraction; legacy DataWeaver data migration and table retirement; provider performance learning; autonomous next‑objective activation.

## Verification per phase (standing)

1. `HELIX_MOCK_MODE=1 pytest tests/unit -q` — no new failures vs baseline (876 pass).
2. `pytest tests/acceptance/test_orch_001_golden_path.py` — skip list only shrinks.
3. Targeted integration for phases touching `/execute`; `tests/workflows` + `tests/evals` when planner/routing/approval behaviour changes (P1, P2, P4, P8).
4. DataWeaver: `cd backend && pytest tests/` (SQLite) and Postgres compose before merging migrations.
5. Artifacts: `artifacts/test_results/<phase>.md`, `artifacts/failure_summaries/` on failure, worklog entry with root cause for any defect fixed.

## Known baseline issues to carry

- 5 pre‑existing unit failures (3 approval‑classifier mock‑mode, 1 S3‑network, 1 sandbox import) — P1 expected to fix the approval cluster; the other two remain documented.
- Commit or stash unrelated working‑tree changes before P0 so baselines are attributable.

## Open decisions (defaults stated)

1. **Principal source** — dev: `HELIX_AUTH_MODE=dev_header` (`X-Helix-User`), refused when `HELIX_ENV=production`; pilot: OIDC (`Principal.identity_provider=oidc`, stable `subject_id`, verified roles, tenant/project); production adds signed audit metadata. P1 not blocked on OIDC.
2. **DataWeaver service vs. library** — default: separate service + HTTP client + local ledger fallback (D2).
3. **`workflow_phases` tables** — default: leave in place, document as legacy; retirement is LATER.
4. **Reference external provider** — resolved: CarolinaCloud as *reference adapter* (P3B), details in `docs/providers/carolina_cloud.md`; not a platform dependency.
5. **Multi‑cloud path** — default: Nextflow executors selected by the execution profile; hosted job‑service providers via P3B‑style adapters; native GCP/Azure adapters LATER.
6. **Profile granularity** — default: one per deployment via `HELIX_EXECUTION_PROFILE`, narrowed per objective; per‑tenant profiles LATER.

## Review‑feedback mapping (rev 2)

| # | Feedback | Where addressed |
|---|---|---|
| 1 | Structured rationale, not `reasoning` | Principle 4; `rationale.py`; all contracts; P1.2, P4, P8 |
| 2 | Split pre‑routing assessment from provider authorization | Target flow; P2 (`PolicyEnvelope`), P4 (`provider_authorizer.py`) |
| 3 | Bind approval to immutable execution intent | `execution_intent.py`; P1.3–1.4; invariants |
| 4 | Idempotent external execution | `execution_request.py`; `idempotency_support`; P3A fabric; P3B |
| 5 | Dev header auth cannot reach production | `auth_mode.py` startup guard; decision 1 |
| 6 | Split Phase 3 | P3A core / P3B external / P3C discovery; sequencing |
| 7 | External provider = reference adapter | P3B wording; `docs/providers/carolina_cloud.md`; principle 8 |
| 8 | Strangler inside DataWeaver | P5 compat adapters; migration → LATER |
| 9 | IDs generated once | Principle 5; `ids.py`; P5 caller‑supplied PKs; P6 `PUT /{id}`; P7 parity |
| 10 | `trace_id` | Principle 6; P0 cross‑cutting field; `GET /traces/{trace_id}` |
| 11 | Provenance versions | `provenance.py` (agent/prompt template/code_commit/schema_version) |
| 12 | Dual‑use triage only | `DualUseTriageCheck` (P2) |
| 13 | ORCH‑001 | Section above; P0 skeleton → P9 green |
| 14 | Execution‑state invariants | `invariants.py` + table |
| 15 | Updated phase structure | Sequencing |
| 16 | MUST / SHOULD / LATER | Priorities |
| 17 | Schedule guidance | Sequencing → schedule guidance |
| 18 | Pre‑P0 checklist | All of the above are in this revision; nothing deferred |
