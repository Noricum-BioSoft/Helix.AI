# Noricum Platform Gap Analysis — Helix.AI + DataWeaver.AI

**Date:** 2026-10-07
**Status:** Analysis complete; no code changed. Implementation plan lives in
`.cursor/plans/noricum-scientific-orchestration-platform.md`.

**Target:** a secure, vendor‑neutral scientific application and orchestration layer
implementing the cycle

> Scientific Objective → Plan → Security Gate → Infrastructure Selection → Human Approval → Execution → Learning → Decision Making

with **Helix.AI** owning reasoning/orchestration and **DataWeaver.AI** owning the
persistent scientific context (what was known, hypothesized, done, measured,
concluded, and why the next decision was made).

---

## 1. Scope and method

| Repo | Inspected | Size | Baseline |
|---|---|---|---|
| Helix.AI (`/Users/eoberortner/git/Helix.AI`, commit `2a062a3`) | `backend/`, `shared/`, `agents/`, `tools/`, `tests/`, `benchmarks/`, `docs/architecture`, `.cursor/plans` | backend ≈34k LOC Python; `main.py` 7.5k | `HELIX_MOCK_MODE=1 pytest tests/unit`: **876 passed / 5 failed / 1 skipped** (pre‑existing; see `artifacts/test_results/unit_baseline_2026-10-07_platform_gap_analysis.md`) |
| DataWeaver.AI (`github.com/Noricum-BioSoft/DataWeaver.AI`, shallow clone to `/tmp`) | whole repo | backend ≈11k LOC Python; ≈72 tests | not executed (SQLite‑backed; runs without LLM key per `backend/tests/conftest.py`) |

Frontends were inventoried only to the extent needed to assess API contracts.

---

## 2. Current state — Helix.AI

### 2.1 Request lifecycle (production path)

```
POST /execute (backend/main.py ~2924)
  → session + WorkflowCheckpoint load
  → approval/staging gate   (orchestration/approval_classifier.py, approval_policy.py)
  → agent.handle_command   (backend/agent.py; HandoffPolicy enforces agent order)
     or CommandRouter.route_command (LLM router)
  → Plan IR staged in checkpoint (WAITING_FOR_INPUTS / WAITING_FOR_APPROVAL)
  → ExecutionBroker.execute_tool (backend/execution_broker.py)
  → sync tool executor | JobManager local/EMR job | Nextflow
  → history_manager.add_run / register_artifact → response envelope
```

A second, cleaner pipeline exists in `backend/orchestrator.py` (Planner → Infra →
Implementation, with contract hashing) but is not on the HTTP path
(`docs/ORCHESTRATION_DUALITY.md`). This duality is a known P1 risk and matters here
because the target chain must be enforced in **one** place.

### 2.2 What exists that the target needs (assets to preserve)

| Asset | Location | Why it matters |
|---|---|---|
| Single‑owner agent roles + enforced handoff chain `IntentDetector → Planner → Infra → [CodeGen] → Broker → Visualizer` | `agents/agent-responsibilities.md`, `agents/handoff-policy.md`, `backend/agent.py::HandoffPolicy`, `backend/policy_checks.py` | The target chain is an *extension* of this: insert `SecurityGate` and `HumanApproval` roles. Policy‑violation tests already exist (`tests/unit/backend/test_handoff_policy.py`). |
| Typed contracts | `backend/contracts/{workflow_plan,infra_decision,execution_spec,dataset_spec}.py`, `backend/plan_ir.py`, `shared/contracts.py` | Pydantic v2; `InfraDecision` already carries `decision_summary`, `reasoning`, `cost_analysis`, `alternatives`, `warnings`, `confidence_score`. |
| Persisted workflow state machine | `backend/workflow_checkpoint.py` (`IDLE … WAITING_FOR_APPROVAL … EXECUTING … COMPLETED`) saved under `session["__checkpoint__"]` | Natural home for `WAITING_FOR_SECURITY_REVIEW`, `DENIED`, `WAITING_FOR_PROVIDER` states. |
| Approval staging | `orchestration/approval_policy.py::should_stage_for_approval`; `READ_ONLY_ROUTER_TOOLS`; `multi_step_workflow` always staged | Plan → Approve → Execute is already the product contract. |
| Run ledger + artifact registry + lineage | `history_manager.add_run`, `register_artifact`, `get_lineage_edges`; `bundle_generator.py` → `run_manifest.json` | Session‑scoped provenance with `parent_artifact_ids`, `derived_from`, `params`. |
| Infrastructure Expert | `infrastructure_decision_agent.py` (LLM + heuristic fallback), `config/environment_capabilities.yaml`, `config/cost_heuristics.yaml`, `tools/env_catalog.py` | YAML env catalog is a proto‑capability registry (cost class, startup, reproducibility, `best_for/avoid_for`). |
| Executors | `execution_broker.py`, `sandbox_executor.py`, `script_executor.py`, `ec2_executor.py`, `nextflow_executor.py`, `job_manager.py` (status / results / cancel / retry / logs) | All the *behaviour* an Execution Fabric needs exists; it lacks a common interface. |
| Ingress policy gate | `orchestration/upload_intake_policy.py` → `allow | allow_with_approval | block`, `sensitivity_class`, `scan_flags`, `approval_required`; `POST /session/{id}/uploads/approve`; `execute_plan` 409 on pending approvals | The only *executable* policy gate today; its outcome shape is the seed of `SecurityAssessment`. |
| Tool inventory | `tool_schemas.py`, `tool_inventory.py` (AST discovery + `which`), `unsupported_tools.py` | Tool‑level capability descriptions (inputs/outputs/tags). |
| Post‑run interpretation | `tabular_qa/analysis_executor.py` (interpretation step), `ds_pipeline/reviewer.py` (recommended next experiments), `advisory_normalizer.py::HelixAdvisory` | Seeds for Observation / Interpretation / Decision. |
| Release gate machinery | `benchmarks/release_thresholds.yaml`, `benchmarks/scoring/*`, `artifacts/release_readiness.json`, `.github/workflows/benchmark-gate.yml` | New capabilities can be gated the same way. |

### 2.3 Helix gaps against required capabilities

| Cap. | Required | Exists today | Gap | Extension point |
|---|---|---|---|---|
| **A** Scientific Objective | Structured objective (biological system, question, constraints, desired evidence, known inputs, expected outputs, success criteria) | Implicit in the command string; `IntentResult{intent, confidence}`; playbook `matches` in `workflow_planner_agent.py` | **No `ScientificObjective` type.** Planner, gate, and recommender all re‑derive intent from prose. | New contract in `backend/contracts/scientific_objective.py`; produced by Intent Detector/Planner; stored on checkpoint and in ledger. |
| **B** Persisted scientific plan | First‑class, versioned plan artifact distinguishing reasoning / computational steps / experimental steps / data / capabilities / assumptions / security requirements | `plan_ir.Plan{steps[tool_name, arguments]}` staged as a **dict** in `WorkflowCheckpoint.pending_plan`; cleared on approval. Two unrelated `WorkflowPlan` classes (`backend/contracts` vs `shared/contracts`). | Plan has no `plan_id`, version, hash, status history, step kind (computational vs experimental), security requirements, or ledger record. Revisions overwrite. | Wrap `plan_ir.Plan` in `ScientificPlan` (new contract); register as artifact `type="scientific_plan"` via `register_artifact`; keep `pending_plan` as a pointer. |
| **C** Secure Science Gate | Pluggable chain → `ALLOW / ALLOW_WITH_APPROVAL / REQUIRE_REVIEW / DENY` + reasons, policies, external screening, required approvals, audit | `upload_intake_policy` (uploads only, single function, profile forced to `p0`); `policy_checks.py` enforces *agent boundaries*, not science safety; `docs/SAFETY_POLICY.md` is prose; `synthesis_submission.py` validates charset/length only | **No gate on scientific actions**, no plugin interface, no sequence‑screening integration point, audit is log‑only (`_emit_policy_audit_event`). Not represented in the handoff chain. | New `backend/security/` package: `SecurityAssessment` contract, `SecurityCheck` protocol, `SecureScienceGate` chain; wrap upload policy as first check; new `SECURITY_GATE` role in `HandoffPolicy` between Planner and Infra. |
| **D** Infrastructure / execution recommender | Rank arbitrary execution targets (local, HPC, cloud, Nextflow/Seqera, APIs, labs, CROs) on suitability, cost, turnaround, throughput, security, privacy, availability, prior performance, provenance; expose criteria | `InfraDecision.infrastructure: Literal["Local","EC2","EMR","Batch","Lambda"]` (duplicated in `execution_spec.py`, `shared/contracts.py`, `environment_capabilities.yaml`); criteria = size/locality/cost/compute; explanation fields exist | Target set is **closed and compute‑only**; no security/privacy/turnaround/prior‑performance criteria; criteria weights opaque (LLM prose + heuristics); Nextflow executes but isn't a decision literal. | Generalise `infrastructure` → `provider_id` validated against the Capability Registry; add `criteria_scores` vector to `InfraDecision`; keep legacy literal as alias for back‑compat. |
| **E** Human approval | Persisted, auditable (who/when/what/hash), action‑specific policies | `WAITING_FOR_APPROVAL` state; approval recognised by **LLM classifier only** (3 unit tests fail under mock mode because of this); upload approvals store note+timestamp but **no principal**; plan approval leaves no record once `pending_plan` is cleared | No `HumanApproval` record, no identity, no plan hash binding, no explicit API. | `HumanApproval` contract; `POST /session/{id}/plans/{plan_id}/approve|reject` endpoint; ledger entry; checkpoint keeps `approval_id`. Keep NL approval as a convenience that *calls* the same endpoint. |
| **F** Execution Fabric | Adapter contract (`describe_capabilities, validate_request, estimate, prepare_execution, execute, get_status, retrieve_results, cancel, provenance`); local + external API + mock wet‑lab providers | `ExecutionBroker` branches on `tool_name` and sync/async; `JobManager` has status/results/cancel; `nextflow_executor`, `ec2_executor`, `sandbox_executor` each have bespoke APIs; no wet‑lab or mock provider | **No common interface**; providers can't be added without editing the broker; no estimate/validate phase; provenance only via bundle. | `backend/execution/providers/` with `ExecutionProvider` Protocol; `LocalComputeProvider` (wraps tool executor + sandbox), `NextflowProvider` (wraps `nextflow_executor`), `EmrProvider` (wraps JobManager EMR), `MockExperimentalProvider`; broker becomes dispatcher behind flag. |
| **G** Learning & decision | Ingest evidence via DataWeaver; assess objective; supported/unsupported hypotheses; uncertainties; next decision; Observation ≠ Interpretation ≠ Decision | `analysis_executor` interpretation step; `ds_pipeline/reviewer.py` next‑experiment text; `HelixAdvisory{summary,next_steps}`; "Decision Packet" is *planned* in `.cursor/plans/system-structural-fixes.md` | No structured Observation/Interpretation/Decision model; no hypothesis linkage; nothing persisted beyond prose; next step never becomes a new objective. | `backend/learning/` with `EvidenceAssessment` contract and `NextDecision` enum; reuse reviewer/interpretation prompts; writes to DataWeaver. |
| **Registry** Scientific Capability Registry | Shared descriptors: provider, capability_id, category, inputs/outputs, systems, constraints, cost, turnaround, integration, security | `environment_capabilities.yaml` (compute envs), `tool_schemas.py` (local tools), `EC2_EXPECTED_CLI_TOOLS` | Two partial registries, neither covers experimental providers, cost models, turnaround, or security flags; not consumed by gate or planner. | `shared/capability_registry.py` (Pydantic schema) + `backend/config/capabilities/*.yaml`; env catalog becomes one provider family. |

### 2.4 Cross‑cutting Helix issues that the plan must respect

1. **Orchestration duality** — new roles must be enforced in the production `agent.py`/`main.py` path, not only in `orchestrator.py`.
2. **Optional persistence** (Slice A in `system-structural-fixes.md`) — new plan/approval/assessment records must go through a single `record_turn`‑style write path or they will exhibit the same "empty bundle" bugs.
3. **Multiple intent layers** — do not add a fourth classifier; the objective should be produced *once* and passed down.
4. **`main.py` size (7.5k LOC)** — new endpoints should live in routers (`backend/api/*.py`) and be mounted, not appended.
5. **Session JSON as the only store** — adequate for MVP ledger; cross‑session scientific memory must live in DataWeaver, not in `sessions/*.json`.

---

## 3. Current state — DataWeaver.AI

### 3.1 Stack

FastAPI + SQLAlchemy (declarative) + Alembic; SQLite default, PostgreSQL via
`DATABASE_URL` (UUID columns already use `sqlalchemy.dialects.postgresql.UUID`);
React 18/CRA + Tailwind + Plotly; OpenAI SDK direct (`gpt-3.5-turbo` hard‑coded in
`data_qa_service.py`, `general_chat.py`); `celery`/`redis`/`python-jose`/`passlib`
declared but unused; `plotly`/`openpyxl` used but undeclared. No auth, no tenancy,
no CI, no Dockerfile in tree.

### 3.2 Data model today

Three partially overlapping model sets:

| Set | Tables | Notes |
|---|---|---|
| `app/models/{workflow,file,dataset}.py` | `workflows`, `workflow_steps` (`external_provider`, `external_config`), `files` (`parent_file_id`), `file_metadata`, `file_relationships` (`relationship_type`, `confidence_score`), `datasets`, `dataset_matches` (`match_type`, `confidence_score`, `is_confirmed`, `confirmed_by/at`) | Integer PKs. **`FileRelationship` is a nascent edge table**; `DatasetMatch` has human‑confirmation fields. |
| `app/models/workflow_phases.py` | `workflow_projects`, `design_phases`, `build_phases`, `test_phases`, `workflow_files`, `workflow_correlations` | UUID PKs; DBTL‑specific columns (vendor, Km/Vmax…); `DesignPhase.project_id` is a bare string, not an FK. |
| `models/bio_entities.py` | `designs`, `builds`, `tests` | UUID PKs; `lineage_hash`, `parent_*_id`, `generation`; `Test.match_confidence/match_method/match_score`. `GET /api/bio/lineage/{design_id}`. |

Ephemeral: `services/workflow_state.py::WorkflowStateManager` and
`app/services/data_context.py::DataContextManager` (`DataContext{type, parent_ids}`)
hold uploaded/merged/visualisation/analysis state **in process memory**, 24 h TTL.

Alembic chain: `002` (bio entities, `down_revision=None`) → `f7c123bd3c74`
(workflows/files/datasets/phases). `database.create_tables()` also calls
`create_all`, and imports from `app.models` only, so the two model trees are not
registered on one metadata consistently.

### 3.3 API today (prefix `/api`)

`/files/*` (upload, metadata, relationships), `/workflows/*` (CRUD + steps),
`/datasets/*` (CRUD, process, match/auto‑match, confirm/reject),
`/bio/*` from **two routers with colliding names** (`app/api/bio_matcher.py`
1,723 LOC session/merge/viz/QA; `api/bio_entities.py` Design/Build/Test CRUD,
`upload-test-results`, `match-preview`, `lineage`), `/intelligent-merge/*`,
`/data-qa/*`, `/general-chat/chat`. Frontend `api.ts` calls `/workflows/{id}/execute|lineage|merge|history`, which are **not implemented**.

### 3.4 Services

`IntelligentMerger` (join/concat strategy suggestion + execution),
`MatchingService` (exact + fuzzy identifier matching; ML TODO),
`BioEntityMatcher` (sequence / mutation / alias → Design/Build with score),
`DataAnalyzer` (quality, stats, correlations), `simple_visualizer` + Plotly JSON,
`DataQAService` (dataframe summary → LLM; rule‑based fallback). No code
generation, no sandbox, no tool calling.

### 3.5 DataWeaver gaps against required capabilities

| Cap. | Required | Exists | Gap | Build on |
|---|---|---|---|---|
| **A** Heterogeneous integration | CSV/Excel/PDF/reports/images/DB/API/analysis outputs; link across modalities, experiments, entities, samples, providers; partial data is normal | CSV (+Excel/JSON/Parquet in merger only); identifier/fuzzy/sequence matching; column‑level join suggestions | No PDF/report/image ingestion; no DB/API connectors (UI connectors are hard‑coded demo); linking is file/column‑level, not entity‑level | `MatchingService`, `BioEntityMatcher`, `IntelligentMerger`, `DatasetMatch` |
| **B** Canonical domain model | Objective, Hypothesis, Evidence, Dataset, BiologicalEntity, Sample, Experiment, ExperimentPlan, Condition, Capability, ExecutionProvider, ExecutionRequest, SecurityAssessment, HumanApproval, ExecutionRun, Observation, Measurement, AnalysisResult, Artifact, Conclusion, Decision, ProvenanceEvent | Dataset, File(≈Artifact), Workflow/Step, Design/Build/Test (≈BiologicalEntity/Experiment/Measurement, narrow), WorkflowProject | **~14 of 21 entities absent**; two DBTL schemas duplicate each other; mixed Integer/UUID keys; no `schema_version` | Promote `Design/Build/Test` → `BiologicalEntity`/`Experiment`/`Measurement`; `Dataset`/`File` → `Dataset`/`Artifact` |
| **C** First‑class relationships | Typed, queryable edges (motivated, produced, analyzed_by, supports, contradicts, based_on, generated) across all entity types | `FileRelationship` (file↔file), `DatasetMatch` (dataset↔file), FK lineage on Design/Build/Test, in‑memory `parent_ids` | Edges exist only per table pair; no polymorphic edge table; no traversal API beyond `/bio/lineage/{design_id}` | `FileRelationship` shape (type + confidence) generalised to a polymorphic `relationships` table |
| **D** Provenance & versioning | Every derived object: source, timestamp, originating experiment, tool, versions, parameters, actor, backend, parents, transformations; traceable to original evidence | `lineage_hash`, `parent_file_id`, `created_at`, `DataContext.parent_ids` (RAM) | No `ProvenanceEvent`, no actor/tool/version/parameters capture, no immutability/supersedes, session trail lost on restart | `lineage_hash` idea → content hashes on Artifact; `DataContext` → persisted events |
| **E** Scientific Context API | `get_objective_context`, `get_related_hypotheses`, `get_available_evidence`, `get_prior_experiments`, `get_measurements`, `get_provenance`, `get_prior_decisions`, `store_plan/execution/observation/decision`; stable service contract, no raw DB access | Session‑scoped `/bio/data-context/{session_id}`, `/bio/lineage/{design_id}`; OpenAPI spec in `docs/openapi.yaml` | None of the named operations; context is per in‑memory session; no auth/API key for a service caller | OpenAPI tooling (`scripts/validate_openapi.py`), existing Pydantic schemas |
| **F** Analysis → structured results | Keep dataframe ops/stats/joins/QC/Plotly; emit `AnalysisResult` + `Artifact` with provenance | All analysis types exist; results returned as JSON and optionally cached in `DataContextManager` | Results are not persisted as entities nor linked to input datasets/experiments; no parameters/version capture | `DataAnalyzer`, `simple_visualizer`, merge endpoints |

### 3.6 Cross‑cutting DataWeaver issues

1. `bio_matcher.py` god‑router (1.7k LOC) mixes session state, merge, viz, QA; should not be extended further — new knowledge endpoints go in a new router package.
2. Two `/api/bio` routers with overlapping paths (`/designs`, `/upload-test-results`) — order‑dependent behaviour.
3. In‑memory session state is incompatible with being Helix's system of record; a restart erases context.
4. No authentication — a service contract consumed by Helix needs at least API‑key/service identity and per‑project scoping.
5. Dependency hygiene (undeclared `plotly`/`openpyxl`, unused Celery/Redis) and missing CI/Docker.

---

## 4. Mapping the target workflow onto existing components

| Stage | Question | Helix owner today | DataWeaver record today | Target owner / record |
|---|---|---|---|---|
| Scientific Objective | What do we need to learn? | prose command; `IntentResult` | — | Helix `ScientificObjective` contract → DW `scientific_objectives` |
| Plan | What should we do? | `plan_ir.Plan` in checkpoint; `WorkflowPlan` | `Workflow/WorkflowStep` (unused by Helix) | Helix `ScientificPlan` artifact → DW `experiment_plans` |
| Security Gate | Should we execute it? | upload intake only | — | Helix `SecureScienceGate` → `SecurityAssessment` → DW `security_assessments` |
| Infrastructure Selection | Where should we execute it? | `InfraDecision` (compute only) | `WorkflowStep.external_provider` (string) | Helix recommender over Capability Registry → DW `execution_providers`, `capabilities` |
| Human Approval | (gate) | checkpoint state + LLM intent | `DatasetMatch.confirmed_by/at` (different purpose) | Helix `HumanApproval` ledger entry → DW `human_approvals` |
| Execution | (do it) | `ExecutionBroker` + executors + `JobManager` | `Workflow.status`, `File` | Helix Execution Fabric providers → DW `execution_requests`, `execution_runs`, `artifacts` |
| Learning | What did we learn? | interpretation prompts, reviewer | `Test.result_value`, `DataAnalyzer` | Helix `EvidenceAssessment` → DW `observations`, `measurements`, `analysis_results`, `conclusions` |
| Decision Making | What should we do next? | `HelixAdvisory.next_steps` (prose) | — | Helix `NextDecision` → DW `decisions` → new `scientific_objectives` (never auto‑executed) |

---

## 5. Architectural decisions (recommended; revisit after Phase 1)

| # | Decision | Rationale |
|---|---|---|
| D1 | **Persistence: PostgreSQL + relational tables + JSONB attrs + one polymorphic `relationships` edge table.** No graph database for MVP. | DataWeaver already runs SQLAlchemy/Alembic/Postgres; edge counts for MVP are small; recursive CTEs cover `get_provenance` depth queries; keeps SQLite for tests. Revisit if traversal queries dominate. |
| D2 | **DataWeaver owns the canonical domain model; Helix consumes it through a typed HTTP client behind a `KnowledgeStore` Protocol with a local‑ledger fallback.** | Avoids a third shared package for MVP; Helix unit tests stay offline; DataWeaver's OpenAPI remains the contract. A shared `noricum-schemas` package is a later extraction once schemas stabilise. |
| D3 | **Capability Registry schema lives in Helix `shared/capability_registry.py` with YAML descriptors; DataWeaver stores `Capability`/`ExecutionProvider` rows mirroring the descriptors (`spec` JSONB).** | Helix's recommender and gate are the primary consumers; DW needs the rows only for provenance of `ExecutionRequest`. One sync script keeps them aligned. |
| D4 | **Security gate = ordered chain of `SecurityCheck` plugins; most‑restrictive outcome wins; every outcome persisted.** External screening is an adapter interface with a mock + "not configured → REQUIRE_REVIEW" default. | Honours "integrate, don't recreate"; never bypasses provider controls; safe failure. |
| D5 | **Execution Fabric via strangler pattern:** `ExecutionBroker` keeps its API; a flag routes to the new provider dispatcher; existing paths wrapped as providers, not rewritten. | Preserves ~880 passing tests and production behaviour; parity measured by existing suites before cutover. |
| D6 | **Human approval becomes an explicit endpoint + ledger record; NL approval detection remains as a UX shortcut that invokes it.** | Fixes the mock‑mode‑fragile LLM‑only approval; gives who/when/what/hash audit. |
| D7 | **Decision never auto‑executes.** A `Decision` may create a new `ScientificObjective` in `proposed` state; the cycle restarts at Plan → Gate → Approval. | Explicit requirement; enforced in `HandoffPolicy` (Decision role is terminal). |
| D8 | **New Helix HTTP surface goes in `backend/api/*` routers**, mounted from `main.py`. | `main.py` is already 7.5k LOC. |

---

## 6. Risks and non‑goals

**Risks**
- Scope creep on the ontology: mitigated by a 21‑entity v1 with JSONB `attrs` and `schema_version`, not a full ontology.
- Contract drift between two repos: mitigated by DW OpenAPI snapshot checked into Helix tests (`tests/fixtures/dataweaver_openapi.json`) and a contract test.
- Behaviour regressions in `/execute`: all new behaviour flag‑gated (`HELIX_SCIENCE_GATE_V1`, `HELIX_EXECUTION_FABRIC_V1`, `HELIX_KNOWLEDGE_STORE=local|dataweaver`), default off until parity.
- Security theatre: a mock screening adapter must be labelled as such in every assessment (`screening_results[].provider="mock"`), and the default when no real adapter is configured is `REQUIRE_REVIEW`, not `ALLOW`.

**Non‑goals (this cycle)**
- Rewriting `agent.py`/`main.py` routing (covered by `.cursor/plans/system-structural-fixes.md` and `eliminate-keyword-routing.md`; this plan depends on Slice A persistence contract but does not re‑do it).
- Paid wet‑lab execution; any real vendor adapter.
- Replacing DataWeaver's frontend.
- A graph database.
