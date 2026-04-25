# Changelog

All notable changes to Helix.AI are documented here.

---

## Unreleased (feature/p0-benchmark-gates)

### No-fallback LLM routing — strict classification throughout
- **Removed all keyword-based matching** from intent, routing, approval, and staging
  classification. Every decision is now LLM-driven with no silent fallbacks.
- `approval_classifier.py` — removed keyword fallback in `classify_approval`;
  LLM failure now raises `ApprovalClassificationError` (no silent `False`).
  Added interrogative-prefix early-rejection (`_INTERROGATIVE_PREFIXES`) and
  expanded `_ANALYTICAL_PATTERNS` with modification/iteration verbs so common
  non-approval commands never reach the LLM.
- `intent_classifier.py` — removed `_classify_intent_with_heuristics` fallback;
  LLM failure raises `IntentClassificationError`.
- `approval_policy.py` — replaced `APPROVAL_COMMANDS` keyword set and all keyword-based
  helpers with an LLM staging classifier (`_classify_staging_intent`).
- `command_router.py` — `HELIX_LLM_ROUTER_FIRST` defaults to `"1"`; raises
  `RoutingError` when LLM returns `None`.

### Test infrastructure
- Added `tests/unit/backend/conftest.py` autouse fixtures (`_mock_staging_classifier`,
  `_mock_intent_classifier`) that patch LLM-dependent classifiers in `HELIX_MOCK_MODE=1`
  so all ~700 unit tests run without a live API key.
- Updated `test_approval_classifier.py`, `test_intent_classifier.py`,
  `test_orchestration_modules.py`, and all routing/dispatch tests to assert that LLM
  failure raises rather than silently falling back.

---

## 3.0.0 — Tabular MVP + Plan → Approve → Execute loop

### Plan → Approve → Execute loop
- Every analysis request now produces a structured plan that must be approved by the
  user before execution. Approval is classified by an LLM (`ApprovalClassifier`) —
  no keyword matching.
- `session_param_extractor.py` — LLM-based parameter and intent extraction.
- `history_manager.py` — plan state tracked per session with approval gate.

### Tabular data MVP (CSV / Excel / TSV)
- Upload-time file profiling for tabular formats.
- `tabular_qa/` — Code Interpreter: LLM generates a sandboxed Python analysis
  script, executes it, and returns structured results.
- Schema preview UI: column types, row counts, and sample values shown in the
  frontend after upload.

### P0 benchmark gates
- `benchmarks/release_thresholds.yaml` — machine-readable pass/fail thresholds for
  unit, integration, workflow, frontend regression, and eval task rates.
- `tests/unit/benchmarks/` — scorer and gate-case fixtures.
- `tests/integration/test_deployed_demo_scenarios.py` — integration gate against
  live backend.

### Frontend
- File upload panel with drag-and-drop; schema metadata displayed after upload.
- Plan display component: structured plan steps shown before execution.
- Approval chip shortcuts ("Proceed", "Revise plan", "Reject") in the chat input.

---

## 2.0.0 — Agent orchestration + AWS deployment

- Multi-agent orchestration: `WorkflowPlannerAgent`, `ImplementationAgent`,
  `InfrastructureDecisionAgent`, `ToolGeneratorAgent`, `ExecutionBroker`.
- Sync/async execution routing; job management for long-running analyses.
- AWS deployment tooling: CDK stack (EC2, EMR, ECS, S3, CloudFront, ALB, ECR),
  GitHub Actions pipeline.
- Frontend UX improvements: demo scenarios, long-running job progress display.

---

## 1.0.0 — Initial release

- FastAPI backend with LLM-driven bioinformatics tool routing.
- 16+ built-in tools: scRNA-seq, bulk RNA-seq, FastQC, alignment, phylogenetics,
  NCBI/UniProt, GO enrichment, plasmid visualization, variant selection, read trimming.
- React + Vite frontend.
- Session management with `parent_run_id` iteration chains.
- Reproducibility bundles: `analysis.py` + `bundle.zip`.
