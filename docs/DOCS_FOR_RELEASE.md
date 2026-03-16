# Documentation for Release — Evaluation

This document evaluates the `docs/` folder and defines what is **needed for release** versus what can be archived or kept as reference. The goal is to keep clear documentation on **system architecture**, **agents**, **user scenarios**, and **tests**.

---

## What we keep for release

### 1. System architecture

| Keep | Path | Reason |
|------|------|--------|
| ✅ | **SYSTEM_OVERVIEW.md** | Single best high-level description: multi-agent flow, PRIMARY vs EXPERIMENTAL, concepts. |
| ✅ | **AWS_USAGE.md** | How AWS is utilized: hosting (ECS, S3, CloudFront, ALB, ECR), data (S3), compute (Local/EC2/EMR), sync vs async, IAM. |
| ✅ | **architecture/ARCHITECTURE_EXPLANATION.md** | Frontend/backend deployment (S3, ECS, CloudFront), local vs production. |
| ✅ | **architecture/BACKEND_DATAFLOW.md** | Request → intent → tool → infra → execution → result; useful for onboarding and debugging. |
| ✅ | **ORCHESTRATION_DUALITY.md** | Clarifies which orchestrator is production (`agent.py`) vs experimental (`orchestrator.py`). |
| ✅ | **architecture/TARGET_ARCHITECTURE_CHECKLIST.md** | Implementation checklist for broker, async routing, EMR. |

### 2. Agents

| Keep | Path | Reason |
|------|------|--------|
| ✅ | **agent.md** | BioAgent system prompt: principles, micro/macroflow, session awareness, output contract. |
| ✅ | **HELIX_VS_CLAWBIO.md** | Positioning and comparison; informs architecture and product decisions. |

### 3. User scenarios

| Keep | Path | Reason |
|------|------|--------|
| ✅ | **USER_SCENARIOS.md** | Canonical scenario set (A–N): personas, triggers, success criteria, roadmap. Drives QA and product. |

### 4. Tests

| Keep | Path | Reason |
|------|------|--------|
| ✅ | **TESTBED.md** | Maps every scenario to test IDs, layers, assertions; points to `tests/testbed/`. |

### 5. Getting started and operations

| Keep | Path | Reason |
|------|------|--------|
| ✅ | **getting-started/ENVIRONMENT_SETUP.md** | Env vars, API keys. |
| ✅ | **deployment/AWS_DEPLOYMENT_GUIDE.md** | Primary deployment guide. |
| ✅ | **development/DEVELOPMENT_GUIDE.md** | Dev setup and contributing. |
| ✅ | **troubleshooting/** | EMR, S3 troubleshooting. |

### 6. Index

| Keep | Path | Reason |
|------|------|--------|
| ✅ | **README.md** | Entry point; lists only the release set (architecture, agents, scenarios, tests, getting started, deployment, development). |

---

## Optional / reference (keep but not required for “release set”)

- **FASTQC_EXECUTION_FLOW.md** — Deep-dive for FastQC routing and execution.
- **guides/** — Feature guides (DNA vendor, EMR, phylogenetic tree, single-cell, etc.).
- **demos/** — Demo walkthroughs.
- **reference/** — EC2, API keys, microservices, UV, etc.
- **api/** — ALB/API URL setup.
- **SAFETY_POLICY.md**, **STATISTICAL_GUIDELINES.md**, **SESSION_MANAGEMENT.md** — Policy and behavior.
- **OUTPUT_SCHEMA.md**, **TRACE_FORMAT_GUIDE.md** — For integrators and debugging.

---

## Archive or don’t promote (one-off / historical)

These have been moved to `docs/archive/` (subdirs: root, analysis, deployment, development, fixes, refactor). They are kept for history only.

**Agent / refactor notes:**  
AGENT_MD_*.md, AGENT_NAMING_MIGRATION.md, AGENT_VERSION_CLEANUP_SUMMARY.md, P1_ORCHESTRATION_CLARITY.md.

**Fix and investigation summaries:**  
CODE_DEBUGGING_COMPLETE.md, CODE_REVIEW_*.md, DATA_FLOW_FIX.md, DEBUGGING_SUCCESS.md, EMR_*_FIX*.md, EMR_*_ANALYSIS.md, FRONTEND_*_VERIFICATION.md, FRONTEND_RENDERING_ISSUES.md, FASTQC_ROUTING_FIX.md, HANDOFF_POLICY_IMPLEMENTATION.md, IN_MEMORY_JOB_TRACKING_ISSUE.md, LOCAL_FASTQC_IMPLEMENTATION.md, P0_FIXES_SUMMARY.md, PROVENANCE_HASHING_FIX.md, SANDBOX_*.md, TIMEOUT_FIX_SUMMARY.md, and similar one-off fix/analysis docs.

**Refactor phase docs:**  
refactor/PHASE_0_RECONNAISSANCE.md through PHASE_5_SUMMARY.md, refactor/FINAL_DELIVERABLES.md.

**Deployment iteration notes:**  
deployment/COPILOT_*.md, CDK_*.md, DEPLOYMENT_ALTERNATIVES*.md, DEPLOYMENT_STATUS_REPORT.md, STACK_ROLLBACK_RECOVERY.md, etc. (keep only AWS_DEPLOYMENT_GUIDE.md as the main deployment doc).

**Analysis and investigations:**  
analysis/*.md (EMR, job validation, PySpark, performance, etc.), CODEBASE_ORGANIZATION_ANALYSIS.md, COMPLETE_SUMMARY.md, COMPLETE_TRACE_CATALOG.md, DEMO_EVAL_SYSTEM_COMPLETE.md, END_TO_END_VERIFICATION.md, INTEGRATION_SUCCESS.md, JOB_STATUS_TRACKING.md, SEMANTIC_INVARIANTS.md, TRACE_*.md, TASK_ARCHITECTURE.md.

**Strategy / GTM / pitch (optional):**  
HELIX_AI_GTM_STRATEGY_INPUT.md, PITCH_DECK_OUTLINE.md, TRL_ASSESSMENT.md, TDS_HELIX_AI_DRAFT.md — keep only if you want them in the repo for release.

**Duplicate or moved content:**  
Various fixes in `docs/fixes/` that duplicate or supersede root-level fix docs; `docs/agent.md` vs `docs/reference/agent.md` (keep one; README points to `docs/agent.md`).

---

## Suggested layout for release (cleanup applied)

One-off and historical docs have been moved to `docs/archive/` (root, analysis, deployment, development, fixes, refactor). The layout below reflects the current release set.

```
docs/
├── README.md                    ← Entry point (release set only)
├── DOCS_FOR_RELEASE.md          ← This evaluation
├── SYSTEM_OVERVIEW.md          ← Architecture (primary)
├── AWS_USAGE.md                ← How AWS is used (hosting, data, compute)
├── ORCHESTRATION_DUALITY.md
├── agent.md                    ← Agents (BioAgent prompt)
├── HELIX_VS_CLAWBIO.md
├── USER_SCENARIOS.md           ← User scenarios
├── TESTBED.md                  ← Tests
├── FASTQC_EXECUTION_FLOW.md    ← Optional
├── architecture/
│   ├── ARCHITECTURE_EXPLANATION.md
│   ├── BACKEND_DATAFLOW.md
│   └── TARGET_ARCHITECTURE_CHECKLIST.md
├── getting-started/
│   └── ENVIRONMENT_SETUP.md
├── deployment/
│   └── AWS_DEPLOYMENT_GUIDE.md
├── development/
│   └── DEVELOPMENT_GUIDE.md
├── troubleshooting/
│   └── (existing)
├── guides/                      ← Optional
├── demos/                       ← Optional
├── reference/                   ← Optional
└── archive/                     ← Optional: move one-off/historical docs here
```

Historical docs are in **docs/archive/**; the **release set** is what **README.md** links to above.

---

## Summary

- **Architecture:** SYSTEM_OVERVIEW.md + architecture/* + ORCHESTRATION_DUALITY.md.
- **Agents:** agent.md + HELIX_VS_CLAWBIO.md.
- **User scenarios:** USER_SCENARIOS.md.
- **Tests:** TESTBED.md (+ `tests/testbed/` in code).
- **Operations:** getting-started, deployment, development, troubleshooting.

**README.md** is updated to list only these; **DOCS_FOR_RELEASE.md** records the full evaluation and what to archive if you want a leaner tree.
