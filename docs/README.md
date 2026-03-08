# Helix.AI Documentation

This directory contains the documentation needed for **release** and ongoing development. Start here for architecture, agents, user scenarios, and tests.

---

## Release documentation (core)

### 1. System architecture

| Document | Description |
|----------|-------------|
| [**SYSTEM_OVERVIEW.md**](SYSTEM_OVERVIEW.md) | High-level system overview: what Helix.AI is, multi-agent flow (Intent → Planner → Infrastructure → Implementation → ExecutionBroker → Visualizer), separation of concerns, foundational concepts. **Primary architecture reference.** |
| [**AWS_USAGE.md**](AWS_USAGE.md) | **How the cloud (AWS) is used**: hosting (ECS, S3, CloudFront, ALB, ECR), data (S3 inputs/outputs, size estimation), compute (Local vs EC2 vs EMR, 3-factor model, sync vs async). |
| [**architecture/ARCHITECTURE_EXPLANATION.md**](architecture/ARCHITECTURE_EXPLANATION.md) | Frontend vs backend deployment (local vs production, S3/CloudFront, ECS), how frontend connects to backend. |
| [**architecture/BACKEND_DATAFLOW.md**](architecture/BACKEND_DATAFLOW.md) | End-to-end backend flow: request reception, preflight, intent, tool selection, infrastructure decision, execution routing, result return. |
| [**ORCHESTRATION_DUALITY.md**](ORCHESTRATION_DUALITY.md) | Explains the two orchestration systems (primary in `agent.py` vs experimental in `orchestrator.py`) and which is used in production. |
| [**architecture/TARGET_ARCHITECTURE_CHECKLIST.md**](architecture/TARGET_ARCHITECTURE_CHECKLIST.md) | Rollout checklist for execution broker, sync/async routing, job management, EMR. |

### 2. Agents

| Document | Description |
|----------|-------------|
| [**agent.md**](agent.md) | BioAgent system prompt: role, principles, task types (micro/macroflow), session awareness, safety, structured output. Defines how the primary agent behaves. |
| [**HELIX_VS_CLAWBIO.md**](HELIX_VS_CLAWBIO.md) | Comparison with ClawBio: philosophy, architecture, feature matrix, strengths/weaknesses. Useful for positioning and design decisions. |

### 3. User scenarios

| Document | Description |
|----------|-------------|
| [**USER_SCENARIOS.md**](USER_SCENARIOS.md) | Scenario definitions for scientists and bioinformaticians: domain context, personas, exploration (A), first run (B), iteration (C), validation (D), reproducibility (E), sequence/plasmid (F–I), roadmap (J–N). Drives product and QA. |

### 4. Tests

| Document | Description |
|----------|-------------|
| [**TESTBED.md**](TESTBED.md) | Testbed definition: every scenario mapped to test IDs, layers (backend unit/integration, frontend, E2E), assertions, and implementation layout. Backend tests live in `tests/testbed/`. |

---

## Getting started and operations

| Document | Description |
|----------|-------------|
| [**getting-started/ENVIRONMENT_SETUP.md**](getting-started/ENVIRONMENT_SETUP.md) | Environment variables, API keys, `.env` setup. |
| [**deployment/AWS_DEPLOYMENT_GUIDE.md**](deployment/AWS_DEPLOYMENT_GUIDE.md) | AWS deployment (CDK, ECS, etc.). |
| [**development/DEVELOPMENT_GUIDE.md**](development/DEVELOPMENT_GUIDE.md) | Development setup, running locally, contributing. |
| [**troubleshooting/**](troubleshooting/) | EMR, S3, and other operational troubleshooting. |

---

## Optional / reference

- **FASTQC_EXECUTION_FLOW.md** — FastQC-specific execution and routing flow.
- **guides/** — Feature guides (DNA vendor, EMR, phylogenetic tree, single-cell, etc.).
- **demos/** — Demo walkthroughs and quick commands.
- **reference/** — Setup and reference notes (EC2, API keys, microservices, etc.).
- **api/** — ALB and API URL setup.
- **SAFETY_POLICY.md**, **STATISTICAL_GUIDELINES.md**, **SESSION_MANAGEMENT.md** — Policy and behavior.
- **OUTPUT_SCHEMA.md**, **TRACE_FORMAT_GUIDE.md** — For integrators and debugging.

---

## Archive

One-off and historical docs (fix summaries, refactor phases, deployment iterations, analysis) have been moved to [**archive/**](archive/). They are kept for history and are not part of the release set. See [**DOCS_FOR_RELEASE.md**](DOCS_FOR_RELEASE.md) for the evaluation and layout.
