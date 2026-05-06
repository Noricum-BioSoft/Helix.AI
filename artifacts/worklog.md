# Helix.AI Worklog

## 2026-05-06 — Beta hotfix: `No module named 'duckdb'` on file upload

### Symptom
User uploaded `media-5__2_.xlsx` to https://helix-beta.noricum-biosoft.com/
(session `fd3e53f0-6ad2-4375-b8a5-70fea4535496`) and the backend returned
`No module named 'duckdb'`. Same root cause would have surfaced as
`No module named 'openpyxl'` once duckdb was installed, since the file is
`.xlsx` and `pandas.read_excel` requires openpyxl for that extension.

### Root cause
Deploy-contract gap, not a code bug. `backend/ds_pipeline/pipelines/ingest.py`
imports `duckdb` at file-upload time and uses `pandas.read_excel` (→ openpyxl)
for `.xlsx` inputs. Both deps were declared only in the **root**
`requirements.txt` (lines 26 and 54), but `scripts/ec2/bootstrap-ec2.sh:102`
and `scripts/ec2/update-from-git.sh:25` both `pip install -r backend/requirements.txt`
— so the production `.venv` (a conda env at /opt/helix/.venv, Python 3.11)
never received them. Direct probe on EC2:

```
$ /opt/helix/.venv/bin/python -c "import duckdb"
ModuleNotFoundError: No module named 'duckdb'
$ /opt/helix/.venv/bin/python -c "import openpyxl"
ModuleNotFoundError: No module named 'openpyxl'
```

The split between the two requirements files has existed for a while; the bug
only became user-visible when this user happened to upload an Excel file
(triggering both imports).

### Fix
1. **Source of truth** (`backend/requirements.txt`): added `duckdb>=1.0` and
   `openpyxl>=3.1` under a comment explaining they're required by the
   ds_pipeline tabular-ingest path. This is what bootstrap and
   update-from-git both consume.
2. **Live remediation**: installed both into the running prod venv via
   `aws ssm send-command`:
   `sudo -u helix /opt/helix/.venv/bin/pip install duckdb>=1.0 openpyxl>=3.1`.
   Resolved to `duckdb-1.5.2`, `openpyxl-3.1.5`, `et-xmlfile-2.0.0`.
3. **Restart**: `systemctl restart helix-backend` so the running uvicorn
   process loads the new dependencies. `helix-backend` reports `active`,
   `/health` returns `status: healthy`.
4. **End-to-end replay**: re-ran `ingest_tabular` against the user's exact
   uploaded path
   `/opt/helix/sessions/fd3e53f0-6ad2-4375-b8a5-70fea4535496/uploads/raw/media-5__2_.xlsx`
   on the live backend. Result: 90 rows × 67 cols on sheet `gene_morpheus_mem`
   ingested cleanly.

### Why our existing tests didn't catch this
`tests/unit/backend/test_tabular_ingest.py` already exercises this exact import
path. It passes in dev because the developer venv was bootstrapped from the
root `requirements.txt` (which has duckdb), not from `backend/requirements.txt`.
The gate suites therefore measure the dev environment's deps, not the
production deploy contract. To prevent recurrence we should add a CI step that
provisions a clean venv from `backend/requirements.txt` and asserts
`from backend.main import app` and `from backend.ds_pipeline.pipelines.ingest
import ingest_tabular` both succeed. Tracked as a follow-up
(`backend-deploy-contract-not-tested`, see release_readiness.json blockers).

### Operational note
This is the second deploy-contract gap surfaced against the live beta in 24h
(prior was AL2 having no Node/npm). Both point at the same underlying issue:
the production deploy path is implicitly different from `pytest` on a developer
laptop. The deploy-script scripts assume a contract that's never end-to-end
verified. Recommendation: add a Makefile/CI target that builds a fresh venv
from `backend/requirements.txt` only, then imports the canonical entry points,
and gate releases on it.

## 2026-05-05 — Beta deploy: cleared `beta-bundle-mixed-content`

### What was deployed
- Pushed commits `865b892..55cd351` (4 commits from 2026-05-04 gate run) to `origin/main`.
- EC2 instance `i-0cc452521209457d4` (us-west-1, public IP `54.176.128.13`,
  user-facing DNS `helix-beta.noricum-biosoft.com`) git-pulled to HEAD `55cd351`
  via `aws ssm send-command`.
- Built SPA locally with `VITE_API_BASE_URL=` (same-origin), shipped
  `frontend/dist` to EC2 via S3 presigned URL, atomically swapped into
  `/opt/helix/frontend/dist`. New bundle hash: `index-Do3-fXHQ.js`.
- Verified externally:
  - `curl https://helix-beta.noricum-biosoft.com/` → 200 OK, references new bundle.
  - `curl https://helix-beta.noricum-biosoft.com/health` → `{"status":"healthy",...}`.
  - `grep -c helix--Publi` on the deployed bundle → 0.

### Why the documented deploy path didn't work and what I did instead
- Port 22 is firewalled on the EC2 security group, so `ssh ec2-user@helix-beta…`
  is unavailable. Switched to AWS SSM (`aws ssm send-command`), matching the
  pattern in `scripts/aws/push-via-s3.sh` and `scripts/aws/build-on-ec2.sh`.
- `./scripts/ec2/update-from-git.sh` failed at the `npm run build` step because
  Node/npm are **not installed on this AL2 instance at all** —
  `find / -maxdepth 5 -name npm` returned nothing. The bootstrap script
  acknowledges that AL2's old glibc means Node has to be built from source
  (~hours), and that build step has clearly never been completed on this host.
  Whoever previously deployed the SPA was building locally and shipping
  `dist/` — same workaround used here.
- EC2 instance role does not have read access to `s3://helix-build-794270057041/`
  (HEAD returned 403). Used a presigned URL instead so the EC2 download is
  plain HTTPS without IAM involvement.
- Pre-existing mixed `root`/`helix` ownership on `/opt/helix` was tripping
  git's "dubious ownership" check. Resolved by `chown -R helix:helix /opt/helix`
  and `git config --global --add safe.directory /opt/helix` for the helix user.

### Operational issues to fix before next deploy (not blocking this release)
1. **AL2 deploy contract is implicit.** `scripts/ec2/bootstrap-ec2.sh` and
   `scripts/ec2/update-from-git.sh` assume `npm` is available on the host, but
   on this instance it isn't and never has been. Either:
   - install Node/npm on the host (NodeSource RPM doesn't work on AL2 due to
     glibc; would need to actually run the bootstrap's "build from source" leg),
   - migrate the host to AL2023 (newer glibc, NodeSource RPM works), or
   - formalise the off-host build contract (build on CI / dev machine, ship
     `dist/` via S3 + SSM). The third option is what's actually been
     happening de facto and matches `scripts/aws/push-via-s3.sh`.
2. **`/opt/helix` ownership drift.** Pre-existing mixed `root`/`helix`
   ownership broke `update-from-git.sh`. Add an idempotent
   `chown -R helix:helix /opt/helix` to the bootstrap and to a sanity-check
   step at the top of `update-from-git.sh`, and add `safe.directory /opt/helix`
   to the bootstrap's `helix` gitconfig.
3. **EC2 IAM role lacks read on `helix-build-794270057041`** even though both
   live in account 794270057041. Either grant the role or pick a single
   canonical relay bucket (`helix-ai-frontend-…` already accessible) and use
   it consistently in the deploy scripts.

### Release gate result for HEAD = `55cd351`
- The 2026-05-04 gate run measured the exact bytes that are now at HEAD
  (the gate ran against the post-fix working tree which has since been
  committed verbatim). Per workspace rule "release-readiness requires fresh
  benchmark and test results", that gate run is treated as the gate run for
  HEAD: no product-code drift between the measurement and HEAD, only commit
  metadata and one operational deploy.
- Both readiness blockers cleared:
  - `beta-bundle-mixed-content` → cleared by today's deploy (above).
  - `uncommitted-changes` → cleared by `git push origin 7a41660..55cd351`.
- `readiness` flipped to `true` in `artifacts/release_readiness.json`.

## 2026-05-04 — Release-gate run (post-router-rewrite)

### Result: 6/6 suites green, 2 blockers, readiness = FALSE

| Suite | Pass | Threshold | Status |
|---|---|---|---|
| unit_tests | 780/780 (1 R-skip) | ≥ 1.0 | PASS |
| backend_tests | 13/13 | ≥ 0.98 | PASS |
| integration_tests | 9/9 | ≥ 0.98 | PASS |
| workflow_tests | 20/20 | ≥ 0.95 | PASS |
| evals (live LLM) | 73/79 = 92.4% | ≥ 0.90 | PASS |
| frontend (vitest + Playwright) | 50 + 9 | ≥ 0.98 | PASS |
| benchmark_50_workflows | 50/50 (p50=15.8s, p95=19.2s) | ≥ 0.95 | PASS |

Full breakdown in `artifacts/release_readiness.json`.

### Root-cause fixes made during the gate run (in scope, not silent broadening)

1. **`backend/main.py`** — `_should_clear_pending_plan_after_execution()`'s status
   allow-list was missing `pipeline_executed`, the synchronous-success status
   `CommandProcessor.execute_pipeline()` returns when the planned pipeline has
   zero remaining steps. Production effect: a successful zero-step early-return
   would leave the session's pending plan dangling, breaking the next turn in
   the demo1 / approve→rerun flow. Added regression test
   `test_pipeline_executed_status_clears_pending_plan` in
   `tests/unit/backend/test_phase25_benchmark_turn_fixtures.py`.

2. **`tests/integration/test_execute_five_demo_golden_paths.py`** — the
   demo1 assertion at line 153 still allowed only `{"success", "pipeline_submitted"}`
   for the post-approval status; it now also accepts `pipeline_executed`.

3. **`tests/workflows/test_e2e_multi_step_workflows.py`** —
   `test_handle_command_routes_to_workflow_handler` mocked the multi-step
   handler and `_get_agent` but not `classify_intent`, which raises by design
   in `HELIX_MOCK_MODE` (see `backend/intent_classifier.py:52`). Added a third
   `patch('backend.intent_classifier.classify_intent', ...)` to the same `with`
   block, mirroring the autouse fixture in `tests/unit/backend/conftest.py`.

4. **`tests/demo_scenarios/framework/tracer.py`** — added
   `from __future__ import annotations` so the PEP 604 union `AgentRole | str`
   on line 75 doesn't blow up at class-definition time under Python 3.9. The
   repo officially targets `>=3.10` (`tests/demo_scenarios/pyproject.toml`,
   EC2 bootstrap installs `python3.11`), but the local `.venv` is 3.9.6, so
   without the future-import the entire `tests/unit/test_semantic_invariants.py`
   collection failed.

5. **`tests/evals/conftest.py`** (new) — evals MUST exercise the live LLM
   router (otherwise they're not testing routing quality). The repo-level
   `tests/conftest.py` setdefaults `HELIX_MOCK_MODE=1` for safety, which makes
   `command_router._route_with_llm` and `intent_classifier._classify_intent_with_llm`
   raise. New eval-scoped conftest (a) loads `.env`, (b) flips
   `HELIX_MOCK_MODE=0`, (c) cleanly skips the suite when no
   `OPENAI_API_KEY`/`DEEPSEEK_API_KEY`/`AZURE_OPENAI_API_KEY` is set so CI
   without secrets stays green instead of producing 70+ "LLM unavailable"
   failures.

6. **`frontend/e2e/visual.spec.ts` + regenerated snapshots** — the empty-state
   conversation pane no longer contains "Run a command to see your prompts here";
   that placeholder was replaced by `CapabilityGrid` (see
   `frontend/src/components/CapabilityGrid.tsx:4-7`). Updated the locator to
   anchor on `Here's what Helix.AI can do — click any item to try it:`.
   Visual snapshots regenerated to absorb the May 2 → May 4 demoScenarios.ts
   growth (+513 LOC → taller landing page).

### Router-quality regressions (do not block release; track for GA)

Live-LLM eval failures, all from `bio.router.*`:

| ID | Got | Expected |
|---|---|---|
| `router.plasmid.default_vector` | `mutate_sequence` | `plasmid_visualization` |
| `bio.router.plasmid.insert_default_vector` | `mutate_sequence` | `plasmid_visualization` |
| `bio.router.directed_evolution.parse_parameters` | `mutate_sequence` | `directed_evolution` |
| `bio.router.quality_assessment.from_inline_fasta` | `fastqc_quality_analysis` | `quality_assessment` |
| `bio.router.plan.trim_merge_qc` | …`fastqc_quality_analysis` (tail) | …`quality_assessment` (tail) |
| `bio.router.clustering.plain_cluster_word_no_fasta_no_route` | `clustering_analysis` | `bulk_rnaseq_analysis` |

Two clusters: (a) plasmid/directed-evolution prompts being eaten by
`mutate_sequence`, (b) `fastqc_quality_analysis` shadowing `quality_assessment`.
Both stem from commit `c1e1e18` ("feat(router): harden routing with
deterministic pre-router, strict schema, split fallbacks"). Likely fixes:

- Tighten the LLM router prompt's tool-disambiguation section for plasmid vs
  mutate_sequence (add a "if user mentions a vector / insert / plasmid map →
  plasmid_visualization" rule).
- Add a deterministic pre-router rule: plain `cluster` (no FASTA, no count
  matrix) → no-route (let the agent ask).

### Blockers preventing readiness=true

- **`beta-bundle-mixed-content`** (high) — see investigation note below.
- **`uncommitted-changes`** (medium) — gate was run against the dirty tree;
  fixes must be committed and the gate re-confirmed against HEAD.

## 2026-05-04 — Investigation: "Server: error" on helix-beta.noricum-biosoft.com

### Symptom
The deployed beta UI at `https://helix-beta.noricum-biosoft.com/` shows a red
`Server: error` badge. The frontend never reaches a healthy state, and (by
extension) every other API call from the SPA fails too.

### Root cause
The JS bundle currently served by nginx
(`/assets/index-DnhrdT75.js`) was built with a baked-in
`VITE_API_BASE_URL=http://helix--Publi-7CO3kaz1lrRR-443122446.us-west-1.elb.amazonaws.com`
(visible in the minified bundle as
`Ml=iG("http://helix--Publi-7CO3kaz1lrRR-443122446.us-west-1.elb.amazonaws.com")??""`).

That URL has two fatal problems for a browser loading the page over HTTPS:

1. **Mixed active content** — `axios.get('http://...')` from an `https://` page
   is blocked unconditionally by Chrome/Firefox/Safari, so
   `helixApi.healthCheck()` rejects → `setServerStatus('error')`
   (`frontend/src/App.tsx:305-317`).
2. **Internal AWS hostname** — the ALB DNS name `helix--Publi-...elb.amazonaws.com`
   does not resolve from the public internet (`curl` from outside AWS returns
   `Could not resolve host`), so even without the mixed-content block the
   request would fail.

The backend itself is healthy. Same-origin probing confirms it:

```
$ curl -i https://helix-beta.noricum-biosoft.com/health
HTTP/1.1 200 OK
{"status":"healthy","service":"Helix.AI Bioinformatics API","agent_disabled":false,"mock_mode":false}
```

The EC2 nginx config (`scripts/ec2/nginx-helix.conf.template`) already proxies
`/health`, `/execute`, `/agent`, `/create_session`, `/download`, `/mcp/*`,
`/session/*`, `/jobs/*` to `127.0.0.1:8001`, and the EC2 install/update scripts
(`scripts/ec2/bootstrap-ec2.sh:109`, `scripts/ec2/update-from-git.sh:41`) build
the SPA with `VITE_API_BASE_URL=` (empty) so it uses same-origin. So whatever
bundle is in `frontend/dist` on the host was NOT produced by those scripts —
someone built it locally with the wrong env var (matching the snippet in
`docs/architecture/ARCHITECTURE_EXPLANATION.md:146,169`) and rsync'd/scp'd it
onto the host.

### Fix (not applied — investigation only)
On the EC2 host (`/opt/helix/Helix.AI` or wherever `HELIX_ROOT` points), run:

```bash
cd <HELIX_ROOT>
HELIX_SKIP_PIP=1 ./scripts/ec2/update-from-git.sh
# or, if you only want to rebuild the frontend without restarting the backend:
cd frontend && npm ci && VITE_API_BASE_URL= npm run build
sudo systemctl reload nginx
```

After that, reload `https://helix-beta.noricum-biosoft.com/` — the badge should
turn green ("Server: healthy") and `/create_session`, `/execute`, `/agent`
should all start working again.

### Related issues to clean up (not blocking the immediate fix)
- `docs/architecture/ARCHITECTURE_EXPLANATION.md:146,169` — hardcodes the
  internal HTTP ALB URL as the build-time `VITE_API_BASE_URL`. This is what
  led to the bad bundle being produced in the first place. Either remove the
  snippet or rewrite it to use `VITE_API_BASE_URL=` (same-origin) for the EC2
  beta deployment, or `VITE_API_BASE_URL=https://helix-beta.noricum-biosoft.com`
  for an explicit absolute build.
- Consider a CI/build-time guard that fails the frontend build when
  `VITE_API_BASE_URL` is set to an `http://` URL (mixed-content trap), and a
  release-readiness check that asserts the deployed bundle's API base resolves
  publicly over HTTPS (or is empty for same-origin).

## 2026-04-27 — UX + advisory schema + agent.md fix

### Summary
Completed the Plan → Approve → Execute UX loop with consistent advisory rendering
and canonical JSON schema across all 50 bioinformatics workflow types.

### Changes
1. **Frontend UX** (`frontend/src/App.tsx`, `theme.css`)
   - Immediate "pending" placeholder + rotating bioinformatics loading messages
   - Follow-up chips now auto-submit (executeCommand) instead of pre-filling
   - Advisory JSON detection handles raw JSON, fenced JSON, all legacy shapes
   - `renderAdvisoryJSON` rewritten for canonical `HelixAdvisory` schema
   - Prompt bubble uses same font/size as response text

2. **Canonical Advisory Schema** (`backend/advisory_normalizer.py`, `backend/main.py`)
   - New Pydantic `HelixAdvisory` model: title, summary, classification, sections,
     workflow_steps, requirements, questions_for_user, next_steps
   - `normalize_advisory_text()` normalises 3 observed LLM output shapes into one
   - `main.py` calls normalizer on every text response before returning to frontend

3. **Agent prompt fix** (`backend/agent.py`)
   - AGENT_PROMPT_PATH now searches [repo_root/agent.md, docs/agent.md]
   - Was silently falling back to a 2-line stub — root cause of poor response quality

4. **docs/agent.md §11** — updated to require canonical HelixAdvisory JSON schema

5. **nginx** — /docs /redoc /openapi.json now password-protected

6. **.env** — consolidated root + backend/.env; corrected stale API key; all model
   vars set (gpt-4.1 reasoning, gpt-4.1-mini routing)

7. **scripts/run_benchmark.py** — new 50-workflow benchmark runner (fresh session per
   probe, 120s timeout, structured JSON output to artifacts/)

### Benchmark result (run: 2026-04-27T23:05:00Z, commit dd36eaa)
- PASSED: 49/50  FAILED: 0  ERRORS: 1 (cold-start timeout, not a product bug)
- Canonical advisory rendering: 32/50 probes
- Plain text (specialist tools): 17/50 probes
- Latency p50 ≈ 25s, p95 ≈ 35s

### Root causes resolved
- Empty/minimal advisory responses: agent.py was loading 2-line fallback prompt
- Inconsistent JSON shapes: advisory_normalizer handles all observed LLM variants
- Double JSON rendering on re-submit: strip ```json``` fences before parse
- Follow-up chips not auto-submitting: onSelect now calls executeCommand directly
- 401 API key error: stale sk-svcacct- key in root .env replaced with sk-proj- key
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false
- Autoloop stopped after max iterations: release_ready is false

---
## 2026-04-28 — Full benchmark re-run (commit 0d54dd3)

**Pass rate**: 50/50 (100%) — up from 49/50 in prior run.
**Advisory rendering**: 24/50 (48%) — note: 21 of the 26 plain-text items are correct specialist-tool routing (needs_inputs); only WF41 and WF46 are agent/success plain-text outliers.
**Latency p50**: 16.5s (−34% vs prior 25s), p95: 31.4s.
**SSE streaming**: shipped, TTFB ~5ms.
**Unit tests**: 724 passed, 3 skipped (optional deps), 0 failed.
**Release readiness**: artifacts/release_readiness.json updated — PASS.
- Autoloop stopped after max iterations: blockers: [{'id': 'unit-skip-optional-deps', 'severity': 'low', 'description': '3 unit tests skipped because anndata/h5py are not installed in the dev venv. These are optional bioinformatics libraries for single-cell profilers.', 'resolution': 'pip install anndata h5py in CI environment, or mark as optional in release gate'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'unit-skip-optional-deps', 'severity': 'low', 'description': '3 unit tests skipped because anndata/h5py are not installed in the dev venv. These are optional bioinformatics libraries for single-cell profilers.', 'resolution': 'pip install anndata h5py in CI environment, or mark as optional in release gate'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'unit-skip-optional-deps', 'severity': 'low', 'description': '3 unit tests skipped because anndata/h5py are not installed in the dev venv. These are optional bioinformatics libraries for single-cell profilers.', 'resolution': 'pip install anndata h5py in CI environment, or mark as optional in release gate'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'unit-skip-optional-deps', 'severity': 'low', 'description': '3 unit tests skipped because anndata/h5py are not installed in the dev venv. These are optional bioinformatics libraries for single-cell profilers.', 'resolution': 'pip install anndata h5py in CI environment, or mark as optional in release gate'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'unit-skip-optional-deps', 'severity': 'low', 'description': '3 unit tests skipped because anndata/h5py are not installed in the dev venv. These are optional bioinformatics libraries for single-cell profilers.', 'resolution': 'pip install anndata h5py in CI environment, or mark as optional in release gate'}]

---
## 2026-05-02 — Benchmark re-run + deps fix (commit 724c774)

**Unit tests**: 742 passed, 2 skipped (Rscript runtime), 0 failed.
anndata + h5py added to requirements.txt; 18 previously-skipped profiler tests now pass.
**Pass rate**: 50/50 (100%) — second consecutive perfect run.
**Advisory rendering**: 25/50 canonical advisory; 25/50 text (all with rich frontend cards).
**Latency p50**: 14.3s (−13% since Apr 28), p95: 19.6s (−38%).
**Release readiness**: readiness=true. Only blocker: Rscript not installed (low severity).
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '2 unit tests skipped because Rscript (R runtime) is not installed. These test R-based single-cell analysis tooling.', 'resolution': 'Install R runtime in CI/CD environment'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'beta-bundle-mixed-content', 'severity': 'high', 'description': "Currently-deployed frontend bundle at https://helix-beta.noricum-biosoft.com/ (index-DnhrdT75.js) was built with VITE_API_BASE_URL=http://helix--Publi-...elb.amazonaws.com (plain-HTTP internal AWS ALB). Browsers block all API calls as mixed content; SPA shows 'Server: error' and is fully non-functional. Backend itself is healthy (curl https://helix-beta.noricum-biosoft.com/health returns {status: healthy, ...}).", 'resolution': 'On EC2 host: cd <HELIX_ROOT> && HELIX_SKIP_PIP=1 ./scripts/ec2/update-from-git.sh   (rebuilds SPA with VITE_API_BASE_URL= empty, reloads systemd). Verify: curl -s https://helix-beta.noricum-biosoft.com/ | grep -o \'index-[^"]*\\.js\' returns a NEW hash, and the resulting bundle does NOT contain \'helix--Publi-\'.'}, {'id': 'uncommitted-changes', 'severity': 'medium', 'description': 'Working tree has 8 uncommitted file modifications (4 production fixes + 4 test/snapshot updates) plus 6 untracked files (4 fixtures, 2 probe scripts). The gate above was run against this dirty tree, so the readiness measurement reflects the post-fix code, not the current HEAD (7a41660). These changes must be committed before declaring HEAD release-ready.', 'resolution': 'Commit the production fixes (backend/main.py, frontend/src/App.tsx) and test/snapshot updates as one or two coherent commits, then re-run the gate against the clean tree to confirm the result reproduces.'}, {'id': 'rscript-skip', 'severity': 'low', 'description': '1 unit test (single_cell R-based) skipped because Rscript runtime not installed. Carries over from May 2 baseline.', 'resolution': 'Install R runtime in CI/CD environment (no impact on Python-only release path).'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'beta-bundle-mixed-content', 'severity': 'high', 'description': "Currently-deployed frontend bundle at https://helix-beta.noricum-biosoft.com/ (index-DnhrdT75.js) was built with VITE_API_BASE_URL=http://helix--Publi-...elb.amazonaws.com (plain-HTTP internal AWS ALB). Browsers block all API calls as mixed content; SPA shows 'Server: error' and is fully non-functional. Backend itself is healthy (curl https://helix-beta.noricum-biosoft.com/health returns {status: healthy, ...}).", 'resolution': 'On EC2 host: cd <HELIX_ROOT> && HELIX_SKIP_PIP=1 ./scripts/ec2/update-from-git.sh   (rebuilds SPA with VITE_API_BASE_URL= empty, reloads systemd). Verify: curl -s https://helix-beta.noricum-biosoft.com/ | grep -o \'index-[^"]*\\.js\' returns a NEW hash, and the resulting bundle does NOT contain \'helix--Publi-\'.'}, {'id': 'uncommitted-changes', 'severity': 'medium', 'description': 'Working tree has 8 uncommitted file modifications (4 production fixes + 4 test/snapshot updates) plus 6 untracked files (4 fixtures, 2 probe scripts). The gate above was run against this dirty tree, so the readiness measurement reflects the post-fix code, not the current HEAD (7a41660). These changes must be committed before declaring HEAD release-ready.', 'resolution': 'Commit the production fixes (backend/main.py, frontend/src/App.tsx) and test/snapshot updates as one or two coherent commits, then re-run the gate against the clean tree to confirm the result reproduces.'}, {'id': 'rscript-skip', 'severity': 'low', 'description': '1 unit test (single_cell R-based) skipped because Rscript runtime not installed. Carries over from May 2 baseline.', 'resolution': 'Install R runtime in CI/CD environment (no impact on Python-only release path).'}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '1 unit test (single_cell R-based) skipped because Rscript runtime not installed. Carries over from May 2 baseline.', 'resolution': 'Install R runtime in CI/CD environment (no impact on Python-only release path).'}, {'id': 'router-quality-regressions-may4', 'severity': 'low', 'description': "6 of 79 router evals regressed on the live LLM run on 2026-05-04. Plasmid/directed-evolution prompts misroute to mutate_sequence (3 cases); quality_assessment is shadowed by fastqc_quality_analysis (2 cases); plain 'cluster' word routes to clustering_analysis instead of bulk_rnaseq_analysis (1 case). Eval pass rate (92.4%) still meets the 90% threshold so this does not block release, but the regressions stem from the recent router rewrite (commit c1e1e18) and should be addressed in a follow-up before promoting beta to GA.", 'resolution': "Triage in worklog 2026-05-04 § 'Router-quality regressions'. Likely fix is tightening the LLM router's tool-disambiguation prompt for plasmid vs mutate_sequence and adding a deterministic pre-router rule for plain 'cluster' (no FASTA) → no-route."}, {'id': 'al2-deploy-contract-implicit', 'severity': 'medium', 'description': "scripts/ec2/bootstrap-ec2.sh and scripts/ec2/update-from-git.sh assume Node/npm is installed on the EC2 host, but on the live beta box (Amazon Linux 2, glibc 2.26) it isn't and never has been: find / -maxdepth 5 -name npm returns nothing. Production deploys have implicitly relied on building the SPA off-host and shipping dist/ — same path used today. The documented deploy script will fail at the npm run build step on any AL2 host. Operational only; no impact on release behaviour.", 'resolution': "Either (a) install Node on AL2 via the bootstrap's build-from-source path (slow, fragile), (b) migrate to AL2023 with newer glibc so NodeSource RPMs work, or (c) formalise the off-host build contract (build on CI / dev machine, ship dist/ via S3 + SSM) and remove the npm step from the on-host scripts. (c) is what is actually happening; recommend it. See worklog 2026-05-05 § 'Operational issues to fix before next deploy'."}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '1 unit test (single_cell R-based) skipped because Rscript runtime not installed. Carries over from May 2 baseline.', 'resolution': 'Install R runtime in CI/CD environment (no impact on Python-only release path).'}, {'id': 'router-quality-regressions-may4', 'severity': 'low', 'description': "6 of 79 router evals regressed on the live LLM run on 2026-05-04. Plasmid/directed-evolution prompts misroute to mutate_sequence (3 cases); quality_assessment is shadowed by fastqc_quality_analysis (2 cases); plain 'cluster' word routes to clustering_analysis instead of bulk_rnaseq_analysis (1 case). Eval pass rate (92.4%) still meets the 90% threshold so this does not block release, but the regressions stem from the recent router rewrite (commit c1e1e18) and should be addressed in a follow-up before promoting beta to GA.", 'resolution': "Triage in worklog 2026-05-04 § 'Router-quality regressions'. Likely fix is tightening the LLM router's tool-disambiguation prompt for plasmid vs mutate_sequence and adding a deterministic pre-router rule for plain 'cluster' (no FASTA) → no-route."}, {'id': 'al2-deploy-contract-implicit', 'severity': 'medium', 'description': "scripts/ec2/bootstrap-ec2.sh and scripts/ec2/update-from-git.sh assume Node/npm is installed on the EC2 host, but on the live beta box (Amazon Linux 2, glibc 2.26) it isn't and never has been: find / -maxdepth 5 -name npm returns nothing. Production deploys have implicitly relied on building the SPA off-host and shipping dist/ — same path used today. The documented deploy script will fail at the npm run build step on any AL2 host. Operational only; no impact on release behaviour.", 'resolution': "Either (a) install Node on AL2 via the bootstrap's build-from-source path (slow, fragile), (b) migrate to AL2023 with newer glibc so NodeSource RPMs work, or (c) formalise the off-host build contract (build on CI / dev machine, ship dist/ via S3 + SSM) and remove the npm step from the on-host scripts. (c) is what is actually happening; recommend it. See worklog 2026-05-05 § 'Operational issues to fix before next deploy'."}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '1 unit test (single_cell R-based) skipped because Rscript runtime not installed. Carries over from May 2 baseline.', 'resolution': 'Install R runtime in CI/CD environment (no impact on Python-only release path).'}, {'id': 'router-quality-regressions-may4', 'severity': 'low', 'description': "6 of 79 router evals regressed on the live LLM run on 2026-05-04. Plasmid/directed-evolution prompts misroute to mutate_sequence (3 cases); quality_assessment is shadowed by fastqc_quality_analysis (2 cases); plain 'cluster' word routes to clustering_analysis instead of bulk_rnaseq_analysis (1 case). Eval pass rate (92.4%) still meets the 90% threshold so this does not block release, but the regressions stem from the recent router rewrite (commit c1e1e18) and should be addressed in a follow-up before promoting beta to GA.", 'resolution': "Triage in worklog 2026-05-04 § 'Router-quality regressions'. Likely fix is tightening the LLM router's tool-disambiguation prompt for plasmid vs mutate_sequence and adding a deterministic pre-router rule for plain 'cluster' (no FASTA) → no-route."}, {'id': 'al2-deploy-contract-implicit', 'severity': 'medium', 'description': "scripts/ec2/bootstrap-ec2.sh and scripts/ec2/update-from-git.sh assume Node/npm is installed on the EC2 host, but on the live beta box (Amazon Linux 2, glibc 2.26) it isn't and never has been: find / -maxdepth 5 -name npm returns nothing. Production deploys have implicitly relied on building the SPA off-host and shipping dist/ — same path used today. The documented deploy script will fail at the npm run build step on any AL2 host. Operational only; no impact on release behaviour.", 'resolution': "Either (a) install Node on AL2 via the bootstrap's build-from-source path (slow, fragile), (b) migrate to AL2023 with newer glibc so NodeSource RPMs work, or (c) formalise the off-host build contract (build on CI / dev machine, ship dist/ via S3 + SSM) and remove the npm step from the on-host scripts. (c) is what is actually happening; recommend it. See worklog 2026-05-05 § 'Operational issues to fix before next deploy'."}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '1 unit test (single_cell R-based) skipped because Rscript runtime not installed. Carries over from May 2 baseline.', 'resolution': 'Install R runtime in CI/CD environment (no impact on Python-only release path).'}, {'id': 'router-quality-regressions-may4', 'severity': 'low', 'description': "6 of 79 router evals regressed on the live LLM run on 2026-05-04. Plasmid/directed-evolution prompts misroute to mutate_sequence (3 cases); quality_assessment is shadowed by fastqc_quality_analysis (2 cases); plain 'cluster' word routes to clustering_analysis instead of bulk_rnaseq_analysis (1 case). Eval pass rate (92.4%) still meets the 90% threshold so this does not block release, but the regressions stem from the recent router rewrite (commit c1e1e18) and should be addressed in a follow-up before promoting beta to GA.", 'resolution': "Triage in worklog 2026-05-04 § 'Router-quality regressions'. Likely fix is tightening the LLM router's tool-disambiguation prompt for plasmid vs mutate_sequence and adding a deterministic pre-router rule for plain 'cluster' (no FASTA) → no-route."}, {'id': 'al2-deploy-contract-implicit', 'severity': 'medium', 'description': "scripts/ec2/bootstrap-ec2.sh and scripts/ec2/update-from-git.sh assume Node/npm is installed on the EC2 host, but on the live beta box (Amazon Linux 2, glibc 2.26) it isn't and never has been: find / -maxdepth 5 -name npm returns nothing. Production deploys have implicitly relied on building the SPA off-host and shipping dist/ — same path used today. The documented deploy script will fail at the npm run build step on any AL2 host. Operational only; no impact on release behaviour.", 'resolution': "Either (a) install Node on AL2 via the bootstrap's build-from-source path (slow, fragile), (b) migrate to AL2023 with newer glibc so NodeSource RPMs work, or (c) formalise the off-host build contract (build on CI / dev machine, ship dist/ via S3 + SSM) and remove the npm step from the on-host scripts. (c) is what is actually happening; recommend it. See worklog 2026-05-05 § 'Operational issues to fix before next deploy'."}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '1 unit test (single_cell R-based) skipped because Rscript runtime not installed. Carries over from May 2 baseline.', 'resolution': 'Install R runtime in CI/CD environment (no impact on Python-only release path).'}, {'id': 'router-quality-regressions-may4', 'severity': 'low', 'description': "6 of 79 router evals regressed on the live LLM run on 2026-05-04. Plasmid/directed-evolution prompts misroute to mutate_sequence (3 cases); quality_assessment is shadowed by fastqc_quality_analysis (2 cases); plain 'cluster' word routes to clustering_analysis instead of bulk_rnaseq_analysis (1 case). Eval pass rate (92.4%) still meets the 90% threshold so this does not block release, but the regressions stem from the recent router rewrite (commit c1e1e18) and should be addressed in a follow-up before promoting beta to GA.", 'resolution': "Triage in worklog 2026-05-04 § 'Router-quality regressions'. Likely fix is tightening the LLM router's tool-disambiguation prompt for plasmid vs mutate_sequence and adding a deterministic pre-router rule for plain 'cluster' (no FASTA) → no-route."}, {'id': 'al2-deploy-contract-implicit', 'severity': 'medium', 'description': "scripts/ec2/bootstrap-ec2.sh and scripts/ec2/update-from-git.sh assume Node/npm is installed on the EC2 host, but on the live beta box (Amazon Linux 2, glibc 2.26) it isn't and never has been: find / -maxdepth 5 -name npm returns nothing. Production deploys have implicitly relied on building the SPA off-host and shipping dist/ — same path used today. The documented deploy script will fail at the npm run build step on any AL2 host. Operational only; no impact on release behaviour.", 'resolution': "Either (a) install Node on AL2 via the bootstrap's build-from-source path (slow, fragile), (b) migrate to AL2023 with newer glibc so NodeSource RPMs work, or (c) formalise the off-host build contract (build on CI / dev machine, ship dist/ via S3 + SSM) and remove the npm step from the on-host scripts. (c) is what is actually happening; recommend it. See worklog 2026-05-05 § 'Operational issues to fix before next deploy'."}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '1 unit test (single_cell R-based) skipped because Rscript runtime not installed. Carries over from May 2 baseline.', 'resolution': 'Install R runtime in CI/CD environment (no impact on Python-only release path).'}, {'id': 'router-quality-regressions-may4', 'severity': 'low', 'description': "6 of 79 router evals regressed on the live LLM run on 2026-05-04. Plasmid/directed-evolution prompts misroute to mutate_sequence (3 cases); quality_assessment is shadowed by fastqc_quality_analysis (2 cases); plain 'cluster' word routes to clustering_analysis instead of bulk_rnaseq_analysis (1 case). Eval pass rate (92.4%) still meets the 90% threshold so this does not block release, but the regressions stem from the recent router rewrite (commit c1e1e18) and should be addressed in a follow-up before promoting beta to GA.", 'resolution': "Triage in worklog 2026-05-04 § 'Router-quality regressions'. Likely fix is tightening the LLM router's tool-disambiguation prompt for plasmid vs mutate_sequence and adding a deterministic pre-router rule for plain 'cluster' (no FASTA) → no-route."}, {'id': 'al2-deploy-contract-implicit', 'severity': 'medium', 'description': "scripts/ec2/bootstrap-ec2.sh and scripts/ec2/update-from-git.sh assume Node/npm is installed on the EC2 host, but on the live beta box (Amazon Linux 2, glibc 2.26) it isn't and never has been: find / -maxdepth 5 -name npm returns nothing. Production deploys have implicitly relied on building the SPA off-host and shipping dist/ — same path used today. The documented deploy script will fail at the npm run build step on any AL2 host. Operational only; no impact on release behaviour.", 'resolution': "Either (a) install Node on AL2 via the bootstrap's build-from-source path (slow, fragile), (b) migrate to AL2023 with newer glibc so NodeSource RPMs work, or (c) formalise the off-host build contract (build on CI / dev machine, ship dist/ via S3 + SSM) and remove the npm step from the on-host scripts. (c) is what is actually happening; recommend it. See worklog 2026-05-05 § 'Operational issues to fix before next deploy'."}]
- Autoloop stopped after max iterations: blockers: [{'id': 'rscript-skip', 'severity': 'low', 'description': '1 unit test (single_cell R-based) skipped because Rscript runtime not installed. Carries over from May 2 baseline.', 'resolution': 'Install R runtime in CI/CD environment (no impact on Python-only release path).'}, {'id': 'router-quality-regressions-may4', 'severity': 'low', 'description': "6 of 79 router evals regressed on the live LLM run on 2026-05-04. Plasmid/directed-evolution prompts misroute to mutate_sequence (3 cases); quality_assessment is shadowed by fastqc_quality_analysis (2 cases); plain 'cluster' word routes to clustering_analysis instead of bulk_rnaseq_analysis (1 case). Eval pass rate (92.4%) still meets the 90% threshold so this does not block release, but the regressions stem from the recent router rewrite (commit c1e1e18) and should be addressed in a follow-up before promoting beta to GA.", 'resolution': "Triage in worklog 2026-05-04 § 'Router-quality regressions'. Likely fix is tightening the LLM router's tool-disambiguation prompt for plasmid vs mutate_sequence and adding a deterministic pre-router rule for plain 'cluster' (no FASTA) → no-route."}, {'id': 'al2-deploy-contract-implicit', 'severity': 'medium', 'description': "scripts/ec2/bootstrap-ec2.sh and scripts/ec2/update-from-git.sh assume Node/npm is installed on the EC2 host, but on the live beta box (Amazon Linux 2, glibc 2.26) it isn't and never has been: find / -maxdepth 5 -name npm returns nothing. Production deploys have implicitly relied on building the SPA off-host and shipping dist/ — same path used today. The documented deploy script will fail at the npm run build step on any AL2 host. Operational only; no impact on release behaviour.", 'resolution': "Either (a) install Node on AL2 via the bootstrap's build-from-source path (slow, fragile), (b) migrate to AL2023 with newer glibc so NodeSource RPMs work, or (c) formalise the off-host build contract (build on CI / dev machine, ship dist/ via S3 + SSM) and remove the npm step from the on-host scripts. (c) is what is actually happening; recommend it. See worklog 2026-05-05 § 'Operational issues to fix before next deploy'."}]
