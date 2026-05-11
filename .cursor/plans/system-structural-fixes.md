# Plan: System Structural Fixes — Stop the "Whack-a-Mole"

**Status:** Proposed
**Goal:** Address the recurring "fix one thing, break another" pattern by tightening structural contracts rather than adding more heuristics. Complements `.cursor/plans/eliminate-keyword-routing.md` (which removes keyword routing); this plan removes the *layered architecture* that makes routing fragile in the first place.

---

## Patterns observed (evidence)

| # | Pattern | Evidence (session IDs) |
|---|---------|------------------------|
| 1 | **Too many parallel intent layers.** Deterministic router → planning gate → LLM router → agent BioinformaticsExecutor — all classifying the same message. | `d7ea20b1`: 3 routing calls, 19s of redundant agent overhead, same answer. |
| 2 | **Persistence is optional, not contractual.** Runs, history entries, artifacts, bundles each have separate write paths that may or may not fire. | `0d362d04`, `71034a68`: empty bundles. `2af19038`: missing conversational turns in `history[]`. Multiple sessions: `produced_artifacts: []` despite artifacts on disk. |
| 3 | **Response shape doesn't match user intent.** Advisory turns return blank text; meta questions return execution plans; conversational requests return `unknown_intent`. | `d7ea20b1` turn 1: `helix_type:advisory` AI message, `result.text=""`. `2af19038`: meta question → tabular_analysis plan. `584fce12`: "profile the workbook" → `unknown_intent`. |
| 4 | **Heuristics accumulate; boundaries don't tighten.** Each failure mode adds a substring or cue word. The real boundary ("code vs prose vs unknown") is never decided in one place. | `is_analytical_request` had 25+ cue words + `len>=8` fallback (now removed). Router has its own keyword exclusion list. Frontend has its own JSON unwrap heuristics. |
| 5 | **Old sessions can't migrate forward.** Each fix only helps fresh sessions; we have no schema version, no Reset UX. | Repeated "still broken" feedback testing in pre-fix sessions. |
| 6 | **Operational fragility hides product bugs.** Frontend build on EC2 fails (GLIBC, npm perms); deploys go via S3 tarball; backend deploys have no behavioral health check. | Earlier deploy thread: NodeSource/AL2 failure; only `git pull` for backend, no verification. |
| 7 | **Per-function tests; no per-flow tests.** Unit tests don't catch user-visible regressions across a multi-turn flow. | Every fix above had passing unit tests *before* the bug was reported. |

---

## Slices (small, shippable)

### Slice A — Persistence contract for a turn (highest impact)

**Single rule:** every `/execute` request that reaches a tool produces exactly one history entry, optionally one run record, with `produced_artifacts` always populated when artifacts were written to disk. No optional paths.

- Centralize history+run+artifact write in a single `record_turn(session_id, command, tool, result, artifacts) -> turn_record` helper called from exactly one place in `main.py` per execution branch (`/execute` and tabular-approval).
- Bundle generator and `/download/*` endpoints read **only** from `history_manager.storage_dir`. Remove any `Path(__file__).parent / "sessions"` references.
- Advisory / `tool=agent` turns must also write history entries (lightweight, no `produced_artifacts`).
- **Test:** `tests/unit/backend/test_turn_persistence_contract.py` — for each tool category (tabular, advisory, bio orchestrator), assert that after a successful call:
  - `len(session["history"])` increased by 1
  - `len(session["runs"])` increased by ≤1
  - if artifacts were produced, `produced_artifacts` is non-empty and every URI exists
  - bundle for that run is non-trivial (manifest + script ∨ at least one artifact)

### Slice B — Collapse intent layers

Decision tree becomes:
1. **Deterministic router** (regex-only, high-confidence patterns) → tool
2. **LLM router** → tool ∈ {known tools, `handle_natural_command`, `unknown_intent`, `multi_step_workflow`}
3. **Stop.**

Concrete changes:

- Delete the `is_analytical_request` planning gate in `main.py` (lines ~3273–3318). Move the tabular-plan staging into the post-router approval path: if the LLM router returns `tabular_analysis` for a multi-step command, *then* generate a plan.
- Remove the agent's `BioinformaticsExecutor` tool-mapping loop. Agent's job becomes "given tool X, produce arguments / code"; never "decide which tool."
- **Test:** `tests/integration/test_routing_layers_single.py` — for each archetypal command, assert the router is called exactly once per request (not 2–3 times as today).

### Slice C — One advisory response shape

- Define `AdvisoryResponse` Pydantic model: `{ title, summary, classification, options, next_steps }`.
- Backend: when an agent / `handle_natural_command` turn produces an advisory, the `result.text` field is the markdown rendering of `AdvisoryResponse`. The structured form lives under `result.advisory`.
- Frontend: render `result.advisory` when present (typed), fall back to `result.text` otherwise. Delete all the "is it in messages[ai].content as a JSON string?" heuristics in `App.tsx`.
- **Test:** `tests/integration/test_advisory_round_trip.py` — backend returns advisory; frontend test asserts the rendered card shows title/summary/options.

### Slice D — End-to-end flow tests, one per archetype

Each test is a multi-turn `/execute` sequence with assertions on user-visible response shape:

| Archetype | Flow |
|-----------|------|
| Tabular planning + approval | upload → "compare X across groups" → plan returned → approve → run → bundle includes inputs+script+plot |
| Tabular profile | upload → "profile the workbook" → tabular_analysis returns rows/cols summary in `text` |
| Meta question after run | (above) → "what's the difference between your next step and the previous analysis?" → response is advisory text, NOT a new plan |
| Advisory | "what can I do with this dataset?" → response is `AdvisoryResponse`, `result.text` non-empty |
| Unknown | "potato" → `unknown_intent` clarification card |

These tests live in `tests/integration/test_user_flows.py`. Each archetype that's broken in a real session adds a corresponding test, so regressions are caught at flow level.

### Slice E — Session schema version + Reset UX

- Add `schema_version` field on session creation. Bump on incompatible changes.
- Frontend "Reset session" button (already exists?) — verify it actually clears storage and starts fresh.
- Optional: when loading an older-version session, surface a banner: "This session was created on an older Helix version; some features (e.g. bundle downloads) may behave differently."

### Slice F — Deployment health check (low effort, high value)

- After `update-from-git.sh` restarts `helix-backend`, hit `/health` and `/execute` with a known-good fixture; fail loudly if response shape doesn't match.
- Print the deployed commit hash in `/health` JSON so we can confirm "is what I'm hitting actually `ea02f2d`?"

---

## Sequencing & estimates

| Slice | Why first | Rough size |
|-------|-----------|------------|
| A | Kills 4 bug classes (downloads, history, repro, "next-step" context) with one change | 1–2 days |
| B | Removes 19s+ latency, removes "three routers, three answers" | 1 day |
| C | User-visible "blank reply" goes away | 0.5 day |
| D | Catches future regressions at flow level | ongoing, start with 4 tests |
| E | Lets users self-recover; removes "test in stale session" confusion | 0.5 day |
| F | Cheap insurance against "did we actually deploy?" | 1 hour |

Suggested order: **A → B → C → F → D → E.** A and B are the structural fixes; everything else is enablement.

---

## Non-goals

- Adding more keyword cues to `is_analytical_request`.
- Adding more router rules.
- Rewriting the agent system prompt.
- Migrating old session JSONs in-place.

---

## Release-gate signal

After Slices A–D, `artifacts/release_readiness.json` should include a `flows_passing` count from `test_user_flows.py`. Release-readiness flips to `false` when any archetype regresses.

---

## Productization: Decision-Ready Workflow

### Positioning

Helix should be positioned as a **decision-ready bioinformatics workflow system**, not
as a generic coding/chat assistant.

Product promise:

> A lab lead can use Helix outputs directly in a decision/review meeting without
> manually reconstructing methods, assumptions, or reproducibility context.

### Decision packet (new primary artifact)

For every completed run/session, Helix should generate a compact "Decision Packet":

- Executive summary (3-6 bullets)
- Inputs used (with names/versions where available)
- Methods + key parameters
- Key plots/tables
- What changed vs previous run
- Caveats/limitations
- Re-run instructions + reproducibility bundle link(s)

This packet becomes the default output for review workflows.

### Decision-readiness checklist (must-pass)

Before a run can be marked "Review Ready", require:

- Non-empty run ledger entry
- Valid artifact registry (`produced_artifacts`) for generated outputs
- Reproducibility bundle generated successfully
- Caveats/assumptions present
- No execution errors

If any item fails: state remains "Draft" with explicit missing items.

### Review states (UI + API)

Add a small lifecycle for each run:

- `draft`
- `review_ready`
- `approved`
- `superseded`

Store state in run metadata so review status is queryable and auditable.

### Run-to-run diff (first-class)

Expose a built-in comparison for consecutive or selected runs:

- Parameter changes
- Input/data changes
- Output deltas (key metrics)
- Interpretation deltas

This is the core review aid for "what changed and why?" decisions.

### Metrics for this goal

Track and publish:

- `% runs marked review_ready`
- `% review_ready runs reopened due to missing context`
- median time from run completion -> review decision
- `% successful reruns from bundle`
- `% sessions with explicit caveats`

### Execution order (product-facing)

1. Decision Packet + checklist gate (visible trust gain)
2. Review states
3. Run-to-run diff view
4. Metrics instrumentation + dashboarding

This stream should run in parallel with Slices A-D above because those slices
provide the data contract needed for decision-ready outputs.
