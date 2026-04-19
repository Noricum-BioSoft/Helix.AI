# System audit for Helix.AI

## Purpose
Map the Helix.AI architecture, identify critical workflows, locate weak testing areas, and produce a practical release-hardening plan.

## When to use
Use this skill when:
- starting a major stabilization effort,
- the codebase has unclear workflow boundaries,
- test coverage needs to be expanded,
- release readiness is unknown.

## Inputs
- repo structure
- current failing tests
- benchmark thresholds
- current docs

## Procedure
1. Inspect major backend runtime files:
   - `backend/main.py`
   - `backend/agent.py`
   - `backend/command_router.py`
   - `backend/execution_broker.py`
2. Inspect workflow/session files:
   - `backend/history_manager.py`
   - `backend/workflow_executor.py`
   - `backend/workflow_checkpoint.py`
3. Inspect domain pipeline directories:
   - `backend/orchestration/`
   - `backend/bio_pipeline/`
   - `backend/ds_pipeline/`
4. Inspect frontend workflow entry points:
   - `frontend/src/App.tsx`
   - `frontend/src/services/helixApi.ts`
   - key components under `frontend/src/components/`
5. Inspect tests under:
   - `tests/unit/`
   - `tests/integration/`
   - `tests/backend/`
   - `tests/frontend/`
   - `tests/evals/`
   - `tests/demo_scenarios/`
   - `tests/workflows/`
6. Write:
   - architecture summary,
   - critical user journeys,
   - missing test coverage,
   - benchmark gaps,
   - prioritized stabilization plan.

## Output
Save results to:
- `.cursor/plans/system-audit-plan.md`
- `artifacts/worklog.md`

## Success criteria
- critical runtime paths identified,
- test gaps identified,
- benchmark candidates proposed,
- next steps prioritized by release impact.
