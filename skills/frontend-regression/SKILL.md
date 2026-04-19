# Frontend regression for Helix.AI

## Purpose
Validate the React frontend for the main iterative workflow user journeys.

## When to use
Use this skill when:
- React components change,
- API integration changes,
- workflow state display changes,
- file upload or plot UX changes.

## Focus areas
- analysis request entry
- plan display or editing
- workflow launch
- progress/state visibility
- output rendering
- plot update flow
- file upload and additional file attachment
- error/retry/resume flows

## Procedure
1. Run frontend tests, lint, build, and type checks.
2. Identify the affected route and components.
3. Validate browser behavior for the changed flow.
4. Ensure UI state matches backend reality.
5. Add or update frontend tests.

## Success criteria
- no broken workflow path,
- state and errors remain visible,
- changed behavior is covered by tests,
- build and type checks pass.
