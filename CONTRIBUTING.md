# Contributing

## Development setup

- **Prereqs**: Python 3.10+, Node 18+
- **Install + run**:

```bash
./start.sh
```

## Project layout

- `backend/`: FastAPI backend
- `frontend/`: React (Vite) UI
- `tools/`: bioinformatics tool implementations
- `infrastructure/`: AWS CDK
- `scripts/`: ops / deployment helpers
- `tests/`: pytest suite

## Quality bar

- **Keep changes small and readable**
- **No secrets in git** (`.env`, keys, tokens)
- **Tests**:

```bash
pytest
cd frontend && npm run build
```

## Pull requests

- **Commit messages**: follow existing `feat:`, `fix:`, `docs:`, `refactor:`, `ci:` style
- **Include**: what changed, why, and how you tested

