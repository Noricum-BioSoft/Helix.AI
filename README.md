# Helix.AI

Helix.AI is a web app for running bioinformatics workflows from natural language prompts.

## Quick start (local dev)

**Prereqs**: Python 3.10+, Node 18+

```bash
./start.sh
```

Frontend runs on `http://localhost:5173`, backend on `http://localhost:8001` (OpenAPI docs at `/docs`).

## Configuration

Copy `backend/.env.example` to `.env` (repo root) and fill in at least one API key:

- `OPENAI_API_KEY` or `DEEPSEEK_API_KEY`

## Repo structure

- `backend/`: FastAPI backend (MCP + orchestration)
- `frontend/`: React UI (Vite)
- `tools/`: tool implementations (bioinformatics)
- `scripts/`: utility + deployment scripts
- `infrastructure/`: AWS CDK stack
- `docs/`: documentation
- `tests/`: pytest suite

## Tests

```bash
pytest
cd frontend && npm run build
```

## Deployment

See `docs/deployment/AWS_DEPLOYMENT_GUIDE.md`.

