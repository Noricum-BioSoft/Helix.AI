#!/usr/bin/env bash
# Update Helix.AI from GitHub and restart the backend.
# Run from the repo root on the EC2 instance, e.g. /opt/helix/Helix.AI
#
# Usage:
#   ./scripts/ec2/update-from-git.sh
#
# Optional env:
#   HELIX_SKIP_PIP=1   skip pip install (faster when only code changed)
#   HELIX_SKIP_FRONTEND=1  skip npm ci + build (when frontend unchanged)

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

echo "[helix] git pull..."
git pull --ff-only

if [[ "${HELIX_SKIP_PIP:-}" != "1" ]]; then
  echo "[helix] pip install..."
  if [[ -d "${ROOT}/.venv" ]]; then
    # shellcheck source=/dev/null
    source "${ROOT}/.venv/bin/activate"
    pip install -r backend/requirements.txt
  else
    echo "No .venv found; run scripts/ec2/bootstrap-ec2.sh first." >&2
    exit 1
  fi
fi

if [[ "${HELIX_SKIP_FRONTEND:-}" != "1" ]] && [[ -d "${ROOT}/frontend" ]]; then
  echo "[helix] frontend build..."
  if ! command -v npm >/dev/null 2>&1; then
    echo "npm not found; set HELIX_SKIP_FRONTEND=1 or install Node.js" >&2
    exit 1
  fi
  (
    cd "${ROOT}/frontend"
    npm ci
    VITE_API_BASE_URL= npm run build
  )
fi

echo "[helix] restart systemd service..."
if systemctl is-active --quiet helix-backend 2>/dev/null; then
  sudo systemctl restart helix-backend
  echo "[helix] helix-backend restarted."
else
  echo "Service helix-backend not active (install scripts/ec2/helix-backend.service first)." >&2
fi

echo "[helix] done."
