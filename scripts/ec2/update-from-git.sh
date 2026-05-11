#!/usr/bin/env bash
# Update Helix.AI from GitHub and restart the backend.
# Run from the repo root on the EC2 instance, e.g. /opt/helix (or /opt/helix/Helix.AI).
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

# Conda-created envs under .venv often have bin/python but no bin/activate.
_helix_use_venv() {
  if [[ -f "${ROOT}/.venv/bin/activate" ]]; then
    # shellcheck source=/dev/null
    source "${ROOT}/.venv/bin/activate"
  elif [[ -x "${ROOT}/.venv/bin/python" ]]; then
    export PATH="${ROOT}/.venv/bin:${PATH}"
  else
    echo "No usable Python under ${ROOT}/.venv (expected bin/activate or bin/python)." >&2
    return 1
  fi
}

echo "[helix] git pull..."
git pull --ff-only

if [[ "${HELIX_SKIP_PIP:-}" != "1" ]]; then
  echo "[helix] pip install..."
  if [[ -d "${ROOT}/.venv" ]]; then
    _helix_use_venv
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


# ---------------------------------------------------------------------------
# Pre-deploy gate: data integrity checks (fast, no backend required)
# ---------------------------------------------------------------------------
echo "[helix] running pre-deploy data integrity checks..."
if [[ -d "${ROOT}/.venv" ]]; then
  _helix_use_venv
fi
if python -m pytest tests/unit/backend/test_demo_data_integrity.py \
    -m "not s3" \
    -q --tb=short 2>&1; then
  echo "[helix] data integrity checks passed."
else
  echo "ERROR: data integrity checks failed — aborting deploy." >&2
  echo "Run:  pytest tests/unit/backend/test_demo_data_integrity.py -v  for details." >&2
  exit 1
fi

echo "[helix] restart systemd service..."
if systemctl is-active --quiet helix-backend 2>/dev/null; then
  # Root can restart directly; deploy user (helix) typically needs passwordless sudo
  # (see /etc/sudoers.d/helix-backend-restart on managed hosts).
  if systemctl restart helix-backend 2>/dev/null; then
    echo "[helix] helix-backend restarted."
  elif sudo -n systemctl restart helix-backend 2>/dev/null; then
    echo "[helix] helix-backend restarted (sudo -n)."
  else
    echo "Service helix-backend is active but restart failed. Run: sudo systemctl restart helix-backend" >&2
    exit 1
  fi
else
  echo "Service helix-backend not active (install scripts/ec2/helix-backend.service first)." >&2
fi

# ---------------------------------------------------------------------------
# Post-deploy behavioral health check
# ---------------------------------------------------------------------------
_BACKEND_URL="${HELIX_BACKEND_URL:-http://localhost:8001}"
_HEALTH_RETRIES=10
_HEALTH_INTERVAL=3

echo "[helix] waiting for backend to become healthy at ${_BACKEND_URL} ..."
for i in $(seq 1 "${_HEALTH_RETRIES}"); do
  if curl -sf "${_BACKEND_URL}/health" -o /tmp/helix_health.json 2>/dev/null; then
    echo "[helix] /health OK:"
    python3 -c "import json,sys; d=json.load(open('/tmp/helix_health.json')); print('  status:', d.get('status')); print('  git_commit:', d.get('git_commit','unknown'))"
    break
  fi
  echo "[helix] health check attempt ${i}/${_HEALTH_RETRIES} failed — retrying in ${_HEALTH_INTERVAL}s..."
  sleep "${_HEALTH_INTERVAL}"
done

# Verify response shape with a known-good fixture request.
_FIXTURE_CMD="List available tools"
echo "[helix] probing /execute with fixture: '${_FIXTURE_CMD}' ..."
if curl -sf -X POST "${_BACKEND_URL}/execute" \
    -H "Content-Type: application/json" \
    -d "{\"command\": \"${_FIXTURE_CMD}\"}" \
    -o /tmp/helix_probe.json 2>/dev/null; then
  _PROBE_SUCCESS=$(python3 -c "import json; d=json.load(open('/tmp/helix_probe.json')); print('ok' if d.get('success') or d.get('status') in ('success','ok') else 'fail')" 2>/dev/null || echo "parse_error")
  if [[ "${_PROBE_SUCCESS}" == "ok" ]]; then
    echo "[helix] /execute probe passed — deployment verified."
  else
    echo "WARNING: /execute probe returned non-success (${_PROBE_SUCCESS}). Check logs." >&2
  fi
else
  echo "WARNING: /execute probe failed (curl error). Backend may not be fully ready." >&2
fi

echo "[helix] done."
