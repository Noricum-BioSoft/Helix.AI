#!/usr/bin/env bash
# Run the FastQC/MultiQC follow-up scenario on the localhost frontend.
#
# Prerequisites (start in separate terminals):
#   1. HELIX_MOCK_MODE=1 uvicorn backend.main_with_mcp:app --host 0.0.0.0 --port 8001
#   2. cd frontend && npm run dev
#
# Then run: ./scripts/run_frontend_e2e.sh

set -e
cd "$(dirname "$0")/.."

echo "Checking backend..."
HEALTH=$(curl -s http://localhost:8001/health) || { echo "Backend not running on 8001. Start with: HELIX_MOCK_MODE=1 uvicorn backend.main_with_mcp:app --host 0.0.0.0 --port 8001"; exit 1; }
AGENT_DISABLED=$(echo "$HEALTH" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('agent_disabled', d.get('mock_mode', False)))" 2>/dev/null || echo "False")
if [ "$AGENT_DISABLED" != "True" ]; then
  echo "Backend must have agent disabled for this E2E (deterministic validation responses)."
  echo "Restart with: HELIX_AGENT_DISABLED=1 uvicorn backend.main_with_mcp:app --host 0.0.0.0 --port 8001"
  echo "Or legacy:   HELIX_MOCK_MODE=1 uvicorn backend.main_with_mcp:app --host 0.0.0.0 --port 8001"
  exit 1
fi
echo "Backend agent_disabled=1 ✓"

echo "Checking frontend..."
curl -s http://localhost:5173 >/dev/null || { echo "Frontend not running on 5173. Start with: cd frontend && npm run dev"; exit 1; }

echo ""
python tests/e2e_frontend_fastqc_multiqc.py
