#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# scripts/setup-nextflow.sh
#
# One-time setup for the Helix.AI local Nextflow execution layer (Tier 3a).
#
# Prerequisites: Java 17+, Docker (running), bash
#
# Usage:
#   bash scripts/setup-nextflow.sh [--skip-docker-check]
# ---------------------------------------------------------------------------

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

SKIP_DOCKER=false
for arg in "$@"; do
  [[ "$arg" == "--skip-docker-check" ]] && SKIP_DOCKER=true
done

echo -e "${BLUE}=== Helix.AI Nextflow Setup ===${NC}"
echo ""

# ── 1. Java ──────────────────────────────────────────────────────────────
echo -n "Checking Java... "
if ! java -version 2>/dev/null; then
  echo -e "${RED}❌ Java not found.${NC}"
  echo "  Install Java 17+: brew install openjdk@17  (macOS) or apt install openjdk-17-jdk (Debian/Ubuntu)"
  exit 1
fi
echo -e "${GREEN}✅ Java found${NC}"

# ── 2. Nextflow ───────────────────────────────────────────────────────────
echo -n "Checking Nextflow... "
NF_TARGET="${HOME}/.local/bin/nextflow"
if command -v nextflow &>/dev/null; then
  echo -e "${GREEN}✅ Nextflow already in PATH: $(command -v nextflow)${NC}"
elif [[ -x "${NF_TARGET}" ]]; then
  echo -e "${GREEN}✅ Nextflow found at ${NF_TARGET}${NC}"
  echo "  Add ~/.local/bin to your PATH to use it from anywhere."
else
  echo -e "${YELLOW}⬇  Installing Nextflow to ${NF_TARGET}...${NC}"
  mkdir -p "${HOME}/.local/bin"
  curl -s https://get.nextflow.io | bash
  mv nextflow "${NF_TARGET}"
  chmod +x "${NF_TARGET}"
  echo -e "${GREEN}✅ Nextflow installed at ${NF_TARGET}${NC}"
  echo "  Run: export PATH=\"\${HOME}/.local/bin:\${PATH}\""
  echo "  Or add the above line to your ~/.zshrc / ~/.bashrc."
fi

# Verify version
NF_BIN="${NF_TARGET}"
command -v nextflow &>/dev/null && NF_BIN="nextflow"
echo "  Version: $("${NF_BIN}" -version 2>&1 | grep 'version' | head -1)"

# ── 3. Docker ─────────────────────────────────────────────────────────────
if [[ "${SKIP_DOCKER}" == "false" ]]; then
  echo -n "Checking Docker... "
  if ! docker --version &>/dev/null; then
    echo -e "${RED}❌ Docker not found.${NC}"
    echo "  Install Docker Desktop from https://www.docker.com/products/docker-desktop"
    exit 1
  fi
  if ! docker info &>/dev/null; then
    echo -e "${YELLOW}⚠  Docker is installed but not running.${NC}"
    echo "  Start Docker Desktop and re-run this script."
    exit 1
  fi
  echo -e "${GREEN}✅ Docker running${NC}"
else
  echo -e "${YELLOW}⚠  Skipping Docker check (--skip-docker-check)${NC}"
fi

# ── 4. Python: sse-starlette ──────────────────────────────────────────────
echo -n "Checking sse-starlette... "
if python -c "import sse_starlette" 2>/dev/null; then
  echo -e "${GREEN}✅ sse-starlette already installed${NC}"
else
  echo -e "${YELLOW}⬇  Installing sse-starlette...${NC}"
  pip install "sse-starlette>=2.1.3,<3.0" --quiet
  echo -e "${GREEN}✅ sse-starlette installed${NC}"
fi

# ── 5. Verify backend modules ─────────────────────────────────────────────
echo -n "Verifying Helix Nextflow modules... "
cd "$(dirname "$0")/.."
python - <<'EOF'
from backend.nextflow_event_bus import get_event_bus
from backend.nextflow_executor import check_prerequisites, HELIX_TO_PIPELINE
prereqs = check_prerequisites()
missing = [k for k, v in prereqs.items() if not v]
if missing:
    print(f"WARN: missing prerequisites: {missing}")
else:
    print(f"OK — all prerequisites present, {len(HELIX_TO_PIPELINE)} pipeline tools registered")
EOF
echo -e "${GREEN}✅ Helix modules OK${NC}"

echo ""
echo -e "${GREEN}=== Setup complete ===${NC}"
echo ""
echo "Nextflow execution tier is ready."
echo "Pipeline tools will launch workflows locally using Docker containers."
echo ""
echo "Environment variable overrides:"
echo "  HELIX_NEXTFLOW_BIN       — path to nextflow binary (default: ~/.local/bin/nextflow)"
echo "  HELIX_NEXTFLOW_WORK_DIR  — Nextflow work directory (default: /tmp/helix_nextflow)"
echo "  HELIX_NEXTFLOW_RESULTS_DIR — results directory (default: /tmp/helix_results)"
echo "  HELIX_NF_LOCAL_MAX_SAMPLES — max samples for local execution (default: 4)"
echo "  HELIX_NF_LOCAL_MAX_GB    — max dataset size for local execution (default: 10.0)"
