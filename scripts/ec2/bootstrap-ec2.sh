#!/usr/bin/env bash
# One-time EC2 setup: packages, user `helix`, venv, frontend build, nginx, systemd.
#
# Prerequisites:
#   - Clone this repo first, e.g. sudo mkdir -p /opt/helix && sudo chown $USER /opt/helix
#     cd /opt/helix && git clone <your-repo-url> Helix.AI && cd Helix.AI
#   - Run this script FROM THE REPO ROOT with sudo:
#       sudo ./scripts/ec2/bootstrap-ec2.sh
#
# After:
#   1. sudo nano /etc/helix/helix.env   # add OPENAI_API_KEY (or DEEPSEEK)
#   2. sudo systemctl start helix-backend nginx
#   3. Open port 80 in the security group; optional: certbot for HTTPS

set -euo pipefail

if [[ "$(id -u)" -ne 0 ]]; then
  echo "Run with sudo: sudo $0" >&2
  exit 1
fi

HELIX_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$HELIX_ROOT"

echo "[bootstrap] HELIX_ROOT=$HELIX_ROOT"

if [[ ! -f "$HELIX_ROOT/backend/main_with_mcp.py" ]]; then
  echo "Run this from the Helix.AI repository root." >&2
  exit 1
fi

install_packages() {
  if command -v dnf >/dev/null 2>&1; then
    dnf install -y python3.11 python3.11-pip nginx git nodejs npm || {
      dnf install -y python3 python3-pip nginx git
      echo "Install Node.js 18+ manually if npm is missing: https://github.com/nodesource/distributions" >&2
    }
  elif [[ -f /etc/os-release ]] && grep -q 'Amazon Linux 2' /etc/os-release && command -v yum >/dev/null 2>&1; then
    # Amazon Linux 2 (yum; common on EC2) — Node via NodeSource
    yum install -y python3 python3-pip nginx git curl gcc python3-devel
    curl -fsSL https://rpm.nodesource.com/setup_20.x | bash -
    yum install -y nodejs
  elif command -v apt-get >/dev/null 2>&1; then
    apt-get update -y
    apt-get install -y python3.11 python3.11-venv python3-pip nginx git curl
    if ! command -v npm >/dev/null 2>&1; then
      curl -fsSL https://deb.nodesource.com/setup_20.x | bash -
      apt-get install -y nodejs
    fi
  else
    echo "Unsupported OS. Install Python 3.11, nginx, git, nodejs, npm manually." >&2
    exit 1
  fi
}

install_packages

if ! id helix >/dev/null 2>&1; then
  useradd -m -s /bin/bash helix
fi

mkdir -p /var/lib/helix/sessions /etc/helix
chown -R helix:helix /var/lib/helix
chown -R helix:helix "$HELIX_ROOT"

if [[ ! -f /etc/helix/helix.env ]]; then
  install -m 600 -o root -g root "$HELIX_ROOT/scripts/ec2/helix.env.example" /etc/helix/helix.env
  echo "[bootstrap] Created /etc/helix/helix.env — EDIT IT and add at least OPENAI_API_KEY"
else
  echo "[bootstrap] Keeping existing /etc/helix/helix.env"
fi

echo "[bootstrap] Python venv + pip install (as user helix)..."
sudo -u helix bash -c "
  cd '$HELIX_ROOT'
  if command -v python3.11 >/dev/null 2>&1; then PY=python3.11; else PY=python3; fi
  \"\$PY\" -m venv .venv
  source .venv/bin/activate
  pip install --upgrade pip
  pip install -r backend/requirements.txt
"

echo "[bootstrap] Frontend production build (same-origin API)..."
sudo -u helix bash -c "
  cd '$HELIX_ROOT/frontend'
  npm ci
  VITE_API_BASE_URL= npm run build
"

escape_sed() { printf '%s\n' "$1" | sed -e 's/[\/&]/\\&/g'; }
ER="$(escape_sed "$HELIX_ROOT")"

echo "[bootstrap] systemd unit..."
sed "s|__HELIX_ROOT__|$ER|g" "$HELIX_ROOT/scripts/ec2/helix-backend.service.template" > /etc/systemd/system/helix-backend.service
systemctl daemon-reload
systemctl enable helix-backend

echo "[bootstrap] nginx site..."
sed "s|__HELIX_ROOT__|$ER|g" "$HELIX_ROOT/scripts/ec2/nginx-helix.conf.template" > /etc/nginx/conf.d/helix.conf
# Default server may conflict; disable default on Amazon Linux if needed
rm -f /etc/nginx/conf.d/default.conf 2>/dev/null || true
nginx -t
systemctl enable nginx

echo ""
echo "=========================================="
echo "Bootstrap complete."
echo "1. Edit secrets:   sudo nano /etc/helix/helix.env"
echo "2. Start services: sudo systemctl start helix-backend nginx"
echo "3. Status:         sudo systemctl status helix-backend"
echo "4. Open EC2 security group inbound TCP 80 (and 443 after certbot)."
echo "5. Updates:        git pull && sudo ./scripts/ec2/update-from-git.sh"
echo "=========================================="
