# EC2 beta deployment (GitHub clone + pull)

This describes a **beta** setup where Helix runs on an EC2 instance **like a developer machine**: clone the repo from GitHub, configure `.env`, run the backend (and optionally **built** static frontend), and **update with `git pull`**.

> **We cannot run AWS or SSH for you.** Everything is automated **on the server** via scripts in **`scripts/ec2/`** after you launch EC2 and clone this repo.

---

## Quick start (on the EC2 instance)

```bash
sudo mkdir -p /opt/helix && sudo chown "$USER:$USER" /opt/helix
cd /opt/helix
git clone https://github.com/<your-org>/Helix.AI.git && cd Helix.AI

sudo ./scripts/ec2/bootstrap-ec2.sh
sudo nano /etc/helix/helix.env    # add OPENAI_API_KEY (or DEEPSEEK_API_KEY)
sudo systemctl start helix-backend nginx
```

Open **port 80** in the EC2 security group. Optional HTTPS: `sudo certbot --nginx -d your.domain` (after DNS → Elastic IP).

**Updates** (after you `git push`):

```bash
cd /opt/helix/Helix.AI && git pull && sudo ./scripts/ec2/update-from-git.sh
```

See **`scripts/ec2/README.md`** for file descriptions.

---

## Is this viable for beta?

**Yes**, with a few caveats:

| Aspect | Notes |
|--------|--------|
| **Workflow** | Push → `git pull` on EC2 → reinstall deps if needed → restart service. Common for small teams. |
| **Parity with “local”** | Same repo layout, `python -m backend.main_with_mcp` (or `uvicorn`), `sessions/` on disk (use `HELIX_SESSIONS_DIR` or default). **Do not** run `npm run dev` for beta testers on the internet; use **`npm run build`** + a web server. |
| **Secrets** | Never commit `.env`. Put API keys on the server only (`/opt/helix/.env` or systemd `EnvironmentFile`). |
| **HTTPS** | Browsers expect HTTPS for real users. Use **Caddy** or **Nginx + Let’s Encrypt** (or **ACM + load balancer** if you add one later). |
| **Updates** | `git pull` can leave a broken venv if dependencies change; run `./scripts/ec2/update-from-git.sh` (or install deps after every pull). |
| **Backups** | Sessions are on the instance disk. **EBS snapshots** or periodic `tar` to a safe place. |
| **Scale / HA** | Single instance; fine for beta. |

---

## Recommended layout on the EC2 instance

```
/opt/helix/Helix.AI/          # git clone
  .env                        # NOT in git; created on server
  .venv/                      # Python virtualenv (optional but recommended)
  sessions/                   # or use HELIX_SESSIONS_DIR=/var/lib/helix/sessions
  frontend/dist/              # after npm run build (if serving UI from this host)
```

Environment (local-only sessions, no S3 for session markers):

```bash
HELIX_LOCAL_SESSIONS_ONLY=1
HELIX_SESSIONS_DIR=/var/lib/helix/sessions
# plus OPENAI_API_KEY / DEEPSEEK_API_KEY from .env
```

---

## One-time setup (outline)

1. **Launch EC2**  
   - Amazon Linux 2023 or Ubuntu, **t4g.medium** or similar, **EBS** for root + enough space for sessions.  
   - Security group: **22** (your IP only), **80** / **443** (testers or `0.0.0.0/0` for beta).

2. **Install dependencies**  
   - Python 3.11+, `git`, `nginx` or Caddy, Node 18+ only if you **build the frontend on the server**.

3. **Clone**  
   ```bash
   sudo mkdir -p /opt/helix && sudo chown $USER:$USER /opt/helix
   cd /opt/helix
   git clone https://github.com/<org>/Helix.AI.git
   cd Helix.AI
   ```

   For a **private** repo, use a [deploy key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/managing-deploy-keys) or HTTPS with a **fine-grained PAT** (stored only on the server).

4. **Python venv + backend**  
   ```bash
   python3.11 -m venv .venv
   source .venv/bin/activate
   pip install -r backend/requirements.txt
   ```

5. **`.env`**  
   Copy from `backend/.env.example` or project docs; set at least one LLM API key. **Do not** commit.

6. **Automated install (recommended)**  
   From the repo root on the server: **`sudo ./scripts/ec2/bootstrap-ec2.sh`**  
   This creates `/var/lib/helix/sessions`, `/etc/helix/helix.env` from the example, Python **venv**, **frontend production build** (`VITE_API_BASE_URL` empty for same-origin), **systemd** `helix-backend`, and **nginx** to serve `frontend/dist` and proxy API routes to `127.0.0.1:8001`.

7. **HTTPS**  
   After DNS points to the instance: **`sudo certbot --nginx -d your.domain`**

---

## Updates (push → pull)

On the server (after you `git push`):

```bash
cd /opt/helix/Helix.AI
./scripts/ec2/update-from-git.sh
```

Or manually:

```bash
git pull
source .venv/bin/activate
pip install -r backend/requirements.txt   # when requirements change
# If you serve frontend from this host:
#   cd frontend && npm ci && VITE_API_BASE_URL= npm run build
sudo systemctl restart helix-backend
sudo systemctl reload nginx   # if you use nginx
```

---

## Differences from your laptop

| Local | EC2 beta |
|-------|-----------|
| `start.sh` runs both Vite dev + backend | Prefer **systemd + uvicorn** + **static** frontend (or backend-only + testers use a build). |
| `localhost` | Use **HTTPS** + DNS or Elastic IP. |
| Casual restarts | Use **`systemctl restart`** after pull. |

---

## Related

- [../scripts/ec2/README.md](../../scripts/ec2/README.md) — script index and copy-paste commands  
- [BETA_DEPLOYMENT_CHEAP.md](BETA_DEPLOYMENT_CHEAP.md) — cost and local vs S3 sessions  
- [ENVIRONMENT_SETUP.md](../getting-started/ENVIRONMENT_SETUP.md) — API keys and `HELIX_*` variables  
