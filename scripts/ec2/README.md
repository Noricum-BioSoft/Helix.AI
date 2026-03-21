# EC2 beta helpers

Run **[docs/deployment/EC2_GIT_BETA.md](../../docs/deployment/EC2_GIT_BETA.md)** on the server — we cannot run AWS/SSH steps from this repo.

## Files

| File | Purpose |
|------|--------|
| **`bootstrap-ec2.sh`** | **Run once** on the instance (with `sudo`) after `git clone`: installs nginx, node, creates `helix` user, venv, builds frontend, installs systemd + nginx configs. |
| **`update-from-git.sh`** | After each `git push`: `git pull`, `pip install`, optional frontend rebuild, `systemctl restart helix-backend`. |
| **`helix-backend.service.template`** | Systemd unit; `bootstrap-ec2.sh` writes `/etc/systemd/system/helix-backend.service`. |
| **`nginx-helix.conf.template`** | Nginx: static `frontend/dist` + proxy to `127.0.0.1:8001`. |
| **`helix.env.example`** | Copy to `/etc/helix/helix.env` — API keys + `HELIX_LOCAL_SESSIONS_ONLY`. |

## Quick sequence (on EC2)

```bash
# 1. Clone (use deploy key or HTTPS if private)
sudo mkdir -p /opt/helix && sudo chown "$USER:$USER" /opt/helix
cd /opt/helix
git clone https://github.com/<org>/Helix.AI.git
cd Helix.AI

# 2. Bootstrap (packages, venv, build, systemd, nginx)
sudo ./scripts/ec2/bootstrap-ec2.sh

# 3. Secrets
sudo nano /etc/helix/helix.env   # OPENAI_API_KEY=...

# 4. Start
sudo systemctl start helix-backend nginx
curl -s http://127.0.0.1/health

# 5. Security group: open TCP 80 (and 443 after HTTPS)
```

HTTPS: `sudo certbot --nginx -d your.domain` (after DNS points to the instance).

Updates after you push to GitHub:

```bash
cd /opt/helix/Helix.AI
git pull
sudo ./scripts/ec2/update-from-git.sh
```
