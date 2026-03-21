# EC2 beta helpers

Run **[docs/deployment/EC2_GIT_BETA.md](../../docs/deployment/EC2_GIT_BETA.md)** on the server — we cannot run AWS/SSH steps from this repo.

## Deploy from S3 via SSM (no git clone on the instance)

1. Create a clean tarball (avoid `git archive ... | gzip` in one pipe on some hosts — use `git archive` then `gzip` the `.tar` file), upload to `s3://…/bootstrap/helix-ai-source.tar.gz` and `ssm-deploy-from-s3.sh`.
2. **Amazon Linux 2:** Node 18 must be **built from source** (official linux-x64 binaries need glibc 2.28+; AL2 has 2.26). Builds can exceed **1 hour**.
3. **`AWS-RunShellScript` uses `executionTimeout` (default 3600 s).** You must pass it in **parameters**, not only `--timeout-seconds` on `send-command`, or the shell step stops after 1 hour.

```bash
aws ssm send-command \
  --region us-west-1 \
  --instance-ids i-xxxxxxxx \
  --document-name AWS-RunShellScript \
  --timeout-seconds 28800 \
  --parameters '{
    "commands":["aws s3 cp s3://YOUR_BUCKET/bootstrap/ssm-deploy-from-s3.sh /tmp/ssm-deploy.sh && chmod +x /tmp/ssm-deploy.sh && bash /tmp/ssm-deploy.sh YOUR_BUCKET bootstrap/helix-ai-source.tar.gz"],
    "executionTimeout":["28800"]
  }'
```

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
