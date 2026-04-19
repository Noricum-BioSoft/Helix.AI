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

### HTTPS — `helix-beta.noricum-biosoft.com`

1. **DNS (where `noricum-biosoft.com` is hosted)**  
   Create an **A record**: **name** `helix-beta` → **value** the instance’s **public IPv4** (today often `54.176.128.13`; if you stop/start the instance or use Elastic IP, update the record).

2. **Security group**  
   Inbound: **TCP 443** from `0.0.0.0/0` (and **80** if not already, for HTTP redirect / Let’s Encrypt).

3. **Nginx `server_name`** (so Certbot matches the vhost)  
   On the server, edit `/etc/nginx/conf.d/helix.conf` and set:
   ```nginx
   server_name helix-beta.noricum-biosoft.com;
   ```
   (Replace `server_name _;` if present.) Then:
   ```bash
   sudo nginx -t && sudo systemctl reload nginx
   ```

4. **Install Certbot** (not installed by default). Pick **one** block that matches the OS:

   **Amazon Linux 2023** (`dnf` available; `/etc/os-release` says “Amazon Linux 2023”):
   ```bash
   sudo dnf install -y certbot python3-certbot-nginx
   ```

   **Amazon Linux 2** (older; uses `yum` / `amazon-linux-extras`): **EPEL is required** — without it, `yum` reports `No package certbot available`.
   ```bash
   sudo amazon-linux-extras install epel -y
   sudo yum install -y certbot python3-certbot-nginx
   ```
   If install still fails, try: `sudo yum --enablerepo=epel install -y certbot python3-certbot-nginx`

   **If `python3-certbot-nginx` is not found** (common on AL2 when EPEL is hidden by repo priorities), install Certbot + the Nginx plugin via **pip**:
   ```bash
   sudo yum install -y python3-pip
   sudo python3 -m pip install --upgrade pip
   sudo python3 -m pip install certbot certbot-nginx
   ```
   **Amazon Linux 2 + pip:** system Python uses **OpenSSL 1.0.2**; **`urllib3` 2.x refuses to load** (`ImportError: urllib3 v2.0 only supports OpenSSL 1.1.1+`). Pin urllib3 to 1.x, then run Certbot:
   ```bash
   sudo python3 -m pip install 'urllib3<2'
   sudo /usr/local/bin/certbot plugins   # should list nginx
   sudo /usr/local/bin/certbot --nginx -d helix-beta.noricum-biosoft.com
   ```

   **Debian/Ubuntu:**
   ```bash
   sudo apt-get update && sudo apt-get install -y certbot python3-certbot-nginx
   ```

5. **Obtain certificate and wire Nginx**
   ```bash
   sudo certbot --nginx -d helix-beta.noricum-biosoft.com
   ```
   Follow prompts (email, agree to terms). Certbot will add TLS and usually redirect HTTP→HTTPS.

6. **Check**  
   Open `https://helix-beta.noricum-biosoft.com` — the browser should show a normal secure connection (no “Not secure” for plain HTTP to the IP).

**Note:** Use the **hostname** in the browser, not the raw IP, so the certificate matches.

---

HTTPS (generic): `sudo certbot --nginx -d your.domain` (after DNS points to the instance and `server_name` matches).

Updates after you push to GitHub:

```bash
cd /opt/helix/Helix.AI
git pull
sudo ./scripts/ec2/update-from-git.sh
```
