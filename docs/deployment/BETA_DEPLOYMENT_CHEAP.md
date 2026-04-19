# Cheapest AWS deployment for beta testers

This guide compares options and recommends a **low-cost** layout for Helix.AI beta. It also explains why **“EC2 + S3 mounted for sessions”** is usually **not** the best fit for this codebase.

---

## How Helix actually uses storage

From `backend/history_manager.py` and related code:

| What | Where | Pattern |
|------|--------|--------|
| Session JSON (`history`, `runs`, metadata) | **`sessions/{session_id}/`** on **local disk** | Frequent **read/write** of small files (`_save_session`) |
| Large job outputs / EMR / uploads | **S3 via boto3** (prefix per session) | API calls, not a POSIX mount |
| Job registry | `sessions/jobs.json` | Local file |

The app is built around **normal filesystem semantics** (mkdir, atomic-ish JSON writes, optional `fcntl` locking). **Mounting S3** (s3fs, goofys, **Mountpoint for Amazon S3**) looks like a disk but:

- Many small writes to session JSON are **slow and expensive** (S3 request pricing + latency).
- **Consistency and locking** are not like a local disk; some tools break under heavy metadata churn.
- You still pay for **S3 requests**; for hot session data, **EBS is cheaper and simpler**.

**Recommendation:** Use **local disk (EBS)** for `sessions/`. For **S3-free beta** (sessions only on the instance), set:

```bash
HELIX_LOCAL_SESSIONS_ONLY=1
HELIX_SESSIONS_DIR=/var/lib/helix/sessions   # optional; default is ./sessions relative to cwd
```

The backend will **not** create S3 marker objects or session prefixes. Session JSON still lives under `HELIX_SESSIONS_DIR` as today.

If you later add S3 for backups or large files, omit `HELIX_LOCAL_SESSIONS_ONLY` and configure `S3_BUCKET_NAME` / IAM as needed. Optional: **`aws s3 sync`** of the sessions directory on a schedule.

---

## Cost drivers to avoid (if “cheapest” is the goal)

| Component | Rough monthly cost | Beta tactic |
|-----------|-------------------|-------------|
| **NAT Gateway** | ~$32+ | Put a **single EC2** in a **public subnet** with a **public IP / Elastic IP**, or use **VPC endpoints** for S3/ECR only if you must stay private |
| **Application Load Balancer** | ~$16–22+ | For beta, **one EC2 + Nginx + Let’s Encrypt** (or **Caddy**) on **443**, or a tiny **reverse proxy** you already run |
| **ECS Fargate** | ~$35+ for 1 vCPU / 2 GB | **One EC2** running Docker (or systemd + uvicorn) is often **cheaper** at low steady load |
| **Multi-AZ NAT** | 2× NAT | Use **one AZ** for beta |

Your **current** stack (VPC + NAT + ALB + Fargate + CloudFront) is operationally nice but **not** the cheapest monthly floor.

---

## Option A — Cheapest practical beta (recommended)

**Single EC2 + EBS + Docker (or native Python)**

1. **Instance:** **Graviton** **`t4g.small`** or **`t4g.medium`** (often best $/perf) or **`t3.small`** if you need x86 for some wheels.  
2. **Disk:** **gp3** EBS **20–30 GB** for OS + `sessions/` (a few $/month).  
3. **Networking:** Elastic IP, security group **80/443** (and **22** only from your IP).  
4. **HTTPS:** **Caddy** or **Nginx + certbot** terminating TLS on the instance.  
5. **App:** Same container as today (`backend/Dockerfile`) or `uvicorn` with the same env as production.  
6. **S3:** Set **`S3_BUCKET_NAME`** (and IAM instance role with `s3:GetObject`/`PutObject` on your bucket). **Do not** rely on an S3 mount for `sessions/`.  
7. **Frontend:** Build once, serve from **S3 + CloudFront** (very cheap) pointing **`VITE_API_BASE_URL`** at `https://your-beta-domain` (your EC2), **or** serve static files from Nginx on the same box (even cheaper, one host only).

**Rough order of magnitude (us-west-1, one region, low traffic):** EC2 + EBS + data transfer **often ~$15–35/month** before domain, excluding API (OpenAI) costs.

---

## Option B — EC2 + “S3 for sessions” without mounting

If the goal is **durability** rather than POSIX:

- Keep **`sessions/` on EBS** for performance.
- Nightly **`aws s3 sync`** of `sessions/` to **`s3://your-backup-bucket/beta-sessions/`** (Instance IAM role, cron).
- Or enable **EBS snapshots** (pay for snapshot size).

This matches how ops teams do “cheap + recoverable” without FUSE.

---

## Option C — Stay on current ECS stack but shrink it

If you prefer **not** to operate EC2:

- **Remove NAT** only if you refactor networking (hard) or accept **public subnets** for tasks (unusual).  
- **Scale Fargate** to **0.5 vCPU / 1 GB** if load allows.  
- **One task**, no second EC2 (already optional in CDK).

Still typically **more $/month** than one small EC2 for quiet beta traffic.

---

## When an S3 mount *can* help

- **Read-mostly** large files (e.g. reference genomes) under a dedicated prefix — **Mountpoint for Amazon S3** or a sync from S3 on boot.  
- **Not** as a drop-in replacement for **`sessions/*.json`** hot paths unless you change the app to batch writes or use a DB.

---

## Minimal checklist (Option A)

- [ ] VPC + public subnet + **1× EC2** + **Elastic IP**  
- [ ] Security groups: 443 from `0.0.0.0/0` (beta) or restrict to tester IPs  
- [ ] IAM role: S3 access for app bucket; Secrets Manager or SSM for API keys  
- [ ] Install Docker **or** Python venv; run Helix backend on **8001** (or behind reverse proxy)  
- [ ] **`S3_BUCKET_NAME`**, **`HELIX_MAX_UPLOAD_MB`** (e.g. 10 for beta), LLM keys  
- [ ] Frontend: static hosting + CORS / same-origin as needed  
- [ ] **Backups:** cron `s3 sync` or EBS snapshots  

---

## Summary

| Idea | Verdict |
|------|--------|
| Cheapest beta | **One small EC2 + EBS for `sessions/`** + optional **S3 sync** for backup |
| S3 **mounted** as primary session storage | **Not recommended** — poor fit for frequent small JSON writes; use **EBS + boto3 S3** for large data |
| Current Fargate + ALB + NAT | **Easier ops**, **higher** baseline cost |

For a **beta**, start with **Option A**; move back to **managed ECS + ALB** when you need HA, autoscaling, and stricter networking.
