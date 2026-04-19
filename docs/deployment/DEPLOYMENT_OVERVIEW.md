# Helix.AI AWS deployment — overview

This document is the **single place to revisit** how deployment works: architecture, two deployment layers, and how to run each step.

---

## 1. Architecture (what runs where)

```
User → CloudFront (HTTPS) → S3 (static React app)
                           → ALB (HTTP) → ECS Fargate (FastAPI backend on port 8001)
```

- **Frontend:** S3 bucket + CloudFront distribution (same origin can route `/health`, `/agent`, `/execute`, etc. to the ALB — see CDK stack).
- **Backend:** Docker image in **ECR**; **ECS Fargate** service pulls `latest` (or your tag) and runs behind the **ALB**.
- **Networking:** VPC with public/private subnets, **one NAT gateway** (cost tradeoff).
- **Secrets:** API keys from **Secrets Manager** (referenced in the ECS task definition).

**Single backend (cost):** By default the CDK stack does **not** create an optional EC2 instance. Only Fargate runs the API. To add EC2 (e.g. build host): `cdk deploy -c createEc2Instance=true`. See [AWS_COST_INVESTIGATION.md](AWS_COST_INVESTIGATION.md).

**Hosted limits:** The Fargate task sets `HELIX_MAX_UPLOAD_MB=10` (configurable). Upload guards in the backend enforce small file uploads; EMR/large-data code remains in the repo but is not used for large uploads on this deployment.

---

## 2. Two layers (do them in order)

| Layer | What | Tool | When |
|-------|------|------|------|
| **A. Infrastructure** | VPC, ECS, ALB, ECR ref, S3 ref, CloudFront, IAM, task definition (env vars) | **AWS CDK** | When you change infra or ECS task env (e.g. new env var in `helix_stack.py`) |
| **B. Application** | Build backend image → push ECR; build frontend → sync S3; invalidate CloudFront; roll ECS to new image | **`scripts/aws/deploy.sh`** or **GitHub Actions** | When you change app code |

Infrastructure does **not** automatically rebuild your Docker image. After CDK changes that only affect task definition env, ECS may roll tasks; you still need a **new image** in ECR when application code changes.

---

## 3. Layer A — CDK (infrastructure)

```bash
cd infrastructure
pip install -r requirements.txt
# If `cdk deploy` fails with "schema version mismatch", use a matching CLI:
npx aws-cdk@latest deploy --require-approval never
```

- **Region/account:** Set `CDK_DEFAULT_ACCOUNT` and `CDK_DEFAULT_REGION` (e.g. `us-west-1`) to match your AWS profile.
- **Outputs:** After deploy, note **ECR repo name**, **ECS cluster/service**, **task definition family**, **S3 bucket**, **CloudFront ID**, **ALB DNS** — use these in `deploy.config` and GitHub secrets.

---

## 4. Layer B — Application deploy

### Option 1: Local script (`deploy.sh`)

1. Copy and edit config:
   ```bash
   cp scripts/aws/deploy.config.example scripts/aws/deploy.config
   ```
2. Fill in values from CDK stack outputs (or AWS Console): `AWS_REGION`, `AWS_ACCOUNT_ID`, `ECR_REPOSITORY`, `S3_BUCKET`, `ECS_*`, `CLOUDFRONT_DISTRIBUTION_ID`, `VITE_API_BASE_URL` (often your CloudFront URL or ALB URL, depending on how the frontend calls the API).

3. Run from repo root:
   ```bash
   ./scripts/aws/deploy.sh
   ```

**Flags in `deploy.config`:** `SKIP_BACKEND=true`, `SKIP_FRONTEND=true`, etc. to run only part of the pipeline.

**Apple Silicon / proxy issues:** CI builds with `linux/amd64` for Fargate. Local `docker build` may produce `arm64` unless you use buildx with `--platform linux/amd64`. If `docker push` to ECR fails with proxy/closed connection, disable Docker proxy for `*.amazonaws.com` or unset `HTTP_PROXY`/`HTTPS_PROXY` for the push. Prefer **GitHub Actions** for a reliable amd64 push.

### Option 2: GitHub Actions (`.github/workflows/deploy.yml`)

Triggers on push to **`main`** or **`version_2.0`**, or **workflow_dispatch**.

**Required repository secrets** (must be set or backend deploy is skipped):

| Secret | Purpose |
|--------|---------|
| `AWS_ACCESS_KEY_ID` | IAM user/role key for deploy |
| `AWS_SECRET_ACCESS_KEY` | IAM secret |
| `AWS_REGION` | e.g. `us-west-1` |
| `ECR_REPOSITORY` | e.g. `helix-ai-backend` |

**For full deploy (frontend too):**

| Secret | Purpose |
|--------|---------|
| `S3_BUCKET` | Frontend static bucket |
| `VITE_API_BASE_URL` | API base URL baked into the frontend build |

**Optional (recommended for automatic ECS roll):**

| Secret | Purpose |
|--------|---------|
| `ECS_CLUSTER_NAME` | Cluster name |
| `ECS_SERVICE_NAME` | Service name |
| `ECS_TASK_DEFINITION_FAMILY` | Task definition **family** name |
| `CLOUDFRONT_DISTRIBUTION_ID` | Cache invalidation after S3 sync |

`AWS_ACCOUNT_ID` is in workflow `env` but not required for the steps if credentials are configured.

---

## 5. Quick verification

- **Health:** `curl -s https://<cloudfront-or-alb>/health` (or `http` if no TLS on ALB).
- **ECS:** AWS Console → ECS → cluster → service → **Events** and **Tasks** (running, new deployment).
- **Logs:** CloudWatch log group `/ecs/helix-ai` (or as defined in CDK).

---

## 6. Related docs

| Doc | Topic |
|-----|--------|
| [AWS_DEPLOYMENT_GUIDE.md](AWS_DEPLOYMENT_GUIDE.md) | Long-form steps, troubleshooting, production hardening |
| [EC2_GIT_BETA.md](EC2_GIT_BETA.md) | EC2 beta: GitHub clone + `git pull`, systemd, local sessions |
| [AWS_COST_INVESTIGATION.md](AWS_COST_INVESTIGATION.md) | Cost drivers, optional EC2, Cost Explorer tips |
| [../../infrastructure/README.md](../../infrastructure/README.md) | CDK commands and stack customization |

---

## 7. Checklist — “I changed code, what do I run?”

- **Only Python/backend:** Push image + ECS roll — `deploy.sh` (backend only) or GitHub Actions backend job; or `SKIP_FRONTEND=true ./scripts/aws/deploy.sh`.
- **Only frontend:** `SKIP_BACKEND=true ./scripts/aws/deploy.sh` or rely on CI frontend job.
- **Infra / env vars on the task:** `cdk deploy` then, if needed, redeploy app so image matches (if code changed too).
- **Secrets (API keys):** Update in Secrets Manager; ECS task already references them — usually **restart/redeploy** service to pick up new secret values.
