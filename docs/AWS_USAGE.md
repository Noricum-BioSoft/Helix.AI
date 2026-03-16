# How Helix.AI Uses AWS

This document describes how the **cloud**, specifically **AWS**, is used by Helix.AI: for hosting, data, and compute. It complements the [AWS Deployment Guide](deployment/AWS_DEPLOYMENT_GUIDE.md) (how to deploy) by focusing on **what** runs where and **why**.

---

## Overview

Helix.AI uses AWS for:

| Purpose | AWS services | Role |
|--------|---------------|------|
| **Hosting the application** | ECS (Fargate), ALB, ECR, S3, CloudFront | Run the backend API and serve the frontend |
| **Data (inputs and outputs)** | S3 | User data (e.g. `s3://bucket/counts.csv`), job results, session artifacts |
| **Compute (optional)** | EC2, EMR | Run tools on larger or distributed compute when data size or tool needs require it |

The system can run **fully local** (backend and frontend on your machine, no AWS). When deployed or when users point to S3 data, AWS is used as above.

---

## 1. Hosting and delivery

### Backend (API)

- **ECS Fargate**: The FastAPI backend runs in a container on ECS. No servers to manage; Fargate handles scaling and availability.
- **ALB**: The Application Load Balancer receives HTTPS traffic and forwards it to the ECS service. Health checks use `/health`.
- **ECR**: Docker images for the backend are stored in Elastic Container Registry and pulled by ECS.

**Relevant env (deployment):** `AWS_REGION`, ECS task definition (CPU/memory), ALB listener.

### Frontend (UI)

- **S3**: The built React app (static files: HTML, JS, CSS) is uploaded to an S3 bucket.
- **CloudFront**: A CloudFront distribution serves the frontend with HTTPS and caching. Users hit the CloudFront URL, which serves from S3 and (for API calls) can proxy to the ALB or use a separate backend URL.

**Relevant env:** `VITE_API_BASE_URL` (backend URL at build time).

### Summary diagram

```
User → CloudFront (frontend) → Browser
         S3 (static assets)

User → ALB → ECS Fargate (backend FastAPI)
              └── ExecutionBroker → Local | EC2 | EMR (see below)
```

See [architecture/ARCHITECTURE_EXPLANATION.md](architecture/ARCHITECTURE_EXPLANATION.md) and [deployment/AWS_DEPLOYMENT_GUIDE.md](deployment/AWS_DEPLOYMENT_GUIDE.md) for deployment details.

---

## 2. Data: S3

### How S3 is used

- **User inputs**: Many tools accept **S3 URIs** as input (e.g. `count_matrix: s3://bucket/counts.csv`, `input_r1: s3://bucket/sample_R1.fastq.gz`). The backend reads from S3 when executing the tool (or streams/copies as needed).
- **Size estimation**: The ExecutionBroker uses **S3 `head_object`** (via boto3) to get object size for routing: e.g. total input size &gt; 100MB may route to EMR instead of local. If credentials or region are missing, size is unknown and fallback rules apply.
- **Job outputs (async)**: For async jobs (e.g. FastQC on EMR), results (reports, logs) are written to S3; the backend stores the result path and returns a `job_id`. The frontend or client can poll for completion and then fetch result URLs.
- **Session/artifacts**: Run artifacts (e.g. `analysis.py`, plots) can be stored under session/run directories. In a deployed setup these may be on ECS-backed storage or, if configured, in S3.

### Credentials and region

- The backend uses **default AWS credentials** (instance role when on ECS, or env/CLI when local).
- **Region**: `AWS_REGION` / `AWS_DEFAULT_REGION` (e.g. `us-east-1`) is used for S3 and other clients. Cross-region buckets may require the bucket’s region for `head_object`; the code attempts to detect and retry with the bucket region when needed.

---

## 3. Compute: Local vs EC2 vs EMR

Helix does not run every tool in the cloud. It chooses **where** to run based on data and tool requirements.

### Who decides

The **Infrastructure Decision Agent** (see [agents/infrastructure-decision-agent.md](../agents/infrastructure-decision-agent.md)) recommends an execution environment. The **ExecutionBroker** uses that plus a **bytes threshold** (e.g. 100MB) to decide **sync** (local or EC2) vs **async** (e.g. EMR).

### The 3-factor model (summary)

1. **Data locality and size**  
   - Inputs on **S3** and **total size &gt; ~100MB** → prefer **EMR** (process where data lives, avoid large egress).  
   - S3 &lt; 100MB or local data → **Local** or **EC2** is fine.  
   - Very large local data → upload to S3 and use EMR (or similar) in practice.

2. **Tool availability**  
   - **EC2**: Pre-installed tools (e.g. BBTools, samtools, FastQC) → good for classic bioinformatics.  
   - **EMR**: Spark/Hadoop; tools via bootstrap or containers.  
   - **Local**: Whatever is installed on the machine running the backend.

3. **Compute requirements**  
   - CPU, memory, runtime. EMR and EC2 support long and heavy jobs; Lambda is not used for heavy bioinformatics (time/memory limits).

### Sync vs async

- **Sync**: Request runs to completion; response includes the result. Used for **Local** and **EC2** (and when the broker forces sync for small/unknown size).
- **Async**: Backend submits a job (e.g. to EMR), returns a **job_id** immediately; client polls for status and fetches results (e.g. from S3) when done. Used for **EMR** (and optionally Batch) when the infrastructure agent recommends it.

### When EC2 is used

- **EC2** is only used if **`HELIX_USE_EC2=true`** (and optionally EC2 instance/keys configured). If not set, EC2 is treated as unavailable and the agent chooses Local or EMR.
- Typical use: tools that benefit from a fixed, tool-rich VM (e.g. FastQC, aligners) on larger-than-local inputs that still fit in one machine.

### When EMR is used

- **EMR** is used for **async** jobs when the infrastructure agent recommends it (e.g. large S3 inputs). The broker submits a job to the cluster; results are written to S3 and linked via `job_id`.
- See [scripts/emr/README.md](../scripts/emr/README.md) and deployment docs for EMR cluster and IAM setup.

---

## 4. Executing a data analysis pipeline on AWS

A **data analysis pipeline** here means either (1) a **single tool run** (e.g. one FastQC or one bulk RNA-seq run), or (2) a **multi-step workflow** (e.g. QC → trim → merge → downstream analysis) planned by the Workflow Planner and executed by the ExecutionBroker.

### Single-tool execution

- The flow in **§5 (End-to-end flow)** applies: one tool, one infra decision (Local / EC2 / EMR), then sync execution or async job submission.
- Inputs are often **S3 URIs**; the broker may estimate total input size via S3 `head_object` to decide sync vs async.
- Outputs (plots, tables, reports) are produced by the tool; for async EMR jobs, outputs are written to S3 and linked via `job_id`.

### Multi-step pipeline execution

- When the user request is turned into a **plan** (Plan IR with multiple steps), the broker can execute it in two ways:
  - **Sync (sequential)**: Steps run one after another on the **same** execution environment (e.g. Local or EC2). The broker resolves references between steps (e.g. step B uses `$ref: steps.A.result.output_path`). Each step's output is passed to the next; no per-step cloud choice today for this path.
  - **Async (plan as a job)**: For large inputs, the **whole plan** can be submitted as one async job (e.g. to EMR). The job payload describes the steps; the cluster runs them (e.g. via a runner that reads from S3 and writes results to S3). The API returns a `job_id`; the client polls for completion and then fetches result URLs (e.g. S3).
- **Where the pipeline runs**: Today, pipeline execution on AWS means either (a) **Local/EC2** (backend or EC2 instance runs the steps sequentially, pulling from S3 as needed), or (b) **EMR** (the plan is submitted as an async job; the EMR step reads inputs from S3 and writes outputs to S3). So AWS is used for **data (S3)** and, when chosen, for **compute (EC2 or EMR)** that runs the pipeline steps.
- **Session and iteration**: Each run (single-tool or pipeline) is recorded in the session with a `run_id`. Iterative "change parameter and re-run" (e.g. patch_and_rerun) creates a new run, often still on the same infra (Local/EC2); artifacts and scripts are stored per run (and in a deployed setup may be on ECS-backed storage or S3 if configured).

### Summary

- **Single-tool**: One infra decision → one sync run (Local/EC2) or one async job (EMR); inputs/outputs on S3 as described above.
- **Pipeline**: Either sequential execution (Local/EC2) with step refs, or one async job (EMR) for the full plan; in both cases AWS is used for data (S3) and, when selected, for compute (EC2/EMR) that executes the pipeline.

---

## 5. End-to-end flow (with AWS)

1. User sends a command (e.g. “Run FastQC on s3://my-bucket/R1.fastq.gz and R2.fastq.gz”).
2. Backend (on ECS or local) receives the request; may call **S3 `head_object`** to get file sizes.
3. **Infrastructure Decision Agent** recommends an environment (Local / EC2 / EMR) using the 3-factor model.
4. **ExecutionBroker**:
   - If **Local/EC2** → runs the tool (or pipeline steps) **synchronously** (download from S3 if needed, then execute).
   - If **EMR** → submits an **async** job (single tool or full plan), returns `job_id`; job reads from S3, writes results to S3; client polls and retrieves result URLs.
5. Response is returned to the frontend (or client); for async, the UI can show job status and a link to the result (e.g. S3 or a backend endpoint that serves it).

---

## 6. Other cloud providers (e.g. GCP) — integration path

Today Helix uses **AWS only** for cloud hosting, storage, and compute. The design does not yet describe GCP, Azure, or other providers, but the same **concepts** (where data lives, where compute runs, sync vs async) can be extended.

### Current scope

- **Hosting**: ECS, S3, CloudFront, ALB, ECR (all AWS).
- **Data**: S3 URIs (`s3://bucket/key`) for inputs and outputs; size estimation and job results assume S3.
- **Compute**: Local (no cloud), EC2, or EMR (all AWS when cloud is used).

### How other clouds could be integrated

- **Abstraction points** in the codebase:
  - **Storage URIs**: Today only `s3://` is handled for size estimation and access. A **storage abstraction** (e.g. a small layer that resolves `gs://`, `s3://`, `abfs://`) would allow **GCS** (GCP), **Azure Blob**, etc. for inputs and outputs. The ExecutionBroker would call this layer instead of S3-only logic.
  - **Compute backend**: The broker today maps "run here" to Local, EC2, or EMR. A **compute backend interface** (e.g. "submit job with payload, return job_id; poll status; get result URLs") would allow plugging in **GCP** (e.g. Cloud Run, Cloud Batch, or Dataflow for distributed steps) or **Azure** (Batch, AKS, etc.) alongside or instead of EMR/EC2.
- **GCP example**: To use **Google Cloud** for pipeline execution you would:
  - **Data**: Accept **GCS** URIs (`gs://bucket/key`) for inputs and outputs; add a GCS client and size estimation (e.g. `get_object` metadata or equivalent); store job results in GCS when the compute runs on GCP.
  - **Compute**: Implement a GCP compute backend (e.g. submit a job to **Cloud Run**, **Cloud Batch**, or **Dataflow** with a payload that describes the tool or plan); return a `job_id` and poll GCP for status; resolve result URLs to GCS.
  - **Config**: Add env or config (e.g. `HELIX_CLOUD_PROVIDER=gcp`, `GCP_PROJECT`, `GCS_BUCKET`) and wire the broker to the GCP backend when chosen.
- **Multi-cloud**: In principle, data could stay on one provider (e.g. S3) while compute runs on another (e.g. GCP Batch) if you handle cross-cloud transfer and IAM. In practice, "data and compute on same provider" (S3 + EMR, or GCS + Cloud Batch) is simpler and cheaper; the docs assume that pattern for AWS today.

### Status and roadmap

- **Today**: Only AWS is implemented and documented. No GCP or Azure code paths.
- **Roadmap**: Adding a storage abstraction and a compute-backend interface would allow integrating GCP (or others) without rewriting the whole broker; the Infrastructure Decision Agent would need to be extended to recommend "GCP" when appropriate (e.g. when inputs are on GCS or when org policy prefers GCP). This is not committed for a specific release; it is the intended **integration path** when other clouds are required.

---

## 7. IAM and permissions (summary)

For the backend to use AWS as above:

| Service | Typical permission need |
|--------|---------------------------|
| **S3** | `s3:GetObject` (read user inputs, job outputs); `s3:PutObject` / `s3:ListBucket` if writing results or listing prefixes. |
| **EMR** | Create and run steps (e.g. `elasticmapreduce:*` or minimal step submission roles). |
| **EC2** | If using EC2: describe instances, run commands or use SSM (depending on implementation). |
| **ECR** | Pull image (for ECS task role). |
| **CloudWatch** | Write logs (for ECS task role). |

Exact policies depend on your deployment; the deployment guide and CDK/infrastructure code define the roles.

---

## 8. Configuration (env) relevant to AWS

| Variable | Purpose |
|----------|---------|
| `AWS_REGION` / `AWS_DEFAULT_REGION` | Region for S3, EMR, and other clients. |
| `HELIX_USE_EC2` | If `true`, EC2 is available as an execution target. |
| `HELIX_ASYNC_BYTES_THRESHOLD` | Optional; override default (e.g. 100MB) for sync vs async. |
| (Deployment) ECS task role, EMR job role, S3 bucket names | Set in infrastructure (CDK/CloudFormation) or deployment config. |

---

## 9. References

- **Deployment**: [deployment/AWS_DEPLOYMENT_GUIDE.md](deployment/AWS_DEPLOYMENT_GUIDE.md) — how to deploy backend and frontend to AWS.
- **Infrastructure as code**: [infrastructure/README.md](../infrastructure/README.md) — CDK stack (ECS, S3, CloudFront, ALB, ECR).
- **Infrastructure decisions**: [agents/infrastructure-decision-agent.md](../agents/infrastructure-decision-agent.md) — 3-factor model and thresholds.
- **Execution flow**: [architecture/BACKEND_DATAFLOW.md](architecture/BACKEND_DATAFLOW.md) — request path and routing.
- **Broker behavior**: `backend/execution_broker.py` — sync/async, S3 size estimation, infra decision usage.
