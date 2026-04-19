# AWS Cost Investigation – Helix.AI

This document summarizes where Helix.AI’s AWS costs come from and how to reduce them.

## 1. Where costs come from

The CDK stack (`infrastructure/helix_infrastructure/helix_stack.py`) and CI/CD (`.github/workflows/deploy.yml`) deploy the following. Approximate monthly costs assume **us-west-1** and 24/7 operation unless noted.

| Resource | Purpose | Approx. monthly cost (baseline) |
|----------|--------|---------------------------------|
| **ECS Fargate** | Backend API (1 task: 1 vCPU, 2 GB) | ~\$35–45 |
| **EC2 t3.medium** | Alternative backend + build host (always created by stack) | ~\$30 |
| **NAT Gateway** (1) | Egress for private subnets (Fargate/EC2) | ~\$33 + \$/GB data |
| **Application Load Balancer** | Routes traffic to backend | ~\$16–22 + LCU usage |
| **CloudFront** | Frontend CDN (PRICE_CLASS_100) | Low (pay per request/GB) |
| **S3** | Frontend bucket + build artifacts | Low |
| **ECR** | Backend Docker images | Low |
| **Secrets Manager** | API keys (2 secrets) | ~\$1 |
| **CloudWatch Logs** | ECS logs, 7-day retention | Low–moderate |

**Single backend (cost-reduced):** The stack is configured for **one** backend only:

- **Fargate** is the only backend by default (used by CI/CD and production).
- **EC2** is **not** created unless you explicitly pass `-c createEc2Instance=true` when deploying. This avoids paying for two backends (~\$30/month savings).

**Optional / on-demand:**

- **EMR**: Not in the CDK stack. Used when workflows run on EMR (see `backend/config/cost_heuristics.yaml` and EMR scripts). Costs apply only when clusters are running.

---

## 2. Recommendations (in order of impact)

### 2.1 Remove or make EC2 optional (high impact)

- **If you only use Fargate:**  
  Deploy with EC2 disabled so the stack does not create the EC2 instance or its target group.  
  **Saves:** ~\$30/month.

- **If you need EC2 sometimes:**  
  Use the new `createEc2Instance` (or equivalent) context/flag so EC2 is only created when needed. When not needed, deploy with the flag off.

**How:** The EC2 instance is **optional**. By default it is **not** created (saves ~\$30/month).

- **Fargate-only (default):**  
  `cdk deploy`  
  No EC2 instance or EC2 target group is created.

- **Include EC2 (e.g. for build-on-EC2 or EC2-backed backend):**  
  `cdk deploy -c createEc2Instance=true`  
  or  
  `HELIX_CREATE_EC2=1 cdk deploy`

**If you already have the stack deployed with EC2:** The next `cdk deploy` without the flag will remove the EC2 instance and related resources (target group, security group, IAM role). To keep EC2, always deploy with `-c createEc2Instance=true`. To remove EC2 and save cost, deploy once without the flag; CloudFormation will delete those resources.

---

### 2.2 Right-size Fargate (medium impact)

- Current: **1 vCPU, 2 GB** per task.
- If the backend is lightly used, try **0.5 vCPU, 1 GB** (min for Fargate).  
  **Saves:** Roughly 40–50% of Fargate compute (~\$15–20/month).

Change in `helix_stack.py`:

```python
# In FargateTaskDefinition
memory_limit_mib=1024,   # 1 GB
cpu=512,                 # 0.5 vCPU
```

Validate with load tests; increase again if you see timeouts or OOMs.

---

### 2.3 NAT Gateway (medium impact)

- **Current:** 1 NAT Gateway (~\$33/month + \$/GB).
- **Options:**  
  - Use a **NAT instance** (e.g. t4g.nano) for much lower fixed cost (few dollars/month), with more ops and lower throughput.  
  - For dev/staging, consider **VPC endpoints** for AWS APIs (e.g. ECR, S3) so Fargate doesn’t need NAT for those; then you can remove NAT if all egress goes via endpoints or a single small NAT instance.

---

### 2.4 Fargate Spot (medium impact, for non-critical workloads)

- Use **Fargate Spot** for the backend if short interruptions are acceptable.  
- **Saves:** Up to ~70% on Fargate compute.

In the stack, set `capacity_provider_strategies` on the Fargate service to use `FARGATE_SPOT` (e.g. 100% Spot or a mix with `FARGATE`).

---

### 2.5 CloudFront and S3 (already lean)

- CloudFront uses **PRICE_CLASS_100** (US/Europe only).  
- Frontend is static on S3.  
No change needed unless you add more regions or logging.

---

### 2.6 Logs and retention (small impact)

- CloudWatch Logs: 7-day retention is already short.  
- You can reduce to 3 days or use a smaller retention for non-production to trim a few dollars.

---

### 2.7 EMR (when used)

- EMR is not in the CDK stack; costs are per cluster and runtime.  
- Use the existing **cost heuristics** in `backend/config/cost_heuristics.yaml` when choosing infrastructure (Local / EC2 / EMR / Batch).  
- For small jobs, prefer **EC2 or Batch** instead of EMR to avoid cluster startup and minimum cluster cost.

---

## 3. Quick wins checklist

- [ ] **Deploy without EC2** (use `createEc2Instance: false` or equivalent) if you only use Fargate.
- [ ] **Right-size Fargate** to 0.5 vCPU / 1 GB if load allows.
- [ ] **Consider Fargate Spot** for the backend service.
- [ ] **Review NAT**: NAT instance or VPC endpoints to cut NAT Gateway cost.
- [ ] **Confirm in AWS Cost Explorer** which services actually drive the bill (ECS, EC2, NAT, ALB, data transfer).

---

## 4. Verifying your actual bill

1. **AWS Cost Explorer** → Group by **Service** to see ECS, EC2, NAT Gateway, CloudFront, S3, etc.
2. **Tags:** If resources are tagged (e.g. `Project=helix-ai`), filter by tag in Cost Explorer.
3. **ECS:** Billing → ECS → Cluster → Service to see Fargate cost per service.
4. **EC2:** Check whether the Helix.AI EC2 instance is running 24/7; if it’s unused, stop it or remove it from the stack.

After applying the optional EC2 and Fargate changes, redeploy and recheck Cost Explorer in the next billing period to confirm savings.

---

## 5. Why your Cost Explorer CSV might show almost no cost

If you export a report from Cost Explorer and the **Total costs** column is near zero (or negative) while your actual AWS bill is high, it’s usually due to one of the following.

### What your export is actually showing

A typical export has:

- **Rows:** Dates (e.g. daily or monthly).
- **Columns:** Cost per service (S3, ECS, EC2, ELB, Data Transfer, etc.).
- **“Service total” row:** Sum of each column over the date range.
- **“Total costs” column:** Often the **net** after credits/refunds.

In your CSV:

- **Service totals** (row 2) show real usage: e.g. S3 ~\$5.80, EC2-Other ~\$0.66, ECS ~\$0.64, EC2-Instances ~\$0.01, ELB very small.
- **Data Transfer** is **negative** (e.g. -\$7.11) — that’s a **credit** (refund or allocation).
- **Total costs** is essentially **zero** because the Data Transfer credit (and possibly other credits) offsets the positive service costs in that view.

So the CSV **does** show where the reported usage is (S3, ECS, EC2-Other, etc.), but the **total** in the export is “after credits” and can look tiny even when your invoice is high.

### Why your bill can still be high

1. **Credits hide cost in the “Total” column**  
   Free tier, promotional credits, or support credits reduce **net** cost in Cost Explorer. The **Billing dashboard** and **invoice (PDF)** show **charges before credits** and **credits** separately. So Cost Explorer’s “Total” can be ~\$0 while the invoice still shows e.g. \$80 charges and \$80 in credits.

2. **Wrong cost metric**  
   Cost Explorer can show **Unblended** (what you’re actually charged), **Blended** (amortized across the org), or **Amortized** (Reserved Instance amortization). If you picked Blended/Amortized or a view that includes heavy credits, the export total won’t match “what I owe.”

3. **Date range or granularity**  
   The export might be for a period when usage was low, or monthly/daily granularity can spread or hide a spike. Compare the export date range to the invoice period.

4. **Filter or grouping**  
   A filter (e.g. by tag, by linked account, by service) can exclude the resources that drive the bill. Untagged resources may appear under “No tag” or in an account you didn’t export.

5. **NAT Gateway and other “EC2-Other”**  
   NAT Gateway and some other networking costs often appear under **EC2-Other** (or **VPC**). Your CSV shows **EC2-Other** (~\$0.66 in the total row) and **VPC** (\$0). If NAT is running 24/7, its share might be in a different time window or in another report (e.g. by **Usage type**).

### How to get numbers that match your bill

1. **Billing dashboard (source of truth)**  
   - **Billing → Bills**: Open the **PDF invoice** for the month you care about.  
   - Check **Charges by service** and **Credits**; the “Total” on the invoice is what you pay.

2. **Cost Explorer report that matches the invoice**  
   - **Cost Explorer → Create report** (or open an existing one).  
   - **Time range:** Same as the invoice (e.g. last month).  
   - **Granularity:** Monthly (or Daily if you want detail).  
   - **Metric:** **Unblended cost** (this is the charge that matches the invoice).  
   - **Group by:** **Service**.  
   - **Do not** add a filter (tag, account, etc.) unless you intentionally want a subset.  
   - Export this report to CSV. The **Total** column and the sum of service columns should be close to the invoice total (before credits).

3. **Find the expensive services**  
   - In the same report, sort by cost. You’ll see which service (e.g. ECS, EC2, S3, Data Transfer, NAT under EC2-Other) drives the bill.  
   - For **NAT Gateway**, group by **Usage type** and look for “NatGateway” or “DataProcessing”.

4. **If you use multiple accounts**  
   - In **Cost Explorer**, ensure you’re in the **payer (management) account** or that you’ve selected the right **linked account**.  
   - Run the same Unblended, by Service, no-filter report for the account that runs Helix.AI.

Once you have an **Unblended, by Service, no-filter** export for the same period as your invoice, the CSV total and the invoice total should align (aside from rounding). The difference from your current export is mainly: use **Unblended** and ignore the “Total” column when heavy **credits** (like Data Transfer) are present; rely on the **invoice** and the **service-level totals** in the CSV to see where the high costs are.
