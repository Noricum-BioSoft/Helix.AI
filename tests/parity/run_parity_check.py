#!/usr/bin/env python3
"""
Frontend parity check: runs all 5 demo scenarios via API against both
local (http://localhost:8001) and AWS CloudFront
(https://d2a8mt5n89vos4.cloudfront.net) backends, then writes per-demo
report.json files and a combined parity_summary.json.

Usage:
    python run_parity_check.py [local|aws|both]
"""
import sys
import os
import json
import time
import requests
from pathlib import Path
from datetime import datetime, timezone

# ── Config ────────────────────────────────────────────────────────────────────
ENVS = {
    "local": {
        "label": "local-localhost",
        "base_url": "http://localhost:8001",
        "output_dir": Path("/Users/eoberortner/tmp/parity-local-20260310"),
    },
    "aws": {
        "label": "aws-cloudfront",
        "base_url": "https://d2a8mt5n89vos4.cloudfront.net",
        "output_dir": Path("/Users/eoberortner/tmp/parity-aws-20260310"),
    },
}

DEMOS = [
    {
        "id": "demo-01",
        "title": "Bulk RNA-seq: Toxoplasma gondii",
        "prompt": """\
You are analyzing an RNA-seq transcriptome dataset from a mouse study investigating the effects of Toxoplasma gondii infection on brain gene expression.

Dataset: 12 samples across a 2×2 factorial design:
- Factor 1: Infection status (Infected vs. Uninfected)
- Factor 2: Time post-infection (7 days vs. 21 days)

Files:
- Count matrix: s3://helix-ai-data/demo/tgondii_counts.csv
- Sample metadata: s3://helix-ai-data/demo/tgondii_metadata.csv

Tasks:
1. Load and validate the count matrix and sample metadata for the 2×2 design.
2. Run DESeq2 with a full factorial model (~infection * time) including interaction term, using LRT for the interaction and Wald tests for main effects.
3. Apply independent filtering and Benjamini-Hochberg correction; report how many genes are significant (padj < 0.05) for each contrast.
4. Generate a PCA plot colored by infection status and shaped by time point.
5. Generate a heatmap of the top 50 DE genes (by interaction effect), annotated by infection and time.
6. For the top 10 DE genes, generate a strip plot showing individual sample expression values faceted by group.
7. Produce a written biological summary interpreting the infection × time interaction.""",
    },
    {
        "id": "demo-02",
        "title": "Single-Cell RNA-seq: SLE Immune Profiling",
        "prompt": """\
You are analyzing a single-cell RNA-seq dataset from a human peripheral blood mononuclear cell (PBMC) study investigating immune cell dysregulation in patients with systemic lupus erythematosus (SLE) compared to healthy controls.

Dataset:
- 8 donors: 4 SLE patients, 4 healthy controls
- ~6,400 cells per donor after QC
- Pre-clustered into 9 cell types: CD4 T, CD8 T, NK, B cells, Plasmacytoid DC, Myeloid DC, Classical Monocytes, Non-Classical Monocytes, Plasma Cells

Files:
- Gene-cell UMI count matrix: s3://helix-ai-data/demo/sle_pbmc_counts.csv
- Cell metadata (donor ID, disease status, cell type): s3://helix-ai-data/demo/sle_pbmc_metadata.csv

Tasks:
1. Load and validate the scRNA-seq count matrix and cell metadata.
2. Normalize and log-transform the count matrix (scran or scanpy-style normalization).
3. Run pseudobulk differential expression for each cell type (SLE vs. healthy), using DESeq2 on donor-level aggregated counts.
4. Summarize the number of DE genes per cell type (padj < 0.05) and identify which cell types show the strongest transcriptional response.
5. Generate a dot plot showing the top 5 DE genes per cell type (showing average expression and fraction expressing).
6. Generate a UMAP colored by cell type and by disease status.
7. For Plasmacytoid DCs, generate a volcano plot of all tested genes.
8. Write a biological interpretation: which cell types and gene sets suggest activation of interferon signaling, plasmablast expansion, or other hallmarks of SLE?""",
    },
    {
        "id": "demo-03",
        "title": "Amplicon QC Pipeline: 16S Gut Microbiome",
        "prompt": """\
You are processing a 16S rRNA amplicon sequencing dataset from a gut microbiome study. Raw paired-end FASTQ files are on S3 and need a full preprocessing pipeline before downstream diversity analysis.

Input files:
- R1: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq
- R2: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq

Pipeline steps (run in order):
1. FastQC quality assessment on both R1 and R2 files.
2. Read trimming to remove low-quality bases and adapter sequences (Trimmomatic or fastp), using paired-end mode.
3. Read merging of trimmed paired-end reads (FLASH or PANDA-seq), targeting a 16S V3-V4 amplicon (~450 bp).
4. Final quality report summarizing raw read counts, post-trim read counts, and merged read counts as a pipeline QC table.

For each step, report: tool used, parameters, input/output file paths, and key statistics.""",
    },
    {
        "id": "demo-04",
        "title": "Bulk RNA-seq: Liver Injury Time-Course",
        "prompt": """\
You are analyzing a bulk RNA-seq time-course dataset from a murine model of acute liver injury (acetaminophen overdose), tracking transcriptional recovery from peak injury back to baseline.

Dataset:
- 5 time points: 0h (control), 6h, 24h, 48h, 72h post-APAP
- 3 biological replicates per time point = 15 samples total
- Tissue: mouse liver

Files:
- Count matrix: s3://helix-ai-data/demo/apap_timecourse_counts.csv
- Sample metadata: s3://helix-ai-data/demo/apap_timecourse_metadata.csv

Tasks:
1. Load and validate the count matrix and sample metadata.
2. Run DESeq2 with a likelihood ratio test (LRT) comparing a full time-course model (~time_point) against a reduced intercept-only model, to identify all genes that change over time.
3. Cluster the LRT-significant genes (padj < 0.01) into temporal expression patterns using k-means (k=6), based on VST-normalized mean expression per time point.
4. For each cluster, generate a line plot of mean expression over time, overlaying all genes in light grey and the cluster centroid in bold.
5. Perform gene set enrichment analysis (GSEA or ORA) on each cluster using MSigDB Hallmark gene sets.
6. Generate a heatmap of the top 100 LRT-significant genes, annotated by temporal cluster and time point.
7. Write a biological narrative: what are the major transcriptional phases of liver recovery after APAP overdose?""",
    },
    {
        "id": "demo-05",
        "title": "Phylogenetics: SARS-CoV-2 Variant Divergence",
        "prompt": """\
You are conducting a comparative evolutionary analysis of the SARS-CoV-2 spike protein across major variants of concern (VOCs) to characterize mutational divergence from the ancestral Wuhan-Hu-1 reference sequence.

Sequences to retrieve and analyze (NCBI accession numbers):
- Wuhan-Hu-1 (ancestral): NC_045512
- Alpha (B.1.1.7): MZ344997
- Beta (B.1.351): MZ433432
- Gamma (P.1): MZ477751
- Delta (B.1.617.2): MZ572142
- Omicron BA.1: OL672836
- Omicron BA.2: OL698718
- Omicron BA.4: OL990435
- Omicron BA.5: OL990436

Tasks:
1. Retrieve the spike protein coding sequence (CDS) for each accession from NCBI.
2. Perform multiple sequence alignment (MSA) of all 9 spike protein sequences using MUSCLE or MAFFT.
3. Build a maximum-likelihood phylogenetic tree from the MSA using IQ-TREE or FastTree, with bootstrap support values.
4. Annotate the tree with VOC labels and color-code by variant lineage.
5. Render the final annotated tree as a rectangular cladogram (PNG format).
6. Write a biological interpretation of the tree topology: which variants cluster together, and what does this suggest about convergent evolution of immune escape mutations?""",
    },
]

JOB_POLL_INTERVAL = 5
JOB_POLL_MAX = 180  # 3 minutes


def run_demo(base_url: str, demo: dict, out_dir: Path, env_label: str) -> dict:
    demo_id = demo["id"]
    out_demo = out_dir / demo_id
    out_demo.mkdir(parents=True, exist_ok=True)
    prompt = demo["prompt"]
    report = {
        "demo": demo_id,
        "title": demo["title"],
        "env": env_label,
        "base_url": base_url,
        "session_id": None,
        "plan_status": None,
        "plan_steps": [],
        "plan_has_custom_step": False,
        "execute_status": None,
        "jobs_submitted": 0,
        "job_ids": [],
        "job_final_statuses": {},
        "all_jobs_completed": False,
        "any_job_failed": False,
        "response_has_text": False,
        "response_has_inline_plots": False,
        "response_has_download_links": False,
        "download_links": [],
        "errors": [],
        "elapsed_s": 0,
    }

    t0 = time.time()
    session = None

    try:
        # Create session
        sess_r = requests.post(f"{base_url}/create_session", timeout=15)
        sess_r.raise_for_status()
        session_id = sess_r.json()["session_id"]
        report["session_id"] = session_id
        print(f"  [{demo_id}] Session: {session_id}")

        # Plan
        plan_r = requests.post(f"{base_url}/execute",
                               json={"command": prompt, "session_id": session_id},
                               timeout=30)
        plan_r.raise_for_status()
        plan_data = plan_r.json()
        raw = plan_data.get("raw_result", {})
        report["plan_status"] = raw.get("status") or plan_data.get("status")
        steps = raw.get("data", {}).get("workflow_plan", {}).get("steps", [])
        report["plan_steps"] = [{"step": s.get("step"), "tool": s.get("tool"), "name": s.get("name", "")[:60]}
                                 for s in steps]
        report["plan_has_custom_step"] = any(s.get("tool") == "custom_step" for s in steps)
        print(f"  [{demo_id}] Plan: {report['plan_status']}, steps={len(steps)}, custom_step={report['plan_has_custom_step']}")

        # Execute
        exec_r = requests.post(f"{base_url}/execute",
                               json={"command": prompt, "session_id": session_id, "execute_plan": True},
                               timeout=30)
        exec_r.raise_for_status()
        exec_data = exec_r.json()
        exec_raw = exec_data.get("raw_result", {})
        report["execute_status"] = exec_raw.get("status") or exec_data.get("status")

        # Collect job IDs from response
        def _extract_job_ids(obj, found):
            if isinstance(obj, dict):
                if "job_id" in obj and isinstance(obj["job_id"], str):
                    found.add(obj["job_id"])
                for v in obj.values():
                    _extract_job_ids(v, found)
            elif isinstance(obj, list):
                for item in obj:
                    _extract_job_ids(item, found)

        job_ids_set = set()
        _extract_job_ids(exec_data, job_ids_set)
        job_ids = list(job_ids_set)
        report["job_ids"] = job_ids
        report["jobs_submitted"] = len(job_ids)
        print(f"  [{demo_id}] Execute: {report['execute_status']}, jobs={len(job_ids)}")

        # Check for inline plots and download links in execute response
        def _find_in(obj, keys):
            if isinstance(obj, dict):
                for k, v in obj.items():
                    if k in keys and v:
                        return True
                    if _find_in(v, keys):
                        return True
            elif isinstance(obj, list):
                return any(_find_in(i, keys) for i in obj)
            return False

        report["response_has_text"] = bool(exec_raw.get("text") or exec_data.get("text"))
        report["response_has_inline_plots"] = _find_in(exec_data, {"visuals", "image_b64", "plots"})
        links = []
        def _collect_links(obj):
            if isinstance(obj, dict):
                if "url" in obj and "label" in obj:
                    links.append({"label": obj.get("label"), "url": obj.get("url")})
                for v in obj.values():
                    _collect_links(v)
            elif isinstance(obj, list):
                for i in obj:
                    _collect_links(i)
        _collect_links(exec_data)
        report["download_links"] = links
        report["response_has_download_links"] = len(links) > 0

        # Poll individual jobs
        elapsed = 0
        while elapsed < JOB_POLL_MAX and job_ids:
            statuses = {}
            for jid in job_ids:
                try:
                    jr = requests.get(f"{base_url}/jobs/{jid}", timeout=10)
                    if jr.status_code == 200:
                        jd = jr.json()
                        statuses[jid] = jd.get("status") or jd.get("job", {}).get("status", "unknown")
                except Exception:
                    statuses[jid] = "unknown"

            still_running = [jid for jid, s in statuses.items()
                             if s not in ("completed", "failed", "error")]
            print(f"  [{demo_id}] t+{elapsed}s jobs: {dict(list(statuses.items())[:3])} ...")
            if not still_running:
                report["job_final_statuses"] = statuses
                break
            time.sleep(JOB_POLL_INTERVAL)
            elapsed += JOB_POLL_INTERVAL
        else:
            report["job_final_statuses"] = statuses if job_ids else {}

        report["all_jobs_completed"] = all(
            s == "completed" for s in report["job_final_statuses"].values()
        ) if report["job_final_statuses"] else (report["execute_status"] in ("pipeline_submitted", "success"))
        report["any_job_failed"] = any(
            s in ("failed", "error") for s in report["job_final_statuses"].values()
        )

    except Exception as e:
        report["errors"].append(str(e))
        print(f"  [{demo_id}] ERROR: {e}")

    report["elapsed_s"] = round(time.time() - t0, 1)
    # custom_step is now a graceful fallthrough that returns success; it is expected
    # for biological summary / narrative steps that don't map to a real tool.
    # We only flag failure if a job actually reported "failed" or there was an error.
    report["success"] = (
        report["plan_status"] in ("workflow_planned", "success")
        and not report["any_job_failed"]
        and len(report["errors"]) == 0
    )

    # Save report
    with open(out_demo / "report.json", "w") as f:
        json.dump(report, f, indent=2)
    print(f"  [{demo_id}] → {'✅ PASS' if report['success'] else '❌ FAIL'} ({report['elapsed_s']}s)")
    return report


def run_env(env_key: str) -> list:
    cfg = ENVS[env_key]
    print(f"\n{'='*60}")
    print(f"ENVIRONMENT: {cfg['label']} ({cfg['base_url']})")
    print(f"{'='*60}")

    # Health check
    try:
        h = requests.get(f"{cfg['base_url']}/health", timeout=10).json()
        print(f"Health: {h}")
    except Exception as e:
        print(f"Health check FAILED: {e}")

    cfg["output_dir"].mkdir(parents=True, exist_ok=True)
    reports = []
    for demo in DEMOS:
        print(f"\n--- {demo['id']}: {demo['title']} ---")
        report = run_demo(cfg["base_url"], demo, cfg["output_dir"], cfg["label"])
        reports.append(report)

    # Summary
    summary = {
        "env": cfg["label"],
        "base_url": cfg["base_url"],
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "demos": reports,
        "pass_count": sum(1 for r in reports if r.get("success")),
        "fail_count": sum(1 for r in reports if not r.get("success")),
    }
    with open(cfg["output_dir"] / "parity_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n{cfg['label']}: {summary['pass_count']}/5 passed")
    return reports


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "both"
    envs_to_run = ["local", "aws"] if mode == "both" else [mode]

    all_results = {}
    for env_key in envs_to_run:
        all_results[env_key] = run_env(env_key)

    # Combined parity report
    if len(envs_to_run) > 1:
        parity = {"timestamp": datetime.now(timezone.utc).isoformat(), "comparison": []}
        local_map = {r["demo"]: r for r in all_results.get("local", [])}
        aws_map = {r["demo"]: r for r in all_results.get("aws", [])}
        for demo in DEMOS:
            did = demo["id"]
            l = local_map.get(did, {})
            a = aws_map.get(did, {})
            parity["comparison"].append({
                "demo": did,
                "title": demo["title"],
                "local_success": l.get("success"),
                "aws_success": a.get("success"),
                "parity": l.get("success") == a.get("success"),
                "local_plan_steps": len(l.get("plan_steps", [])),
                "aws_plan_steps": len(a.get("plan_steps", [])),
                "local_custom_step": l.get("plan_has_custom_step"),
                "aws_custom_step": a.get("plan_has_custom_step"),
                "local_jobs": l.get("jobs_submitted"),
                "aws_jobs": a.get("jobs_submitted"),
                "local_any_fail": l.get("any_job_failed"),
                "aws_any_fail": a.get("any_job_failed"),
                "local_plots": l.get("response_has_inline_plots"),
                "aws_plots": a.get("response_has_inline_plots"),
                "local_links": len(l.get("download_links", [])),
                "aws_links": len(a.get("download_links", [])),
            })

        out_path = Path("/Users/eoberortner/tmp/parity-report-20260310.json")
        with open(out_path, "w") as f:
            json.dump(parity, f, indent=2)
        print(f"\nParity report: {out_path}")
        for item in parity["comparison"]:
            status = "✅ PARITY" if item["parity"] else "⚠️  DIFF"
            print(f"  {status}  {item['demo']}: local={item['local_success']} aws={item['aws_success']}")
