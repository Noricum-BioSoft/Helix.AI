#!/usr/bin/env python3
"""
Strict end-to-end rerun of all 5 Helix demos against AWS backend.

What this script does:
- Creates fresh session per demo
- Executes plan/input flow and execution flow
- Polls async jobs until completion
- Retrieves visualization outputs for pipeline demos
- Downloads available artifacts and bundles
- Writes normalized output structure per demo:
  01_plan_input, 02_execute, 03_jobs, 04_artifacts, 05_bundles, 06_verification
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import pathlib
import re
import sys
import time
import urllib.parse
import zipfile
from typing import Any

import requests


DEFAULT_AWS_URL = "http://HelixA-ALBAE-D7BksiQIynZb-1051248867.us-west-1.elb.amazonaws.com"


DEMOS = [
    {
        "id": "demo-01",
        "kind": "needs_inputs",
        "prompt": (
            "You are analyzing an RNA-seq transcriptome dataset from a mouse study "
            "investigating the effects of Toxoplasma gondii infection on brain gene expression."
        ),
        "followup": (
            "count_matrix: s3://noricum-ngs-data/demo/rnaseq/tgondii_counts.csv\n"
            "sample_metadata: s3://noricum-ngs-data/demo/rnaseq/tgondii_metadata.csv\n"
            "design_formula: ~infection_status + time_point + infection_status:time_point"
        ),
    },
    {
        "id": "demo-02",
        "kind": "needs_inputs",
        "prompt": (
            "You are analyzing a single-cell RNA-seq dataset from a human peripheral blood "
            "mononuclear cell (PBMC) study investigating immune cell dysregulation in patients "
            "with systemic lupus erythematosus (SLE) compared to healthy controls."
        ),
        "followup": (
            "data_file: s3://noricum-ngs-data/demo/scrna/sle_pbmc_filtered_feature_bc_matrix.h5\n"
            "data_format: 10x\n"
            "resolution: 0.5\n"
            "steps: all"
        ),
    },
    {
        "id": "demo-03",
        "kind": "pipeline",
        "prompt": (
            "You are processing a 16S rRNA amplicon sequencing dataset from a gut microbiome "
            "study. Raw paired-end FASTQ files are on S3 and need a full preprocessing pipeline "
            "before downstream diversity analysis.\n\n"
            "Dataset\n"
            "  Illumina MiSeq 2x250 bp paired-end reads; V3-V4 hypervariable region.\n"
            "  Forward reads: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq\n"
            "  Reverse reads: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq\n"
            "  Output prefix:  s3://noricum-ngs-data/test-output/amplicon-demo/\n\n"
            "Pipeline Steps\n"
            "  1. Run FastQC quality assessment on both raw R1 and R2 files.\n"
            "  2. Trim adapter sequences (CTGTCTCTTATACACATCT) and low-quality bases (Phred < 20).\n"
            "  3. Merge overlapping paired-end reads with minimum overlap of 20 bp.\n"
            "  4. Generate a quality report summarizing read counts before and after each step."
        ),
    },
    {
        "id": "demo-04",
        "kind": "needs_inputs",
        "prompt": (
            "You are analyzing a bulk RNA-seq time-course dataset from a murine model of acute "
            "liver injury (acetaminophen overdose), tracking transcriptional recovery from peak "
            "injury back to baseline."
        ),
        "followup": (
            "count_matrix: s3://noricum-ngs-data/demo/rnaseq/apap_timecourse_counts.csv\n"
            "sample_metadata: s3://noricum-ngs-data/demo/rnaseq/apap_timecourse_metadata.csv\n"
            "design_formula: ~time_point"
        ),
    },
    {
        "id": "demo-05",
        "kind": "pipeline",
        "prompt": (
            "You are conducting a comparative evolutionary analysis of the SARS-CoV-2 spike protein "
            "across major variants of concern (VOCs) to characterize mutational divergence from "
            "the ancestral Wuhan-Hu-1 reference sequence.\n\n"
            "Dataset\n"
            "  Use accessions: MN908947.3,OQ898928.1,OR353131.1,MW642250.1,OR323381.1,PP847536.1,PP848071.1,PP405604.1\n\n"
            "Objectives\n"
            "  1. Retrieve spike sequences from NCBI.\n"
            "  2. Perform multiple sequence alignment.\n"
            "  3. Build phylogenetic tree.\n"
            "  4. Annotate mutation sites.\n"
            "  5. Compute pairwise identity matrix."
        ),
    },
]


def _sanitize(name: str) -> str:
    val = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")
    return val or "artifact"


def _save_json(path: pathlib.Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str))


def _extract_links(payload: dict[str, Any]) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    if isinstance(payload.get("data"), dict):
        data = payload["data"]
        for key in ("downloadable_artifacts", "links"):
            arr = data.get(key)
            if isinstance(arr, list):
                for item in arr:
                    if isinstance(item, dict) and isinstance(item.get("url"), str):
                        candidates.append(item)
    if isinstance(payload.get("raw_result"), dict):
        raw = payload["raw_result"]
        arr = raw.get("links")
        if isinstance(arr, list):
            for item in arr:
                if isinstance(item, dict) and isinstance(item.get("url"), str):
                    candidates.append(item)
    seen: set[str] = set()
    out: list[dict[str, Any]] = []
    for c in candidates:
        url = c.get("url")
        if url in seen:
            continue
        seen.add(url)
        out.append(c)
    return out


def _to_absolute_url(base_url: str, url: str) -> str:
    if url.startswith("http://") or url.startswith("https://"):
        return url
    if url.startswith("/"):
        return base_url + url
    return f"{base_url}/{url}"


def _discover_job_ids(stage_payload: dict[str, Any]) -> list[str]:
    ids: list[str] = []
    for src in (
        stage_payload.get("jobs"),
        (stage_payload.get("raw_result") or {}).get("jobs") if isinstance(stage_payload.get("raw_result"), dict) else None,
        (stage_payload.get("data") or {}).get("results", {}).get("jobs")
        if isinstance(stage_payload.get("data"), dict)
        and isinstance(stage_payload["data"].get("results"), dict)
        else None,
    ):
        if isinstance(src, list):
            for item in src:
                if isinstance(item, dict) and item.get("job_id"):
                    ids.append(str(item["job_id"]))
    for one in (
        (stage_payload.get("raw_result") or {}).get("result") if isinstance(stage_payload.get("raw_result"), dict) else None,
        (stage_payload.get("data") or {}).get("results", {}).get("result")
        if isinstance(stage_payload.get("data"), dict)
        and isinstance(stage_payload["data"].get("results"), dict)
        else None,
    ):
        if isinstance(one, dict) and one.get("job_id"):
            ids.append(str(one["job_id"]))
    return list(dict.fromkeys(ids))


def _parse_run_id_from_links(links: list[dict[str, Any]]) -> str | None:
    for link in links:
        url = str(link.get("url", ""))
        if "run_id=" in url:
            query = urllib.parse.urlparse(url).query
            parsed = urllib.parse.parse_qs(query)
            run_ids = parsed.get("run_id")
            if run_ids:
                return run_ids[0]
    return None


def run(base_url: str, out_root: pathlib.Path, timeout_jobs_s: int) -> int:
    session = requests.Session()
    health = session.get(f"{base_url}/health", timeout=30)
    health.raise_for_status()

    timestamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    root = out_root / f"helix-aws-strict-rerun-{timestamp}"
    root.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []

    for demo in DEMOS:
        demo_id = demo["id"]
        print(f"[{demo_id}] starting")
        demo_root = root / demo_id
        d01 = demo_root / "01_plan_input"
        d02 = demo_root / "02_execute"
        d03 = demo_root / "03_jobs"
        d04 = demo_root / "04_artifacts"
        d05 = demo_root / "05_bundles"
        d06 = demo_root / "06_verification"
        for d in (d01, d02, d03, d04, d05, d06):
            d.mkdir(parents=True, exist_ok=True)

        rec: dict[str, Any] = {"demo": demo_id, "errors": []}

        create_resp = session.post(f"{base_url}/create_session", timeout=30)
        create_resp.raise_for_status()
        session_id = create_resp.json()["session_id"]
        rec["session_id"] = session_id

        stage1_payload = {
            "command": demo["prompt"],
            "session_id": session_id,
            "execute_plan": False,
        }
        stage1 = session.post(f"{base_url}/execute", json=stage1_payload, timeout=1200).json()
        _save_json(d01 / "stage1.json", stage1)

        if demo["kind"] == "needs_inputs":
            stage2_payload = {
                "command": demo["followup"],
                "session_id": session_id,
                "execute_plan": False,
            }
        else:
            stage2_payload = {
                "command": demo["prompt"],
                "session_id": session_id,
                "execute_plan": True,
            }
        stage2 = session.post(f"{base_url}/execute", json=stage2_payload, timeout=1800).json()
        _save_json(d02 / "stage2.json", stage2)

        rec["stage2_status"] = stage2.get("status")
        rec["stage2_success"] = stage2.get("success")
        rec["run_id"] = stage2.get("run_id")

        job_ids = _discover_job_ids(stage2)
        rec["jobs_discovered_count"] = len(job_ids)
        completed = 0
        failed: list[dict[str, Any]] = []
        completed_job_ids: list[str] = []

        for jid in job_ids:
            t0 = time.time()
            final_job: dict[str, Any] | None = None
            while time.time() - t0 < timeout_jobs_s:
                job_resp = session.get(f"{base_url}/jobs/{jid}", timeout=60)
                if job_resp.status_code == 200:
                    payload = job_resp.json()
                    job = payload.get("job", payload) if isinstance(payload, dict) else {}
                    final_job = job if isinstance(job, dict) else {"raw": payload}
                    st = str((final_job or {}).get("status", "")).lower()
                    if st in {"completed", "failed", "error", "cancelled"}:
                        break
                time.sleep(3)
            if final_job is None:
                final_job = {"job_id": jid, "status": "timeout", "error": "Polling timeout"}
            _save_json(d03 / f"job_{jid}.json", final_job)
            st = str(final_job.get("status", "")).lower()
            if st == "completed":
                completed += 1
                completed_job_ids.append(jid)
            else:
                failed.append(
                    {
                        "job_id": jid,
                        "status": final_job.get("status"),
                        "error": final_job.get("error"),
                    }
                )
        rec["jobs_completed_count"] = completed
        rec["jobs_failed"] = failed

        downloads_dir = d04 / "downloads"
        downloads_dir.mkdir(parents=True, exist_ok=True)
        downloaded_files: list[str] = []

        # Download artifacts from stage2 response.
        links = _extract_links(stage2)
        for idx, link in enumerate(links, 1):
            label = _sanitize(str(link.get("label", f"artifact_{idx}")))
            fmt = str(link.get("format", "")).strip()
            if "." not in label and fmt:
                label = f"{label}.{fmt}"
            filename = f"{idx:02d}_{label}"
            abs_url = _to_absolute_url(base_url, str(link["url"]))
            try:
                response = session.get(abs_url, timeout=600)
                if response.status_code == 200:
                    dst = downloads_dir / filename
                    dst.write_bytes(response.content)
                    downloaded_files.append(str(dst.relative_to(demo_root)))
            except Exception as exc:
                rec["errors"].append(f"download_failed:{filename}:{exc}")

        # For pipeline demos, ask to visualize results of last completed job.
        if demo["kind"] == "pipeline" and completed_job_ids:
            last_job = completed_job_ids[-1]
            try:
                copy_payload = {"session_id": session_id}
                session.post(
                    f"{base_url}/jobs/{last_job}/copy-to-session",
                    json=copy_payload,
                    timeout=600,
                )
            except Exception:
                # Non-fatal; visualization may still work.
                pass
            viz_payload = {
                "command": f"visualize the results of job {last_job}",
                "session_id": session_id,
                "execute_plan": False,
            }
            viz = session.post(f"{base_url}/execute", json=viz_payload, timeout=1800).json()
            _save_json(d02 / f"stage3_visualize_{last_job}.json", viz)
            viz_links = _extract_links(viz)
            for idx, link in enumerate(viz_links, 1):
                label = _sanitize(str(link.get("label", f"viz_artifact_{idx}")))
                fmt = str(link.get("format", "")).strip()
                if "." not in label and fmt:
                    label = f"{label}.{fmt}"
                filename = f"viz_{idx:02d}_{label}"
                abs_url = _to_absolute_url(base_url, str(link["url"]))
                try:
                    response = session.get(abs_url, timeout=600)
                    if response.status_code == 200:
                        dst = downloads_dir / filename
                        dst.write_bytes(response.content)
                        downloaded_files.append(str(dst.relative_to(demo_root)))
                except Exception as exc:
                    rec["errors"].append(f"viz_download_failed:{filename}:{exc}")

            # Parse run_id from any bundle links after visualization.
            if rec.get("run_id") is None:
                rec["run_id"] = _parse_run_id_from_links(viz_links)

        # Download explicit bundle when run_id available.
        if rec.get("run_id"):
            run_id = str(rec["run_id"])
            bundle_url = (
                f"{base_url}/download/bundle?"
                f"session_id={urllib.parse.quote(str(session_id))}&run_id={urllib.parse.quote(run_id)}"
            )
            try:
                bundle_resp = session.get(bundle_url, timeout=900)
                if bundle_resp.status_code == 200:
                    bundle_path = d05 / f"bundle_{run_id}.zip"
                    bundle_path.write_bytes(bundle_resp.content)
                    downloaded_files.append(str(bundle_path.relative_to(demo_root)))
                    try:
                        unzip_dir = d05 / f"bundle_{run_id}_unzipped"
                        unzip_dir.mkdir(parents=True, exist_ok=True)
                        with zipfile.ZipFile(bundle_path, "r") as zf:
                            zf.extractall(unzip_dir)
                    except Exception as exc:
                        rec["errors"].append(f"bundle_unzip_failed:{exc}")
            except Exception as exc:
                rec["errors"].append(f"bundle_download_failed:{exc}")

        rec["downloaded_files"] = sorted(set(downloaded_files))

        # Strict success criteria for AWS validation.
        strict = bool(stage2.get("success"))
        if demo["kind"] == "pipeline":
            strict = strict and len(job_ids) > 0 and len(failed) == 0 and completed == len(job_ids)
        strict = strict and len(rec["downloaded_files"]) > 0
        rec["strict_success"] = strict
        if not strict:
            rec["errors"].append("strict_validation_failed")

        _save_json(d06 / "verification.json", rec)
        records.append(rec)
        print(
            f"[{demo_id}] strict={rec['strict_success']} "
            f"stage2={rec['stage2_status']} jobs={len(job_ids)} downloads={len(rec['downloaded_files'])}"
        )

    summary = {
        "base_url": base_url,
        "root": str(root),
        "all_strict_success": all(bool(r.get("strict_success")) for r in records),
        "records": records,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    _save_json(root / "summary.json", summary)
    (root / "README_NAVIGATION.md").write_text(
        "# AWS Strict Rerun\n\n"
        "Each demo contains:\n"
        "- 01_plan_input\n"
        "- 02_execute\n"
        "- 03_jobs\n"
        "- 04_artifacts\n"
        "- 05_bundles\n"
        "- 06_verification\n"
    )
    print(f"ROOT {root}")
    print(f"ALL_STRICT_SUCCESS {summary['all_strict_success']}")
    return 0 if summary["all_strict_success"] else 1


def main() -> int:
    parser = argparse.ArgumentParser(description="Run strict 5-demo rerun against AWS backend.")
    parser.add_argument(
        "--base-url",
        default=os.getenv("BACKEND_URL", DEFAULT_AWS_URL),
        help="Backend base URL (default: BACKEND_URL env or known AWS ALB).",
    )
    parser.add_argument(
        "--out-root",
        default=str(pathlib.Path.home() / "tmp"),
        help="Root directory for outputs (default: ~/tmp).",
    )
    parser.add_argument(
        "--job-timeout-s",
        type=int,
        default=1800,
        help="Per-job poll timeout in seconds (default: 1800).",
    )
    args = parser.parse_args()
    try:
        return run(
            base_url=args.base_url.rstrip("/"),
            out_root=pathlib.Path(args.out_root),
            timeout_jobs_s=args.job_timeout_s,
        )
    except requests.HTTPError as exc:
        print(f"HTTP error: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:
        print(f"Unexpected error: {exc}", file=sys.stderr)
        return 3


if __name__ == "__main__":
    raise SystemExit(main())

