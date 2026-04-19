#!/usr/bin/env python3
"""
Cloud smoke test for Helix.AI backend.

Runs a small suite of end-to-end checks against a running backend URL:
- /health
- /tools/list
- /execute: Q&A
- /execute: sequence alignment
- /execute: phylogenetic tree
- /execute: FastQC (S3 inputs) → local async job → poll /jobs/{id} → verify S3 outputs exist

Usage:
  python tests/demo_scenarios/cloud_smoke_test.py --backend-url https://<cloudfront-domain>
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from typing import Any, Dict, Optional, Tuple


def _http_json(method: str, url: str, *, timeout: int = 30, payload: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    import requests

    if method.upper() == "GET":
        r = requests.get(url, timeout=timeout)
    elif method.upper() == "POST":
        r = requests.post(url, json=payload or {}, timeout=timeout)
    else:
        raise ValueError(f"Unsupported method: {method}")
    r.raise_for_status()
    return r.json()


def _extract_job_id(execute_response: Dict[str, Any]) -> Optional[str]:
    # Standard response shape from backend/main.py
    data = execute_response.get("data") if isinstance(execute_response, dict) else None
    if isinstance(data, dict):
        results = data.get("results")
        if isinstance(results, dict):
            result = results.get("result")
            if isinstance(result, dict) and isinstance(result.get("job_id"), str):
                return result["job_id"]
    # Fallbacks (older shapes)
    for path in [
        ("job_id",),
        ("result", "job_id"),
        ("raw_result", "result", "job_id"),
    ]:
        cur: Any = execute_response
        ok = True
        for k in path:
            if not isinstance(cur, dict) or k not in cur:
                ok = False
                break
            cur = cur[k]
        if ok and isinstance(cur, str):
            return cur
    return None


def _poll_job(backend_url: str, job_id: str, *, max_seconds: int = 240) -> Dict[str, Any]:
    deadline = time.time() + max_seconds
    last_status = None
    while time.time() < deadline:
        payload = _http_json("GET", f"{backend_url}/jobs/{job_id}", timeout=15)
        job = payload.get("job") if isinstance(payload, dict) else None
        if not isinstance(job, dict):
            raise RuntimeError(f"Unexpected job payload: {payload!r}")
        status = (job.get("status") or "unknown").lower()
        if status != last_status:
            print(f"  job {job_id}: {status}")
            last_status = status
        if status in ("completed", "failed", "cancelled", "error"):
            return job
        time.sleep(10)
    raise TimeoutError(f"Job {job_id} did not finish within {max_seconds}s")


def _s3_split(uri: str) -> Tuple[str, str]:
    if not uri.startswith("s3://"):
        raise ValueError(f"Not an S3 URI: {uri}")
    rest = uri[len("s3://") :]
    bucket, _, key = rest.partition("/")
    return bucket, key


def _s3_list(prefix_uri: str) -> list[str]:
    import boto3

    bucket, prefix = _s3_split(prefix_uri)
    s3 = boto3.client("s3")
    keys: list[str] = []
    paginator = s3.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix.rstrip("/") + "/"):
        for obj in page.get("Contents", []) or []:
            keys.append(obj["Key"])
    return keys


def main() -> int:
    ap = argparse.ArgumentParser(description="Helix.AI cloud backend smoke test")
    ap.add_argument(
        "--backend-url",
        default="http://localhost:8001",
        help="Backend base URL (e.g. https://<cloudfront-domain>)",
    )
    ap.add_argument(
        "--fastqc-output",
        default="s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/fastqc",
        help="S3 output prefix for FastQC artifacts",
    )
    args = ap.parse_args()
    backend = args.backend_url.rstrip("/")

    print(f"== Smoke test against {backend} ==")

    health = _http_json("GET", f"{backend}/health", timeout=15)
    assert health.get("status") == "healthy", health
    print("✓ /health")

    tools = _http_json("GET", f"{backend}/tools/list", timeout=30)
    assert isinstance(tools.get("tools"), list), tools
    print(f"✓ /tools/list ({len(tools['tools'])} tools)")

    # Q&A
    qa = _http_json("POST", f"{backend}/execute", timeout=60, payload={"command": "What is a FASTQ file?"})
    assert qa.get("success") is True, qa
    print("✓ /execute (ask)")

    # Alignment
    aln_cmd = "Align these DNA sequences in FASTA format and return the alignment:\n>seq1\nACGTACGTACGT\n>seq2\nACGTTCGTACGT\n>seq3\nACGTACGTTCGT"
    aln = _http_json("POST", f"{backend}/execute", timeout=90, payload={"command": aln_cmd})
    assert aln.get("success") is True, aln
    print("✓ /execute (alignment)")

    # Tree
    tree_cmd = "Visualize the phylogenetic tree for these sequences:\n>seq1\nACGTACGTACGT\n>seq2\nACGTTCGTACGT\n>seq3\nACGTACGTTCGT"
    tree = _http_json("POST", f"{backend}/execute", timeout=120, payload={"command": tree_cmd})
    assert tree.get("success") is True, tree
    print("✓ /execute (phylogenetic_tree)")

    # FastQC (local async job)
    fastqc_cmd = (
        "Perform FastQC analysis the following forward R1 and reverse R2 reads: "
        "R1: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq "
        "R2: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq "
        f"Output the results in the following directory: {args.fastqc_output}"
    )
    fastqc = _http_json("POST", f"{backend}/execute", timeout=90, payload={"command": fastqc_cmd})
    assert fastqc.get("success") is True, fastqc
    job_id = _extract_job_id(fastqc)
    assert job_id, json.dumps(fastqc)[:2000]
    print(f"✓ /execute (fastqc) submitted job_id={job_id}")

    job = _poll_job(backend, job_id, max_seconds=300)
    if job.get("status") != "completed":
        raise RuntimeError(f"FastQC job did not complete successfully: {job}")
    print("✓ /jobs/{id} completed")

    results = _http_json("GET", f"{backend}/jobs/{job_id}/results", timeout=60)
    assert results.get("success") is True, results
    print("✓ /jobs/{id}/results")

    keys = _s3_list(args.fastqc_output)
    # require at least the two zip files
    want1 = args.fastqc_output.replace("s3://", "").split("/", 1)[1].rstrip("/") + "/test_mate_R1_fastqc.zip"
    want2 = args.fastqc_output.replace("s3://", "").split("/", 1)[1].rstrip("/") + "/test_mate_R2_fastqc.zip"
    assert any(k.endswith("test_mate_R1_fastqc.zip") for k in keys), want1
    assert any(k.endswith("test_mate_R2_fastqc.zip") for k in keys), want2
    print(f"✓ S3 outputs present under {args.fastqc_output}")

    print("== OK ==")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as e:
        print(f"FAILED: {e}", file=sys.stderr)
        raise

