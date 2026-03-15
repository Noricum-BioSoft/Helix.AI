#!/usr/bin/env python3
"""
End-to-end test: 16S demo + FastQC/MultiQC follow-ups via API.

Run with backend on localhost:8001:
  HELIX_MOCK_MODE=1 uvicorn backend.main_with_mcp:app --host 0.0.0.0 --port 8001

Then: python tests/e2e_fastqc_multiqc_followups.py

HELIX_MOCK_MODE=1 uses the deterministic router for FastQC/MultiQC (avoids agent timeout).
"""

import json
import os
import sys
import time
from pathlib import Path

import requests

BASE = os.environ.get("HELIX_API_BASE", "http://localhost:8001")


def create_session():
    r = requests.post(f"{BASE}/session/create", json={}, timeout=10)
    r.raise_for_status()
    return r.json().get("session_id")


def execute(session_id: str, prompt: str):
    r = requests.post(
        f"{BASE}/execute",
        json={"command": prompt, "session_id": session_id},
        timeout=120,
    )
    r.raise_for_status()
    return r.json()


def main():
    print("E2E: 16S demo + FastQC/MultiQC follow-ups")
    print("=" * 60)

    sid = create_session()
    print(f"Session: {sid}")

    # 1. Run 16S amplicon demo (may take a while)
    demo = """You are processing a 16S rRNA amplicon sequencing dataset from a gut microbiome study. Raw paired-end FASTQ files are on S3 and need a full preprocessing pipeline before downstream diversity analysis.

Input files:
- R1: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq
- R2: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq

Pipeline steps (run in order):
1. FastQC quality assessment on both R1 and R2 files.
2. Read trimming to remove low-quality bases and adapter sequences.
3. Read merging of trimmed paired-end reads.
4. Final quality report.

Output prefix: s3://noricum-ngs-data/test-output/amplicon-demo/"""

    print("\n1. Running 16S amplicon demo...")
    t0 = time.time()
    r1 = execute(sid, demo)
    print(f"   Done in {time.time()-t0:.1f}s")
    success1 = r1.get("success", False)
    print(f"   success={success1}")

    # 2. "Run FastQC on merged reads" -> expect clear merged-FASTA explanation
    print("\n2. 'Run FastQC quality assessment on the merged reads'")
    r2 = execute(sid, "Run FastQC quality assessment on the merged reads.")
    text2 = (r2.get("result") or {}).get("text") or r2.get("text") or ""
    ok2 = "FastQC cannot run on merged" in text2 or "FASTA" in text2
    print(f"   ok={ok2} | text[:200]={text2[:200]}...")

    # 3. "perform FastQC analysis on merged.fasta"
    print("\n3. 'perform FastQC analysis on s3://.../merged.fasta'")
    r3 = execute(sid, "perform FastQC analysis on s3://noricum-ngs-data/test-output/amplicon-demo/merged.fasta")
    text3 = (r3.get("result") or {}).get("text") or r3.get("text") or ""
    ok3 = "FASTA" in text3 or "cannot run" in text3.lower()
    print(f"   ok={ok3} | text[:200]={text3[:200]}...")

    # 4. "Run FastQC on the raw reads" -> session-aware R1/R2
    print("\n4. 'Run FastQC on the raw reads from the pipeline'")
    r4 = execute(sid, "Run FastQC on the raw reads from the pipeline.")
    # May run FastQC (success) or return needs_inputs - both acceptable for E2E
    success4 = r4.get("success", False)
    err4 = (r4.get("result") or {}).get("text") or ""
    ok4 = success4 or ("test_mate" in err4) or ("input_r1" in str(r4))
    print(f"   ok={ok4} | success={success4}")

    # 5. "run MultiQC"
    print("\n5. 'run MultiQC on s3://.../merged.fasta'")
    r5 = execute(sid, "run MultiQC on s3://noricum-ngs-data/test-output/amplicon-demo/merged.fasta")
    text5 = (r5.get("result") or {}).get("text") or r5.get("text") or ""
    ok5 = "not supported" in text5.lower() or "MultiQC" in text5
    ok5 = ok5 and ("Alternative" in text5 or "FastQC" in text5)
    print(f"   ok={ok5} | text[:200]={text5[:200]}...")

    # Summary
    print("\n" + "=" * 60)
    all_ok = ok2 and ok3 and ok4 and ok5
    print(f"Result: {'PASS' if all_ok else 'FAIL'}")
    print(f"  Test 2 (merged reads): {ok2}")
    print(f"  Test 3 (merged.fasta): {ok3}")
    print(f"  Test 4 (raw reads):    {ok4}")
    print(f"  Test 5 (MultiQC):      {ok5}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
