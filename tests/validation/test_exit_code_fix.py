#!/usr/bin/env python3
"""Test script to verify exit code fix works."""
import sys
import os
import subprocess
import json
import time
from pathlib import Path

# Add project root
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

from backend.job_manager import get_job_manager

print("=" * 70)
print("TESTING EXIT CODE FIX - PROOF")
print("=" * 70)
print()
print("This script will:")
print("  1. Submit a test job")
print("  2. Wait for completion")
print("  3. Check EMR status vs results.json")
print("  4. Check runner.log for exit code message")
print()

# Unset EMR_CLUSTER_ID to test auto-creation
original_cluster_id = os.environ.pop("EMR_CLUSTER_ID", None)

try:
    jm = get_job_manager()
    
    # Submit a small test job
    r1_path = "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq"
    r2_path = "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq"
    output_path = f"s3://noricum-ngs-data/test-output/exit_code_test_{int(time.time())}.fq"
    
    print(f"Submitting test job...")
    print(f"  R1: {r1_path}")
    print(f"  R2: {r2_path}")
    print(f"  Output: {output_path}")
    print()
    
    job_id = jm.submit_universal_emr_job(
        tool_name="read_merging",
        tool_args={
            "forward_reads": r1_path,
            "reverse_reads": r2_path,
            "output": output_path,
            "min_overlap": 12,
        },
    )
    
    print(f"✅ Job submitted: {job_id}")
    job = jm.jobs.get(job_id)
    step_id = job.get("step_id") if job else None
    cluster_id = job.get("cluster_id") if job else None
    
    print(f"Step ID: {step_id}")
    print(f"Cluster ID: {cluster_id}")
    print()
    print("Waiting for job to complete (max 5 minutes)...")
    
    # Wait for completion
    for i in range(30):  # 30 * 10 seconds = 5 minutes
        time.sleep(10)
        job = jm.get_job_status(job_id)
        status = job.get("status")
        print(f"  [{i+1}/30] Status: {status}", end="\r")
        sys.stdout.flush()
        
        if status in ("completed", "failed"):
            break
    
    print()
    print()
    print("=" * 70)
    print("RESULTS")
    print("=" * 70)
    
    job = jm.get_job_status(job_id)
    emr_status = job.get("status")
    print(f"EMR Status: {emr_status}")
    
    results_status = None
    # Check results.json
    if step_id:
        region = os.getenv("AWS_REGION", "us-west-1")
        s3_bucket = os.getenv("S3_DATASET_BUCKET", "noricum-ngs-data")
        results_path = f"s3://{s3_bucket}/emr-results/{job_id}/{step_id}/results.json"
        
        try:
            result = subprocess.run(
                ["aws", "s3", "cp", results_path, "-", "--region", region],
                capture_output=True,
                text=True,
                timeout=30,
            )
            if result.returncode == 0:
                results = json.loads(result.stdout)
                results_status = results.get("status")
                print(f"results.json status: {results_status}")
                
                if results_status == "success":
                    summary = results.get("result", {}).get("summary", {})
                    print(f"  Total pairs: {summary.get('total_pairs', 0)}")
                    print(f"  Merged pairs: {summary.get('merged_pairs', 0)}")
        except Exception as e:
            print(f"Could not read results.json: {e}")
    
    # Check runner.log for exit code message
    exit_code_found = False
    if step_id:
        region = os.getenv("AWS_REGION", "us-west-1")
        s3_bucket = os.getenv("S3_DATASET_BUCKET", "noricum-ngs-data")
        
        # Try runner.log in results directory
        log_paths = [
            f"s3://{s3_bucket}/emr-results/{job_id}/{step_id}/runner.log",
        ]
        
        # Also try wrapper.log
        if cluster_id:
            log_paths.append(f"s3://{s3_bucket}/emr-logs/{cluster_id}/steps/{step_id}/wrapper.log")
        
        for log_path in log_paths:
            try:
                result = subprocess.run(
                    ["aws", "s3", "cp", log_path, "-", "--region", region],
                    capture_output=True,
                    text=True,
                    timeout=30,
                )
                if result.returncode == 0 and result.stdout.strip():
                    print()
                    print(f"Found log: {log_path}")
                    print("-" * 70)
                    lines = result.stdout.strip().split('\n')
                    # Look for exit code message
                    for line in lines:
                        if "exiting with code" in line.lower() or "Job completed successfully" in line:
                            print(f"✅ {line}")
                            exit_code_found = True
                    
                    # Show last 10 lines
                    print("\nLast 10 lines of log:")
                    for line in lines[-10:]:
                        print(f"  {line}")
                    break
            except Exception:
                continue
    
    print()
    print("=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    
    success = False
    if emr_status == "completed":
        print("🎉 SUCCESS: EMR marked job as completed!")
        print("✅ Exit code fix is working correctly!")
        success = True
    elif emr_status == "failed":
        if results_status == "success":
            print("❌ FAILURE: EMR marked job as failed but results.json shows success")
            print("   This means the exit code fix is NOT working")
            if exit_code_found:
                print("   However, exit code message was found in logs")
                print("   This suggests the fix is in place but something else is wrong")
            else:
                print("   Exit code message not found in logs")
        else:
            print("⚠️  Job actually failed (both EMR and results.json show failure)")
    else:
        print(f"⚠️  Job status: {emr_status} (incomplete)")
    
    sys.exit(0 if success else 1)
        
finally:
    if original_cluster_id:
        os.environ["EMR_CLUSTER_ID"] = original_cluster_id

