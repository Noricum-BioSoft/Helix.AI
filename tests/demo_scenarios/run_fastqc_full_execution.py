#!/usr/bin/env python3
"""
Run actual FastQC execution (not just tool mapping).

This script triggers the full execution pipeline:
  1. Tool mapping (via agent)
  2. Infrastructure decision (Local vs EMR)
  3. ExecutionBroker submits job
  4. Wait for results
  5. Display output location

Usage:
    python run_fastqc_full_execution.py --small   # Use test dataset
    python run_fastqc_full_execution.py --large   # Use full dataset
"""

import sys
import os
import asyncio
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "backend"))

async def run_full_fastqc_execution(
    dataset_size: str = "small",
    backend_url: str = "http://localhost:8001",
    *,
    max_wait_seconds: int = 7200,
    poll_interval_seconds: int = 30,
):
    """
    Run complete FastQC workflow including actual execution.
    
    Args:
        dataset_size: "small" (test data) or "large" (full dataset)
        backend_url: URL of the running backend (default: http://localhost:8001)
    """
    import requests
    import json
    import time
    from datetime import datetime
    
    # Dataset configuration
    datasets = {
        "small": {
            "files": [
                "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq",
                "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq"
            ],
            "environment": "Local"
        },
        "large": {
            "files": [
                "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/mate_R1.fq",
                "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/mate_R2.fq"
            ],
            "environment": "EMR"
        }
    }
    
    config = datasets[dataset_size]
    
    print(f"\n{'='*70}")
    print(f"🧬 FULL FASTQC EXECUTION - {dataset_size.upper()} DATASET")
    print(f"{'='*70}\n")
    
    print(f"📁 Input files:")
    for f in config["files"]:
        print(f"   • {f}")
    print(f"\n🖥️  Environment: {config['environment']}")
    
    if config['environment'] == "EMR":
        print(f"\n⏱️  Note: This script will stay running and poll for job status")
        print(f"   until the EMR job completes or times out (max 2 hours).")
    
    print(f"{'='*70}\n")
    
    # Create session
    print(f"📞 Connecting to backend at {backend_url}...")
    try:
        session_response = requests.post(f"{backend_url}/create_session", timeout=10)
        session_response.raise_for_status()
        session_id = session_response.json()["session_id"]
        print(f"✅ Session created: {session_id}\n")
    except Exception as e:
        print(f"❌ Failed to create session: {e}")
        print(f"   Is the backend running? Try: uv run python -m backend.main")
        raise
    
    # Build command
    if len(config["files"]) > 1:
        command = f"Run FastQC quality analysis on paired-end reads: {config['files'][0]} and {config['files'][1]}"
    else:
        command = f"Run FastQC quality analysis on {config['files'][0]}"
    
    print(f"🚀 Starting FastQC execution...")
    print(f"   Command: {command}\n")
    
    try:
        # Execute command via backend API
        execute_response = requests.post(
            f"{backend_url}/execute",
            json={
                "command": command,
                "session_id": session_id
            },
            timeout=60  # Short timeout - should return immediately with job_id for async jobs
        )
        execute_response.raise_for_status()
        result = execute_response.json()
        
        print(f"\n✅ Command submitted!\n")
        
        # Check if this is an async job (returns job_id)
        job_id = None
        if isinstance(result, dict):
            # Check multiple locations for job_id in order of likelihood
            # Direct top-level
            if "job_id" in result:
                job_id = result["job_id"]
            # In result field
            elif "result" in result and isinstance(result["result"], dict):
                if "job_id" in result["result"]:
                    job_id = result["result"]["job_id"]
            # In data field (direct)
            elif "data" in result and isinstance(result["data"], dict):
                if "job_id" in result["data"]:
                    job_id = result["data"]["job_id"]
                # In data.results.result (nested structure from MCP response)
                elif "results" in result["data"] and isinstance(result["data"]["results"], dict):
                    results = result["data"]["results"]
                    if "job_id" in results:
                        job_id = results["job_id"]
                    elif "result" in results and isinstance(results["result"], dict):
                        if "job_id" in results["result"]:
                            job_id = results["result"]["job_id"]
            # In raw_result field
            elif "raw_result" in result and isinstance(result["raw_result"], dict):
                raw = result["raw_result"]
                if "job_id" in raw:
                    job_id = raw["job_id"]
                elif "result" in raw and isinstance(raw["result"], dict):
                    if "job_id" in raw["result"]:
                        job_id = raw["result"]["job_id"]
            
            # For large datasets, we expect an async job with job_id
            if not job_id and dataset_size == "large":
                print(f"⚠️  Warning: Expected job_id for large dataset but didn't find one")
                print(f"   Response structure: {list(result.keys())}")
                print(f"\n   Full response:\n{json.dumps(result, indent=2)}")
        
        if job_id:
            # Async execution - poll for status
            print(f"\n🔄 Job submitted to EMR")
            print(f"   Job ID: {job_id}")
            print(f"{'='*70}\n")
            print(f"⏳ Polling job status (checking every 30 seconds)...")
            print(f"   This script will keep running until the job completes.\n")
            
            max_wait_time = max_wait_seconds
            poll_interval = poll_interval_seconds
            elapsed = 0
            poll_count = 0
            last_status = None
            
            while elapsed < max_wait_time:
                poll_count += 1
                timestamp = datetime.now().strftime("%H:%M:%S")
                
                try:
                    status_response = requests.get(
                        f"{backend_url}/jobs/{job_id}",
                        timeout=10
                    )
                    status_response.raise_for_status()
                    status_payload = status_response.json()
                    status = status_payload.get("job") if isinstance(status_payload, dict) and "job" in status_payload else status_payload
                    if not isinstance(status, dict):
                        status = {}
                    
                    current_status = status.get("status", "unknown")
                    
                    # Only print if status changed or every 5 polls
                    if current_status != last_status or poll_count % 5 == 0:
                        print(f"   [{timestamp}] Poll #{poll_count}: {current_status}", end="")
                        
                        if "progress" in status:
                            print(f" ({status['progress']}%)", end="")
                        
                        # Show additional details if available
                        if "step" in status:
                            print(f" - {status['step']}", end="")
                        
                        print()
                        last_status = current_status
                    
                    if current_status in ["completed", "success", "succeeded"]:
                        print(f"\n{'='*70}")
                        print(f"✅ JOB COMPLETED SUCCESSFULLY!")
                        print(f"{'='*70}\n")
                        
                        # Get detailed results
                        print(f"📥 Fetching job results...\n")
                        
                        try:
                            results_response = requests.get(
                                f"{backend_url}/jobs/{job_id}/results",
                                timeout=30
                            )
                            results_response.raise_for_status()
                            results = results_response.json()
                            
                            print(f"{'='*70}")
                            print(f"JOB RESULTS (JSON)")
                            print(f"{'='*70}\n")
                            print(json.dumps(results, indent=2))
                            print(f"\n{'='*70}\n")
                            
                            # Extract key information if available
                            if isinstance(results, dict):
                                if "output_location" in results:
                                    print(f"📍 Output Location: {results['output_location']}")
                                if "artifacts" in results:
                                    print(f"📦 Artifacts:")
                                    for artifact in results.get("artifacts", []):
                                        print(f"   • {artifact}")
                                if "summary" in results:
                                    print(f"\n📊 Summary:")
                                    print(f"   {results['summary']}")
                            
                        except requests.RequestException as e:
                            print(f"⚠️  Could not fetch results: {e}")
                            print(f"   Status code: {getattr(e.response, 'status_code', 'N/A')}")
                            print(f"   Check manually: {backend_url}/jobs/{job_id}/results")
                        except Exception as e:
                            print(f"⚠️  Error processing results: {e}")
                            print(f"   Check manually: {backend_url}/jobs/{job_id}/results")
                        
                        return result
                    
                    elif current_status in ["failed", "error", "cancelled"]:
                        print(f"\n{'='*70}")
                        print(f"❌ JOB FAILED: {current_status}")
                        print(f"{'='*70}\n")
                        
                        print(f"ERROR DETAILS (JSON):")
                        print(json.dumps(status, indent=2))
                        
                        # Try to get detailed logs
                        print(f"\n{'='*70}")
                        print(f"FETCHING JOB LOGS...")
                        print(f"{'='*70}\n")
                        
                        try:
                            logs_response = requests.get(
                                f"{backend_url}/jobs/{job_id}/logs",
                                timeout=30
                            )
                            logs_response.raise_for_status()
                            logs = logs_response.json()
                            print(json.dumps(logs, indent=2))
                        except Exception as e:
                            print(f"⚠️  Could not fetch logs: {e}")
                            print(f"   Check manually: {backend_url}/jobs/{job_id}/logs")
                        
                        print(f"\n{'='*70}\n")
                        raise RuntimeError(f"Job failed with status: {current_status}")
                    
                    # Still running - wait and poll again
                    else:
                        # Job is still in progress
                        if poll_count == 1:
                            print(f"   Job is {current_status}, continuing to poll...")
                    
                    time.sleep(poll_interval)
                    elapsed += poll_interval
                    
                except requests.RequestException as e:
                    print(f"   [{timestamp}] ⚠️  Network error checking status: {e}")
                    print(f"   Will retry in {poll_interval}s...")
                    time.sleep(poll_interval)
                    elapsed += poll_interval
                except Exception as e:
                    print(f"   [{timestamp}] ⚠️  Unexpected error: {e}")
                    print(f"   Will retry in {poll_interval}s...")
                    time.sleep(poll_interval)
                    elapsed += poll_interval
            
            if elapsed >= max_wait_time:
                print(f"\n{'='*70}")
                print(f"⚠️  TIMEOUT REACHED")
                print(f"{'='*70}\n")
                print(f"Job still running after {max_wait_time}s ({max_wait_time//3600}h)")
                print(f"\nManual commands:")
                print(f"  • Check status: curl {backend_url}/jobs/{job_id}")
                print(f"  • Get results:  curl {backend_url}/jobs/{job_id}/results")
                print(f"  • View logs:    curl {backend_url}/jobs/{job_id}/logs")
                print(f"\n{'='*70}\n")
        else:
            # Sync execution - result is already available
            print(f"{'='*70}")
            print(f"✅ EXECUTION COMPLETED!")
            print(f"{'='*70}\n")
            
            print(f"EXECUTION RESPONSE (JSON):")
            print(f"{'='*70}\n")
            print(json.dumps(result, indent=2))
            print(f"\n{'='*70}\n")
        
        print(f"\n{'='*70}")
        print(f"NEXT STEPS")
        print(f"{'='*70}\n")
        
        if job_id:
            print(f"✅ Job ID: {job_id}")
            if config["environment"] == "EMR":
                print(f"\n1. View results:")
                print(f"   curl {backend_url}/jobs/{job_id}/results")
                print(f"\n2. Check logs:")
                print(f"   curl {backend_url}/jobs/{job_id}/logs")
                print(f"\n3. Monitor EMR (optional):")
                print(f"   aws emr list-clusters --active")
        else:
            print(f"📝 Response received from backend")
            print(f"   Check the response above for details")
        
        print(f"\n{'='*70}\n")
        
        return result
        
    except Exception as e:
        print(f"\n❌ Execution failed: {e}")
        print(f"\nTroubleshooting:")
        print(f"  • Check AWS credentials are configured")
        print(f"  • Verify S3 bucket access permissions")
        print(f"  • Ensure EMR cluster permissions (for large dataset)")
        print(f"  • Check backend logs for details")
        raise


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Run full FastQC execution with results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_fastqc_full_execution.py --small
  python run_fastqc_full_execution.py --large

This runs the COMPLETE workflow:
  1. ✅ Tool mapping (identify fastqc_quality_analysis)
  2. ✅ Infrastructure decision (Local vs EMR)
  3. ✅ Job submission via ExecutionBroker
  4. ✅ Actual FastQC analysis
  5. ✅ Results saved to S3
        """
    )
    
    parser.add_argument(
        "--small",
        action="store_const",
        const="small",
        dest="size",
        help="Use small test dataset (fast, runs locally)"
    )
    parser.add_argument(
        "--large",
        action="store_const",
        const="large",
        dest="size",
        help="Use large dataset (slower, runs on EMR)"
    )
    parser.add_argument(
        "--backend-url",
        default="http://localhost:8001",
        help="Backend URL (default: http://localhost:8001)"
    )
    parser.add_argument(
        "--max-wait-seconds",
        type=int,
        default=7200,
        help="Max seconds to wait/poll for async jobs (default: 7200)"
    )
    parser.add_argument(
        "--poll-interval-seconds",
        type=int,
        default=30,
        help="Polling interval in seconds for async jobs (default: 30)"
    )
    
    args = parser.parse_args()
    
    if not args.size:
        parser.print_help()
        sys.exit(1)
    
    # Check dependencies
    try:
        import boto3
    except ImportError:
        print("❌ Error: boto3 not installed")
        print("   Install with: pip install boto3")
        sys.exit(1)
    
    # Check AWS credentials
    try:
        import boto3
        sts = boto3.client('sts')
        sts.get_caller_identity()
    except Exception as e:
        print(f"❌ Error: AWS credentials not configured")
        print(f"   {e}")
        print(f"\n   Configure with: aws configure")
        sys.exit(1)
    
    # Run execution
    asyncio.run(
        run_full_fastqc_execution(
            args.size,
            args.backend_url,
            max_wait_seconds=args.max_wait_seconds,
            poll_interval_seconds=args.poll_interval_seconds,
        )
    )


if __name__ == "__main__":
    main()
