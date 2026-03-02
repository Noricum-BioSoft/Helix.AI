#!/usr/bin/env python3
"""
Check the status of a running FastQC job.

Usage:
    # Check job status once
    python3 scripts/check_job_status.py <job_id>
    
    # Watch job status (poll every 30 seconds)
    python3 scripts/check_job_status.py <job_id> --watch
    
    # Check job status with custom backend URL
    python3 scripts/check_job_status.py <job_id> --backend-url http://localhost:8001

Examples:
    python3 scripts/check_job_status.py 2632bfbd-5e36-45a3-bae5-9d0a7dcab5a3
    python3 scripts/check_job_status.py 2632bfbd-5e36-45a3-bae5-9d0a7dcab5a3 --watch
"""

import argparse
import json
import sys
import time
from datetime import datetime
from typing import Dict, Any
import requests

# ANSI color codes
GREEN = "\033[92m"
YELLOW = "\033[93m"
RED = "\033[91m"
BLUE = "\033[94m"
CYAN = "\033[96m"
RESET = "\033[0m"
BOLD = "\033[1m"

# Status to color mapping
STATUS_COLORS = {
    "submitted": CYAN,
    "pending": CYAN,
    "running": BLUE,
    "completed": GREEN,
    "failed": RED,
    "cancelled": YELLOW,
}

def colorize(text: str, color: str) -> str:
    """Add color to text."""
    return f"{color}{text}{RESET}"

def format_duration(seconds: float) -> str:
    """Format duration in human-readable format."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        return f"{seconds/60:.1f}m"
    else:
        return f"{seconds/3600:.1f}h"

def get_job_status(job_id: str, backend_url: str = "http://localhost:8001") -> Dict[str, Any]:
    """Fetch job status from backend API."""
    try:
        response = requests.get(f"{backend_url}/jobs/{job_id}", timeout=10)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"{RED}❌ Failed to fetch job status: {e}{RESET}", file=sys.stderr)
        sys.exit(1)

def get_job_logs(job_id: str, backend_url: str = "http://localhost:8001") -> Dict[str, Any]:
    """Fetch job logs from backend API."""
    try:
        response = requests.get(f"{backend_url}/jobs/{job_id}/logs", timeout=10)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException:
        return {}

def print_separator(char="=", length=70):
    """Print a separator line."""
    print(char * length)

def display_job_status(job_data: Dict[str, Any], show_full: bool = True):
    """Display formatted job status."""
    if not job_data.get("success"):
        print(f"{RED}❌ Error: {job_data.get('error', 'Unknown error')}{RESET}")
        return
    
    job = job_data.get("job", {})
    
    # Header
    print_separator()
    print(f"{BOLD}JOB STATUS{RESET}")
    print_separator()
    
    # Basic info
    status = job.get("status", "unknown")
    status_color = STATUS_COLORS.get(status, RESET)
    
    print(f"Job ID:     {BOLD}{job.get('job_id', 'N/A')}{RESET}")
    print(f"Tool:       {job.get('tool_name', 'N/A')}")
    print(f"Status:     {colorize(status.upper(), status_color)}")
    
    # Timing info
    created_at = job.get("created_at")
    if created_at:
        print(f"Created:    {created_at}")
    
    started_at = job.get("started_at")
    if started_at:
        print(f"Started:    {started_at}")
    
    completed_at = job.get("completed_at")
    if completed_at:
        print(f"Completed:  {completed_at}")
    
    # Duration
    if status == "completed" and "duration_s" in job:
        duration = format_duration(job["duration_s"])
        print(f"Duration:   {duration}")
    elif started_at and not completed_at:
        # Calculate elapsed time for running jobs
        try:
            start_time = datetime.fromisoformat(started_at.replace('Z', '+00:00'))
            elapsed = (datetime.now(start_time.tzinfo) - start_time).total_seconds()
            print(f"Elapsed:    {format_duration(elapsed)}")
        except Exception:
            pass
    
    print()
    
    # Infrastructure info
    if show_full:
        print_separator("-")
        print(f"{BOLD}INFRASTRUCTURE{RESET}")
        print_separator("-")
        
        cluster_id = job.get("cluster_id") or job.get("emr_cluster_id")
        if cluster_id:
            print(f"Cluster ID: {cluster_id}")
        
        step_id = job.get("step_id") or job.get("emr_step_id")
        if step_id:
            print(f"Step ID:    {step_id}")
        
        print()
    
    # Input/Output info
    if show_full:
        print_separator("-")
        print(f"{BOLD}FILES{RESET}")
        print_separator("-")
        
        # Inputs
        inputs = job.get("inputs", {})
        if isinstance(inputs, dict):
            r1 = inputs.get("r1_path") or inputs.get("input_r1")
            r2 = inputs.get("r2_path") or inputs.get("input_r2")
            if r1:
                print(f"Input R1:   {r1}")
            if r2:
                print(f"Input R2:   {r2}")
        
        # Output
        output = job.get("output_path")
        if output:
            print(f"Output:     {output}")
        
        print()
    
    # Error info (if failed)
    if status == "failed":
        print_separator("-")
        print(f"{BOLD}{RED}ERROR DETAILS{RESET}")
        print_separator("-")
        
        error = job.get("error", "No error details available")
        print(f"{error}")
        print()
    
    # Progress indicator for running jobs
    if status in ("pending", "running", "submitted"):
        print_separator("-")
        print(f"{BOLD}PROGRESS{RESET}")
        print_separator("-")
        
        if status == "pending":
            print(f"{CYAN}⏳ Job is pending... (EMR cluster may be starting up){RESET}")
            print(f"{CYAN}   This typically takes 5-10 minutes for a new cluster.{RESET}")
        elif status == "running":
            print(f"{BLUE}🔄 Job is running...{RESET}")
        elif status == "submitted":
            print(f"{CYAN}📤 Job submitted, waiting for cluster...{RESET}")
        
        print()
    
    # Success message
    if status == "completed":
        print_separator("-")
        print(f"{GREEN}✅ Job completed successfully!{RESET}")
        print_separator("-")
        
        # Check for results
        results_path = job.get("results_path")
        if results_path:
            print(f"\n{BOLD}View results:{RESET}")
            print(f"  python3 scripts/check_job_output.py {job.get('job_id')}")
        
        print()

def watch_job_status(job_id: str, backend_url: str = "http://localhost:8001", interval: int = 30):
    """Watch job status with periodic updates."""
    print(f"\n{BOLD}Watching job: {job_id}{RESET}")
    print(f"(Polling every {interval} seconds. Press Ctrl+C to stop)\n")
    
    try:
        while True:
            # Clear screen (optional - comment out if you want history)
            # print("\033[2J\033[H")
            
            # Get and display status
            job_data = get_job_status(job_id, backend_url)
            display_job_status(job_data, show_full=False)
            
            # Check if job is terminal
            job = job_data.get("job", {})
            status = job.get("status", "unknown")
            
            if status in ("completed", "failed", "cancelled"):
                print(f"\n{BOLD}Job reached terminal state: {status}{RESET}")
                print(f"\nView full details:")
                print(f"  python3 scripts/check_job_status.py {job_id}")
                if status == "completed":
                    print(f"\nView results:")
                    print(f"  python3 scripts/check_job_output.py {job_id}")
                break
            
            # Wait before next check
            print(f"\n{colorize(f'Next check in {interval}s...', CYAN)}", end="", flush=True)
            time.sleep(interval)
            print("\r" + " " * 50 + "\r", end="")  # Clear the line
            
    except KeyboardInterrupt:
        print(f"\n\n{YELLOW}Stopped watching job.{RESET}")
        print(f"\nCheck status anytime:")
        print(f"  python3 scripts/check_job_status.py {job_id}")

def main():
    parser = argparse.ArgumentParser(
        description="Check FastQC job status",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Check job status once
  python3 scripts/check_job_status.py 2632bfbd-5e36-45a3-bae5-9d0a7dcab5a3
  
  # Watch job status (auto-update every 30s)
  python3 scripts/check_job_status.py 2632bfbd-5e36-45a3-bae5-9d0a7dcab5a3 --watch
  
  # Check with custom backend URL
  python3 scripts/check_job_status.py <job_id> --backend-url http://localhost:8001
        """
    )
    
    parser.add_argument("job_id", help="Job ID to check")
    parser.add_argument(
        "--backend-url",
        default="http://localhost:8001",
        help="Backend URL (default: http://localhost:8001)"
    )
    parser.add_argument(
        "--watch",
        action="store_true",
        help="Watch job status with periodic updates"
    )
    parser.add_argument(
        "--interval",
        type=int,
        default=30,
        help="Polling interval in seconds for --watch mode (default: 30)"
    )
    
    args = parser.parse_args()
    
    if args.watch:
        watch_job_status(args.job_id, args.backend_url, args.interval)
    else:
        job_data = get_job_status(args.job_id, args.backend_url)
        display_job_status(job_data)

if __name__ == "__main__":
    main()
