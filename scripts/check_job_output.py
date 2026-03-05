#!/usr/bin/env python3
"""
Helper script to check job results and verify output file existence.

Usage:
    python3 scripts/check_job_output.py <job_id> [results_path]
"""

import json
import subprocess
import sys
import os
from typing import Optional, Dict, Any, List, Tuple

def download_s3_file(s3_path: str) -> Optional[str]:
    """Download S3 file and return contents as string."""
    try:
        result = subprocess.run(
            ['aws', 's3', 'cp', s3_path, '-'],
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to download {s3_path}: {e.stderr}", file=sys.stderr)
        return None
    except FileNotFoundError:
        print("❌ AWS CLI not found. Please install AWS CLI.", file=sys.stderr)
        return None

def check_s3_file_exists(s3_path: str) -> Tuple[bool, Optional[str]]:
    """Check if S3 file exists."""
    try:
        result = subprocess.run(
            ['aws', 's3', 'ls', s3_path],
            capture_output=True,
            text=True,
            check=True
        )
        return True, None
    except subprocess.CalledProcessError as e:
        return False, e.stderr

def get_s3_file_size(s3_path: str) -> Optional[int]:
    """Get S3 file size in bytes."""
    try:
        result = subprocess.run(
            ['aws', 's3', 'ls', s3_path],
            capture_output=True,
            text=True,
            check=True
        )
        # Parse output: "2024-01-01 12:00:00     123456 file.fq"
        parts = result.stdout.strip().split()
        if len(parts) >= 4:
            return int(parts[2])
        return None
    except Exception:
        return None

def find_output_paths(results: Dict[str, Any]) -> List[Tuple[str, str]]:
    """Extract all potential output paths from results."""
    paths = []
    
    # Check tool result
    if 'result' in results:
        tool_result = results['result']
        
        # Direct output_path
        if 'output_path' in tool_result:
            paths.append(('result.output_path', tool_result['output_path']))
        
        # In summary
        if 'summary' in tool_result and isinstance(tool_result['summary'], dict):
            if 'output_path' in tool_result['summary']:
                paths.append(('result.summary.output_path', tool_result['summary']['output_path']))
        
        # In output dict
        if 'output' in tool_result and isinstance(tool_result['output'], dict):
            if 'output_path' in tool_result['output']:
                paths.append(('result.output.output_path', tool_result['output']['output_path']))
    
    # Check payload arguments
    if 'payload' in results:
        payload = results['payload']
        if 'arguments' in payload:
            args = payload['arguments']
            if 'output' in args:
                paths.append(('payload.arguments.output', args['output']))
    
    return paths

def check_job_output(job_id: str, results_path: Optional[str] = None):
    """Check job results and verify output file existence."""
    print("=" * 70)
    print(f"Checking Job Output: {job_id}")
    print("=" * 70)
    
    # If results_path not provided, try to get from job manager
    if not results_path:
        try:
            from backend.job_manager import get_job_manager
            jm = get_job_manager()
            job = jm.get_job_status(job_id)
            if 'results_path' in job:
                results_path = job['results_path']
            elif 'output_path' in job:
                output_path = job['output_path']
                results_path = output_path.rstrip('/') + '/results.json'
        except Exception as e:
            print(f"⚠️  Could not get results_path from job manager: {e}")
            print("   Please provide results_path manually.")
            return
    
    if not results_path:
        print("❌ No results_path provided or found")
        return
    
    if not results_path.startswith('s3://'):
        print(f"❌ Invalid results_path (must be S3 URI): {results_path}")
        return
    
    print(f"\n📥 Downloading results.json from: {results_path}")
    results_json = download_s3_file(results_path)
    
    if not results_json:
        print("❌ Failed to download results.json")
        return
    
    try:
        results = json.loads(results_json)
    except json.JSONDecodeError as e:
        print(f"❌ Failed to parse results.json: {e}")
        return
    
    print("✅ Successfully loaded results.json\n")
    
    # Display job status
    print("=" * 70)
    print("JOB STATUS")
    print("=" * 70)
    print(f"Status: {results.get('status', 'unknown')}")
    if 'error' in results and results['error']:
        print(f"Error: {results['error']}")
    if 'duration_s' in results:
        print(f"Duration: {results['duration_s']:.2f}s")
    print()
    
    # Display payload info
    if 'payload' in results:
        payload = results['payload']
        print("=" * 70)
        print("PAYLOAD INFO")
        print("=" * 70)
        print(f"Tool: {payload.get('tool_name', 'unknown')}")
        if 'arguments' in payload:
            args = payload['arguments']
            print("\nArguments:")
            for key, value in args.items():
                if isinstance(value, str) and len(value) > 100:
                    print(f"  {key}: {value[:100]}... (truncated)")
                else:
                    print(f"  {key}: {value}")
        print()
    
    # Display tool result
    if 'result' in results:
        tool_result = results['result']
        print("=" * 70)
        print("TOOL RESULT")
        print("=" * 70)
        print(f"Status: {tool_result.get('status', 'unknown')}")
        if 'text' in tool_result:
            print(f"Text: {tool_result['text']}")
        if 'error' in tool_result:
            print(f"Error: {tool_result['error']}")
        
        if 'summary' in tool_result:
            print("\nSummary:")
            summary = tool_result['summary']
            if isinstance(summary, dict):
                for key, value in summary.items():
                    print(f"  {key}: {value}")
        print()
    
    # Find and check output paths
    output_paths = find_output_paths(results)
    
    print("=" * 70)
    print("OUTPUT FILE VERIFICATION")
    print("=" * 70)
    
    if not output_paths:
        print("⚠️  No output_path found in results!")
        print("\nSearched in:")
        print("  - result.output_path")
        print("  - result.summary.output_path")
        print("  - result.output.output_path")
        print("  - payload.arguments.output")
        
        # Show full result structure for debugging
        print("\n📋 Full result structure (first 2000 chars):")
        result_str = json.dumps(results.get('result', {}), indent=2)
        print(result_str[:2000])
        if len(result_str) > 2000:
            print("... (truncated)")
    else:
        all_exist = True
        for location, path in output_paths:
            print(f"\n📍 {location}")
            print(f"   Path: {path}")
            
            exists, error = check_s3_file_exists(path)
            if exists:
                size = get_s3_file_size(path)
                size_str = f" ({size:,} bytes)" if size else ""
                print(f"   ✅ EXISTS{size_str}")
            else:
                all_exist = False
                print(f"   ❌ DOES NOT EXIST")
                if error:
                    print(f"   Error: {error.strip()}")
        
        print()
        if all_exist:
            print("✅ All output files exist!")
        else:
            print("❌ Some output files are missing!")
    
    print("\n" + "=" * 70)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python3 scripts/check_job_output.py <job_id> [results_path]")
        sys.exit(1)
    
    job_id = sys.argv[1]
    results_path = sys.argv[2] if len(sys.argv) > 2 else None
    
    check_job_output(job_id, results_path)



