#!/usr/bin/env python3
"""
Test validation with a REAL job (checks actual S3 files).

Usage:
    python3 scripts/test_validation_real_job.py <job_id>
"""

import sys
import os
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from backend.job_manager import get_job_manager
import json

def run_real_job_validation(job_id: str):
    """Run validation with a real job ID (CLI helper; not a pytest test)."""
    print("=" * 70)
    print(f"Testing Validation for REAL Job: {job_id}")
    print("=" * 70)
    print()
    
    jm = get_job_manager()
    
    try:
        # Get job status (this will trigger validation if job is completed)
        print("1. Getting job status...")
        job = jm.get_job_status(job_id)
        print(f"   Status: {job.get('status')}")
        print(f"   Tool: {job.get('tool_name')}")
        print(f"   Output path: {job.get('output_path')}")
        print(f"   Step ID: {job.get('step_id')}")
        print()
        
        if job.get('status') != 'completed':
            print(f"⚠️  Job is not completed (status: {job.get('status')})")
            print("   Validation only runs for completed jobs.")
            return
        
        # Get job results (this includes validation)
        print("2. Getting job results with validation...")
        try:
            results = jm.get_job_results(job_id)
            print(f"   Results path: {results.get('results_path')}")
            print()
            
            # Show validation results
            validation = results.get('output_validation', {})
            print("3. Validation Results:")
            print("-" * 70)
            print(f"   results_json_exists: {validation.get('results_json_exists')}")
            print(f"   all_files_exist: {validation.get('all_files_exist')}")
            print(f"   actual_results_path: {validation.get('actual_results_path')}")
            print()
            
            if validation.get('validation_warnings'):
                print("   Warnings:")
                for i, warning in enumerate(validation.get('validation_warnings', []), 1):
                    print(f"     {i}. {warning}")
                print()
            
            if validation.get('output_files_found'):
                print("   Output Files Checked:")
                for i, file_info in enumerate(validation.get('output_files_found', []), 1):
                    status = "✅ EXISTS" if file_info.get('exists') else "❌ MISSING"
                    print(f"     {i}. {status}: {file_info.get('location')}")
                    print(f"        Path: {file_info.get('path')}")
                    if file_info.get('error'):
                        print(f"        Error: {file_info.get('error')}")
                print()
            
            # Summary
            print("=" * 70)
            print("SUMMARY")
            print("=" * 70)
            if validation.get('results_json_exists'):
                print("✅ results.json exists")
                if validation.get('all_files_exist'):
                    print("✅ All output files exist")
                else:
                    print("❌ Some output files are missing")
                    print("   Check the output_files_found list above for details")
            else:
                print("❌ results.json NOT FOUND")
                print("   Checked multiple locations (see warnings above)")
                print("   Job may not have completed successfully")
            
        except ValueError as e:
            print(f"❌ Error: {e}")
            print("   Job may not be completed yet")
        
    except Exception as e:
        print(f"❌ Error testing job: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python3 scripts/test_validation_real_job.py <job_id>")
        print()
        print("Example:")
        print("  python3 scripts/test_validation_real_job.py 92bd8bec-0575-4681-a6a2-a6f1cedf48ed")
        sys.exit(1)
    
    job_id = sys.argv[1]
    run_real_job_validation(job_id)

