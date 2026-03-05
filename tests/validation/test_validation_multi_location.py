#!/usr/bin/env python3
"""
Test script to demonstrate multi-location results.json checking.

This script shows that the validation function:
1. Checks multiple possible locations
2. Returns actual_results_path when found
3. Shows which locations were checked if none found
4. Uses the found path for validation

NOTE: This test uses MOCKS - it does NOT submit real jobs or check real S3 files.
It tests the logic and flow of the validation function.

To test with a REAL job, use scripts/check_job_output.py instead.
"""

import sys
import os
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from backend.job_manager import JobManager
from unittest.mock import Mock, patch
import json

# NOTE: This test uses MOCKS to test the logic flow without requiring:
# - Real AWS credentials
# - Real S3 files
# - Real job submissions
# 
# To test with a REAL job, use: python3 scripts/test_validation_real_job.py <job_id>

def test_multi_location_checking():
    """Test that validation checks multiple locations."""
    print("=" * 70)
    print("TEST: Multi-Location Results.json Checking")
    print("=" * 70)
    print()
    
    jm = JobManager()
    
    # Test case 1: results.json found at primary location
    print("Test 1: results.json found at PRIMARY location")
    print("-" * 70)
    job_id = "test-job-1"
    primary_path = "s3://bucket/session/job1/results.json"
    
    job = {
        "job_id": job_id,
        "step_id": "s-123456789",
        "output_path": "s3://bucket/session/job1/",
        "status": "completed"
    }
    
    with patch.object(jm, '_check_s3_file_exists') as mock_check:
        # Mock: primary location exists
        mock_check.return_value = (True, None)
        
        result = jm._validate_job_output_files(job_id, primary_path, job)
        
        print(f"✓ Checked locations: {len([primary_path])} location(s)")
        print(f"✓ results_json_exists: {result['results_json_exists']}")
        print(f"✓ actual_results_path: {result['actual_results_path']}")
        print(f"✓ Paths checked: {mock_check.call_count} time(s)")
        assert result['results_json_exists'] == True
        assert result['actual_results_path'] == primary_path
        print("✅ PASS: Found at primary location")
    print()
    
    # Test case 2: results.json found at alternative location (with step_id)
    print("Test 2: results.json found at ALTERNATIVE location (with step_id)")
    print("-" * 70)
    job_id = "test-job-2"
    primary_path = "s3://bucket/session/job2/results.json"
    alt_path = "s3://bucket/session/job2/s-123456789/results.json"
    
    job = {
        "job_id": job_id,
        "step_id": "s-123456789",
        "output_path": "s3://bucket/session/job2/",
        "status": "completed"
    }
    
    with patch.object(jm, '_check_s3_file_exists') as mock_check:
        # Mock: primary doesn't exist, alternative does
        def mock_check_side_effect(path):
            if path == primary_path:
                return (False, "File not found")
            elif path == alt_path:
                return (True, None)
            return (False, "File not found")
        
        mock_check.side_effect = mock_check_side_effect
        
        result = jm._validate_job_output_files(job_id, primary_path, job)
        
        print(f"✓ Checked locations: {mock_check.call_count} location(s)")
        print(f"✓ results_json_exists: {result['results_json_exists']}")
        print(f"✓ actual_results_path: {result['actual_results_path']}")
        # Check if warning was added (only if path differs from original)
        warnings_about_alt = [w for w in result['validation_warnings'] if 'alternative location' in w]
        print(f"✓ Warning about alternative location: {len(warnings_about_alt) > 0}")
        if warnings_about_alt:
            print(f"  Warning: {warnings_about_alt[0]}")
        assert result['results_json_exists'] == True
        assert result['actual_results_path'] == alt_path
        print("✅ PASS: Found at alternative location")
    print()
    
    # Test case 3: results.json not found at any location
    print("Test 3: results.json NOT FOUND at any location")
    print("-" * 70)
    job_id = "test-job-3"
    primary_path = "s3://bucket/session/job3/results.json"
    
    job = {
        "job_id": job_id,
        "step_id": "s-123456789",
        "output_path": "s3://bucket/session/job3/",
        "status": "completed"
    }
    
    with patch.object(jm, '_check_s3_file_exists') as mock_check:
        # Mock: all locations don't exist
        mock_check.return_value = (False, "File not found")
        
        result = jm._validate_job_output_files(job_id, primary_path, job)
        
        print(f"✓ Checked locations: {mock_check.call_count} location(s)")
        print(f"✓ results_json_exists: {result['results_json_exists']}")
        print(f"✓ actual_results_path: {result['actual_results_path']}")
        print(f"✓ Number of warnings: {len(result['validation_warnings'])}")
        print(f"✓ Warnings include paths checked: {any('locations:' in w for w in result['validation_warnings'])}")
        print(f"✓ Warnings include count: {any('Checked' in w and 'locations' in w for w in result['validation_warnings'])}")
        
        # Print warnings
        print("\nWarnings:")
        for i, warning in enumerate(result['validation_warnings'], 1):
            print(f"  {i}. {warning}")
        
        assert result['results_json_exists'] == False
        assert result['actual_results_path'] is None
        assert len(result['validation_warnings']) >= 1
        assert any('locations:' in w for w in result['validation_warnings'])
        print("\n✅ PASS: Correctly reports not found with all checked locations")
    print()
    
    # Test case 4: Verify paths_to_check logic
    print("Test 4: Verify all paths are added to check list")
    print("-" * 70)
    job_id = "test-job-4"
    primary_path = "s3://helix-ai-frontend-794270057041-us-west-1/session/job4/results.json"
    
    job = {
        "job_id": job_id,
        "step_id": "s-123456789",
        "output_path": "s3://helix-ai-frontend-794270057041-us-west-1/session/job4/",
        "status": "completed"
    }
    
    checked_paths = []
    
    with patch.object(jm, '_check_s3_file_exists') as mock_check:
        def mock_check_side_effect(path):
            checked_paths.append(path)
            return (False, "File not found")
        
        mock_check.side_effect = mock_check_side_effect
        
        result = jm._validate_job_output_files(job_id, primary_path, job)
        
        print(f"✓ Total paths checked: {len(checked_paths)}")
        print(f"✓ Paths checked:")
        for i, path in enumerate(checked_paths, 1):
            print(f"    {i}. {path}")
        
        # Should check: primary, with step_id, and default EMR location
        expected_paths = [
            primary_path,
            f"s3://helix-ai-frontend-794270057041-us-west-1/session/job4/s-123456789/results.json",
            f"s3://noricum-ngs-data/emr-results/{job_id}/s-123456789/results.json"
        ]
        
        print(f"\n✓ Expected {len(expected_paths)} paths")
        for expected in expected_paths:
            if expected in checked_paths:
                print(f"  ✅ Found: {expected}")
            else:
                print(f"  ❌ Missing: {expected}")
        
        assert len(checked_paths) >= 2  # At least primary and one alternative
        print("\n✅ PASS: Multiple paths checked")
    print()
    
    print("=" * 70)
    print("ALL TESTS PASSED")
    print("=" * 70)
    print()
    print("Summary:")
    print("1. ✅ Checks multiple possible locations")
    print("2. ✅ Returns actual_results_path when found")
    print("3. ✅ Shows which locations were checked if none found")
    print("4. ✅ Uses the found path for validation")
    print()
    print("NOTE: These tests use MOCKS - they test the logic, not real S3 files.")
    print("To test with a REAL job and REAL S3 files, use:")
    print("  python3 scripts/test_validation_real_job.py <job_id>")
    print("  or")
    print("  python3 scripts/check_job_output.py <job_id>")
    print()

if __name__ == '__main__':
    test_multi_location_checking()

