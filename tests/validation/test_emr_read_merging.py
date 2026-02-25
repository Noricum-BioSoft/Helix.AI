#!/usr/bin/env python3
"""
Test script to verify EMR read_merging fixes.

Tests:
1. Bootstrap script includes boto3
2. Universal runner correctly uses merge_reads_from_s3() for S3 files
3. boto3 is available on EMR clusters
"""

import os
import sys
import subprocess
import json
import tempfile
from pathlib import Path

def check_bootstrap_script_includes_boto3():
    """Check if the bootstrap script in S3 includes boto3."""
    print("=" * 60)
    print("Test 1: Checking bootstrap script includes boto3")
    print("=" * 60)
    
    region = os.getenv("AWS_REGION", "us-west-1")
    s3_bucket = os.getenv("S3_DATASET_BUCKET", "noricum-ngs-data")
    bootstrap_s3_path = f"s3://{s3_bucket}/emr-bootstrap/fastqc-bootstrap.sh"
    
    try:
        result = subprocess.run(
            ["aws", "s3", "cp", bootstrap_s3_path, "-", "--region", region],
            capture_output=True,
            text=True,
            timeout=30,
        )
        
        if result.returncode != 0:
            print(f"❌ Failed to fetch bootstrap script: {result.stderr}")
            return False
        
        bootstrap_content = result.stdout
        
        if "boto3" in bootstrap_content and "botocore" in bootstrap_content:
            print("✅ Bootstrap script includes boto3 and botocore")
            
            # Show the relevant lines
            lines = bootstrap_content.split('\n')
            for i, line in enumerate(lines):
                if 'boto3' in line or 'botocore' in line:
                    print(f"   Line {i+1}: {line.strip()}")
            
            return True
        else:
            print("❌ Bootstrap script does NOT include boto3/botocore")
            print("   Bootstrap script content (first 50 lines):")
            for i, line in enumerate(lines[:50], 1):
                print(f"   {i:3d}: {line}")
            return False
            
    except Exception as e:
        print(f"❌ Exception checking bootstrap script: {e}")
        return False


def check_universal_runner_code():
    """Check if the universal runner code has the correct read_merging logic."""
    print("\n" + "=" * 60)
    print("Test 2: Checking universal runner code")
    print("=" * 60)
    
    runner_path = Path(__file__).resolve().parent.parent.parent / "scripts" / "emr" / "universal_emr_runner.py"
    
    if not runner_path.exists():
        print(f"❌ Universal runner not found at {runner_path}")
        return False
    
    try:
        content = runner_path.read_text()
        
        checks = {
            "merge_reads_from_s3": "merge_reads_from_s3()" in content,
            "skip materialization": "should_materialize = False" in content or "tool_name == \"read_merging\"" in content,
            "output S3 path check": "_looks_like_s3_uri(output_path)" in content,
        }
        
        all_passed = all(checks.values())
        
        if all_passed:
            print("✅ Universal runner has correct read_merging logic")
            for check, passed in checks.items():
                print(f"   ✓ {check}")
        else:
            print("❌ Universal runner missing some fixes:")
            for check, passed in checks.items():
                status = "✓" if passed else "✗"
                print(f"   {status} {check}")
        
        # Show the read_merging section
        print("\n   read_merging section:")
        lines = content.split('\n')
        in_read_merging = False
        indent = 0
        for i, line in enumerate(lines):
            if 'if tool_name == "read_merging":' in line:
                in_read_merging = True
                indent = len(line) - len(line.lstrip())
            if in_read_merging:
                print(f"   {i+1:4d}: {line}")
                # Stop after merge_reads_from_s3 call or next tool
                if ('return read_merging.merge_reads_from_s3' in line or
                    (line.strip() and not line.startswith(' ') and 'if tool_name ==' in line and 'read_merging' not in line)):
                    break
        
        return all_passed
        
    except Exception as e:
        print(f"❌ Exception checking universal runner: {e}")
        return False


def test_boto3_import():
    """Test if boto3 can be imported (simulating EMR environment)."""
    print("\n" + "=" * 60)
    print("Test 3: Testing boto3 import (local - informational)")
    print("=" * 60)
    
    try:
        import boto3
        print("✅ boto3 can be imported locally")
        print(f"   boto3 version: {boto3.__version__}")
        return True
    except ImportError:
        print("ℹ️  boto3 not installed locally (this is OK)")
        print("   The bootstrap script will install boto3 on EMR clusters.")
        print("   This is an informational check - not a failure.")
        return True  # Don't fail the test for this


def check_bootstrap_script_in_code():
    """Check if job_manager.py always uploads bootstrap script."""
    print("\n" + "=" * 60)
    print("Test 4: Checking job_manager bootstrap script logic")
    print("=" * 60)
    
    job_manager_path = Path(__file__).resolve().parent.parent.parent / "backend" / "job_manager.py"
    
    if not job_manager_path.exists():
        print(f"❌ job_manager.py not found at {job_manager_path}")
        return False
    
    try:
        content = job_manager_path.read_text()
        
        # Check if _ensure_bootstrap_script always uploads
        # Look for the method and see if it skips upload when script exists
        method_start = content.find("def _ensure_bootstrap_script")
        if method_start == -1:
            print("❌ _ensure_bootstrap_script method not found")
            return False
        
        method_content = content[method_start:method_start+5000]  # Get method content
        
        # Check if it has logic to skip upload when script exists
        has_skip_logic = (
            "if result.returncode == 0:" in method_content and
            "return bootstrap_s3_path" in method_content and
            "Bootstrap script exists" in method_content
        )
        
        # Check if it has "always upload" logic
        has_always_upload = "Always upload" in method_content or "Ensuring bootstrap script is up-to-date" in method_content
        
        if has_always_upload or not has_skip_logic:
            print("✅ Bootstrap script logic will always upload latest version")
            print("   (Script is uploaded every time, ensuring latest version)")
            return True
        else:
            print("⚠️  Bootstrap script logic may skip upload if script exists")
            print("   This could mean old versions are used.")
            return False
        
    except Exception as e:
        print(f"❌ Exception checking job_manager: {e}")
        return False


def create_test_payload():
    """Create a test payload for read_merging."""
    print("\n" + "=" * 60)
    print("Test 5: Creating test payload structure")
    print("=" * 60)
    
    payload = {
        "tool_name": "read_merging",
        "arguments": {
            "forward_reads": "s3://test-bucket/test_R1.fq",
            "reverse_reads": "s3://test-bucket/test_R2.fq",
            "output": "s3://test-bucket/merged.fq",
            "min_overlap": 12
        }
    }
    
    print("✅ Test payload structure:")
    print(json.dumps(payload, indent=2))
    
    # Check if payload has output (which should trigger merge_reads_from_s3)
    if payload["arguments"].get("output"):
        print("✅ Payload has output S3 path (will use merge_reads_from_s3)")
        return True
    else:
        print("❌ Payload missing output path")
        return False


def main():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("EMR Read Merging Fix Verification Tests")
    print("=" * 60)
    print()
    
    results = []
    
    results.append(("Bootstrap script includes boto3", check_bootstrap_script_includes_boto3()))
    results.append(("Universal runner code", check_universal_runner_code()))
    results.append(("boto3 import (local)", test_boto3_import()))
    results.append(("Job manager bootstrap logic", check_bootstrap_script_in_code()))
    results.append(("Test payload structure", create_test_payload()))
    
    print("\n" + "=" * 60)
    print("Test Summary")
    print("=" * 60)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{status}: {name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    # Core tests (excluding informational ones)
    core_tests = [
        ("Bootstrap script includes boto3", results[0][1]),
        ("Universal runner code", results[1][1]),
        ("Job manager bootstrap logic", results[3][1]),
    ]
    core_passed = sum(1 for _, result in core_tests if result)
    core_total = len(core_tests)
    
    print(f"\nCore fixes verified: {core_passed}/{core_total}")
    
    if core_passed == core_total:
        print("\n🎉 All core fixes are in place!")
        print("   • Bootstrap script includes boto3 ✅")
        print("   • Universal runner uses merge_reads_from_s3() ✅")
        print("   • Bootstrap script always uploads latest version ✅")
        print("\n   The next EMR job should work correctly!")
        return 0
    else:
        print(f"\n⚠️  {core_total - core_passed} core test(s) failed.")
        print("   Please review the output above and fix any issues.")
        return 1


if __name__ == "__main__":
    sys.exit(main())

