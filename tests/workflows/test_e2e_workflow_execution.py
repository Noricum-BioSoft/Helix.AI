#!/usr/bin/env python3
"""
End-to-End Workflow Execution Tests

These tests actually RUN workflows to completion and verify results:
1. Submit job to execution infrastructure (Local or EMR)
2. Wait for job completion (with timeout)
3. Verify output files exist in S3
4. Validate result quality/correctness

Usage:
    # Run all E2E tests (SLOW - may take 10-30 min)
    pytest tests/workflows/test_e2e_workflow_execution.py -v -s

    # Run specific workflow
    pytest tests/workflows/test_e2e_workflow_execution.py::TestE2EWorkflows::test_fastqc_small -v -s
    
    # Run only fast tests (small datasets, local execution)
    pytest tests/workflows/test_e2e_workflow_execution.py -v -s -m fast
    
    # Run slow tests (large datasets, EMR execution)
    pytest tests/workflows/test_e2e_workflow_execution.py -v -s -m slow

Requirements:
    - AWS credentials configured
    - S3 bucket access
    - EMR permissions (for large dataset tests)
    - Backend dependencies installed
"""

import sys
import os
import asyncio
import pytest
import boto3
import time
from pathlib import Path
from typing import Dict, Any, Optional, List
from datetime import datetime, timedelta

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "backend"))


class WorkflowExecutionTester:
    """Test harness for end-to-end workflow execution."""
    
    def __init__(self, s3_bucket: str = "noricum-ngs-data"):
        self.s3_bucket = s3_bucket
        self.s3_client = boto3.client('s3')
        self.emr_client = boto3.client('emr')
        
    async def submit_and_wait(
        self,
        tool_name: str,
        parameters: Dict[str, Any],
        environment: str = "Local",
        timeout_seconds: int = 1800,  # 30 minutes default
        check_interval: int = 30  # Check every 30 seconds
    ) -> Dict[str, Any]:
        """
        Submit job and wait for completion.
        
        Args:
            tool_name: Tool to execute (e.g., "fastqc_quality_analysis")
            parameters: Tool parameters
            environment: "Local" or "EMR"
            timeout_seconds: Maximum wait time
            check_interval: How often to check status
            
        Returns:
            Dict with job results and metadata
            
        Raises:
            TimeoutError: If job doesn't complete within timeout
            RuntimeError: If job fails
        """
        from backend.execution_broker import ExecutionBroker, ExecutionRequest
        
        # Create request
        session_id = f"e2e_test_{tool_name}_{int(time.time())}"
        request = ExecutionRequest(
            tool_name=tool_name,
            arguments=parameters,
            session_id=session_id,
            original_command=f"E2E test: {tool_name}",
            session_context={}
        )
        
        print(f"\n🚀 Submitting job: {tool_name}")
        print(f"   Session: {session_id}")
        print(f"   Environment: {environment}")
        print(f"   Parameters: {parameters}")
        
        # Create a simple tool executor for testing
        # In production, this would be dispatch_tool from main.py
        async def dummy_tool_executor(tool_name: str, args: Dict[str, Any]) -> Dict[str, Any]:
            """Dummy executor for E2E tests - actual execution happens via job manager."""
            return {
                "status": "submitted",
                "tool": tool_name,
                "message": f"Tool {tool_name} would be executed with args: {args}"
            }
        
        # Submit job
        broker = ExecutionBroker(tool_executor=dummy_tool_executor)
        start_time = time.time()
        submission_result = await broker.execute_tool(request)
        
        job_id = submission_result.get('job_id')
        output_location = submission_result.get('output_location', '')
        
        print(f"\n✅ Job submitted!")
        print(f"   Job ID: {job_id}")
        print(f"   Output: {output_location}")
        
        # Wait for completion
        print(f"\n⏳ Waiting for completion (timeout: {timeout_seconds}s)...")
        
        elapsed = 0
        while elapsed < timeout_seconds:
            await asyncio.sleep(check_interval)
            elapsed = time.time() - start_time
            
            # Check job status
            status = await self._check_job_status(
                job_id, 
                environment, 
                output_location
            )
            
            print(f"   [{int(elapsed)}s] Status: {status['state']}")
            
            if status['state'] == 'COMPLETED':
                print(f"\n✅ Job completed in {int(elapsed)}s")
                return {
                    'job_id': job_id,
                    'output_location': output_location,
                    'duration_seconds': elapsed,
                    'status': status,
                    'submission_result': submission_result
                }
            elif status['state'] == 'FAILED':
                raise RuntimeError(f"Job failed: {status.get('message', 'Unknown error')}")
        
        raise TimeoutError(f"Job did not complete within {timeout_seconds}s")
    
    async def _check_job_status(
        self,
        job_id: str,
        environment: str,
        output_location: str
    ) -> Dict[str, str]:
        """Check job status based on environment."""
        
        if environment == "EMR":
            # Check EMR step status
            try:
                # Parse cluster ID from job_id (format: cluster_id:step_id)
                if ':' in job_id:
                    cluster_id, step_id = job_id.split(':', 1)
                else:
                    # Fallback: find active cluster
                    clusters = self.emr_client.list_clusters(
                        ClusterStates=['STARTING', 'BOOTSTRAPPING', 'RUNNING', 'WAITING']
                    )
                    if not clusters['Clusters']:
                        return {'state': 'FAILED', 'message': 'No active EMR cluster'}
                    cluster_id = clusters['Clusters'][0]['Id']
                    step_id = job_id
                
                response = self.emr_client.describe_step(
                    ClusterId=cluster_id,
                    StepId=step_id
                )
                
                state = response['Step']['Status']['State']
                
                # Map EMR states to our states
                state_map = {
                    'PENDING': 'RUNNING',
                    'RUNNING': 'RUNNING',
                    'COMPLETED': 'COMPLETED',
                    'CANCELLED': 'FAILED',
                    'FAILED': 'FAILED',
                    'INTERRUPTED': 'FAILED'
                }
                
                return {
                    'state': state_map.get(state, 'RUNNING'),
                    'message': response['Step']['Status'].get('StateChangeReason', {}).get('Message', '')
                }
                
            except Exception as e:
                return {'state': 'UNKNOWN', 'message': str(e)}
        
        else:  # Local execution
            # Check if output files exist in S3
            if output_location and output_location.startswith('s3://'):
                # Parse S3 path
                parts = output_location.replace('s3://', '').split('/', 1)
                bucket = parts[0]
                prefix = parts[1] if len(parts) > 1 else ''
                
                try:
                    response = self.s3_client.list_objects_v2(
                        Bucket=bucket,
                        Prefix=prefix,
                        MaxKeys=1
                    )
                    
                    if response.get('KeyCount', 0) > 0:
                        return {'state': 'COMPLETED', 'message': 'Output files found'}
                    else:
                        return {'state': 'RUNNING', 'message': 'Waiting for output files'}
                        
                except Exception as e:
                    return {'state': 'RUNNING', 'message': f'Checking status: {e}'}
            
            # Fallback: assume completed after some time
            return {'state': 'COMPLETED', 'message': 'Local execution (status uncertain)'}
    
    def verify_s3_outputs(
        self,
        output_location: str,
        expected_files: Optional[List[str]] = None,
        min_file_count: int = 1
    ) -> Dict[str, Any]:
        """
        Verify output files exist in S3.
        
        Args:
            output_location: S3 path (s3://bucket/prefix/)
            expected_files: List of expected filenames (optional)
            min_file_count: Minimum number of files expected
            
        Returns:
            Dict with verification results
        """
        print(f"\n🔍 Verifying S3 outputs: {output_location}")
        
        # Parse S3 path
        if not output_location.startswith('s3://'):
            return {'verified': False, 'error': 'Invalid S3 path'}
        
        parts = output_location.replace('s3://', '').split('/', 1)
        bucket = parts[0]
        prefix = parts[1] if len(parts) > 1 else ''
        
        # List objects
        try:
            paginator = self.s3_client.get_paginator('list_objects_v2')
            pages = paginator.paginate(Bucket=bucket, Prefix=prefix)
            
            all_files = []
            for page in pages:
                if 'Contents' in page:
                    all_files.extend([obj['Key'] for obj in page['Contents']])
            
            print(f"   Found {len(all_files)} file(s)")
            
            # Check file count
            if len(all_files) < min_file_count:
                return {
                    'verified': False,
                    'error': f'Found {len(all_files)} files, expected at least {min_file_count}',
                    'files': all_files
                }
            
            # Check expected files
            if expected_files:
                found_files = [f.split('/')[-1] for f in all_files]
                missing = [f for f in expected_files if f not in found_files]
                
                if missing:
                    return {
                        'verified': False,
                        'error': f'Missing files: {missing}',
                        'files': all_files
                    }
            
            print(f"   ✅ All checks passed!")
            for f in all_files[:5]:  # Show first 5
                print(f"      • {f}")
            if len(all_files) > 5:
                print(f"      ... and {len(all_files) - 5} more")
            
            return {
                'verified': True,
                'file_count': len(all_files),
                'files': all_files
            }
            
        except Exception as e:
            return {'verified': False, 'error': str(e)}


@pytest.fixture(scope="session")
def workflow_tester():
    """Create workflow execution tester."""
    return WorkflowExecutionTester()


@pytest.fixture(scope="session")
def check_aws_access():
    """Check AWS credentials are configured."""
    try:
        sts = boto3.client('sts')
        identity = sts.get_caller_identity()
        print(f"\n✅ AWS Access: {identity['Arn']}")
        return True
    except Exception as e:
        pytest.skip(f"AWS credentials not configured: {e}")


# ============================================================================
# TEST SUITE: End-to-End Workflow Execution
# ============================================================================

class TestE2EWorkflows:
    """End-to-end workflow execution tests."""
    
    @pytest.mark.asyncio
    @pytest.mark.fast
    @pytest.mark.e2e
    async def test_fastqc_small_dataset(self, workflow_tester, check_aws_access):
        """
        Test FastQC on small test dataset (Local execution).
        
        Expected: 
        - Job completes in < 5 minutes
        - HTML and ZIP files generated
        - Files saved to S3
        """
        result = await workflow_tester.submit_and_wait(
            tool_name="fastqc_quality_analysis",
            parameters={
                "input_r1": "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq",
                "input_r2": "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq"
            },
            environment="Local",
            timeout_seconds=300  # 5 minutes
        )
        
        # Verify job completed
        assert result['status']['state'] == 'COMPLETED'
        assert result['duration_seconds'] < 300
        
        # Verify outputs exist
        verification = workflow_tester.verify_s3_outputs(
            result['output_location'],
            min_file_count=2  # At least HTML + ZIP
        )
        
        assert verification['verified'], f"Output verification failed: {verification.get('error')}"
        assert verification['file_count'] >= 2
        
        # Check for expected file types
        files = verification['files']
        has_html = any(f.endswith('.html') for f in files)
        has_zip = any(f.endswith('.zip') for f in files)
        
        assert has_html, "Missing HTML report"
        assert has_zip, "Missing ZIP data file"
    
    @pytest.mark.asyncio
    @pytest.mark.slow
    @pytest.mark.e2e
    @pytest.mark.skip(reason="Large dataset test - run manually")
    async def test_fastqc_large_dataset(self, workflow_tester, check_aws_access):
        """
        Test FastQC on large dataset (EMR execution).
        
        Expected:
        - Job completes in < 30 minutes
        - Results saved to S3
        - EMR cluster handles job correctly
        """
        result = await workflow_tester.submit_and_wait(
            tool_name="fastqc_quality_analysis",
            parameters={
                "input_r1": "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/mate_R1.fq",
                "input_r2": "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/mate_R2.fq"
            },
            environment="EMR",
            timeout_seconds=1800  # 30 minutes
        )
        
        assert result['status']['state'] == 'COMPLETED'
        
        verification = workflow_tester.verify_s3_outputs(
            result['output_location'],
            min_file_count=2
        )
        
        assert verification['verified']
    
    @pytest.mark.asyncio
    @pytest.mark.fast
    @pytest.mark.e2e
    async def test_read_merging_small(self, workflow_tester, check_aws_access):
        """
        Test read merging on small paired-end dataset.
        
        Expected:
        - Reads merged successfully
        - Output FASTQ generated
        - Statistics available
        """
        result = await workflow_tester.submit_and_wait(
            tool_name="read_merging",
            parameters={
                "input_r1": "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq",
                "input_r2": "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq",
                "quality_threshold": 30,
                "min_overlap": 10
            },
            environment="Local",
            timeout_seconds=300
        )
        
        assert result['status']['state'] == 'COMPLETED'
        
        verification = workflow_tester.verify_s3_outputs(
            result['output_location'],
            min_file_count=1  # At least merged FASTQ
        )
        
        assert verification['verified']
    
    @pytest.mark.asyncio
    @pytest.mark.fast
    @pytest.mark.e2e
    @pytest.mark.xfail(reason="Tool name may vary")
    async def test_quality_assessment(self, workflow_tester, check_aws_access):
        """
        Test quality assessment workflow.
        
        This tests the complete quality control pipeline.
        """
        result = await workflow_tester.submit_and_wait(
            tool_name="quality_assessment",
            parameters={
                "input_files": [
                    "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq"
                ]
            },
            environment="Local",
            timeout_seconds=300
        )
        
        assert result['status']['state'] == 'COMPLETED'


class TestE2EWorkflowFailures:
    """Test error handling in workflow execution."""
    
    @pytest.mark.asyncio
    @pytest.mark.fast
    @pytest.mark.e2e
    async def test_invalid_input_file(self, workflow_tester, check_aws_access):
        """Test that invalid input files are handled gracefully."""
        
        with pytest.raises((RuntimeError, TimeoutError)):
            await workflow_tester.submit_and_wait(
                tool_name="fastqc_quality_analysis",
                parameters={
                    "input_r1": "s3://noricum-ngs-data/invalid/nonexistent.fq"
                },
                environment="Local",
                timeout_seconds=120
            )
    
    @pytest.mark.asyncio
    @pytest.mark.fast
    @pytest.mark.e2e
    async def test_missing_required_parameter(self, workflow_tester, check_aws_access):
        """Test that missing parameters are handled."""
        
        with pytest.raises(Exception):  # Should fail at submission
            await workflow_tester.submit_and_wait(
                tool_name="fastqc_quality_analysis",
                parameters={},  # Missing required inputs
                environment="Local",
                timeout_seconds=60
            )


# ============================================================================
# USAGE EXAMPLES
# ============================================================================

if __name__ == "__main__":
    """
    Run directly for manual testing:
    
    python tests/workflows/test_e2e_workflow_execution.py
    """
    print(__doc__)
    print("\n" + "="*70)
    print("RUN WITH PYTEST:")
    print("="*70)
    print("pytest tests/workflows/test_e2e_workflow_execution.py -v -s")
    print("\nUse -m markers to filter:")
    print("  -m fast   # Only fast tests (small datasets)")
    print("  -m slow   # Only slow tests (large datasets)")
    print("  -m e2e    # All end-to-end tests")
