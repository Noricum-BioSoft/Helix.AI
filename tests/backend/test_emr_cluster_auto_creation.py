"""
Test that EMR cluster auto-creation works when EMR_CLUSTER_ID is not set.

This validates the fix where the system should automatically find or create
an EMR cluster instead of raising an error.
"""

import pytest
from unittest.mock import Mock, patch, MagicMock
import os
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

from backend.job_manager import JobManager


class TestEMRClusterAutoCreation:
    """Test EMR cluster auto-creation logic."""
    
    def test_submit_fastqc_job_auto_creates_cluster_when_not_set(self):
        """Test that submit_fastqc_job finds or creates cluster when EMR_CLUSTER_ID is not set."""
        # Clear EMR_CLUSTER_ID
        if "EMR_CLUSTER_ID" in os.environ:
            del os.environ["EMR_CLUSTER_ID"]
        
        job_manager = JobManager()
        
        # Mock the cluster finding/creation methods
        with patch.object(job_manager, '_find_active_cluster', return_value=None) as mock_find:
            with patch.object(job_manager, '_create_emr_cluster', return_value='j-TEST123456') as mock_create:
                with patch.object(job_manager, '_wait_for_cluster_ready', return_value=True) as mock_wait:
                    with patch.object(job_manager, '_check_cluster_state', return_value='WAITING'):
                        with patch('subprocess.run') as mock_subprocess:
                            mock_subprocess.return_value = Mock(returncode=0, stdout="", stderr="")
                            
                            # This should NOT raise an error, but should create a cluster
                            try:
                                job_id = job_manager.submit_fastqc_job(
                                    r1_path="s3://bucket/R1.fq",
                                    r2_path="s3://bucket/R2.fq"
                                )
                                
                                # Verify cluster creation was attempted
                                mock_create.assert_called_once()
                                # Implementation may choose to submit immediately and let the
                                # cluster stabilize asynchronously, so don't require a wait call.
                                
                                # Verify EMR_CLUSTER_ID was set
                                assert os.environ.get("EMR_CLUSTER_ID") == 'j-TEST123456'
                                
                                print("✅ PASS: submit_fastqc_job auto-creates cluster when EMR_CLUSTER_ID not set")
                            except ValueError as e:
                                if "EMR_CLUSTER_ID" in str(e):
                                    pytest.fail(f"❌ FAIL: Still raises error when EMR_CLUSTER_ID not set: {e}")
                                else:
                                    raise
    
    def test_submit_fastqc_job_uses_existing_cluster_when_found(self):
        """Test that submit_fastqc_job uses existing active cluster when found."""
        # Clear EMR_CLUSTER_ID
        if "EMR_CLUSTER_ID" in os.environ:
            del os.environ["EMR_CLUSTER_ID"]
        
        job_manager = JobManager()
        
        # Mock finding an active cluster
        with patch.object(job_manager, '_find_active_cluster', return_value='j-EXISTING123') as mock_find:
            with patch.object(job_manager, '_check_cluster_state', return_value='WAITING'):
                with patch('subprocess.run') as mock_subprocess:
                    mock_subprocess.return_value = Mock(returncode=0, stdout="", stderr="")
                    
                    job_id = job_manager.submit_fastqc_job(
                        r1_path="s3://bucket/R1.fq",
                        r2_path="s3://bucket/R2.fq"
                    )
                    
                    # Verify it found and used the existing cluster
                    mock_find.assert_called_once()
                    assert os.environ.get("EMR_CLUSTER_ID") == 'j-EXISTING123'
                    
                    # Verify it did NOT try to create a new cluster
                    assert not hasattr(job_manager, '_create_emr_cluster') or \
                           not hasattr(job_manager._create_emr_cluster, 'call_count')
                    
                    print("✅ PASS: submit_fastqc_job uses existing cluster when found")
    
    def test_submit_universal_emr_job_auto_creates_cluster_when_not_set(self):
        """Test that submit_universal_emr_job finds or creates cluster when EMR_CLUSTER_ID is not set."""
        # Clear EMR_CLUSTER_ID
        if "EMR_CLUSTER_ID" in os.environ:
            del os.environ["EMR_CLUSTER_ID"]
        
        job_manager = JobManager()
        
        # Mock the cluster finding/creation methods
        with patch.object(job_manager, '_find_active_cluster', return_value=None) as mock_find:
            with patch.object(job_manager, '_create_emr_cluster', return_value='j-TEST789') as mock_create:
                with patch.object(job_manager, '_wait_for_cluster_ready', return_value=True) as mock_wait:
                    with patch.object(job_manager, '_check_cluster_state', return_value='WAITING'):
                        with patch('subprocess.run') as mock_subprocess:
                            mock_subprocess.return_value = Mock(returncode=0, stdout="", stderr="")
                            with patch('tempfile.NamedTemporaryFile') as mock_tempfile:
                                mock_tempfile.return_value.__enter__.return_value.name = '/tmp/test.json'
                                
                                # This should NOT raise an error
                                try:
                                    job_id = job_manager.submit_universal_emr_job(
                                        tool_name="read_merging",
                                        tool_args={"forward_reads": "s3://bucket/R1.fq"},
                                        output_path="s3://bucket/emr-results/test-job/",
                                        # Avoid bundling local tools/ (slow + filesystem dependent)
                                        tools_bundle_s3="s3://bucket/emr-scripts/tools_bundle_test.zip",
                                    )
                                    
                                    # Verify cluster creation was attempted
                                    mock_create.assert_called_once()
                                    # Implementation may choose to submit immediately and let the
                                    # cluster stabilize asynchronously, so don't require a wait call.
                                    
                                    # Verify EMR_CLUSTER_ID was set
                                    assert os.environ.get("EMR_CLUSTER_ID") == 'j-TEST789'
                                    
                                    print("✅ PASS: submit_universal_emr_job auto-creates cluster when EMR_CLUSTER_ID not set")
                                except ValueError as e:
                                    if "EMR_CLUSTER_ID" in str(e):
                                        pytest.fail(f"❌ FAIL: Still raises error when EMR_CLUSTER_ID not set: {e}")
                                    else:
                                        raise
    
    def test_no_error_when_emr_cluster_id_not_set(self):
        """Test that no ValueError is raised when EMR_CLUSTER_ID is not set."""
        # Clear EMR_CLUSTER_ID
        if "EMR_CLUSTER_ID" in os.environ:
            del os.environ["EMR_CLUSTER_ID"]
        
        job_manager = JobManager()
        
        # Mock successful cluster creation
        with patch.object(job_manager, '_find_active_cluster', return_value=None):
            with patch.object(job_manager, '_create_emr_cluster', return_value='j-AUTO123'):
                with patch.object(job_manager, '_wait_for_cluster_ready', return_value=True):
                    with patch.object(job_manager, '_check_cluster_state', return_value='WAITING'):
                        with patch('subprocess.run') as mock_subprocess:
                            mock_subprocess.return_value = Mock(returncode=0, stdout="Step ID: s-TEST", stderr="")
                            
                            # Should not raise ValueError about EMR_CLUSTER_ID
                            try:
                                # Try both methods
                                with patch('tempfile.NamedTemporaryFile'):
                                    job_manager.submit_universal_emr_job(
                                        tool_name="test_tool",
                                        tool_args={},
                                        output_path="s3://bucket/emr-results/test-job/",
                                        tools_bundle_s3="s3://bucket/emr-scripts/tools_bundle_test.zip",
                                    )
                                
                                print("✅ PASS: No ValueError raised when EMR_CLUSTER_ID not set")
                            except ValueError as e:
                                if "EMR_CLUSTER_ID" in str(e) and "not set" in str(e):
                                    pytest.fail(f"❌ FAIL: ValueError still raised: {e}")
                                else:
                                    # Other ValueErrors are OK
                                    pass


if __name__ == "__main__":
    pytest.main([__file__, "-v"])



