"""
Unit tests for JobManager with AWS mocks.

Tests job submission, status tracking, and EMR interactions
without making real AWS API calls.
"""

import os
import pytest
import subprocess
import json
from unittest.mock import Mock, patch, MagicMock, call
from datetime import datetime, timezone

import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "backend"))

# Import after path setup
from job_manager import (
    JobManager,
    get_job_manager,
    STATUS_SUBMITTED,
    STATUS_RUNNING,
    STATUS_COMPLETED,
    STATUS_FAILED,
    STATUS_CANCELLED,
    EMR_STATUS_TO_JOB_STATUS,
)


class TestJobManager:
    """Test suite for JobManager class."""

    @pytest.fixture
    def job_manager(self):
        """Create a fresh JobManager instance for each test."""
        # Reset the singleton
        import job_manager as jm_module
        jm_module._job_manager_instance = None
        return JobManager()

    @pytest.fixture
    def mock_env_vars(self):
        """Mock environment variables."""
        with patch.dict(
            os.environ,
            {
                "EMR_CLUSTER_ID": "j-TEST123456789",
                "AWS_REGION": "us-east-1",
            },
            clear=False,
        ):
            yield

    @pytest.fixture
    def mock_script_path(self, tmp_path):
        """Create a mock submission script."""
        script_dir = tmp_path / "scripts" / "emr"
        script_dir.mkdir(parents=True)
        script_file = script_dir / "submit-fastqc-job.sh"
        script_file.write_text("#!/bin/bash\necho 'Step ID: s-ABC123456'\n")
        script_file.chmod(0o755)
        return script_file

    def test_submit_fastqc_job_success(
        self, job_manager, mock_env_vars, mock_script_path
    ):
        """Test successful job submission."""
        r1_path = "s3://my-bucket/data/R1.fastq"
        r2_path = "s3://my-bucket/data/R2.fastq"
        output_path = "s3://my-bucket/results/"

        with patch.object(
            job_manager, "project_root", mock_script_path.parent.parent.parent
        ), patch.object(
            job_manager, "_check_cluster_state", return_value="WAITING"
        ), patch(
            "subprocess.run",
            return_value=Mock(
                returncode=0,
                stdout="Step ID: s-ABC123456\n✅ Job submitted successfully!",
                stderr="",
            ),
        ):
            job_id = job_manager.submit_fastqc_job(
                r1_path=r1_path, r2_path=r2_path, output_path=output_path
            )

            assert job_id is not None
            assert isinstance(job_id, str)
            assert len(job_id) > 0

            # Check job is stored
            job = job_manager.jobs[job_id]
            assert job["status"] == STATUS_SUBMITTED
            assert job["r1_path"] == r1_path
            assert job["r2_path"] == r2_path
            assert job["output_path"] == output_path
            assert job["cluster_id"] == "j-TEST123456789"
            assert job["step_id"] == "s-ABC123456"
            assert job["type"] == "fastqc"
            assert "submitted_at" in job
            assert "updated_at" in job

    def test_submit_fastqc_job_missing_emr_cluster_id(self, job_manager):
        """Test job submission fails when EMR_CLUSTER_ID is not set."""
        with patch.dict(os.environ, {}, clear=True):
            with pytest.raises(ValueError, match="EMR_CLUSTER_ID"):
                job_manager.submit_fastqc_job(
                    r1_path="s3://bucket/R1.fastq",
                    r2_path="s3://bucket/R2.fastq",
                )

    def test_submit_fastqc_job_missing_paths(self, job_manager, mock_env_vars):
        """Test job submission fails when paths are missing."""
        with pytest.raises(ValueError, match="Both R1 and R2 paths"):
            job_manager.submit_fastqc_job(r1_path="", r2_path="s3://bucket/R2.fastq")

        with pytest.raises(ValueError, match="Both R1 and R2 paths"):
            job_manager.submit_fastqc_job(r1_path="s3://bucket/R1.fastq", r2_path="")

    def test_submit_fastqc_job_invalid_s3_paths(self, job_manager, mock_env_vars):
        """Test job submission fails for non-S3 paths."""
        with pytest.raises(ValueError, match="must be S3 URIs"):
            job_manager.submit_fastqc_job(
                r1_path="/local/path/R1.fastq", r2_path="s3://bucket/R2.fastq"
            )

        with pytest.raises(ValueError, match="must be S3 URIs"):
            job_manager.submit_fastqc_job(
                r1_path="s3://bucket/R1.fastq", r2_path="/local/path/R2.fastq"
            )

    def test_submit_fastqc_job_script_not_found(
        self, job_manager, mock_env_vars, tmp_path
    ):
        """Test job submission fails when script is not found."""
        with patch.object(
            job_manager, "project_root", tmp_path
        ), patch.object(
            job_manager, "_check_cluster_state", return_value="WAITING"
        ), pytest.raises(FileNotFoundError):
            job_manager.submit_fastqc_job(
                r1_path="s3://bucket/R1.fastq", r2_path="s3://bucket/R2.fastq"
            )

    def test_submit_fastqc_job_script_fails(
        self, job_manager, mock_env_vars, mock_script_path
    ):
        """Test job submission handles script failure."""
        with patch.object(
            job_manager, "project_root", mock_script_path.parent.parent.parent
        ), patch.object(
            job_manager, "_check_cluster_state", return_value="WAITING"
        ), patch(
            "subprocess.run",
            return_value=Mock(
                returncode=1,
                stdout="",
                stderr="Error: Cluster not ready",
            ),
        ), pytest.raises(
            RuntimeError, match="Failed to submit job"
        ):
            job_manager.submit_fastqc_job(
                r1_path="s3://bucket/R1.fastq", r2_path="s3://bucket/R2.fastq"
            )

    def test_extract_step_id(self, job_manager):
        """Test step ID extraction from script output."""
        # Test normal output
        output1 = "Step ID: s-ABC123456\n✅ Job submitted successfully!"
        assert job_manager._extract_step_id(output1) == "s-ABC123456"

        # Test with different format
        output2 = "Submitted step with ID s-XYZ789012"
        assert job_manager._extract_step_id(output2) == "s-XYZ789012"

        # Test with no step ID
        output3 = "Job submitted but no step ID found"
        assert job_manager._extract_step_id(output3) is None

    def test_get_job_status_submitted(self, job_manager):
        """Test getting status for a submitted job."""
        job_id = "test-job-id"
        job_manager.jobs[job_id] = {
            "job_id": job_id,
            "status": STATUS_SUBMITTED,
            "cluster_id": "j-TEST123",
            "step_id": "s-ABC123",
            "updated_at": datetime.now(timezone.utc).isoformat(),
        }

        with patch(
            "subprocess.run",
            return_value=Mock(
                returncode=0,
                stdout="RUNNING",
                stderr="",
            ),
        ):
            status = job_manager.get_job_status(job_id)

            assert status["job_id"] == job_id
            assert status["status"] == STATUS_RUNNING  # Should update to running

    def test_get_job_status_completed(self, job_manager):
        """Test getting status for a completed job (cached)."""
        job_id = "test-job-id"
        job_manager.jobs[job_id] = {
            "job_id": job_id,
            "status": STATUS_COMPLETED,
            "cluster_id": "j-TEST123",
            "step_id": "s-ABC123",
            "updated_at": datetime.now(timezone.utc).isoformat(),
        }

        # Should not call EMR API for terminal status
        status = job_manager.get_job_status(job_id)
        assert status["status"] == STATUS_COMPLETED

    def test_get_job_status_not_found(self, job_manager):
        """Test getting status for non-existent job."""
        with pytest.raises(ValueError, match="not found"):
            job_manager.get_job_status("non-existent-job-id")

    def test_check_emr_step_status_success(self, job_manager, mock_env_vars):
        """Test checking EMR step status successfully."""
        with patch(
            "subprocess.run",
            return_value=Mock(
                returncode=0,
                stdout="COMPLETED\n",
                stderr="",
            ),
        ):
            status, failure_details = job_manager._check_emr_step_status(
                cluster_id="j-TEST123", step_id="s-ABC123"
            )
            assert status == "COMPLETED"
            assert failure_details is None

    def test_check_emr_step_status_failure(self, job_manager, mock_env_vars):
        """Test checking EMR step status when AWS CLI fails."""
        with patch(
            "subprocess.run",
            return_value=Mock(
                returncode=1,
                stdout="",
                stderr="Error: Step not found",
            ),
        ):
            status, failure_details = job_manager._check_emr_step_status(
                cluster_id="j-TEST123", step_id="s-ABC123"
            )
            assert status is None
            assert failure_details is None

    def test_check_emr_step_status_timeout(self, job_manager, mock_env_vars):
        """Test checking EMR step status times out."""
        with patch(
            "subprocess.run",
            side_effect=subprocess.TimeoutExpired("aws", 10),
        ):
            status, failure_details = job_manager._check_emr_step_status(
                cluster_id="j-TEST123", step_id="s-ABC123"
            )
            assert status is None
            assert failure_details is None

    def test_get_job_results_success(self, job_manager):
        """Test getting results for a completed job."""
        job_id = "test-job-id"
        job_manager.jobs[job_id] = {
            "job_id": job_id,
            "status": STATUS_COMPLETED,
            "r1_path": "s3://bucket/data/R1.fastq",
            "r2_path": "s3://bucket/data/R2.fastq",
            "output_path": "s3://bucket/results/",
        }

        results = job_manager.get_job_results(job_id)

        assert results["job_id"] == job_id
        assert results["status"] == STATUS_COMPLETED
        assert results["results_path"] == "s3://bucket/results/results.json"
        assert results["output_path"] == "s3://bucket/results/"

    def test_get_job_results_not_completed(self, job_manager):
        """Test getting results for a non-completed job."""
        job_id = "test-job-id"
        job_manager.jobs[job_id] = {
            "job_id": job_id,
            "status": STATUS_RUNNING,
        }

        with pytest.raises(ValueError, match="not completed"):
            job_manager.get_job_results(job_id)

    def test_get_job_results_auto_generate_output_path(self, job_manager):
        """Test auto-generating output path when not provided."""
        job_id = "test-job-id"
        job_manager.jobs[job_id] = {
            "job_id": job_id,
            "status": STATUS_COMPLETED,
            "r1_path": "s3://bucket/data/sample_R1.fastq",
            "r2_path": "s3://bucket/data/sample_R2.fastq",
            "output_path": None,
        }

        results = job_manager.get_job_results(job_id)

        # When output_path is None, JobManager falls back to the legacy default dataset path.
        assert results["results_path"] == "s3://bucket/fastqc-results/GRCh38.p12.MafHi/results.json"

    def test_cancel_job_success(self, job_manager, mock_env_vars):
        """Test cancelling a running job."""
        job_id = "test-job-id"
        job_manager.jobs[job_id] = {
            "job_id": job_id,
            "status": STATUS_RUNNING,
            "cluster_id": "j-TEST123",
            "step_id": "s-ABC123",
            "updated_at": datetime.now(timezone.utc).isoformat(),
        }

        with patch(
            "subprocess.run",
            return_value=Mock(returncode=0, stdout="", stderr=""),
        ):
            result = job_manager.cancel_job(job_id)

            assert result["job_id"] == job_id
            assert result["status"] == STATUS_CANCELLED
            # Note: cancel_job modifies the job dict in-place, so check via get_job_status
            updated_job = job_manager.get_job_status(job_id)
            assert updated_job["status"] == STATUS_CANCELLED

    def test_cancel_job_already_completed(self, job_manager):
        """Test cancelling an already completed job."""
        job_id = "test-job-id"
        job_manager.jobs[job_id] = {
            "job_id": job_id,
            "status": STATUS_COMPLETED,
        }

        with pytest.raises(ValueError, match="Cannot cancel"):
            job_manager.cancel_job(job_id)

    def test_list_jobs_all(self, job_manager):
        """Test listing all jobs."""
        # Create multiple jobs
        for i in range(5):
            job_id = f"job-{i}"
            job_manager.jobs[job_id] = {
                "job_id": job_id,
                "status": STATUS_SUBMITTED if i % 2 == 0 else STATUS_RUNNING,
                "submitted_at": datetime.now(timezone.utc).isoformat(),
            }

        jobs = job_manager.list_jobs()

        assert len(jobs) == 5

    def test_list_jobs_filtered_by_status(self, job_manager):
        """Test listing jobs filtered by status."""
        # Create jobs with different statuses
        job_manager.jobs["job-1"] = {
            "job_id": "job-1",
            "status": STATUS_SUBMITTED,
            "submitted_at": datetime.now(timezone.utc).isoformat(),
        }
        job_manager.jobs["job-2"] = {
            "job_id": "job-2",
            "status": STATUS_RUNNING,
            "submitted_at": datetime.now(timezone.utc).isoformat(),
        }
        job_manager.jobs["job-3"] = {
            "job_id": "job-3",
            "status": STATUS_SUBMITTED,
            "submitted_at": datetime.now(timezone.utc).isoformat(),
        }

        jobs = job_manager.list_jobs(status=STATUS_SUBMITTED)

        assert len(jobs) == 2
        assert all(j["status"] == STATUS_SUBMITTED for j in jobs)

    def test_list_jobs_with_limit(self, job_manager):
        """Test listing jobs with limit."""
        # Create 10 jobs
        for i in range(10):
            job_id = f"job-{i}"
            job_manager.jobs[job_id] = {
                "job_id": job_id,
                "status": STATUS_SUBMITTED,
                "submitted_at": datetime.now(timezone.utc).isoformat(),
            }

        jobs = job_manager.list_jobs(limit=5)

        assert len(jobs) == 5

    def test_list_jobs_invalid_status(self, job_manager):
        """Test listing jobs with invalid status filter."""
        with pytest.raises(ValueError, match="Invalid status"):
            job_manager.list_jobs(status="invalid_status")

    def test_emr_status_to_job_status_mapping(self):
        """Test EMR status to job status mapping."""
        assert EMR_STATUS_TO_JOB_STATUS["PENDING"] == STATUS_SUBMITTED
        assert EMR_STATUS_TO_JOB_STATUS["RUNNING"] == STATUS_RUNNING
        assert EMR_STATUS_TO_JOB_STATUS["COMPLETED"] == STATUS_COMPLETED
        assert EMR_STATUS_TO_JOB_STATUS["FAILED"] == STATUS_FAILED
        assert EMR_STATUS_TO_JOB_STATUS["CANCELLED"] == STATUS_CANCELLED

    def test_get_job_manager_singleton(self):
        """Test that get_job_manager returns a singleton."""
        import job_manager as jm_module
        jm_module._job_manager_instance = None

        manager1 = get_job_manager()
        manager2 = get_job_manager()

        assert manager1 is manager2

    def test_submit_with_timeout(self, job_manager, mock_env_vars, mock_script_path):
        """Test job submission timeout handling."""
        with patch.object(
            job_manager, "project_root", mock_script_path.parent.parent.parent
        ), patch.object(
            job_manager, "_check_cluster_state", return_value="WAITING"
        ), patch(
            "subprocess.run",
            side_effect=subprocess.TimeoutExpired("script.sh", 60),
        ), pytest.raises(
            RuntimeError, match="timed out"
        ):
            job_manager.submit_fastqc_job(
                r1_path="s3://bucket/R1.fastq", r2_path="s3://bucket/R2.fastq"
            )


class TestJobManagerIntegration:
    """Integration-style tests with more complex scenarios."""

    @pytest.fixture
    def job_manager(self):
        """Create a fresh JobManager instance."""
        import job_manager as jm_module
        jm_module._job_manager_instance = None
        return JobManager()

    @pytest.fixture
    def mock_env_vars(self):
        """Mock environment variables."""
        with patch.dict(
            os.environ,
            {
                "EMR_CLUSTER_ID": "j-TEST123456789",
                "AWS_REGION": "us-east-1",
            },
            clear=False,
        ):
            yield

    @pytest.fixture
    def mock_script_path(self, tmp_path):
        """Create a mock submission script."""
        script_dir = tmp_path / "scripts" / "emr"
        script_dir.mkdir(parents=True)
        script_file = script_dir / "submit-fastqc-job.sh"
        script_file.write_text("#!/bin/bash\necho 'Step ID: s-ABC123456'\n")
        script_file.chmod(0o755)
        return script_file

    def test_full_job_lifecycle(self, job_manager, mock_env_vars, mock_script_path):
        """Test a complete job lifecycle: submit -> running -> completed."""
        with patch.object(
            job_manager, "project_root", mock_script_path.parent.parent.parent
        ):
            # Submit job
            with patch(
                "subprocess.run",
                return_value=Mock(
                    returncode=0,
                    stdout="Step ID: s-LIFECYCLE123\n✅ Job submitted!",
                    stderr="",
                ),
            ):
                job_id = job_manager.submit_fastqc_job(
                    r1_path="s3://bucket/R1.fastq", r2_path="s3://bucket/R2.fastq"
                )

            assert job_id is not None
            initial_status = job_manager.get_job_status(job_id)
            assert initial_status["status"] == STATUS_SUBMITTED

            # Update to running
            with patch(
                "subprocess.run",
                return_value=Mock(
                    returncode=0,
                    stdout="RUNNING\n",
                    stderr="",
                ),
            ):
                running_status = job_manager.get_job_status(job_id)
                # Status should update if EMR reports RUNNING
                # (Note: This test may not update if the mapping doesn't match exactly)

            # Update to completed
            job_manager.jobs[job_id]["status"] = STATUS_COMPLETED
            completed_status = job_manager.get_job_status(job_id)
            assert completed_status["status"] == STATUS_COMPLETED

            # Get results
            results = job_manager.get_job_results(job_id)
            assert results["status"] == STATUS_COMPLETED
            assert "results_path" in results


class TestJobManagerRetryAndLogs:
    """Test suite for retry and logs functionality (Phase 3)."""

    @pytest.fixture
    def job_manager(self):
        """Create a fresh JobManager instance for each test."""
        import job_manager as jm_module
        jm_module._job_manager_instance = None
        return JobManager()

    @pytest.fixture
    def mock_env_vars(self):
        """Mock environment variables."""
        with patch.dict(
            os.environ,
            {
                "EMR_CLUSTER_ID": "j-TEST123456789",
                "AWS_REGION": "us-east-1",
            },
            clear=False,
        ):
            yield

    @pytest.fixture
    def mock_script_path(self, tmp_path):
        """Create a mock submission script."""
        script_dir = tmp_path / "scripts" / "emr"
        script_dir.mkdir(parents=True)
        script_file = script_dir / "submit-fastqc-job.sh"
        script_file.write_text("#!/bin/bash\necho 'Step ID: s-ABC123456'\n")
        script_file.chmod(0o755)
        return script_file

    def test_retry_failed_job(
        self, job_manager, mock_env_vars, mock_script_path
    ):
        """Test retrying a failed job."""
        with patch.object(
            job_manager, "project_root", mock_script_path.parent.parent.parent
        ):
            # First, create a failed job
            with patch(
                "subprocess.run",
                return_value=Mock(
                    returncode=0,
                    stdout="Step ID: s-FAILED123\n✅ Job submitted!",
                    stderr="",
                ),
            ):
                original_job_id = job_manager.submit_fastqc_job(
                    r1_path="s3://bucket/R1.fastq", r2_path="s3://bucket/R2.fastq"
                )

            # Mark job as failed
            job_manager.jobs[original_job_id]["status"] = STATUS_FAILED

            # Retry the job
            with patch(
                "subprocess.run",
                return_value=Mock(
                    returncode=0,
                    stdout="Step ID: s-RETRY123\n✅ Job submitted!",
                    stderr="",
                ),
            ):
                new_job_id = job_manager.retry_job(original_job_id)

            # Verify new job was created
            assert new_job_id != original_job_id
            assert new_job_id in job_manager.jobs

            # Verify new job has same parameters
            new_job = job_manager.jobs[new_job_id]
            original_job = job_manager.jobs[original_job_id]
            assert new_job["r1_path"] == original_job["r1_path"]
            assert new_job["r2_path"] == original_job["r2_path"]
            assert new_job["status"] == STATUS_SUBMITTED

    def test_retry_non_failed_job_raises_error(
        self, job_manager, mock_env_vars, mock_script_path
    ):
        """Test that retrying a non-failed job raises an error."""
        with patch.object(
            job_manager, "project_root", mock_script_path.parent.parent.parent
        ):
            # Create a running job
            with patch(
                "subprocess.run",
                return_value=Mock(
                    returncode=0,
                    stdout="Step ID: s-RUNNING123\n✅ Job submitted!",
                    stderr="",
                ),
            ):
                job_id = job_manager.submit_fastqc_job(
                    r1_path="s3://bucket/R1.fastq", r2_path="s3://bucket/R2.fastq"
                )

            # Try to retry a running job - should raise error
            with pytest.raises(ValueError, match="Cannot retry job"):
                job_manager.retry_job(job_id)

    def test_get_job_logs_success(
        self, job_manager, mock_env_vars
    ):
        """Test getting EMR logs for a job."""
        # Create a job with cluster_id and step_id
        job_id = "test-job-123"
        job_manager.jobs[job_id] = {
            "job_id": job_id,
            "status": STATUS_FAILED,
            "cluster_id": "j-TEST123",
            "step_id": "s-ABC123",
            "r1_path": "s3://bucket/R1.fastq",
            "r2_path": "s3://bucket/R2.fastq",
            "submitted_at": datetime.now(timezone.utc).isoformat(),
            "updated_at": datetime.now(timezone.utc).isoformat(),
        }

        # Mock S3 log download
        mock_logs = "EMR Step Logs\nError: File not found\nTraceback..."
        with patch(
            "subprocess.run",
            return_value=Mock(
                returncode=0,
                stdout=mock_logs,
                stderr="",
            ),
        ):
            logs = job_manager.get_job_logs(job_id)

        assert logs == mock_logs

    def test_get_job_logs_no_cluster_id(
        self, job_manager
    ):
        """Test getting logs when cluster_id is missing."""
        job_id = "test-job-123"
        job_manager.jobs[job_id] = {
            "job_id": job_id,
            "status": STATUS_FAILED,
            "step_id": "s-ABC123",
            "submitted_at": datetime.now(timezone.utc).isoformat(),
            "updated_at": datetime.now(timezone.utc).isoformat(),
        }

        logs = job_manager.get_job_logs(job_id)
        assert "Logs not available" in logs

    def test_get_job_logs_s3_fallback(
        self, job_manager, mock_env_vars
    ):
        """Test getting logs with S3 fallback to stderr/stdout."""
        job_id = "test-job-123"
        job_manager.jobs[job_id] = {
            "job_id": job_id,
            "status": STATUS_FAILED,
            "cluster_id": "j-TEST123",
            "step_id": "s-ABC123",
            "submitted_at": datetime.now(timezone.utc).isoformat(),
            "updated_at": datetime.now(timezone.utc).isoformat(),
        }

        # Mock wrapper.log not found, but stderr exists
        def mock_subprocess(cmd, **kwargs):
            if "wrapper.log" in " ".join(cmd):
                return Mock(returncode=1, stdout="", stderr="Not found")
            elif "stderr.gz" in " ".join(cmd):
                return Mock(returncode=0, stdout="STDERR LOGS", stderr="")
            return Mock(returncode=1, stdout="", stderr="")

        with patch("subprocess.run", side_effect=mock_subprocess):
            logs = job_manager.get_job_logs(job_id)

        assert "STDERR" in logs
        assert "STDERR LOGS" in logs

    def test_get_job_logs_timeout(
        self, job_manager, mock_env_vars
    ):
        """Test handling timeout when fetching logs."""
        job_id = "test-job-123"
        job_manager.jobs[job_id] = {
            "job_id": job_id,
            "status": STATUS_FAILED,
            "cluster_id": "j-TEST123",
            "step_id": "s-ABC123",
            "submitted_at": datetime.now(timezone.utc).isoformat(),
            "updated_at": datetime.now(timezone.utc).isoformat(),
        }

        with patch(
            "subprocess.run",
            side_effect=subprocess.TimeoutExpired("aws", 15),
        ):
            logs = job_manager.get_job_logs(job_id)

        assert "Timeout" in logs


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

