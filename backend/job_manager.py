"""
Job Manager for tracking and managing long-running EMR FastQC jobs.

This module provides a JobManager class that tracks job status, metadata,
and handles communication with AWS EMR for FastQC analysis jobs.
"""

import os
import uuid
import logging
import subprocess
import re
import json
from datetime import datetime, timezone
from typing import Dict, Optional, List, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)

# Job status constants
STATUS_SUBMITTED = "submitted"
STATUS_RUNNING = "running"
STATUS_COMPLETED = "completed"
STATUS_FAILED = "failed"
STATUS_CANCELLED = "cancelled"

VALID_STATUSES = {
    STATUS_SUBMITTED,
    STATUS_RUNNING,
    STATUS_COMPLETED,
    STATUS_FAILED,
    STATUS_CANCELLED,
}

# EMR step status to job status mapping
EMR_STATUS_TO_JOB_STATUS = {
    "PENDING": STATUS_SUBMITTED,
    "RUNNING": STATUS_RUNNING,
    "COMPLETED": STATUS_COMPLETED,
    "CANCELLED": STATUS_CANCELLED,
    "FAILED": STATUS_FAILED,
    "INTERRUPTED": STATUS_FAILED,
    "CANCEL_PENDING": STATUS_CANCELLED,
    "CANCEL_PENDING": STATUS_CANCELLED,
}


class JobManager:
    """
    Manages FastQC job lifecycle and status tracking.
    
    For MVP, uses in-memory storage. Can be upgraded to Redis/DB later
    for persistence across restarts.
    """

    def __init__(self):
        """Initialize the JobManager with empty job store."""
        self.jobs: Dict[str, Dict] = {}
        self.project_root = Path(__file__).resolve().parent.parent

    def submit_fastqc_job(
        self,
        r1_path: str,
        r2_path: str,
        output_path: Optional[str] = None,
        session_id: Optional[str] = None,
    ) -> str:
        """
        Submit a FastQC job to EMR cluster.
        
        Args:
            r1_path: S3 URI for forward/read 1 FASTQ file (e.g., s3://bucket/data/R1.fastq)
            r2_path: S3 URI for reverse/read 2 FASTQ file (e.g., s3://bucket/data/R2.fastq)
            output_path: Optional S3 URI for results directory. If None and session_id provided,
                        will use session's dedicated S3 path.
            session_id: Optional session ID. If provided and output_path is None, results will be
                       saved to the session's dedicated S3 location.
        
        Returns:
            job_id: Unique identifier for tracking this job
        
        Raises:
            ValueError: If required environment variables are missing
            RuntimeError: If job submission fails
        """
        # Validate inputs
        if not r1_path or not r2_path:
            raise ValueError("Both R1 and R2 paths are required")
        
        if not r1_path.startswith("s3://") or not r2_path.startswith("s3://"):
            raise ValueError("Both R1 and R2 paths must be S3 URIs (s3://...)")
        
        # Generate unique job ID first (needed for session-specific folder naming)
        job_id = str(uuid.uuid4())
        
        # Determine output path: use session S3 path if session_id provided and output_path is None
        if not output_path and session_id:
            try:
                from history_manager import history_manager
                session = history_manager.get_session(session_id)
                if session:
                    s3_bucket = session.get("metadata", {}).get("s3_bucket")
                    s3_path = session.get("metadata", {}).get("s3_path")
                    if s3_bucket and s3_path:
                        # Use session's dedicated S3 path with job_id as sub-folder: s3://bucket/session_id/job_id/
                        output_path = f"s3://{s3_bucket}/{s3_path}{job_id}/"
                        logger.info(f"Using session S3 path for job output: {output_path}")
                    else:
                        logger.warning(f"Session {session_id} exists but has no S3 path. Using default output location.")
                else:
                    logger.warning(f"Session {session_id} not found. Using default output location.")
            except Exception as e:
                logger.warning(f"Failed to get session S3 path: {e}. Using default output location.")
        
        # Get EMR cluster ID from environment
        cluster_id = os.getenv("EMR_CLUSTER_ID")
        if not cluster_id:
            raise ValueError(
                "EMR_CLUSTER_ID environment variable not set. "
                "Please set it to your EMR cluster ID (e.g., j-XXXXXXXXXXXXX)"
            )
        
        # Check cluster state before submission
        cluster_state = self._check_cluster_state(cluster_id)
        if not cluster_state or cluster_state not in ("WAITING", "RUNNING"):
            # Cluster doesn't exist, is terminated, or not ready
            if cluster_state:
                logger.info(f"Cluster {cluster_id} is in state '{cluster_state}' (not ready)")
            else:
                logger.info(f"Cluster {cluster_id} not found or inaccessible")
            
            # Try to find an active cluster
            active_cluster_id = self._find_active_cluster()
            if active_cluster_id:
                logger.info(
                    f"Found active cluster {active_cluster_id}, using it instead."
                )
                cluster_id = active_cluster_id
                # Update environment for the script
                os.environ["EMR_CLUSTER_ID"] = cluster_id
            else:
                # No active cluster found - create a new one
                logger.info(
                    f"No active clusters found. Creating a new EMR cluster..."
                )
                new_cluster_id = self._create_emr_cluster()
                if not new_cluster_id:
                    raise RuntimeError(
                        f"Failed to create new EMR cluster. "
                        f"Please create a cluster manually using ./scripts/aws/setup-emr-cluster.sh"
                    )
                
                # Wait for cluster to be ready
                logger.info(f"Waiting for new cluster {new_cluster_id} to be ready...")
                if not self._wait_for_cluster_ready(new_cluster_id, max_wait_minutes=15):
                    raise RuntimeError(
                        f"New EMR cluster {new_cluster_id} did not become ready within 15 minutes. "
                        f"Please check the cluster status manually."
                    )
                
                cluster_id = new_cluster_id
                # Update environment for the script
                os.environ["EMR_CLUSTER_ID"] = cluster_id
                logger.info(f"✅ Using newly created cluster: {cluster_id}")
        
        # Call submission script
        script_path = self.project_root / "scripts" / "emr" / "submit-fastqc-job.sh"
        
        if not script_path.exists():
            raise FileNotFoundError(
                f"Submission script not found at {script_path}"
            )
        
        # Prepare environment for subprocess
        env = os.environ.copy()
        env["EMR_CLUSTER_ID"] = cluster_id
        
        # Build command - use absolute path and set working directory to project root
        cmd = [str(script_path)]
        if r1_path:
            cmd.append(r1_path)
        if r2_path:
            cmd.append(r2_path)
        if output_path:
            cmd.append(output_path)
        
        logger.info(f"Submitting FastQC job {job_id}: R1={r1_path}, R2={r2_path}")
        logger.info(f"Using script: {script_path}")
        logger.info(f"Working directory: {self.project_root}")
        
        try:
            # Run submission script and capture output
            # Set cwd to project root so relative paths in script work correctly
            result = subprocess.run(
                cmd,
                env=env,
                cwd=str(self.project_root),
                capture_output=True,
                text=True,
                timeout=60,  # 60 second timeout for submission
            )
            
            if result.returncode != 0:
                error_msg = result.stderr or result.stdout or "Unknown error"
                logger.error(f"Job submission failed: {error_msg}")
                raise RuntimeError(f"Failed to submit job: {error_msg}")
            
            # Extract step ID from output
            # The script outputs: "Step ID: s-XXXXXXXXXXXXX"
            step_id = self._extract_step_id(result.stdout)
            
            if not step_id:
                logger.warning(
                    f"Could not extract step ID from output: {result.stdout}"
                )
                # Still create job record, but without step_id
                # Status checking will be limited
        
        except subprocess.TimeoutExpired:
            raise RuntimeError("Job submission timed out after 60 seconds")
        except Exception as e:
            logger.error(f"Exception during job submission: {e}")
            raise RuntimeError(f"Job submission failed: {str(e)}")
        
        # Create job directory in session's directory if session_id is provided
        job_local_path = None
        if session_id:
            try:
                from history_manager import history_manager
                # Get session (this will ensure directory exists)
                session = history_manager.get_session(session_id)
                if session:
                    session_local_path = session.get("metadata", {}).get("local_path")
                    if session_local_path:
                        job_dir = Path(session_local_path) / job_id
                        job_dir.mkdir(parents=True, exist_ok=True)
                        job_local_path = str(job_dir)
                        logger.info(f"Created job directory: {job_dir}")
                    else:
                        # Fallback: create directory structure from storage_dir
                        storage_dir = history_manager.storage_dir
                        session_dir = storage_dir / session_id
                        session_dir.mkdir(parents=True, exist_ok=True)  # Ensure session dir exists
                        job_dir = session_dir / job_id
                        job_dir.mkdir(parents=True, exist_ok=True)
                        job_local_path = str(job_dir)
                        logger.info(f"Created job directory (fallback): {job_dir}")
                else:
                    # Session doesn't exist, create it
                    logger.warning(f"Session {session_id} not found, creating directory structure anyway")
                    storage_dir = history_manager.storage_dir
                    session_dir = storage_dir / session_id
                    session_dir.mkdir(parents=True, exist_ok=True)
                    job_dir = session_dir / job_id
                    job_dir.mkdir(parents=True, exist_ok=True)
                    job_local_path = str(job_dir)
                    logger.info(f"Created job directory (session not found): {job_dir}")
            except Exception as e:
                logger.warning(f"Failed to create job directory for session {session_id}: {e}")
                import traceback
                traceback.print_exc()
        
        # Store job metadata
        self.jobs[job_id] = {
            "job_id": job_id,
            "status": STATUS_SUBMITTED,
            "type": "fastqc",
            # Phase 1.4: generalized job fields (keep legacy keys too)
            "job_type": "tool",
            "tool_name": "fastqc_quality_analysis",
            "args": {
                "r1_path": r1_path,
                "r2_path": r2_path,
                "output_path": output_path,
            },
            "infra": {
                "provider": "aws",
                "service": "emr",
                "cluster_id": cluster_id,
                "step_id": step_id,
            },
            "results": {},
            "cluster_id": cluster_id,
            "step_id": step_id,
            "r1_path": r1_path,
            "r2_path": r2_path,
            "output_path": output_path,
            "session_id": session_id,
            "local_path": job_local_path,
            "submitted_at": datetime.now(timezone.utc).isoformat(),
            "updated_at": datetime.now(timezone.utc).isoformat(),
            "error": None,
        }
        
        # Save job metadata to local directory if job_local_path exists
        if job_local_path:
            try:
                job_metadata_file = Path(job_local_path) / "job_metadata.json"
                with open(job_metadata_file, 'w') as f:
                    json.dump(self.jobs[job_id], f, indent=2)
                logger.info(f"Saved job metadata to: {job_metadata_file}")
            except Exception as e:
                logger.warning(f"Failed to save job metadata to local directory: {e}")
        
        logger.info(f"Job {job_id} submitted successfully with step_id {step_id}")
        return job_id

    def _extract_step_id(self, output: str) -> Optional[str]:
        """
        Extract EMR step ID from script output.
        
        Looks for patterns like:
        - "Step ID: s-XXXXXXXXXXXXX"
        - "s-[A-Z0-9]+" in the output
        """
        # Try to find "Step ID: s-XXXXXXXXXXXXX"
        match = re.search(r"Step ID:\s*(s-[A-Z0-9]+)", output, re.IGNORECASE)
        if match:
            return match.group(1)
        
        # Try to find any step ID pattern
        match = re.search(r"(s-[A-Z0-9]+)", output)
        if match:
            return match.group(1)
        
        return None

    def get_job_status(self, job_id: str) -> Dict:
        """
        Get current status of a job, checking EMR if needed.
        
        Args:
            job_id: Unique job identifier
        
        Returns:
            Dict with job status and metadata
        
        Raises:
            ValueError: If job_id not found
        """
        if job_id not in self.jobs:
            raise ValueError(f"Job {job_id} not found")
        
        job = self.jobs[job_id]
        
        # If job is already terminal, return cached status
        if job["status"] in (STATUS_COMPLETED, STATUS_FAILED, STATUS_CANCELLED):
            return job.copy()
        
        # Check EMR step status
        if job.get("step_id") and job.get("cluster_id"):
            try:
                emr_status, failure_details = self._check_emr_step_status(
                    job["cluster_id"], job["step_id"]
                )
                if emr_status:
                    # Update job status based on EMR status
                    old_status = job["status"]
                    job_status = EMR_STATUS_TO_JOB_STATUS.get(
                        emr_status, STATUS_RUNNING
                    )
                    if job_status != old_status:
                        job["status"] = job_status
                        job["updated_at"] = datetime.now(timezone.utc).isoformat()
                        logger.info(
                            f"Job {job_id} status updated: {old_status} -> {job_status}"
                        )
                        
                        # Save updated job metadata to local directory
                        self._save_job_metadata_to_local(job_id)
                        
                        # If job just completed successfully, generate HTML visualizations
                        if job_status == STATUS_COMPLETED and old_status != STATUS_COMPLETED:
                            # Job just completed - trigger HTML generation
                            try:
                                self._generate_and_upload_html_visualizations(job_id)
                            except Exception as e:
                                logger.warning(f"Failed to generate HTML visualizations for job {job_id}: {e}")
                                # Don't fail the job status check if HTML generation fails
                    
                    # If job failed, store failure details
                    if job_status == STATUS_FAILED:
                        # Try to get detailed error from logs
                        error_message = self._extract_user_friendly_error(
                            job["cluster_id"], job["step_id"], failure_details
                        )
                        job["error"] = error_message
                        logger.info(
                            f"Job {job_id} failed: {error_message}"
                        )
            except Exception as e:
                logger.warning(
                    f"Failed to check EMR status for job {job_id}: {e}"
                )
                # Continue with cached status
        
        return job.copy()

    def _check_emr_step_status(self, cluster_id: str, step_id: str) -> Tuple[Optional[str], Optional[Dict]]:
        """
        Check EMR step status and failure details using AWS CLI.
        
        Args:
            cluster_id: EMR cluster ID
            step_id: EMR step ID
        
        Returns:
            Tuple of (EMR step status string or None, failure details dict or None)
        """
        try:
            region = os.getenv("AWS_REGION", "us-east-1")
            
            # Get step status
            status_result = subprocess.run(
                [
                    "aws",
                    "emr",
                    "describe-step",
                    "--cluster-id",
                    cluster_id,
                    "--step-id",
                    step_id,
                    "--region",
                    region,
                    "--query",
                    "Step.Status.State",
                    "--output",
                    "text",
                ],
                capture_output=True,
                text=True,
                timeout=10,
            )
            
            status = None
            failure_details = None
            
            if status_result.returncode == 0:
                status = status_result.stdout.strip()
                
                # If status is FAILED, get failure details
                if status == "FAILED":
                    failure_result = subprocess.run(
                        [
                            "aws",
                            "emr",
                            "describe-step",
                            "--cluster-id",
                            cluster_id,
                            "--step-id",
                            step_id,
                            "--region",
                            region,
                            "--query",
                            "Step.Status.FailureDetails",
                            "--output",
                            "json",
                        ],
                        capture_output=True,
                        text=True,
                        timeout=10,
                    )
                    
                    if failure_result.returncode == 0 and failure_result.stdout.strip():
                        import json
                        try:
                            failure_details = json.loads(failure_result.stdout.strip())
                            # Handle case where AWS returns empty dict or "null"
                            if not failure_details or failure_details == {}:
                                failure_details = None
                        except json.JSONDecodeError:
                            logger.warning(f"Failed to parse failure details JSON: {failure_result.stdout}")
                            failure_details = None
                
                return (status if status else None, failure_details)
            else:
                logger.warning(
                    f"AWS CLI error checking step status: {status_result.stderr}"
                )
                return (None, None)
                
        except subprocess.TimeoutExpired:
            logger.warning("Timeout checking EMR step status")
            return (None, None)
        except Exception as e:
            logger.warning(f"Exception checking EMR step status: {e}")
            return (None, None)

    def _extract_user_friendly_error(
        self, cluster_id: str, step_id: str, failure_details: Optional[Dict]
    ) -> str:
        """
        Extract and format user-friendly error message from EMR failure.
        
        Tries multiple sources:
        1. EMR FailureDetails (if available)
        2. Wrapper log from S3 (most detailed)
        3. Generic helpful message
        
        Returns a user-friendly error message with actionable guidance.
        """
        # First, try to get error from wrapper log in S3
        wrapper_log_error = self._extract_error_from_wrapper_log(cluster_id, step_id)
        if wrapper_log_error:
            return wrapper_log_error
        
        # If we have EMR failure details, use them
        if failure_details:
            reason = failure_details.get("Reason", "")
            message = failure_details.get("Message", "")
            
            if reason or message:
                # Translate common EMR error reasons to user-friendly messages
                error_text = reason or message
                return self._translate_error_to_user_friendly(error_text)
        
        # Fallback: Generic but helpful message
        return (
            "The FastQC analysis job failed on the EMR cluster. "
            "This could be due to:\n"
            "• Issues with the input FASTQ files (file format, permissions, or file not found)\n"
            "• Problems with the analysis script\n"
            "• Insufficient cluster resources\n"
            "• Network or S3 access issues\n\n"
            f"To investigate further, check the detailed logs at:\n"
            f"s3://noricum-ngs-data/emr-logs/{cluster_id}/steps/{step_id}/"
        )
    
    def _extract_error_from_wrapper_log(self, cluster_id: str, step_id: str) -> Optional[str]:
        """
        Try to extract error message from wrapper.log in S3.
        Returns user-friendly error message or None if not available.
        """
        try:
            region = os.getenv("AWS_REGION", "us-east-1")
            s3_bucket = os.getenv("S3_DATASET_BUCKET", "noricum-ngs-data")
            log_path = f"s3://{s3_bucket}/emr-logs/{cluster_id}/steps/{step_id}/wrapper.log"
            
            # Try to download and parse the wrapper log
            result = subprocess.run(
                [
                    "aws",
                    "s3",
                    "cp",
                    log_path,
                    "-",
                    "--region",
                    region,
                ],
                capture_output=True,
                text=True,
                timeout=15,
            )
            
            if result.returncode == 0 and result.stdout:
                # Look for common error patterns in the log
                log_content = result.stdout
                
                # Check for Python syntax errors
                if "IndentationError" in log_content or "SyntaxError" in log_content:
                    return (
                        "The analysis script has a syntax error. "
                        "This is a code issue that needs to be fixed by the development team. "
                        "The error has been logged and will be addressed."
                    )
                
                # Check for file not found errors
                if "FileNotFoundError" in log_content or "No such file" in log_content:
                    return (
                        "A required file was not found. "
                        "This could mean:\n"
                        "• The input FASTQ files don't exist at the specified S3 location\n"
                        "• The analysis script is missing\n"
                        "• A required dependency is not available\n\n"
                        "Please verify that your input files exist and are accessible."
                    )
                
                # Check for permission errors
                if "Permission denied" in log_content or "AccessDenied" in log_content:
                    return (
                        "Permission denied error. "
                        "The EMR cluster doesn't have permission to access the required files. "
                        "Please check:\n"
                        "• S3 bucket permissions for the input files\n"
                        "• IAM role permissions for the EMR cluster\n"
                        "• Output path write permissions"
                    )
                
                # Check for Spark/execution errors
                if "SparkException" in log_content or "ApplicationMaster" in log_content:
                    return (
                        "The Spark job failed during execution. "
                        "This could be due to:\n"
                        "• Insufficient memory or CPU resources on the cluster\n"
                        "• Issues processing the input data\n"
                        "• Problems with the data format\n\n"
                        "Try reducing the input file size or increasing cluster resources."
                    )
                
                # Check for timeout errors
                if "timeout" in log_content.lower() or "Timeout" in log_content:
                    return (
                        "The job timed out. "
                        "The analysis took too long to complete. "
                        "This could mean:\n"
                        "• The input files are very large\n"
                        "• The cluster is under heavy load\n"
                        "• There are performance issues\n\n"
                        "Try using a larger cluster or splitting the analysis into smaller chunks."
                    )
                
                # Extract last few lines of error output
                lines = log_content.split('\n')
                error_lines = []
                for line in reversed(lines[-50:]):  # Check last 50 lines
                    if any(keyword in line.lower() for keyword in ['error', 'failed', 'exception', 'traceback']):
                        error_lines.insert(0, line.strip())
                        if len(error_lines) >= 5:  # Get up to 5 error lines
                            break
                
                if error_lines:
                    error_snippet = '\n'.join(error_lines[-3:])  # Last 3 error lines
                    return (
                        f"The job failed with the following error:\n\n{error_snippet}\n\n"
                        "If this error is unclear, please check the detailed logs for more information."
                    )
            
        except subprocess.TimeoutExpired:
            logger.warning("Timeout fetching wrapper log from S3")
        except Exception as e:
            logger.warning(f"Failed to extract error from wrapper log: {e}")
        
        return None
    
    def _generate_and_upload_html_visualizations(self, job_id: str) -> None:
        """
        Generate HTML visualizations for a completed FastQC job and upload to session S3.
        
        Args:
            job_id: Job identifier
        """
        try:
            job = self.jobs.get(job_id)
            if not job:
                logger.warning(f"Job {job_id} not found for HTML generation")
                return
            
            # Get results path
            results_path = None
            # First, try to get from session_results_path if available
            if job.get("session_results_path"):
                results_path = job.get("session_results_path")
            else:
                # Don't call get_job_results() here as it might trigger HTML generation again
                # Instead, construct the path directly from available job metadata
                output_path = job.get("output_path")
                if output_path:
                    results_path = output_path.rstrip("/") + "/results.json"
                else:
                    # Construct from default pattern (same logic as get_job_results but without recursion)
                    r1_path = job.get("r1_path", "")
                    s3_bucket = os.getenv("S3_DATASET_BUCKET", "noricum-ngs-data")
                    
                    if r1_path:
                        parts = r1_path.replace("s3://", "").split("/", 1)
                        bucket = parts[0] if parts else s3_bucket
                    else:
                        bucket = s3_bucket
                    
                    # Try to extract dataset identifier from R1 path
                    dataset_name = None
                    if r1_path:
                        r1_parts = r1_path.split("/")
                        for i, part in enumerate(r1_parts):
                            if "GRCh38" in part and "MafHi" in part:
                                dataset_name = part
                                break
                            elif "GRCh38" in part and i + 1 < len(r1_parts):
                                if "MafHi" in r1_parts[i + 1]:
                                    dataset_name = f"{part}/{r1_parts[i + 1]}"
                                    break
                    
                    if dataset_name:
                        results_path = f"s3://{bucket}/fastqc-results/{dataset_name}/results.json"
                    else:
                        results_path = f"s3://{bucket}/fastqc-results/GRCh38.p12.MafHi/results.json"
            
            if not results_path:
                logger.warning(f"No results path for job {job_id}")
                return
            
            # Get session S3 path if available
            session_id = job.get("session_id")
            session_s3_path = None
            if session_id:
                try:
                    from history_manager import history_manager
                    session = history_manager.get_session(session_id)
                    if session:
                        s3_bucket = session.get("metadata", {}).get("s3_bucket")
                        s3_path = session.get("metadata", {}).get("s3_path")
                        if s3_bucket and s3_path:
                            # Use job_id as sub-folder name: s3://bucket/session_id/job_id/
                            session_s3_path = f"s3://{s3_bucket}/{s3_path}{job_id}/"
                except Exception as e:
                    logger.warning(f"Failed to get session S3 path: {e}")
            
            # If no session path, use the output_path from job
            if not session_s3_path:
                output_path = job.get("output_path")
                if output_path:
                    session_s3_path = output_path.rstrip("/") + "/"
                else:
                    logger.warning(f"No session S3 path or output_path for job {job_id}, skipping HTML upload")
                    return
            
            logger.info(f"Generating HTML visualizations for job {job_id}, uploading to {session_s3_path}")
            
            # Import visualization functions
            import sys
            import tempfile
            import shutil
            visualization_script = self.project_root / "scripts" / "emr" / "visualize_fastqc_results.py"
            
            if not visualization_script.exists():
                logger.warning(f"Visualization script not found at {visualization_script}")
                return
            
            # Create temporary directory for HTML files
            # Note: This directory will be automatically deleted when the 'with' block exits
            with tempfile.TemporaryDirectory() as temp_dir:
                logger.info(f"Created temporary directory for HTML generation: {temp_dir}")
                # Import the visualization module
                sys.path.insert(0, str(visualization_script.parent))
                try:
                    import importlib.util
                    spec = importlib.util.spec_from_file_location("visualize_fastqc_results", visualization_script)
                    viz_module = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(viz_module)
                    
                    # Load results
                    logger.info(f"Loading results from: {results_path}")
                    results = viz_module.load_results(results_path)
                    logger.info("✅ Results loaded successfully")
                    
                    # Generate visualizations
                    html_dir = os.path.join(temp_dir, "html")
                    os.makedirs(html_dir, exist_ok=True)
                    logger.info(f"Generating HTML visualizations in: {html_dir}")
                    main_html_path = viz_module.create_visualizations(results, html_dir)
                    
                    if not main_html_path:
                        logger.warning(f"Failed to generate HTML visualizations for job {job_id}")
                        return
                    
                    logger.info(f"✅ HTML visualizations generated successfully at: {main_html_path}")
                    
                    # Upload HTML files to session S3 path
                    import boto3
                    s3_client = boto3.client('s3')
                    
                    # Parse S3 path
                    s3_path_parts = session_s3_path.replace("s3://", "").split("/", 1)
                    bucket = s3_path_parts[0]
                    prefix = s3_path_parts[1] if len(s3_path_parts) > 1 else ""
                    
                    # Upload main HTML file
                    main_html_filename = os.path.basename(main_html_path)
                    s3_key = f"{prefix}{main_html_filename}" if prefix else main_html_filename
                    try:
                        s3_client.upload_file(main_html_path, bucket, s3_key)
                        logger.info(f"✅ Uploaded {main_html_filename} to s3://{bucket}/{s3_key}")
                    except Exception as e:
                        logger.error(f"❌ Failed to upload {main_html_filename}: {e}")
                        raise
                    
                    # Upload individual chart HTML files if they exist
                    html_files = [f for f in os.listdir(html_dir) if f.endswith('.html') and f != main_html_filename]
                    for html_file in html_files:
                        local_path = os.path.join(html_dir, html_file)
                        chart_s3_key = f"{prefix}{html_file}" if prefix else html_file
                        try:
                            s3_client.upload_file(local_path, bucket, chart_s3_key)
                            logger.info(f"✅ Uploaded {html_file} to s3://{bucket}/{chart_s3_key}")
                            
                            # Also copy to local job directory if it exists
                            self._copy_file_to_job_directory(job_id, local_path, html_file)
                        except Exception as e:
                            logger.error(f"❌ Failed to upload {html_file}: {e}")
                            # Continue with other files even if one fails
                    
                    # Copy main HTML file to local job directory
                    self._copy_file_to_job_directory(job_id, main_html_path, main_html_filename)
                    
                    # Always copy results.json to session S3 path (even if it's the same path, ensures it exists)
                    if results_path:
                        try:
                            # Parse original location
                            orig_parts = results_path.replace("s3://", "").split("/", 1)
                            orig_bucket = orig_parts[0]
                            orig_key = orig_parts[1] if len(orig_parts) > 1 else ""
                            
                            # Upload to session location
                            session_results_key = f"{prefix}results.json" if prefix else "results.json"
                            
                            # Only copy if source and destination are different
                            if results_path != f"{session_s3_path}results.json":
                                s3_client.copy_object(
                                    CopySource={'Bucket': orig_bucket, 'Key': orig_key},
                                    Bucket=bucket,
                                    Key=session_results_key
                                )
                                logger.info(f"Copied results.json to session S3 path: s3://{bucket}/{session_results_key}")
                            else:
                                logger.info(f"Results.json already at session S3 path: s3://{bucket}/{session_results_key}")
                            
                            # Also download and save to local job directory
                            try:
                                import tempfile
                                with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.json') as tmp_file:
                                    tmp_path = tmp_file.name
                                s3_client.download_file(orig_bucket, orig_key, tmp_path)
                                self._copy_file_to_job_directory(job_id, tmp_path, "results.json")
                                os.unlink(tmp_path)  # Clean up temp file
                            except Exception as e:
                                logger.warning(f"Failed to copy results.json to local job directory: {e}")
                        except Exception as e:
                            logger.warning(f"Failed to copy results.json to session S3 path: {e}")
                            # Don't fail the whole process if copy fails
                    
                    # Update job with session results path
                    job["session_results_path"] = f"{session_s3_path}results.json"
                    job["session_html_path"] = f"{session_s3_path}{main_html_filename}"
                    
                    # Ensure job directory exists and save updated metadata
                    self._ensure_job_directory(job_id, job)
                    self._save_job_metadata_to_local(job_id)
                    
                    logger.info(f"✅ Successfully generated and uploaded HTML visualizations for job {job_id}")
                    logger.info(f"   Results JSON: s3://{bucket}/{session_results_key}")
                    logger.info(f"   Main HTML: s3://{bucket}/{s3_key}")
                    logger.info(f"   Note: Temporary local files have been cleaned up (this is expected)")
                    
                except Exception as e:
                    logger.error(f"❌ Failed to generate HTML visualizations for job {job_id}: {e}")
                    import traceback
                    traceback.print_exc()
                    raise
        
        except Exception as e:
            logger.error(f"Error in _generate_and_upload_html_visualizations for job {job_id}: {e}")
            import traceback
            traceback.print_exc()
    
    def _save_job_metadata_to_local(self, job_id: str) -> None:
        """Save job metadata to local job directory if it exists.
        
        If directory doesn't exist, tries to create it from job metadata.
        """
        try:
            job = self.jobs.get(job_id)
            if not job:
                return
            
            job_local_path = job.get("local_path")
            
            # If no local path, try to create it
            if not job_local_path:
                job_local_path = self._ensure_job_directory(job_id, job)
            
            if not job_local_path:
                return
            
            job_metadata_file = Path(job_local_path) / "job_metadata.json"
            job_metadata_file.parent.mkdir(parents=True, exist_ok=True)
            with open(job_metadata_file, 'w') as f:
                json.dump(job, f, indent=2)
            logger.debug(f"Saved job metadata to: {job_metadata_file}")
        except Exception as e:
            logger.warning(f"Failed to save job metadata to local directory: {e}")
    
    def _ensure_job_directory(self, job_id: str, job: Optional[Dict] = None) -> Optional[str]:
        """Ensure job directory exists, creating it if necessary.
        
        Returns the job_local_path if successful, None otherwise.
        """
        if job is None:
            job = self.jobs.get(job_id)
            if not job:
                return None
        
        session_id = job.get("session_id")
        if not session_id:
            return None
        
        try:
            from history_manager import history_manager
            # Get session (this will ensure session directory exists)
            session = history_manager.get_session(session_id)
            
            session_local_path = None
            if session:
                session_local_path = session.get("metadata", {}).get("local_path")
            
            if not session_local_path:
                # Fallback: use storage_dir
                storage_dir = history_manager.storage_dir
                session_dir = storage_dir / session_id
                session_dir.mkdir(parents=True, exist_ok=True)
                session_local_path = str(session_dir)
            
            job_dir = Path(session_local_path) / job_id
            job_dir.mkdir(parents=True, exist_ok=True)
            job_local_path = str(job_dir)
            
            # Update job with local path
            job["local_path"] = job_local_path
            
            logger.info(f"Created/ensured job directory: {job_dir}")
            return job_local_path
        except Exception as e:
            logger.warning(f"Failed to ensure job directory for job {job_id}: {e}")
            return None
    
    def _copy_file_to_job_directory(self, job_id: str, source_path: str, filename: str) -> None:
        """Copy a file to the job's local directory if it exists.
        
        If job directory doesn't exist, try to create it from job metadata.
        """
        try:
            job = self.jobs.get(job_id)
            if not job:
                return
            
            job_local_path = job.get("local_path")
            
            # If no local path, try to create it
            if not job_local_path:
                job_local_path = self._ensure_job_directory(job_id, job)
            
            if not job_local_path or not os.path.exists(source_path):
                return
            
            dest_path = Path(job_local_path) / filename
            # Ensure parent directory exists
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            import shutil
            shutil.copy2(source_path, dest_path)
            logger.info(f"✅ Copied {filename} to local job directory: {dest_path}")
        except Exception as e:
            logger.warning(f"Failed to copy {filename} to local job directory: {e}")
    
    def _translate_error_to_user_friendly(self, error_text: str) -> str:
        """
        Translate technical error messages to user-friendly language.
        """
        error_lower = error_text.lower()
        
        # Common error translations
        if "indentationerror" in error_lower or "syntaxerror" in error_lower:
            return (
                "The analysis script has a syntax error. "
                "This is a code issue that needs to be fixed by the development team."
            )
        
        if "filenotfound" in error_lower or "no such file" in error_lower:
            return (
                "A required file was not found. "
                "Please verify that your input FASTQ files exist at the specified S3 locations."
            )
        
        if "permission" in error_lower or "access denied" in error_lower:
            return (
                "Permission denied. "
                "The EMR cluster doesn't have permission to access the required files. "
                "Please check S3 bucket and IAM permissions."
            )
        
        if "timeout" in error_lower:
            return (
                "The job timed out. "
                "The analysis took too long to complete. "
                "Try using a larger cluster or smaller input files."
            )
        
        if "memory" in error_lower or "outofmemory" in error_lower:
            return (
                "Out of memory error. "
                "The cluster doesn't have enough memory to process the data. "
                "Try using a larger cluster instance type or reducing the input file size."
            )
        
        # Return the original error if we can't translate it, but format it nicely
        return f"The job failed: {error_text}"

    def get_job_results(self, job_id: str) -> Dict:
        """
        Get results for a completed job.
        
        Args:
            job_id: Unique job identifier
        
        Returns:
            Dict with job results
        
        Raises:
            ValueError: If job_id not found or job not completed
        """
        job = self.get_job_status(job_id)
        
        if job["status"] != STATUS_COMPLETED:
            raise ValueError(
                f"Job {job_id} is not completed (current status: {job['status']})"
            )
        
        # If job has session_id but no session_results_path, try to copy results now
        # But only if we're not already in the process of generating HTML (avoid infinite loop)
        session_id = job.get("session_id")
        if session_id and not job.get("session_results_path") and not job.get("_generating_html"):
            logger.info(f"Job {job_id} has session_id but no session_results_path, attempting to copy results")
            try:
                # Set flag to prevent recursive calls
                job["_generating_html"] = True
                self._generate_and_upload_html_visualizations(job_id)
                # Re-fetch job to get updated session_results_path
                job = self.get_job_status(job_id)
            except Exception as e:
                logger.warning(f"Failed to copy results to session for job {job_id}: {e}")
                # Continue with original results path
            finally:
                # Clear flag
                if job_id in self.jobs:
                    self.jobs[job_id].pop("_generating_html", None)
        
        # Results should be at output_path/results.json
        output_path = job.get("output_path")
        if not output_path:
            # When output_path is None, the script uses its default pattern.
            # The script default is: s3://${S3_BUCKET}/fastqc-results/GRCh38.p12.MafHi/
            # But it may also use a pattern based on the input file location.
            # We'll try multiple strategies to find the results.
            
            r1_path = job.get("r1_path", "")
            s3_bucket = os.getenv("S3_DATASET_BUCKET", "noricum-ngs-data")
            
            # Strategy 1: Try script's default pattern (most common when output_path is None)
            # Extract bucket from R1 path or use env var
            if r1_path:
                parts = r1_path.replace("s3://", "").split("/", 1)
                bucket = parts[0] if parts else s3_bucket
            else:
                bucket = s3_bucket
            
            # Try to extract dataset identifier from R1 path
            dataset_name = None
            if r1_path:
                r1_parts = r1_path.split("/")
                # Look for GRCh38.p12.MafHi pattern
                for i, part in enumerate(r1_parts):
                    if "GRCh38" in part and "MafHi" in part:
                        dataset_name = part
                        break
                    elif "GRCh38" in part and i + 1 < len(r1_parts):
                        # Check if next part contains MafHi
                        if "MafHi" in r1_parts[i + 1]:
                            dataset_name = f"{part}/{r1_parts[i + 1]}"
                            break
            
            # Use script's default pattern
            if dataset_name:
                output_path = f"s3://{bucket}/fastqc-results/{dataset_name}/"
            else:
                # Fallback to script's hardcoded default
                output_path = f"s3://{bucket}/fastqc-results/GRCh38.p12.MafHi/"
        
        if not output_path:
            raise ValueError("Cannot determine output path for results")
        
        results_path = output_path.rstrip("/") + "/results.json"
        
        # If job has session_results_path (from HTML generation), use that instead
        session_results_path = job.get("session_results_path")
        if session_results_path:
            results_path = session_results_path
        
        return {
            "job_id": job_id,
            "status": job["status"],
            "results_path": results_path,
            "output_path": output_path,
            "r1_path": job.get("r1_path"),
            "r2_path": job.get("r2_path"),
            "session_results_path": session_results_path,
            "session_html_path": job.get("session_html_path"),
        }

    def cancel_job(self, job_id: str) -> Dict:
        """
        Cancel a running job.
        
        Args:
            job_id: Unique job identifier
        
        Returns:
            Dict with cancellation result
        
        Raises:
            ValueError: If job_id not found or job cannot be cancelled
        """
        job = self.get_job_status(job_id)
        
        if job["status"] in (STATUS_COMPLETED, STATUS_FAILED, STATUS_CANCELLED):
            raise ValueError(
                f"Cannot cancel job {job_id} - already {job['status']}"
            )
        
        # Cancel EMR step if it exists
        if job.get("step_id") and job.get("cluster_id"):
            try:
                region = os.getenv("AWS_REGION", "us-east-1")
                result = subprocess.run(
                    [
                        "aws",
                        "emr",
                        "cancel-steps",
                        "--cluster-id",
                        job["cluster_id"],
                        "--step-ids",
                        job["step_id"],
                        "--region",
                        region,
                    ],
                    capture_output=True,
                    text=True,
                    timeout=10,
                )
                
                if result.returncode != 0:
                    logger.warning(
                        f"Failed to cancel EMR step: {result.stderr}"
                    )
            except Exception as e:
                logger.warning(f"Exception cancelling EMR step: {e}")
        
        # Update job status in stored job dict
        self.jobs[job_id]["status"] = STATUS_CANCELLED
        self.jobs[job_id]["updated_at"] = datetime.now(timezone.utc).isoformat()
        
        return {
            "job_id": job_id,
            "status": STATUS_CANCELLED,
            "message": f"Job {job_id} cancelled",
        }

    def backfill_job_directories(self, job_ids: Optional[List[str]] = None) -> Dict[str, bool]:
        """Backfill/create local directories for existing jobs.
        
        Args:
            job_ids: Optional list of specific job IDs to backfill. If None, backfills all jobs.
        
        Returns:
            Dict mapping job_id to success status (True if directory created, False otherwise)
        """
        results = {}
        jobs_to_process = job_ids if job_ids else list(self.jobs.keys())
        
        for job_id in jobs_to_process:
            try:
                job = self.jobs.get(job_id)
                if not job:
                    results[job_id] = False
                    continue
                
                job_local_path = self._ensure_job_directory(job_id, job)
                if job_local_path:
                    # Save metadata to the directory
                    self._save_job_metadata_to_local(job_id)
                    results[job_id] = True
                    logger.info(f"✅ Backfilled directory for job {job_id}: {job_local_path}")
                else:
                    results[job_id] = False
                    logger.warning(f"❌ Failed to backfill directory for job {job_id}")
            except Exception as e:
                results[job_id] = False
                logger.error(f"Error backfilling directory for job {job_id}: {e}")
        
        return results
    
    def list_jobs(
        self, status: Optional[str] = None, limit: int = 100
    ) -> List[Dict]:
        """
        List all jobs, optionally filtered by status.
        
        Args:
            status: Optional status filter
            limit: Maximum number of jobs to return
        
        Returns:
            List of job metadata dicts
        """
        jobs = list(self.jobs.values())
        
        if status:
            if status not in VALID_STATUSES:
                raise ValueError(f"Invalid status: {status}")
            jobs = [j for j in jobs if j["status"] == status]
        
        # Sort by submitted_at descending (most recent first)
        jobs.sort(
            key=lambda x: x.get("submitted_at", ""), reverse=True
        )
        
        return jobs[:limit]

    def retry_job(self, job_id: str) -> str:
        """
        Retry a failed job with the same parameters.
        
        Args:
            job_id: Unique job identifier of the failed job
        
        Returns:
            new_job_id: Unique identifier for the new job
        
        Raises:
            ValueError: If job_id not found or job is not failed
        """
        original_job = self.get_job_status(job_id)
        
        if original_job["status"] != STATUS_FAILED:
            raise ValueError(
                f"Cannot retry job {job_id} - current status is {original_job['status']}, must be 'failed'"
            )
        
        # Extract parameters from original job
        r1_path = original_job.get("r1_path")
        r2_path = original_job.get("r2_path")
        output_path = original_job.get("output_path")
        session_id = original_job.get("session_id")
        
        if not r1_path or not r2_path:
            raise ValueError(
                f"Cannot retry job {job_id} - missing R1 or R2 path"
            )
        
        # Submit new job with same parameters
        new_job_id = self.submit_fastqc_job(
            r1_path=r1_path,
            r2_path=r2_path,
            output_path=output_path,
            session_id=session_id,
        )
        
        logger.info(
            f"Retried job {job_id} as new job {new_job_id} "
            f"(R1={r1_path}, R2={r2_path})"
        )
        
        return new_job_id

    def _check_cluster_state(self, cluster_id: str) -> Optional[str]:
        """
        Check the state of an EMR cluster.
        
        Args:
            cluster_id: EMR cluster ID
        
        Returns:
            Cluster state string or None if cluster doesn't exist
        """
        try:
            region = os.getenv("AWS_REGION", "us-east-1")
            result = subprocess.run(
                [
                    "aws",
                    "emr",
                    "describe-cluster",
                    "--cluster-id",
                    cluster_id,
                    "--region",
                    region,
                    "--query",
                    "Cluster.Status.State",
                    "--output",
                    "text",
                ],
                capture_output=True,
                text=True,
                timeout=10,
            )
            
            if result.returncode == 0:
                return result.stdout.strip()
            else:
                logger.warning(f"Failed to check cluster state: {result.stderr}")
                return None
        except Exception as e:
            logger.warning(f"Exception checking cluster state: {e}")
            return None
    
    def _find_active_cluster(self) -> Optional[str]:
        """
        Find an active EMR cluster (WAITING or RUNNING state).
        
        Returns:
            Active cluster ID or None if no active cluster found
        """
        try:
            region = os.getenv("AWS_REGION", "us-east-1")
            result = subprocess.run(
                [
                    "aws",
                    "emr",
                    "list-clusters",
                    "--region",
                    region,
                    "--active",
                    "--query",
                    "Clusters[?Status.State==`WAITING` || Status.State==`RUNNING`].[Id,Name,Status.State]",
                    "--output",
                    "text",
                ],
                capture_output=True,
                text=True,
                timeout=10,
            )
            
            if result.returncode == 0 and result.stdout.strip():
                # Get the first active cluster
                lines = result.stdout.strip().split('\n')
                if lines:
                    # Format: cluster_id cluster_name state
                    cluster_id = lines[0].split()[0] if lines[0].split() else None
                    if cluster_id:
                        logger.info(f"Found active cluster: {cluster_id}")
                        return cluster_id
            
            return None
        except Exception as e:
            logger.warning(f"Exception finding active cluster: {e}")
            return None
    
    def _ensure_iam_roles(self) -> Tuple[Optional[str], Optional[str]]:
        """
        Ensure EMR IAM roles exist, creating them if necessary.
        
        Returns:
            Tuple of (service_role_arn, instance_profile_arn) or (None, None) on error
        """
        try:
            region = os.getenv("AWS_REGION", "us-east-1")
            service_role_name = "EMR_DefaultRole"
            instance_profile_name = "EMR_EC2_DefaultRole"
            ec2_role_name = "EMR_EC2_DefaultRole"
            
            # Check if service role exists
            result = subprocess.run(
                [
                    "aws",
                    "iam",
                    "get-role",
                    "--role-name",
                    service_role_name,
                    "--query",
                    "Role.Arn",
                    "--output",
                    "text",
                ],
                capture_output=True,
                text=True,
                timeout=10,
            )
            
            if result.returncode == 0:
                service_role_arn = result.stdout.strip()
                logger.info(f"Using existing EMR service role: {service_role_arn}")
            else:
                logger.warning(
                    f"EMR service role {service_role_name} does not exist. "
                    f"Please run ./scripts/aws/setup-emr-cluster.sh to create it."
                )
                return (None, None)
            
            # Check if instance profile exists
            result = subprocess.run(
                [
                    "aws",
                    "iam",
                    "get-instance-profile",
                    "--instance-profile-name",
                    instance_profile_name,
                    "--query",
                    "InstanceProfile.Arn",
                    "--output",
                    "text",
                ],
                capture_output=True,
                text=True,
                timeout=10,
            )
            
            if result.returncode == 0:
                instance_profile_arn = result.stdout.strip()
                logger.info(f"Using existing EMR instance profile: {instance_profile_arn}")
            else:
                logger.warning(
                    f"EMR instance profile {instance_profile_name} does not exist. "
                    f"Please run ./scripts/aws/setup-emr-cluster.sh to create it."
                )
                return (None, None)
            
            return (service_role_arn, instance_profile_arn)
        except Exception as e:
            logger.error(f"Exception ensuring IAM roles: {e}")
            return (None, None)
    
    def _ensure_bootstrap_script(self) -> Optional[str]:
        """
        Ensure bootstrap script exists in S3, uploading if necessary.
        
        Returns:
            S3 path to bootstrap script or None on error
        """
        try:
            region = os.getenv("AWS_REGION", "us-east-1")
            s3_bucket = os.getenv("S3_DATASET_BUCKET", "noricum-ngs-data")
            bootstrap_s3_path = f"s3://{s3_bucket}/emr-bootstrap/fastqc-bootstrap.sh"
            
            # Check if bootstrap script exists in S3
            result = subprocess.run(
                [
                    "aws",
                    "s3",
                    "ls",
                    bootstrap_s3_path,
                    "--region",
                    region,
                ],
                capture_output=True,
                text=True,
                timeout=10,
            )
            
            if result.returncode == 0:
                logger.info(f"Bootstrap script exists at {bootstrap_s3_path}")
                return bootstrap_s3_path
            
            # Bootstrap script doesn't exist, create and upload it
            logger.info("Bootstrap script not found, creating and uploading...")
            
            bootstrap_content = """#!/bin/bash
# Bootstrap script for EMR cluster
# Installs Python dependencies for FASTQ processing

# Update system
sudo yum update -y

# Install additional Python packages
sudo pip3 install biopython pandas numpy matplotlib seaborn plotly

# Install FASTQ parsing libraries
sudo pip3 install pyspark[fastq] 2>/dev/null || echo "pyspark already installed"

# Create directories for logs
sudo mkdir -p /mnt/var/log/fastqc
sudo chmod 777 /mnt/var/log/fastqc

echo "Bootstrap script completed successfully"
"""
            
            # Create temporary file
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as tmp_file:
                tmp_file.write(bootstrap_content)
                tmp_file_path = tmp_file.name
            
            try:
                # Upload to S3
                result = subprocess.run(
                    [
                        "aws",
                        "s3",
                        "cp",
                        tmp_file_path,
                        bootstrap_s3_path,
                        "--region",
                        region,
                    ],
                    capture_output=True,
                    text=True,
                    timeout=30,
                )
                
                if result.returncode == 0:
                    logger.info(f"✅ Uploaded bootstrap script to {bootstrap_s3_path}")
                    return bootstrap_s3_path
                else:
                    logger.error(f"Failed to upload bootstrap script: {result.stderr}")
                    return None
            finally:
                # Clean up temp file
                try:
                    os.unlink(tmp_file_path)
                except:
                    pass
                    
        except Exception as e:
            logger.error(f"Exception ensuring bootstrap script: {e}")
            return None
    
    def _create_emr_cluster(self) -> Optional[str]:
        """
        Create a new EMR cluster.
        
        Returns:
            Cluster ID or None on error
        """
        try:
            region = os.getenv("AWS_REGION", "us-east-1")
            s3_bucket = os.getenv("S3_DATASET_BUCKET", "noricum-ngs-data")
            cluster_name = os.getenv("EMR_CLUSTER_NAME", "helix-ai-fastqc")
            instance_type = os.getenv("EMR_INSTANCE_TYPE", "m5.xlarge")
            instance_count = int(os.getenv("EMR_INSTANCE_COUNT", "3"))
            release_label = "emr-6.15.0"
            
            logger.info(f"Creating new EMR cluster: {cluster_name}")
            
            # Ensure IAM roles exist
            service_role_arn, instance_profile_arn = self._ensure_iam_roles()
            if not service_role_arn or not instance_profile_arn:
                raise RuntimeError(
                    "EMR IAM roles not found. Please run ./scripts/aws/setup-emr-cluster.sh "
                    "to create the required IAM roles."
                )
            
            # Ensure bootstrap script exists
            bootstrap_script = self._ensure_bootstrap_script()
            if not bootstrap_script:
                raise RuntimeError("Failed to ensure bootstrap script exists in S3")
            
            # Create configurations JSON
            import tempfile
            import json
            
            configurations = [
                {
                    "Classification": "spark-defaults",
                    "Properties": {
                        "spark.sql.adaptive.enabled": "true",
                        "spark.sql.adaptive.coalescePartitions.enabled": "true",
                        "spark.serializer": "org.apache.spark.serializer.KryoSerializer",
                        "spark.sql.execution.arrow.pyspark.enabled": "true"
                    }
                },
                {
                    "Classification": "spark-env",
                    "Configurations": [
                        {
                            "Classification": "export",
                            "Properties": {
                                "PYSPARK_PYTHON": "/usr/bin/python3",
                                "PYSPARK_DRIVER_PYTHON": "/usr/bin/python3"
                            }
                        }
                    ]
                }
            ]
            
            instance_groups = [
                {
                    "Name": "Master",
                    "InstanceGroupType": "MASTER",
                    "InstanceType": instance_type,
                    "InstanceCount": 1
                },
                {
                    "Name": "Core",
                    "InstanceGroupType": "CORE",
                    "InstanceType": instance_type,
                    "InstanceCount": instance_count - 1,
                    "BidPrice": "0.10"  # Spot instance pricing
                }
            ]
            
            # Create temporary files for JSON configs
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as config_file:
                json.dump(configurations, config_file)
                config_file_path = config_file.name
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as groups_file:
                json.dump(instance_groups, groups_file)
                groups_file_path = groups_file.name
            
            try:
                # Build AWS CLI command
                # Note: --applications needs to be formatted as separate Name=... arguments
                # --bootstrap-actions needs Path= and Name= in the same string
                cmd = [
                    "aws",
                    "emr",
                    "create-cluster",
                    "--name", cluster_name,
                    "--release-label", release_label,
                    "--applications", "Name=Spark", "Name=Hadoop", "Name=Zeppelin",
                    "--instance-groups", f"file://{groups_file_path}",
                    "--bootstrap-actions", f"Path={bootstrap_script},Name=Install Python Dependencies",
                    "--configurations", f"file://{config_file_path}",
                    "--service-role", service_role_arn,
                    "--ec2-attributes", f"InstanceProfile={instance_profile_arn}",
                    "--log-uri", f"s3://{s3_bucket}/emr-logs/",
                    "--region", region,
                    "--tags", "Project=Helix.AI", "Purpose=FASTQ Quality Analysis",
                    "--query", "ClusterId",
                    "--output", "text",
                ]
                
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=60,
                )
                
                if result.returncode == 0:
                    cluster_id = result.stdout.strip()
                    logger.info(f"✅ EMR cluster created: {cluster_id}")
                    return cluster_id
                else:
                    logger.error(f"Failed to create EMR cluster: {result.stderr}")
                    return None
            finally:
                # Clean up temp files
                try:
                    os.unlink(config_file_path)
                    os.unlink(groups_file_path)
                except:
                    pass
                    
        except Exception as e:
            logger.error(f"Exception creating EMR cluster: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def _wait_for_cluster_ready(self, cluster_id: str, max_wait_minutes: int = 15) -> bool:
        """
        Wait for EMR cluster to be in WAITING or RUNNING state.
        
        Args:
            cluster_id: EMR cluster ID
            max_wait_minutes: Maximum time to wait in minutes
        
        Returns:
            True if cluster is ready, False if timeout or error
        """
        import time
        
        region = os.getenv("AWS_REGION", "us-east-1")
        max_wait_seconds = max_wait_minutes * 60
        check_interval = 30  # Check every 30 seconds
        start_time = time.time()
        
        logger.info(f"Waiting for cluster {cluster_id} to be ready (max {max_wait_minutes} minutes)...")
        
        while time.time() - start_time < max_wait_seconds:
            state = self._check_cluster_state(cluster_id)
            
            if state in ("WAITING", "RUNNING"):
                logger.info(f"✅ Cluster {cluster_id} is ready (state: {state})")
                return True
            elif state in ("TERMINATED", "TERMINATED_WITH_ERRORS", "TERMINATING"):
                logger.error(f"Cluster {cluster_id} is in terminal state: {state}")
                return False
            
            # Log progress every 2 minutes
            elapsed_minutes = int((time.time() - start_time) / 60)
            if elapsed_minutes > 0 and elapsed_minutes % 2 == 0:
                logger.info(f"Cluster {cluster_id} still starting... (state: {state}, elapsed: {elapsed_minutes} minutes)")
            
            time.sleep(check_interval)
        
        logger.warning(f"Timeout waiting for cluster {cluster_id} to be ready")
        return False

    def get_job_logs(self, job_id: str) -> str:
        """
        Get EMR logs for a job.
        
        Args:
            job_id: Unique job identifier
        
        Returns:
            Log content as string
        
        Raises:
            ValueError: If job_id not found
        """
        job = self.get_job_status(job_id)
        
        cluster_id = job.get("cluster_id")
        step_id = job.get("step_id")
        
        if not cluster_id or not step_id:
            return "Logs not available - job does not have cluster_id or step_id"
        
        try:
            region = os.getenv("AWS_REGION", "us-east-1")
            s3_bucket = os.getenv("S3_DATASET_BUCKET", "noricum-ngs-data")
            log_path = f"s3://{s3_bucket}/emr-logs/{cluster_id}/steps/{step_id}/wrapper.log"
            
            # Try to download wrapper.log from S3
            result = subprocess.run(
                [
                    "aws",
                    "s3",
                    "cp",
                    log_path,
                    "-",
                    "--region",
                    region,
                ],
                capture_output=True,
                text=True,
                timeout=15,
            )
            
            if result.returncode == 0 and result.stdout:
                return result.stdout
            else:
                # Try to get stderr or stdout logs
                stderr_path = f"s3://{s3_bucket}/emr-logs/{cluster_id}/steps/{step_id}/stderr.gz"
                stdout_path = f"s3://{s3_bucket}/emr-logs/{cluster_id}/steps/{step_id}/stdout.gz"
                
                # Try stderr first
                stderr_result = subprocess.run(
                    [
                        "aws",
                        "s3",
                        "cp",
                        stderr_path,
                        "-",
                        "--region",
                        region,
                    ],
                    capture_output=True,
                    text=True,
                    timeout=15,
                )
                
                if stderr_result.returncode == 0 and stderr_result.stdout:
                    return f"=== STDERR ===\n{stderr_result.stdout}"
                
                # Try stdout
                stdout_result = subprocess.run(
                    [
                        "aws",
                        "s3",
                        "cp",
                        stdout_path,
                        "-",
                        "--region",
                        region,
                    ],
                    capture_output=True,
                    text=True,
                    timeout=15,
                )
                
                if stdout_result.returncode == 0 and stdout_result.stdout:
                    return f"=== STDOUT ===\n{stdout_result.stdout}"
                
                return (
                    f"Logs not yet available in S3.\n"
                    f"Expected location: {log_path}\n"
                    f"Logs may still be uploading to S3."
                )
        except subprocess.TimeoutExpired:
            return "Timeout fetching logs from S3. Please try again later."
        except Exception as e:
            logger.warning(f"Failed to fetch logs for job {job_id}: {e}")
            return f"Error fetching logs: {str(e)}"


# Global singleton instance
_job_manager_instance: Optional[JobManager] = None


def get_job_manager() -> JobManager:
    """Get the global JobManager instance (singleton pattern)."""
    global _job_manager_instance
    if _job_manager_instance is None:
        _job_manager_instance = JobManager()
    return _job_manager_instance

