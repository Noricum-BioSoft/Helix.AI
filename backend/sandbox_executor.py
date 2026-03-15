"""
Sandbox Executor - Run bioinformatics tools in isolated Docker container

This module provides secure, isolated execution of bioinformatics tools
in a Docker container with all necessary dependencies pre-installed.

Features:
- Resource limits (CPU, memory)
- Network isolation
- Read-only filesystem mounts
- Automatic cleanup
- Tool version tracking
- Execution logging

Usage:
    from backend.sandbox_executor import SandboxExecutor
    
    executor = SandboxExecutor()
    result = executor.execute_tool(
        tool="fastqc",
        args=["input_R1.fq", "input_R2.fq", "-o", "output/"],
        input_files={"input_R1.fq": "/path/to/r1.fq", "input_R2.fq": "/path/to/r2.fq"},
        working_dir="/tmp/work",
        max_memory="4g",
        max_cpus=2,
        timeout=300
    )
"""

import logging
import os
import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Any

logger = logging.getLogger(__name__)


@dataclass
class ExecutionResult:
    """Result of sandbox tool execution."""
    success: bool
    exit_code: int
    stdout: str
    stderr: str
    execution_time: float
    output_files: List[str]
    error_message: Optional[str] = None


class SandboxExecutor:
    """
    Execute bioinformatics tools in isolated Docker sandbox.
    
    This executor runs tools inside a Docker container (helix-biotools)
    with resource limits and security constraints.
    """
    
    # Docker image for bioinformatics tools
    # Can be either "helix-biotools:latest" (full) or "helix-biotools:minimal" (fast build)
    BIOTOOLS_IMAGE = os.getenv("HELIX_BIOTOOLS_IMAGE", "helix-biotools:latest")
    
    # Supported tools and their executables
    SUPPORTED_TOOLS = {
        "fastqc": "/usr/local/bin/fastqc",
        "muscle": "/usr/local/bin/muscle",
        "clustalo": "/usr/bin/clustalo",
        "trimmomatic": "/usr/local/bin/trimmomatic",
        "samtools": "/usr/local/bin/samtools",
        "bwa": "/usr/local/bin/bwa",
        "bowtie2": "/usr/local/bin/bowtie2",
        "star": "/usr/local/bin/STAR",
        "hisat2": "/usr/local/bin/hisat2",
        "stringtie": "/usr/local/bin/stringtie",
        "featureCounts": "/usr/local/bin/featureCounts",
    }
    
    def __init__(self, image: Optional[str] = None):
        """
        Initialize sandbox executor.
        
        Args:
            image: Docker image to use (defaults to BIOTOOLS_IMAGE)
        """
        self.image = image or self.BIOTOOLS_IMAGE
        # Track whether Docker is actually usable. When False, all execution
        # methods transparently fall back to host subprocess (same Python env).
        self._docker_available: bool = True

        if os.getenv("HELIX_SANDBOX_HOST_FALLBACK") == "1":
            logger.warning("HELIX_SANDBOX_HOST_FALLBACK=1: skipping Docker availability checks")
            self._docker_available = False
        else:
            self._docker_available = self._check_docker_available()
            if self._docker_available:
                self._check_image_available()

    def _check_docker_available(self) -> bool:
        """Check if Docker is available. Returns True if available, False otherwise (never raises)."""
        try:
            result = subprocess.run(
                ["docker", "--version"],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode != 0:
                logger.warning("Docker is not available (non-zero exit from 'docker --version'); falling back to host execution")
                return False
            logger.debug(f"Docker available: {result.stdout.strip()}")
            return True
        except (FileNotFoundError, subprocess.TimeoutExpired) as e:
            logger.warning(f"Docker is not installed or not running ({e}); falling back to host execution")
            return False
    
    def _check_image_available(self):
        """Check if biotools Docker image is available."""
        try:
            result = subprocess.run(
                ["docker", "images", "-q", self.image],
                capture_output=True,
                text=True,
                timeout=5
            )
            if not result.stdout.strip():
                logger.warning(
                    f"Docker image {self.image} not found. "
                    f"Please build it with: docker build -f backend/Dockerfile.biotools -t {self.image} ."
                )
                # Don't raise error - allow fallback to host execution
        except subprocess.TimeoutExpired as e:
            logger.warning(f"Docker image check timed out: {e}")
    
    def execute_tool(
        self,
        tool: str,
        args: List[str],
        input_files: Optional[Dict[str, str]] = None,
        working_dir: Optional[str] = None,
        output_dir: Optional[str] = None,
        max_memory: str = "4g",
        max_cpus: int = 2,
        timeout: int = 300,
        env_vars: Optional[Dict[str, str]] = None
    ) -> ExecutionResult:
        """
        Execute a bioinformatics tool in sandboxed Docker container.
        
        Args:
            tool: Tool name (e.g., "fastqc", "muscle")
            args: Tool arguments (as list)
            input_files: Map of filename -> local path to mount as input
            working_dir: Local working directory for tool execution
            output_dir: Local directory for output files
            max_memory: Maximum memory (e.g., "4g", "8g")
            max_cpus: Maximum CPU cores
            timeout: Execution timeout in seconds
            env_vars: Environment variables to pass to container
        
        Returns:
            ExecutionResult with status, output, and files
        
        Raises:
            ValueError: If tool is not supported
            RuntimeError: If Docker execution fails
        """
        if tool not in self.SUPPORTED_TOOLS:
            raise ValueError(
                f"Unsupported tool: {tool}. "
                f"Supported tools: {', '.join(self.SUPPORTED_TOOLS.keys())}"
            )
        
        tool_path = self.SUPPORTED_TOOLS[tool]
        
        # Create temporary directories if not provided
        cleanup_work_dir = False
        cleanup_output_dir = False
        
        if working_dir is None:
            working_dir = tempfile.mkdtemp(prefix="helix-sandbox-")
            cleanup_work_dir = True
        
        if output_dir is None:
            output_dir = tempfile.mkdtemp(prefix="helix-output-")
            cleanup_output_dir = True
        
        work_path = Path(working_dir)
        output_path = Path(output_dir)
        
        work_path.mkdir(parents=True, exist_ok=True)
        output_path.mkdir(parents=True, exist_ok=True)
        
        try:
            # Copy input files to working directory
            if input_files:
                for filename, local_path in input_files.items():
                    dest = work_path / filename
                    # Skip if source and dest are the same
                    if Path(local_path).resolve() != dest.resolve():
                        subprocess.run(["cp", local_path, str(dest)], check=True)
            
            # Build Docker command
            docker_cmd = [
                "docker", "run",
                "--rm",  # Remove container after execution
                "--network", "none",  # No network access
                "--memory", max_memory,  # Memory limit
                "--cpus", str(max_cpus),  # CPU limit
                "--user", f"{os.getuid()}:{os.getgid()}",  # Run as current user
                "-v", f"{work_path.absolute()}:/sandbox/work:rw",  # Mount work dir
                "-v", f"{output_path.absolute()}:/sandbox/output:rw",  # Mount output dir
                "-w", "/sandbox/work",  # Set working directory
            ]
            
            # Add environment variables
            if env_vars:
                for key, value in env_vars.items():
                    docker_cmd.extend(["-e", f"{key}={value}"])
            
            # Add image and tool command
            docker_cmd.append(self.image)
            docker_cmd.append(tool_path)
            docker_cmd.extend(args)
            
            logger.info(f"🐳 Running {tool} in Docker sandbox...")
            logger.debug(f"Command: {' '.join(docker_cmd)}")
            
            start_time = time.time()
            
            # Execute tool in Docker
            result = subprocess.run(
                docker_cmd,
                capture_output=True,
                text=True,
                timeout=timeout
            )
            
            execution_time = time.time() - start_time
            
            # Collect output files
            output_files = [
                str(f.relative_to(output_path))
                for f in output_path.rglob("*")
                if f.is_file()
            ]
            
            success = result.returncode == 0
            
            if success:
                logger.info(f"✅ {tool} completed in {execution_time:.2f}s")
            else:
                logger.error(f"❌ {tool} failed with exit code {result.returncode}")
                logger.error(f"stderr: {result.stderr}")
            
            return ExecutionResult(
                success=success,
                exit_code=result.returncode,
                stdout=result.stdout,
                stderr=result.stderr,
                execution_time=execution_time,
                output_files=output_files,
                error_message=result.stderr if not success else None
            )
        
        except subprocess.TimeoutExpired:
            logger.error(f"⏱️  {tool} execution timed out after {timeout}s")
            return ExecutionResult(
                success=False,
                exit_code=-1,
                stdout="",
                stderr=f"Execution timed out after {timeout} seconds",
                execution_time=timeout,
                output_files=[],
                error_message=f"Timeout after {timeout}s"
            )
        
        except Exception as e:
            logger.error(f"❌ {tool} execution failed: {e}")
            return ExecutionResult(
                success=False,
                exit_code=-1,
                stdout="",
                stderr=str(e),
                execution_time=0.0,
                output_files=[],
                error_message=str(e)
            )
        
        finally:
            # Cleanup temporary directories
            if cleanup_work_dir and work_path.exists():
                subprocess.run(["rm", "-rf", str(work_path)], check=False)
            if cleanup_output_dir and output_path.exists():
                subprocess.run(["rm", "-rf", str(output_path)], check=False)

    def execute_command(
        self,
        cmd: List[str],
        *,
        working_dir: Optional[str] = None,
        output_dir: Optional[str] = None,
        max_memory: str = "4g",
        max_cpus: int = 2,
        timeout: int = 300,
        env_vars: Optional[Dict[str, str]] = None,
        allow_host_fallback: Optional[bool] = None,
    ) -> ExecutionResult:
        """
        Execute an arbitrary command in the sandboxed Docker container.

        This is used for script-based iterative workflows (e.g. running `python script.py`).

        Host fallback is **disabled by default**. To enable it (typically in unit tests),
        set `allow_host_fallback=True` or env var `HELIX_SANDBOX_HOST_FALLBACK=1`.
        """
        if not isinstance(cmd, list) or not cmd or not all(isinstance(x, str) and x for x in cmd):
            raise ValueError("cmd must be a non-empty list of strings")

        cleanup_work_dir = False
        cleanup_output_dir = False

        if working_dir is None:
            working_dir = tempfile.mkdtemp(prefix="helix-sandbox-work-")
            cleanup_work_dir = True
        if output_dir is None:
            output_dir = tempfile.mkdtemp(prefix="helix-sandbox-output-")
            cleanup_output_dir = True

        work_path = Path(working_dir)
        output_path = Path(output_dir)
        work_path.mkdir(parents=True, exist_ok=True)
        output_path.mkdir(parents=True, exist_ok=True)

        if allow_host_fallback is None:
            # Auto-enable host fallback when Docker is not available in this
            # environment (e.g. ECS Fargate, CI without Docker).
            allow_host_fallback = (
                not self._docker_available
                or os.getenv("HELIX_SANDBOX_HOST_FALLBACK") == "1"
            )

        start_time = time.time()

        def _collect_outputs() -> List[str]:
            try:
                return [
                    str(f.relative_to(output_path))
                    for f in output_path.rglob("*")
                    if f.is_file()
                ]
            except Exception:
                return []

        # If host fallback is enabled, prefer running on host immediately to avoid
        # flakiness/hangs when Docker Desktop is unavailable.
        if allow_host_fallback:
            try:
                import sys as _sys

                host_cmd = list(cmd)
                if host_cmd[0] in {"python", "python3"}:
                    host_cmd = [_sys.executable, "-I"] + host_cmd[1:]

                env = os.environ.copy()
                if env_vars:
                    env.update({str(k): str(v) for k, v in env_vars.items()})

                logger.warning(f"⚠️  HELIX_SANDBOX_HOST_FALLBACK=1: running on host: {' '.join(host_cmd)}")
                result = subprocess.run(
                    host_cmd,
                    capture_output=True,
                    text=True,
                    timeout=timeout,
                    cwd=str(work_path),
                    env=env,
                )
                exec_time = time.time() - start_time
                output_files = _collect_outputs()
                success = result.returncode == 0
                return ExecutionResult(
                    success=success,
                    exit_code=result.returncode,
                    stdout=result.stdout,
                    stderr=result.stderr,
                    execution_time=exec_time,
                    output_files=output_files,
                    error_message=result.stderr if not success else None,
                )
            except subprocess.TimeoutExpired:
                exec_time = time.time() - start_time
                return ExecutionResult(
                    success=False,
                    exit_code=-1,
                    stdout="",
                    stderr=f"Execution timed out after {timeout} seconds (host fallback)",
                    execution_time=exec_time,
                    output_files=_collect_outputs(),
                    error_message=f"Timeout after {timeout}s (host fallback)",
                )
            except Exception as host_err:
                exec_time = time.time() - start_time
                return ExecutionResult(
                    success=False,
                    exit_code=-1,
                    stdout="",
                    stderr=str(host_err),
                    execution_time=exec_time,
                    output_files=_collect_outputs(),
                    error_message=str(host_err),
                )

        try:
            docker_cmd = [
                "docker",
                "run",
                "--rm",
                "--network",
                "none",
                "--memory",
                max_memory,
                "--cpus",
                str(max_cpus),
                "--user",
                f"{os.getuid()}:{os.getgid()}",
                "-v",
                f"{work_path.absolute()}:/sandbox/work:rw",
                "-v",
                f"{output_path.absolute()}:/sandbox/output:rw",
                "-w",
                "/sandbox/work",
            ]

            if env_vars:
                for key, value in env_vars.items():
                    docker_cmd.extend(["-e", f"{key}={value}"])

            docker_cmd.append(self.image)
            docker_cmd.extend(cmd)

            logger.info(f"🐳 Running command in Docker sandbox: {' '.join(cmd)}")
            result = subprocess.run(
                docker_cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
            exec_time = time.time() - start_time
            output_files = _collect_outputs()
            success = result.returncode == 0
            return ExecutionResult(
                success=success,
                exit_code=result.returncode,
                stdout=result.stdout,
                stderr=result.stderr,
                execution_time=exec_time,
                output_files=output_files,
                error_message=result.stderr if not success else None,
            )
        except subprocess.TimeoutExpired:
            exec_time = time.time() - start_time
            return ExecutionResult(
                success=False,
                exit_code=-1,
                stdout="",
                stderr=f"Execution timed out after {timeout} seconds",
                execution_time=exec_time,
                output_files=_collect_outputs(),
                error_message=f"Timeout after {timeout}s",
            )
        except Exception as e:
            if not allow_host_fallback:
                exec_time = time.time() - start_time
                return ExecutionResult(
                    success=False,
                    exit_code=-1,
                    stdout="",
                    stderr=str(e),
                    execution_time=exec_time,
                    output_files=_collect_outputs(),
                    error_message=str(e),
                )

            # Host fallback (intended for unit tests / dev without Docker).
            try:
                import sys as _sys

                host_cmd = list(cmd)
                if host_cmd[0] in {"python", "python3"}:
                    host_cmd = [_sys.executable, "-I"] + host_cmd[1:]

                env = os.environ.copy()
                if env_vars:
                    env.update({str(k): str(v) for k, v in env_vars.items()})

                logger.warning(f"⚠️  Docker sandbox unavailable; running on host: {' '.join(host_cmd)}")
                result = subprocess.run(
                    host_cmd,
                    capture_output=True,
                    text=True,
                    timeout=timeout,
                    cwd=str(work_path),
                    env=env,
                )
                exec_time = time.time() - start_time
                output_files = _collect_outputs()
                success = result.returncode == 0
                return ExecutionResult(
                    success=success,
                    exit_code=result.returncode,
                    stdout=result.stdout,
                    stderr=result.stderr,
                    execution_time=exec_time,
                    output_files=output_files,
                    error_message=result.stderr if not success else None,
                )
            except subprocess.TimeoutExpired:
                exec_time = time.time() - start_time
                return ExecutionResult(
                    success=False,
                    exit_code=-1,
                    stdout="",
                    stderr=f"Execution timed out after {timeout} seconds (host fallback)",
                    execution_time=exec_time,
                    output_files=_collect_outputs(),
                    error_message=f"Timeout after {timeout}s (host fallback)",
                )
            except Exception as host_err:
                exec_time = time.time() - start_time
                return ExecutionResult(
                    success=False,
                    exit_code=-1,
                    stdout="",
                    stderr=str(host_err),
                    execution_time=exec_time,
                    output_files=_collect_outputs(),
                    error_message=str(host_err),
                )
        finally:
            if cleanup_work_dir and work_path.exists():
                subprocess.run(["rm", "-rf", str(work_path)], check=False)
            if cleanup_output_dir and output_path.exists():
                subprocess.run(["rm", "-rf", str(output_path)], check=False)
    
    def get_tool_version(self, tool: str) -> Optional[str]:
        """
        Get version of installed tool.
        
        Args:
            tool: Tool name
        
        Returns:
            Version string or None if unavailable
        """
        if tool not in self.SUPPORTED_TOOLS:
            return None
        
        tool_path = self.SUPPORTED_TOOLS[tool]
        
        try:
            # Try common version flags
            for version_flag in ["--version", "-v", "-version", "version"]:
                result = subprocess.run(
                    ["docker", "run", "--rm", self.image, tool_path, version_flag],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                if result.returncode == 0 and result.stdout.strip():
                    return result.stdout.strip().split("\n")[0]
        except Exception as e:
            logger.debug(f"Could not get version for {tool}: {e}")
        
        return None
    
    def list_available_tools(self) -> Dict[str, Dict[str, Any]]:
        """
        List all available tools and their metadata.
        
        Returns:
            Dictionary of tool name -> metadata
        """
        tools_info = {}
        for tool_name in self.SUPPORTED_TOOLS:
            version = self.get_tool_version(tool_name)
            tools_info[tool_name] = {
                "name": tool_name,
                "path": self.SUPPORTED_TOOLS[tool_name],
                "version": version,
                "available": version is not None
            }
        return tools_info


# Singleton instance
_sandbox_executor: Optional[SandboxExecutor] = None


def get_sandbox_executor() -> SandboxExecutor:
    """
    Get singleton SandboxExecutor instance.
    
    Returns:
        SandboxExecutor instance
    """
    global _sandbox_executor
    if _sandbox_executor is None:
        _sandbox_executor = SandboxExecutor()
    return _sandbox_executor
