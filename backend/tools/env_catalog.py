"""
Environment Capabilities Catalog - Static configuration of execution environments.

This tool provides factual information about each infrastructure environment
(Local, EC2, EMR, Batch, Lambda) to inform infrastructure decisions.

Key principle: Static, versioned configuration (not LLM-generated).
All values are documented with sources and assumptions.
"""

import os
import yaml
import logging
from pathlib import Path
from typing import Dict, Any, List, Optional, Literal
from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)


class EnvironmentCapability(BaseModel):
    """Capabilities of a single execution environment."""
    
    name: Literal["Local", "EC2", "EMR", "Batch", "Lambda"] = Field(
        ...,
        description="Environment name"
    )
    
    # Startup characteristics
    startup_time_seconds: float = Field(
        ...,
        ge=0.0,
        description="Typical cold start time in seconds"
    )
    
    startup_overhead: Literal["None", "Low", "Medium", "High"] = Field(
        ...,
        description="Startup overhead class"
    )
    
    # Compute capabilities
    max_cpu_cores: Optional[int] = Field(
        default=None,
        description="Max CPU cores (None = unlimited/depends on instance type)"
    )
    
    max_memory_gb: Optional[float] = Field(
        default=None,
        description="Max memory in GB (None = unlimited/depends on instance type)"
    )
    
    max_runtime_minutes: Optional[float] = Field(
        default=None,
        description="Max runtime in minutes (None = unlimited)"
    )
    
    supports_distributed: bool = Field(
        default=False,
        description="Supports distributed/parallel processing"
    )
    
    supports_gpu: bool = Field(
        default=False,
        description="GPU support available"
    )
    
    # Tool availability
    pre_installed_tools: List[str] = Field(
        default_factory=list,
        description="List of pre-installed bioinformatics tools"
    )
    
    container_support: bool = Field(
        default=False,
        description="Supports containerized execution"
    )
    
    # Data access
    s3_native: bool = Field(
        default=False,
        description="Can process S3 data in-place (no download)"
    )
    
    data_transfer_required: bool = Field(
        default=True,
        description="Requires data transfer from S3"
    )
    
    # Operational characteristics
    reproducibility: Literal["Low", "Medium", "High"] = Field(
        ...,
        description="Reproducibility guarantee"
    )
    
    debuggability: Literal["Low", "Medium", "High"] = Field(
        ...,
        description="Ease of debugging"
    )
    
    failure_blast_radius: Literal["Low", "Medium", "High"] = Field(
        ...,
        description="Impact of failure"
    )
    
    # Cost characteristics
    cost_class: Literal["Free", "Low", "Medium", "High"] = Field(
        ...,
        description="Relative cost class"
    )
    
    cost_model: str = Field(
        ...,
        description="Cost model description"
    )
    
    # Availability
    available: bool = Field(
        default=True,
        description="Whether this environment is currently available"
    )
    
    availability_notes: Optional[str] = Field(
        default=None,
        description="Notes on availability (e.g., requires HELIX_USE_EC2=true)"
    )
    
    # Metadata
    description: str = Field(
        ...,
        description="Human-readable description"
    )
    
    best_for: List[str] = Field(
        default_factory=list,
        description="Use cases this environment is best suited for"
    )
    
    avoid_for: List[str] = Field(
        default_factory=list,
        description="Use cases to avoid this environment for"
    )


class EnvironmentCapabilityCatalog:
    """
    Static catalog of execution environment capabilities.
    
    Loads from YAML configuration file and provides lookup/filtering methods.
    
    Usage:
        catalog = EnvironmentCapabilityCatalog()
        ec2_caps = catalog.get_environment("EC2")
        available = catalog.get_available_environments()
    """
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize catalog from YAML config.
        
        Args:
            config_path: Path to YAML config file (uses default if None)
        """
        if config_path is None:
            # Use default config from backend/config/
            project_root = Path(__file__).resolve().parent.parent.parent
            config_path = project_root / "backend" / "config" / "environment_capabilities.yaml"
        
        self.config_path = Path(config_path)
        self.environments: Dict[str, EnvironmentCapability] = {}
        
        # Load config
        self._load_config()
    
    def _load_config(self):
        """Load environment capabilities from YAML config."""
        if not self.config_path.exists():
            logger.warning(f"Environment catalog not found at {self.config_path}, using defaults")
            self._load_defaults()
            return
        
        try:
            with open(self.config_path, "r") as f:
                data = yaml.safe_load(f)
            
            if not data or "environments" not in data:
                logger.warning("Invalid catalog format, using defaults")
                self._load_defaults()
                return
            
            for env_data in data["environments"]:
                env = EnvironmentCapability(**env_data)
                self.environments[env.name] = env
            
            logger.info(f"Loaded {len(self.environments)} environments from {self.config_path}")
        
        except Exception as e:
            logger.error(f"Failed to load environment catalog: {e}")
            self._load_defaults()
    
    def _load_defaults(self):
        """Load default environment capabilities (fallback)."""
        # Check if EC2 is enabled
        use_ec2 = os.getenv("HELIX_USE_EC2", "false").lower() == "true"
        
        self.environments = {
            "Local": EnvironmentCapability(
                name="Local",
                startup_time_seconds=0.0,
                startup_overhead="None",
                max_cpu_cores=None,  # Depends on machine
                max_memory_gb=None,  # Depends on machine
                max_runtime_minutes=None,
                supports_distributed=False,
                supports_gpu=False,
                pre_installed_tools=[],
                container_support=False,
                s3_native=False,
                data_transfer_required=True,
                reproducibility="Low",
                debuggability="High",
                failure_blast_radius="Low",
                cost_class="Free",
                cost_model="No cloud costs",
                available=True,
                description="Local execution on current machine",
                best_for=["Small files (<100MB)", "Quick prototyping", "Development"],
                avoid_for=["Large files (>1GB)", "Distributed operations", "Production workloads"]
            ),
            "EC2": EnvironmentCapability(
                name="EC2",
                startup_time_seconds=2.0,
                startup_overhead="Low",
                max_cpu_cores=None,  # Depends on instance type
                max_memory_gb=None,  # Depends on instance type
                max_runtime_minutes=None,
                supports_distributed=False,
                supports_gpu=True,
                pre_installed_tools=["bbtools", "samtools", "bcftools", "bwa", "bowtie2", "fastqc"],
                container_support=False,
                s3_native=False,
                data_transfer_required=True,
                reproducibility="Medium",
                debuggability="High",
                failure_blast_radius="Low",
                cost_class="Low",
                cost_model="EC2 instance hours (t3.medium ~$0.04/hr, m5.xlarge ~$0.19/hr)",
                available=use_ec2,
                availability_notes="Requires HELIX_USE_EC2=true",
                description="EC2 instance with pre-installed bioinformatics tools",
                best_for=["Medium files (10-100MB)", "Tools requiring installation", "SSH access for debugging"],
                avoid_for=["Very large files (>10GB)", "Distributed operations", "Long-running jobs (>1hr)"]
            ),
            "EMR": EnvironmentCapability(
                name="EMR",
                startup_time_seconds=180.0,  # 3min cluster startup
                startup_overhead="High",
                max_cpu_cores=None,  # Scalable
                max_memory_gb=None,  # Scalable
                max_runtime_minutes=None,
                supports_distributed=True,
                supports_gpu=False,
                pre_installed_tools=[],
                container_support=True,
                s3_native=True,
                data_transfer_required=False,
                reproducibility="High",
                debuggability="Medium",
                failure_blast_radius="Medium",
                cost_class="Medium",
                cost_model="EMR cluster hours (m5.xlarge ~$0.27/hr with EMR markup)",
                available=True,
                description="AWS EMR for distributed processing of large S3 datasets",
                best_for=["Large S3 files (>100MB)", "Distributed operations", "Spark workflows"],
                avoid_for=["Small files (<100MB)", "Local files", "Quick prototyping"]
            ),
            "Batch": EnvironmentCapability(
                name="Batch",
                startup_time_seconds=60.0,  # Container pull + startup
                startup_overhead="Medium",
                max_cpu_cores=None,  # Depends on job definition
                max_memory_gb=None,  # Depends on job definition
                max_runtime_minutes=None,
                supports_distributed=False,
                supports_gpu=True,
                pre_installed_tools=[],
                container_support=True,
                s3_native=False,
                data_transfer_required=True,
                reproducibility="High",
                debuggability="Medium",
                failure_blast_radius="Low",
                cost_class="Low",
                cost_model="EC2 instance hours (pay for underlying compute)",
                available=False,  # Not fully implemented yet
                availability_notes="Batch routing not yet implemented",
                description="AWS Batch for containerized, medium-sized jobs",
                best_for=["Containerized tools", "Medium jobs (10-100MB)", "Reproducible workflows"],
                avoid_for=["Very large files (>10GB)", "Quick prototyping", "Tools requiring custom setup"]
            ),
            "Lambda": EnvironmentCapability(
                name="Lambda",
                startup_time_seconds=1.0,  # Cold start
                startup_overhead="Low",
                max_cpu_cores=6,
                max_memory_gb=10.0,
                max_runtime_minutes=15.0,  # Hard limit
                supports_distributed=False,
                supports_gpu=False,
                pre_installed_tools=[],
                container_support=True,
                s3_native=True,
                data_transfer_required=False,
                reproducibility="High",
                debuggability="Low",
                failure_blast_radius="Low",
                cost_class="Low",
                cost_model="Invocation + compute time (first 1M requests free)",
                available=False,  # Not implemented yet
                availability_notes="Lambda routing not yet implemented",
                description="AWS Lambda for lightweight, event-driven preprocessing",
                best_for=["Lightweight preprocessing (<15min)", "Event-driven workflows", "Serverless"],
                avoid_for=["Complex bioinformatics tools", "Long-running jobs (>15min)", "Large datasets"]
            ),
        }
    
    def get_environment(self, name: str) -> Optional[EnvironmentCapability]:
        """Get capabilities for a specific environment."""
        return self.environments.get(name)
    
    def get_available_environments(self) -> List[EnvironmentCapability]:
        """Get list of currently available environments."""
        return [env for env in self.environments.values() if env.available]
    
    def get_environments_for_use_case(
        self,
        file_size_mb: Optional[float] = None,
        supports_distributed: Optional[bool] = None,
        s3_native: Optional[bool] = None,
        max_runtime_minutes: Optional[float] = None
    ) -> List[EnvironmentCapability]:
        """
        Filter environments by use case requirements.
        
        Args:
            file_size_mb: Required file size support
            supports_distributed: Whether distributed processing is needed
            s3_native: Whether S3-native processing is needed
            max_runtime_minutes: Maximum acceptable runtime
        
        Returns:
            List of suitable environments
        """
        suitable = []
        
        for env in self.get_available_environments():
            # Check distributed requirement
            if supports_distributed is not None and supports_distributed and not env.supports_distributed:
                continue
            
            # Check S3 native requirement
            if s3_native is not None and s3_native and not env.s3_native:
                continue
            
            # Check runtime limit
            if max_runtime_minutes is not None and env.max_runtime_minutes is not None:
                if max_runtime_minutes > env.max_runtime_minutes:
                    continue
            
            suitable.append(env)
        
        return suitable
    
    def to_dict(self) -> Dict[str, Any]:
        """Export catalog to dictionary."""
        return {
            name: env.model_dump()
            for name, env in self.environments.items()
        }
    
    def to_yaml(self, filepath: str):
        """Export catalog to YAML file."""
        data = {
            "environments": [
                env.model_dump()
                for env in self.environments.values()
            ]
        }
        
        with open(filepath, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)
        
        logger.info(f"Exported environment catalog to {filepath}")
