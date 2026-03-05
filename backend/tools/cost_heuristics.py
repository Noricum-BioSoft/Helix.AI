"""
Cost Heuristic Table - Relative cost classes and ranges (not exact dollars).

This tool provides cost guidance to the Infrastructure Decision Agent without
fake precision. All costs are:
- Ranges (min, max) not exact values
- Based on documented assumptions
- Marked with confidence levels
- Relative classes (Free, Low, Medium, High) not absolute

Key principle: Cost estimates should acknowledge uncertainty, not pretend precision.
"""

import os
import yaml
import logging
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple, Literal
from pydantic import BaseModel, Field, field_validator

logger = logging.getLogger(__name__)


class CostHeuristic(BaseModel):
    """Cost heuristic for a specific operation on an environment."""
    
    environment: Literal["Local", "EC2", "EMR", "Batch", "Lambda"] = Field(
        ...,
        description="Execution environment"
    )
    
    operation_class: Literal["Small", "Medium", "Large", "Distributed"] = Field(
        ...,
        description="Operation size class"
    )
    
    # Cost range
    cost_range_usd: Tuple[float, float] = Field(
        ...,
        description="Cost range in USD (min, max)"
    )
    
    cost_class: Literal["Free", "Low", "Medium", "High"] = Field(
        ...,
        description="Relative cost class"
    )
    
    cost_confidence: float = Field(
        ...,
        ge=0.0,
        le=1.0,
        description="Confidence in cost estimate (0-1)"
    )
    
    @field_validator('cost_range_usd')
    @classmethod
    def validate_cost_range(cls, v):
        """Validate cost range is non-negative and min <= max."""
        min_cost, max_cost = v
        if min_cost < 0 or max_cost < 0:
            raise ValueError("Cost range cannot contain negative values")
        if min_cost > max_cost:
            raise ValueError(f"Cost range invalid: min ({min_cost}) > max ({max_cost})")
        return v
    
    # Assumptions
    assumptions: str = Field(
        ...,
        description="Explicit assumptions (region, instance type, runtime, etc.)"
    )
    
    # Breakdown
    compute_cost_usd: Optional[Tuple[float, float]] = Field(
        default=None,
        description="Compute cost range"
    )
    
    data_transfer_cost_usd: Optional[float] = Field(
        default=None,
        description="Data transfer cost (if applicable)"
    )
    
    storage_cost_usd: Optional[Tuple[float, float]] = Field(
        default=None,
        description="Storage cost range (if applicable)"
    )
    
    # Metadata
    description: str = Field(
        ...,
        description="Human-readable description"
    )
    
    notes: Optional[str] = Field(
        default=None,
        description="Additional notes or caveats"
    )


class CostHeuristicTable:
    """
    Static table of cost heuristics for different environment + operation combinations.
    
    Provides cost ranges and relative cost classes, not exact dollar amounts.
    All costs acknowledge uncertainty and document assumptions.
    
    Usage:
        cost_table = CostHeuristicTable()
        heuristic = cost_table.get_cost_heuristic("EMR", "Large")
        cost_range = heuristic.cost_range_usd  # (min, max)
    """
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize cost table from YAML config.
        
        Args:
            config_path: Path to YAML config file (uses default if None)
        """
        if config_path is None:
            # Use default config from backend/config/
            project_root = Path(__file__).resolve().parent.parent.parent
            config_path = project_root / "backend" / "config" / "cost_heuristics.yaml"
        
        self.config_path = Path(config_path)
        self.heuristics: Dict[Tuple[str, str], CostHeuristic] = {}
        
        # Load config
        self._load_config()
    
    def _load_config(self):
        """Load cost heuristics from YAML config."""
        if not self.config_path.exists():
            logger.warning(f"Cost heuristic table not found at {self.config_path}, using defaults")
            self._load_defaults()
            return
        
        try:
            with open(self.config_path, "r") as f:
                data = yaml.safe_load(f)
            
            if not data or "heuristics" not in data:
                logger.warning("Invalid cost table format, using defaults")
                self._load_defaults()
                return
            
            for heuristic_data in data["heuristics"]:
                heuristic = CostHeuristic(**heuristic_data)
                key = (heuristic.environment, heuristic.operation_class)
                self.heuristics[key] = heuristic
            
            logger.info(f"Loaded {len(self.heuristics)} cost heuristics from {self.config_path}")
        
        except Exception as e:
            logger.error(f"Failed to load cost heuristic table: {e}")
            self._load_defaults()
    
    def _load_defaults(self):
        """Load default cost heuristics (fallback)."""
        default_region = os.getenv("AWS_REGION", os.getenv("AWS_DEFAULT_REGION", "us-east-1"))
        
        # Local costs (free)
        self.heuristics[("Local", "Small")] = CostHeuristic(
            environment="Local",
            operation_class="Small",
            cost_range_usd=(0.0, 0.0),
            cost_class="Free",
            cost_confidence=1.0,
            assumptions="No cloud costs, local compute",
            description="Small operation on local machine",
        )
        
        self.heuristics[("Local", "Medium")] = CostHeuristic(
            environment="Local",
            operation_class="Medium",
            cost_range_usd=(0.0, 0.0),
            cost_class="Free",
            cost_confidence=1.0,
            assumptions="No cloud costs, local compute",
            description="Medium operation on local machine",
        )
        
        # EC2 costs
        self.heuristics[("EC2", "Small")] = CostHeuristic(
            environment="EC2",
            operation_class="Small",
            cost_range_usd=(0.1, 0.5),
            cost_class="Low",
            cost_confidence=0.6,
            assumptions=f"{default_region}, t3.medium (~$0.04/hr), 5-10min runtime",
            compute_cost_usd=(0.05, 0.2),
            data_transfer_cost_usd=0.01,
            description="Small operation on EC2 instance",
            notes="Data transfer cost depends on file sizes"
        )
        
        self.heuristics[("EC2", "Medium")] = CostHeuristic(
            environment="EC2",
            operation_class="Medium",
            cost_range_usd=(0.5, 2.0),
            cost_class="Low",
            cost_confidence=0.5,
            assumptions=f"{default_region}, m5.xlarge (~$0.19/hr), 15-30min runtime",
            compute_cost_usd=(0.3, 1.5),
            data_transfer_cost_usd=0.05,
            description="Medium operation on EC2 instance",
            notes="Actual cost depends on instance type and runtime"
        )
        
        # EMR costs
        self.heuristics[("EMR", "Medium")] = CostHeuristic(
            environment="EMR",
            operation_class="Medium",
            cost_range_usd=(1.5, 4.0),
            cost_class="Medium",
            cost_confidence=0.5,
            assumptions=f"{default_region}, m5.xlarge nodes (2-3), 15-30min runtime, EMR markup ~40%",
            compute_cost_usd=(1.0, 3.0),
            storage_cost_usd=(0.1, 0.5),
            description="Medium distributed operation on EMR",
            notes="Cluster startup time adds ~3min overhead"
        )
        
        self.heuristics[("EMR", "Large")] = CostHeuristic(
            environment="EMR",
            operation_class="Large",
            cost_range_usd=(3.0, 15.0),
            cost_class="Medium",
            cost_confidence=0.4,
            assumptions=f"{default_region}, m5.xlarge nodes (4-8), 30-60min runtime, EMR markup ~40%",
            compute_cost_usd=(2.0, 12.0),
            storage_cost_usd=(0.5, 2.0),
            description="Large distributed operation on EMR",
            notes="Cost scales with cluster size and runtime"
        )
        
        self.heuristics[("EMR", "Distributed")] = CostHeuristic(
            environment="EMR",
            operation_class="Distributed",
            cost_range_usd=(5.0, 30.0),
            cost_class="High",
            cost_confidence=0.3,
            assumptions=f"{default_region}, m5.xlarge nodes (8-16), 1-3hr runtime, EMR markup ~40%",
            compute_cost_usd=(4.0, 25.0),
            storage_cost_usd=(1.0, 5.0),
            description="Distributed operation with parallelization on EMR",
            notes="High variability depending on parallelization efficiency"
        )
        
        # Batch costs
        self.heuristics[("Batch", "Small")] = CostHeuristic(
            environment="Batch",
            operation_class="Small",
            cost_range_usd=(0.2, 0.8),
            cost_class="Low",
            cost_confidence=0.5,
            assumptions=f"{default_region}, t3.medium compute, 10-20min runtime",
            compute_cost_usd=(0.1, 0.5),
            data_transfer_cost_usd=0.01,
            description="Small containerized job on Batch",
            notes="Container pull adds startup overhead"
        )
        
        self.heuristics[("Batch", "Medium")] = CostHeuristic(
            environment="Batch",
            operation_class="Medium",
            cost_range_usd=(1.0, 3.0),
            cost_class="Low",
            cost_confidence=0.4,
            assumptions=f"{default_region}, m5.xlarge compute, 20-40min runtime",
            compute_cost_usd=(0.8, 2.5),
            data_transfer_cost_usd=0.05,
            description="Medium containerized job on Batch",
            notes="Good for reproducible, medium-sized jobs"
        )
        
        # Lambda costs
        self.heuristics[("Lambda", "Small")] = CostHeuristic(
            environment="Lambda",
            operation_class="Small",
            cost_range_usd=(0.0, 0.2),
            cost_class="Low",
            cost_confidence=0.7,
            assumptions=f"{default_region}, 2GB memory, <5min runtime",
            compute_cost_usd=(0.0, 0.1),
            description="Small serverless function on Lambda",
            notes="First 1M requests/month free"
        )
    
    def get_cost_heuristic(
        self,
        environment: str,
        operation_class: str
    ) -> Optional[CostHeuristic]:
        """
        Get cost heuristic for environment + operation class.
        
        Args:
            environment: Environment name (Local, EC2, EMR, Batch, Lambda)
            operation_class: Operation class (Small, Medium, Large, Distributed)
        
        Returns:
            CostHeuristic or None if not found
        """
        key = (environment, operation_class)
        return self.heuristics.get(key)
    
    def estimate_cost_for_files(
        self,
        environment: str,
        total_size_mb: float,
        file_count: int = 1
    ) -> Optional[CostHeuristic]:
        """
        Estimate cost based on file sizes.
        
        Classification:
        - Small: <100MB total
        - Medium: 100MB-10GB total
        - Large: >10GB total
        - Distributed: >10GB + parallelizable
        
        Args:
            environment: Environment name
            total_size_mb: Total file size in MB
            file_count: Number of files (for parallelization hint)
        
        Returns:
            CostHeuristic or None if not found
        """
        # Classify operation
        if total_size_mb < 100:
            operation_class = "Small"
        elif total_size_mb < 10 * 1024:  # 10GB
            operation_class = "Medium"
        else:
            # Large files - could be distributed if multiple files
            if file_count > 1 and environment in ["EMR", "Batch"]:
                operation_class = "Distributed"
            else:
                operation_class = "Large"
        
        return self.get_cost_heuristic(environment, operation_class)
    
    def compare_costs(
        self,
        environments: List[str],
        operation_class: str
    ) -> List[Tuple[str, CostHeuristic]]:
        """
        Compare costs across multiple environments for the same operation class.
        
        Args:
            environments: List of environment names
            operation_class: Operation class (Small, Medium, Large, Distributed)
        
        Returns:
            List of (environment, heuristic) tuples, sorted by cost (low to high)
        """
        results = []
        for env in environments:
            heuristic = self.get_cost_heuristic(env, operation_class)
            if heuristic:
                results.append((env, heuristic))
        
        # Sort by average cost
        results.sort(key=lambda x: sum(x[1].cost_range_usd) / 2)
        return results
    
    def to_dict(self) -> Dict[str, Any]:
        """Export cost table to dictionary."""
        return {
            f"{env}_{op_class}": heuristic.model_dump()
            for (env, op_class), heuristic in self.heuristics.items()
        }
    
    def to_yaml(self, filepath: str):
        """Export cost table to YAML file."""
        data = {
            "heuristics": [
                heuristic.model_dump()
                for heuristic in self.heuristics.values()
            ]
        }
        
        with open(filepath, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)
        
        logger.info(f"Exported cost heuristic table to {filepath}")
