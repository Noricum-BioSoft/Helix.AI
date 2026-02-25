"""
Contracts module for Helix.AI multi-agent system.

This module defines strict Pydantic contracts for agent inputs and outputs,
ensuring:
- Type safety
- Automatic validation
- Clear interface boundaries
- Easy serialization/deserialization

All inter-agent communication should use these contracts.
"""

from backend.contracts.infra_decision import (
    InfraDecision,
    FileAnalysis,
    ComputationalRequirements,
    CostAnalysis,
    InfraAlternative,
)

from backend.contracts.workflow_plan import (
    WorkflowPlan,
    DataInput,
    OperationSpec,
    ConstraintSpec,
)

from backend.contracts.dataset_spec import (
    DatasetSpec,
    FileSpec,
)

from backend.contracts.execution_spec import (
    ExecutionToolSpec,
    ContainerSpec,
    CommandSpec,
    RetryPolicy,
    ResourceRequirements,
    OutputSpec,
)

__all__ = [
    # Infrastructure Decision
    "InfraDecision",
    "FileAnalysis",
    "ComputationalRequirements",
    "CostAnalysis",
    "InfraAlternative",
    # Workflow Plan
    "WorkflowPlan",
    "DataInput",
    "OperationSpec",
    "ConstraintSpec",
    # Dataset Spec
    "DatasetSpec",
    "FileSpec",
    # Execution Spec
    "ExecutionToolSpec",
    "ContainerSpec",
    "CommandSpec",
    "RetryPolicy",
    "ResourceRequirements",
    "OutputSpec",
]
