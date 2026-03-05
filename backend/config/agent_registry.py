"""
Canonical Agent Registry - Single source of truth for agent names and roles.

This registry defines the official agent names, roles, and responsibilities.
All code MUST use these constants to ensure consistency across:
- Orchestrator invocations
- Test assertions
- Trace/contract logging
- Analytics and observability

DO NOT hardcode agent names anywhere else in the codebase.
"""

from enum import Enum
from typing import Dict, List, Optional
from dataclasses import dataclass


class AgentName(str, Enum):
    """
    Canonical agent names.
    
    These are the ONLY valid agent names in the system.
    Use these constants everywhere:
    - orchestrator.py: When invoking agents
    - Test files: When asserting agent sequences
    - Tracing: When logging agent invocations
    - Analytics: When analyzing agent performance
    """
    
    # Infrastructure planning
    INFRASTRUCTURE_DECISION = "InfrastructureDecisionAgent"
    
    # Implementation planning
    IMPLEMENTATION = "ImplementationAgent"
    
    # Workflow planning
    WORKFLOW_PLANNER = "WorkflowPlannerAgent"
    
    # Specialized agents (if/when added)
    TOOL_GENERATOR = "ToolGeneratorAgent"
    REFACTOR = "RefactorAgent"


@dataclass
class AgentInfo:
    """Information about an agent."""
    
    name: AgentName
    display_name: str
    role: str
    responsibilities: List[str]
    input_contract: Optional[str] = None
    output_contract: Optional[str] = None
    
    def __str__(self) -> str:
        return f"{self.name.value} ({self.role})"


# Agent registry with detailed information
AGENT_REGISTRY: Dict[AgentName, AgentInfo] = {
    AgentName.INFRASTRUCTURE_DECISION: AgentInfo(
        name=AgentName.INFRASTRUCTURE_DECISION,
        display_name="Infrastructure Decision Agent",
        role="Infrastructure Selector",
        responsibilities=[
            "Analyze data inputs (size, location, format)",
            "Estimate computational requirements (CPU, memory, runtime)",
            "Recommend execution environment (Local/EC2/EMR/Batch/Lambda)",
            "Provide cost estimates with ranges and assumptions",
            "List alternative infrastructure options with tradeoffs"
        ],
        input_contract="WorkflowPlan",
        output_contract="InfraDecision"
    ),
    
    AgentName.IMPLEMENTATION: AgentInfo(
        name=AgentName.IMPLEMENTATION,
        display_name="Implementation Agent",
        role="Execution Planner",
        responsibilities=[
            "Design execution plan (commands, containers, resources)",
            "Select appropriate tools and container images",
            "Define retry policies and timeouts",
            "Specify expected outputs and validation criteria",
            "Generate ExecutionToolSpec for external runner"
        ],
        input_contract="WorkflowPlan + InfraDecision",
        output_contract="ExecutionToolSpec"
    ),
    
    AgentName.WORKFLOW_PLANNER: AgentInfo(
        name=AgentName.WORKFLOW_PLANNER,
        display_name="Workflow Planner Agent",
        role="Workflow Designer",
        responsibilities=[
            "Parse user intent into workflow steps",
            "Identify required tools and data inputs",
            "Define workflow dependencies and ordering",
            "Generate WorkflowPlan with structured steps"
        ],
        input_contract="str (user command)",
        output_contract="WorkflowPlan"
    ),
    
    AgentName.TOOL_GENERATOR: AgentInfo(
        name=AgentName.TOOL_GENERATOR,
        display_name="Tool Generator Agent",
        role="Tool Creator",
        responsibilities=[
            "Generate new tool implementations from specifications",
            "Create tool wrappers for bioinformatics software",
            "Validate tool contracts and interfaces"
        ],
        input_contract="ToolSpec",
        output_contract="Tool implementation"
    ),
    
    AgentName.REFACTOR: AgentInfo(
        name=AgentName.REFACTOR,
        display_name="Refactor Agent",
        role="Code Optimizer",
        responsibilities=[
            "Refactor existing code for maintainability",
            "Optimize performance bottlenecks",
            "Apply best practices and patterns"
        ],
        input_contract="Code + RefactorSpec",
        output_contract="Refactored code"
    ),
}


def get_agent_info(agent_name: AgentName) -> AgentInfo:
    """Get information about an agent."""
    return AGENT_REGISTRY[agent_name]


def get_all_agents() -> List[AgentInfo]:
    """Get all registered agents."""
    return list(AGENT_REGISTRY.values())


def validate_agent_name(name: str) -> bool:
    """Check if a name is a valid agent name."""
    try:
        AgentName(name)
        return True
    except ValueError:
        return False


def get_agent_sequence_display(sequence: List[AgentName]) -> str:
    """
    Format agent sequence for display.
    
    Example: "InfrastructureDecisionAgent → ImplementationAgent"
    """
    return " → ".join([agent.value for agent in sequence])


# Common agent sequences (for validation and testing)
COMMON_SEQUENCES = {
    "infrastructure_to_implementation": [
        AgentName.INFRASTRUCTURE_DECISION,
        AgentName.IMPLEMENTATION
    ],
    "full_workflow": [
        AgentName.WORKFLOW_PLANNER,
        AgentName.INFRASTRUCTURE_DECISION,
        AgentName.IMPLEMENTATION
    ],
}


def get_common_sequence(name: str) -> Optional[List[AgentName]]:
    """Get a predefined agent sequence by name."""
    return COMMON_SEQUENCES.get(name)
