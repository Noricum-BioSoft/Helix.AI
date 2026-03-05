"""
Demo & Eval Framework for Helix.AI Multi-Agent System

This package provides scenario-driven testing and demonstration capabilities.
"""

from .scenario import (
    Scenario,
    ScenarioMetadata,
    ScenarioInput,
    ExpectedBehavior,
    ScenarioLoader,
    ScenarioCategory,
)
from .executor import ScenarioExecutor, ExecutionTrace
from .tracer import AgentCallTracer, AgentCall
from .validator import ContractValidator, ValidationResult
from .reporter import ScenarioReporter

__all__ = [
    "Scenario",
    "ScenarioMetadata",
    "ScenarioInput",
    "ExpectedBehavior",
    "ScenarioLoader",
    "ScenarioCategory",
    "ScenarioExecutor",
    "ExecutionTrace",
    "AgentCallTracer",
    "AgentCall",
    "ContractValidator",
    "ValidationResult",
    "ScenarioReporter",
]
