"""
Configuration module for Helix.AI backend.

Contains:
- agent_registry: Canonical agent names and roles
"""

from backend.config.agent_registry import (
    AgentName,
    AgentInfo,
    AGENT_REGISTRY,
    get_agent_info,
    get_all_agents,
    validate_agent_name,
    get_agent_sequence_display,
    COMMON_SEQUENCES,
    get_common_sequence,
)

__all__ = [
    "AgentName",
    "AgentInfo",
    "AGENT_REGISTRY",
    "get_agent_info",
    "get_all_agents",
    "validate_agent_name",
    "get_agent_sequence_display",
    "COMMON_SEQUENCES",
    "get_common_sequence",
]
