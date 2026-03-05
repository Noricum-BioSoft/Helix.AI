"""
Agent call tracer - captures agent invocations and contracts.

This module provides observability into the multi-agent system by recording:
- Which agents were called, in what order
- What contracts each agent produced
- Timestamps and durations
"""

import time
from datetime import datetime
from typing import Any, Dict, List, Optional
from pydantic import BaseModel, Field
from enum import Enum

# Import AgentRole with fallback for standalone usage
try:
    from backend.agent import AgentRole
except ImportError:
    # Fallback for when running outside of main project (LEGACY names)
    # See backend.config.agent_registry for canonical names
    class AgentRole(str, Enum):
        """Agent roles in the system."""
        INTENT_DETECTOR = "IntentDetector"
        GURU = "BioinformaticsGuru"  # DEPRECATED: Use WorkflowPlannerAgent
        PLANNER = "BioinformaticsExecutor"  # DEPRECATED: Use ImplementationAgent
        INFRA = "InfrastructureExpert"  # DEPRECATED: Use InfrastructureDecisionAgent
        CODEGEN = "CodeGenerator"
        BROKER = "ExecutionBroker"
        VISUALIZER = "DataVisualizer"


class AgentCall(BaseModel):
    """Record of a single agent invocation."""
    agent: str
    timestamp: str
    duration_ms: Optional[int] = None
    input_data: Dict[str, Any] = Field(default_factory=dict)
    output_contract: Optional[Dict[str, Any]] = None
    contract_type: Optional[str] = None
    error: Optional[str] = None
    metadata: Dict[str, Any] = Field(default_factory=dict)


class AgentCallTracer:
    """
    Traces agent calls during scenario execution.
    
    This is the primary observability mechanism for the demo/eval framework.
    It captures the complete call graph of agents, including:
    - Sequence of agent invocations
    - Contracts exchanged between agents
    - Timing information
    - Errors and warnings
    
    Usage:
        tracer = AgentCallTracer()
        
        with tracer.trace_agent(AgentRole.INTENT_DETECTOR) as call:
            result = await run_intent_detector(prompt)
            call.set_output(result, "IntentResult")
        
        # Later, analyze the trace
        trace = tracer.get_trace()
        assert trace.agent_sequence == ["IntentDetector", "BioinformaticsGuru"]
    """
    
    def __init__(self):
        self.calls: List[AgentCall] = []
        self._current_call: Optional[AgentCall] = None
        self._start_time: Optional[float] = None
    
    def start_agent_call(
        self,
        agent: AgentRole | str,
        input_data: Optional[Dict[str, Any]] = None
    ) -> AgentCall:
        """Start tracing an agent call."""
        agent_name = agent.value if hasattr(agent, "value") else str(agent)
        call = AgentCall(
            agent=agent_name,
            timestamp=datetime.utcnow().isoformat(),
            input_data=input_data or {},
        )
        self._current_call = call
        self._start_time = time.time()
        return call

    def record_contract(self, contract: Dict[str, Any]) -> None:
        """
        Record a contract on the current agent call.

        Back-compat helper used by some unit tests that manually construct
        contract dicts like {"agent": "...", "contract_type": "...", "data": {...}}.
        """
        if self._current_call is None:
            raise RuntimeError("No agent call in progress")

        data = contract.get("data", {})
        contract_type = contract.get("contract_type")

        if not isinstance(data, dict):
            # Keep it permissive; validator can decide if it's invalid.
            self._current_call.output_contract = {"value": data}
        else:
            self._current_call.output_contract = data

        if contract_type:
            self._current_call.contract_type = str(contract_type)
    
    def end_agent_call(
        self,
        output_contract: Optional[Any] = None,
        contract_type: Optional[str] = None,
        error: Optional[str] = None
    ) -> None:
        """End the current agent call and record it."""
        if self._current_call is None:
            raise RuntimeError("No agent call in progress")
        
        if self._start_time is not None:
            duration = int((time.time() - self._start_time) * 1000)
            self._current_call.duration_ms = duration
        
        if output_contract is not None:
            # Convert Pydantic models to dict
            if hasattr(output_contract, 'model_dump'):
                self._current_call.output_contract = output_contract.model_dump()
            elif hasattr(output_contract, 'dict'):
                self._current_call.output_contract = output_contract.dict()
            else:
                self._current_call.output_contract = output_contract
            
            self._current_call.contract_type = contract_type or type(output_contract).__name__
        
        if error is not None:
            self._current_call.error = error
        
        self.calls.append(self._current_call)
        self._current_call = None
        self._start_time = None
    
    def add_metadata(self, key: str, value: Any) -> None:
        """Add metadata to the current agent call."""
        if self._current_call is not None:
            self._current_call.metadata[key] = value
    
    def get_agent_sequence(self) -> List[str]:
        """Get the sequence of agents called."""
        return [call.agent for call in self.calls]
    
    def get_contracts(self) -> List[Dict[str, Any]]:
        """Get all contracts produced."""
        contracts = []
        for call in self.calls:
            if call.output_contract is not None:
                contracts.append({
                    "agent": call.agent,
                    "contract_type": call.contract_type,
                    "data": call.output_contract
                })
        return contracts
    
    def get_total_duration_ms(self) -> int:
        """Get total duration of all agent calls."""
        return sum(call.duration_ms or 0 for call in self.calls)
    
    def has_errors(self) -> bool:
        """Check if any agent calls had errors."""
        return any(call.error is not None for call in self.calls)
    
    def get_errors(self) -> List[str]:
        """Get all errors from agent calls."""
        return [call.error for call in self.calls if call.error is not None]
    
    def clear(self) -> None:
        """Clear the trace history."""
        self.calls = []
        self._current_call = None
        self._start_time = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Export trace as a dictionary."""
        return {
            "calls": [call.model_dump() for call in self.calls],
            "agent_sequence": self.get_agent_sequence(),
            "total_duration_ms": self.get_total_duration_ms(),
            "has_errors": self.has_errors(),
        }


class TracedAgentContext:
    """Context manager for tracing an agent call."""
    
    def __init__(self, tracer: AgentCallTracer, agent: AgentRole, input_data: Optional[Dict[str, Any]] = None):
        self.tracer = tracer
        self.agent = agent
        self.input_data = input_data
        self.call: Optional[AgentCall] = None
    
    def __enter__(self) -> 'TracedAgentContext':
        self.call = self.tracer.start_agent_call(self.agent, self.input_data)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is not None:
            self.tracer.end_agent_call(error=str(exc_val))
        else:
            self.tracer.end_agent_call()
        return False  # Don't suppress exceptions
    
    def set_output(self, contract: Any, contract_type: Optional[str] = None) -> None:
        """Set the output contract for this agent call."""
        if self.call is not None:
            if hasattr(contract, 'model_dump'):
                self.call.output_contract = contract.model_dump()
            elif hasattr(contract, 'dict'):
                self.call.output_contract = contract.dict()
            else:
                self.call.output_contract = contract
            
            self.call.contract_type = contract_type or type(contract).__name__
    
    def add_metadata(self, key: str, value: Any) -> None:
        """Add metadata to this agent call."""
        if self.call is not None:
            self.call.metadata[key] = value


# Helper function to create traced agent context
def trace_agent(tracer: AgentCallTracer, agent: AgentRole, input_data: Optional[Dict[str, Any]] = None) -> TracedAgentContext:
    """Create a traced agent context."""
    return TracedAgentContext(tracer, agent, input_data)
