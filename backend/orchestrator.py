"""
EXPERIMENTAL ORCHESTRATOR - PROTOTYPE ONLY 🟡

Orchestrator - Coordinates multi-agent execution with request tracking.

**STATUS:** This is an EXPERIMENTAL clean 2-agent pipeline with Pydantic contracts
and provenance hashing. Currently used only for CLI demos and unit tests.

- Status: 🟡 PROTOTYPE (not used in production)
- Entry Point: Orchestrator.run_pipeline()
- Used By: cli_demo.py, test_phase3_agents.py only
- Architecture: Clean 2-agent pipeline (InfrastructureDecisionAgent → ImplementationAgent)

For the PRIMARY orchestrator used in production, see backend/agent.py
For detailed comparison, see docs/ORCHESTRATION_DUALITY.md

---

The Orchestrator is responsible for:
1. Invoking agents in sequence (InfraDecisionAgent → ImplementationAgent → ...)
2. Tracking request lifecycle and agent handoffs
3. Logging inputs/outputs for each agent
4. Managing errors and fallbacks
5. Providing trace/audit logs

Key principle: Composable, observable multi-agent pipelines.
"""

import logging
import uuid
import time
import json
import hashlib
from typing import Optional, Dict, Any, List
from dataclasses import dataclass, field
from datetime import datetime

from backend.contracts.workflow_plan import WorkflowPlan
from backend.contracts.infra_decision import InfraDecision
from backend.contracts.execution_spec import ExecutionToolSpec
from backend.infrastructure_decision_agent import decide_infrastructure
from backend.implementation_agent import plan_implementation
from backend.config.agent_registry import AgentName

logger = logging.getLogger(__name__)


def compute_contract_hash(contract_json: str) -> str:
    """
    Compute stable SHA256 hash of a contract JSON string.
    
    This enables:
    - Cheap deduplication (detect if same input already processed)
    - Provenance tracking (detect if contract changed between runs)
    - Stronger baseline diffs (quick comparison without deep inspection)
    
    Args:
        contract_json: JSON string of the contract (from model_dump_json())
    
    Returns:
        Hex string of SHA256 hash (64 characters)
    
    Note:
        We hash the raw JSON string (not parsed dict) to ensure stability.
        Pydantic model_dump_json() produces deterministic output (sorted keys).
    """
    if not contract_json:
        return ""
    
    # Compute SHA256 hash of UTF-8 encoded JSON
    hash_obj = hashlib.sha256(contract_json.encode('utf-8'))
    return hash_obj.hexdigest()


@dataclass
class AgentInvocation:
    """Record of a single agent invocation."""
    
    agent_name: str
    request_id: str
    start_time: float
    end_time: Optional[float] = None
    duration_ms: Optional[float] = None
    success: bool = False
    error: Optional[str] = None
    
    # Input/output contracts
    input_contract: Optional[str] = None  # JSON serialized
    output_contract: Optional[str] = None  # JSON serialized
    
    # Metadata
    input_hash: Optional[str] = None
    output_hash: Optional[str] = None
    confidence_score: Optional[float] = None
    warnings: List[str] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "agent_name": self.agent_name,
            "request_id": self.request_id,
            "start_time": self.start_time,
            "end_time": self.end_time,
            "duration_ms": self.duration_ms,
            "success": self.success,
            "error": self.error,
            "input_hash": self.input_hash,
            "output_hash": self.output_hash,
            "confidence_score": self.confidence_score,
            "warnings": self.warnings,
            "timestamp": datetime.fromtimestamp(self.start_time).isoformat()
        }


@dataclass
class OrchestratorTrace:
    """Complete trace of multi-agent execution."""
    
    request_id: str
    user_command: str
    start_time: float
    end_time: Optional[float] = None
    duration_ms: Optional[float] = None
    success: bool = False
    error: Optional[str] = None
    
    # Agent invocations in sequence
    invocations: List[AgentInvocation] = field(default_factory=list)
    
    # Final outputs
    infra_decision: Optional[InfraDecision] = None
    execution_spec: Optional[ExecutionToolSpec] = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "request_id": self.request_id,
            "user_command": self.user_command,
            "start_time": self.start_time,
            "end_time": self.end_time,
            "duration_ms": self.duration_ms,
            "success": self.success,
            "error": self.error,
            "invocations": [inv.to_dict() for inv in self.invocations],
            "infra_decision": self.infra_decision.model_dump() if self.infra_decision else None,
            "execution_spec": self.execution_spec.model_dump() if self.execution_spec else None,
            "timestamp": datetime.fromtimestamp(self.start_time).isoformat()
        }
    
    def to_json(self, indent: int = 2) -> str:
        """Serialize trace to JSON."""
        return json.dumps(self.to_dict(), indent=indent)
    
    def get_invocation_hashes(self) -> Dict[str, tuple[str, str]]:
        """
        Get input/output hashes for all invocations.
        
        Returns:
            Dict mapping agent_name -> (input_hash, output_hash)
        
        Use cases:
        - Detect if same input was already processed (deduplication)
        - Check if contract changed between runs (provenance)
        - Quick baseline comparison without deep inspection
        """
        return {
            inv.agent_name: (inv.input_hash or "", inv.output_hash or "")
            for inv in self.invocations
        }
    
    def has_same_contracts_as(self, other: "OrchestratorTrace") -> bool:
        """
        Check if this trace has same contracts as another trace.
        
        Uses hashes for fast comparison (O(n) vs O(n*m) for deep inspection).
        
        Args:
            other: Another OrchestratorTrace to compare against
        
        Returns:
            True if all agents have matching input/output hashes
        """
        self_hashes = self.get_invocation_hashes()
        other_hashes = other.get_invocation_hashes()
        
        # Check if agent sequences match
        if set(self_hashes.keys()) != set(other_hashes.keys()):
            return False
        
        # Check if hashes match for each agent
        for agent_name in self_hashes:
            if self_hashes[agent_name] != other_hashes[agent_name]:
                return False
        
        return True
    
    def get_contract_diff_summary(self, other: "OrchestratorTrace") -> Dict[str, str]:
        """
        Get summary of contract differences between two traces.
        
        Returns:
            Dict mapping agent_name -> diff status ("added", "removed", "input_changed", 
            "output_changed", "both_changed", "unchanged")
        
        Use cases:
        - Baseline regression testing (did behavior change?)
        - Provenance tracking (what changed between versions?)
        - Debug multi-agent systems (which agent introduced variance?)
        """
        self_hashes = self.get_invocation_hashes()
        other_hashes = other.get_invocation_hashes()
        
        diff_summary = {}
        
        # Check all agents in self
        for agent_name, (self_input, self_output) in self_hashes.items():
            if agent_name not in other_hashes:
                diff_summary[agent_name] = "added"
            else:
                other_input, other_output = other_hashes[agent_name]
                if self_input != other_input and self_output != other_output:
                    diff_summary[agent_name] = "both_changed"
                elif self_input != other_input:
                    diff_summary[agent_name] = "input_changed"
                elif self_output != other_output:
                    diff_summary[agent_name] = "output_changed"
                else:
                    diff_summary[agent_name] = "unchanged"
        
        # Check for removed agents
        for agent_name in other_hashes:
            if agent_name not in self_hashes:
                diff_summary[agent_name] = "removed"
        
        return diff_summary


class Orchestrator:
    """
    Orchestrator for multi-agent workflows.
    
    Coordinates agent execution, tracks request lifecycle, provides observability.
    
    Usage:
        orchestrator = Orchestrator()
        
        # Option 1: Run full pipeline
        trace = await orchestrator.run_pipeline(
            command="Run FastQC on sample.fq",
            workflow_plan=plan
        )
        
        # Option 2: Run agents individually
        request_id = orchestrator.generate_request_id()
        infra_decision = await orchestrator.invoke_infrastructure_agent(
            workflow_plan=plan,
            request_id=request_id
        )
        execution_spec = await orchestrator.invoke_implementation_agent(
            workflow_plan=plan,
            infra_decision=infra_decision,
            request_id=request_id
        )
    """
    
    def __init__(self):
        """Initialize Orchestrator."""
        self.traces: Dict[str, OrchestratorTrace] = {}
    
    @staticmethod
    def generate_request_id() -> str:
        """Generate unique request ID."""
        return f"req_{uuid.uuid4().hex[:12]}"
    
    async def invoke_infrastructure_agent(
        self,
        workflow_plan: WorkflowPlan,
        request_id: str,
        trace: Optional[OrchestratorTrace] = None
    ) -> InfraDecision:
        """
        Invoke Infrastructure Decision Agent.
        
        Args:
            workflow_plan: Workflow to analyze
            request_id: Request ID for tracing
            trace: Optional trace to append to
        
        Returns:
            InfraDecision with infrastructure recommendation
        
        Raises:
            Exception: If agent fails
        """
        invocation = AgentInvocation(
            agent_name=AgentName.INFRASTRUCTURE_DECISION.value,
            request_id=request_id,
            start_time=time.time()
        )
        
        try:
            logger.info(f"[{request_id}] Invoking InfrastructureDecisionAgent")
            
            # Serialize input and compute hash
            invocation.input_contract = workflow_plan.model_dump_json()
            invocation.input_hash = compute_contract_hash(invocation.input_contract)
            
            # Invoke agent with full workflow_plan (contains rich dataset metadata)
            infra_decision = await decide_infrastructure(
                command=workflow_plan.description,
                workflow_plan=workflow_plan,  # Pass full plan with size_bytes, location_type, etc.
                request_id=request_id
            )
            
            # Record success
            invocation.end_time = time.time()
            invocation.duration_ms = (invocation.end_time - invocation.start_time) * 1000
            invocation.success = True
            
            # Serialize output and compute hash
            invocation.output_contract = infra_decision.model_dump_json()
            invocation.output_hash = compute_contract_hash(invocation.output_contract)
            
            invocation.confidence_score = infra_decision.confidence_score
            invocation.warnings = infra_decision.warnings
            
            logger.info(
                f"[{request_id}] InfrastructureDecisionAgent completed: "
                f"{infra_decision.infrastructure} (confidence={infra_decision.confidence_score:.2f}) "
                f"[input_hash={invocation.input_hash[:8]}..., output_hash={invocation.output_hash[:8]}...]"
            )
            
            if trace:
                trace.invocations.append(invocation)
                trace.infra_decision = infra_decision
            
            return infra_decision
        
        except Exception as e:
            # Record failure
            invocation.end_time = time.time()
            invocation.duration_ms = (invocation.end_time - invocation.start_time) * 1000
            invocation.success = False
            invocation.error = f"{type(e).__name__}: {str(e)}"
            
            logger.error(f"[{request_id}] InfrastructureDecisionAgent failed: {invocation.error}")
            
            if trace:
                trace.invocations.append(invocation)
            
            raise
    
    async def invoke_implementation_agent(
        self,
        workflow_plan: WorkflowPlan,
        infra_decision: InfraDecision,
        request_id: str,
        trace: Optional[OrchestratorTrace] = None
    ) -> ExecutionToolSpec:
        """
        Invoke Implementation Agent.
        
        Args:
            workflow_plan: Workflow to execute
            infra_decision: Infrastructure decision from previous agent
            request_id: Request ID for tracing
            trace: Optional trace to append to
        
        Returns:
            ExecutionToolSpec with execution plan
        
        Raises:
            Exception: If agent fails
        """
        invocation = AgentInvocation(
            agent_name=AgentName.IMPLEMENTATION.value,
            request_id=request_id,
            start_time=time.time()
        )
        
        try:
            logger.info(f"[{request_id}] Invoking ImplementationAgent")
            
            # Serialize inputs and compute hash
            invocation.input_contract = json.dumps({
                "workflow_plan": workflow_plan.model_dump(),
                "infra_decision": infra_decision.model_dump()
            }, sort_keys=True)  # sort_keys for deterministic hash
            invocation.input_hash = compute_contract_hash(invocation.input_contract)
            
            # Invoke agent
            execution_spec = await plan_implementation(
                workflow_plan=workflow_plan,
                infra_decision=infra_decision,
                request_id=request_id
            )
            
            # Record success
            invocation.end_time = time.time()
            invocation.duration_ms = (invocation.end_time - invocation.start_time) * 1000
            invocation.success = True
            
            # Serialize output and compute hash
            invocation.output_contract = execution_spec.model_dump_json()
            invocation.output_hash = compute_contract_hash(invocation.output_contract)
            
            invocation.confidence_score = execution_spec.confidence_score
            invocation.warnings = execution_spec.warnings
            
            logger.info(
                f"[{request_id}] ImplementationAgent completed: "
                f"{execution_spec.tool_name} on {execution_spec.infrastructure} "
                f"(confidence={execution_spec.confidence_score:.2f}) "
                f"[input_hash={invocation.input_hash[:8]}..., output_hash={invocation.output_hash[:8]}...]"
            )
            
            if trace:
                trace.invocations.append(invocation)
                trace.execution_spec = execution_spec
            
            return execution_spec
        
        except Exception as e:
            # Record failure
            invocation.end_time = time.time()
            invocation.duration_ms = (invocation.end_time - invocation.start_time) * 1000
            invocation.success = False
            invocation.error = f"{type(e).__name__}: {str(e)}"
            
            logger.error(f"[{request_id}] ImplementationAgent failed: {invocation.error}")
            
            if trace:
                trace.invocations.append(invocation)
            
            raise
    
    async def run_pipeline(
        self,
        command: str,
        workflow_plan: WorkflowPlan,
        request_id: Optional[str] = None
    ) -> OrchestratorTrace:
        """
        Run complete multi-agent pipeline.
        
        Pipeline:
        1. InfrastructureDecisionAgent (WHERE to execute)
        2. ImplementationAgent (HOW to execute)
        
        Args:
            command: User command
            workflow_plan: Parsed workflow plan
            request_id: Optional request ID (generated if None)
        
        Returns:
            OrchestratorTrace with complete execution history
        """
        if request_id is None:
            request_id = self.generate_request_id()
        
        trace = OrchestratorTrace(
            request_id=request_id,
            user_command=command,
            start_time=time.time()
        )
        
        try:
            logger.info(f"[{request_id}] Starting multi-agent pipeline for: {command}")
            
            # Step 1: Infrastructure Decision
            infra_decision = await self.invoke_infrastructure_agent(
                workflow_plan=workflow_plan,
                request_id=request_id,
                trace=trace
            )
            
            # Step 2: Implementation Planning
            execution_spec = await self.invoke_implementation_agent(
                workflow_plan=workflow_plan,
                infra_decision=infra_decision,
                request_id=request_id,
                trace=trace
            )
            
            # Pipeline succeeded
            trace.end_time = time.time()
            trace.duration_ms = (trace.end_time - trace.start_time) * 1000
            trace.success = True
            
            logger.info(
                f"[{request_id}] Pipeline completed successfully in {trace.duration_ms:.0f}ms"
            )
            
            # Store trace
            self.traces[request_id] = trace
            
            return trace
        
        except Exception as e:
            # Pipeline failed
            trace.end_time = time.time()
            trace.duration_ms = (trace.end_time - trace.start_time) * 1000
            trace.success = False
            trace.error = f"{type(e).__name__}: {str(e)}"
            
            logger.error(f"[{request_id}] Pipeline failed: {trace.error}")
            
            # Store trace (even for failures)
            self.traces[request_id] = trace
            
            raise
    
    def get_trace(self, request_id: str) -> Optional[OrchestratorTrace]:
        """Get trace for a request."""
        return self.traces.get(request_id)
    
    def get_all_traces(self) -> List[OrchestratorTrace]:
        """Get all traces."""
        return list(self.traces.values())
    
    def export_trace_json(self, request_id: str, filepath: str):
        """Export trace to JSON file."""
        trace = self.get_trace(request_id)
        if trace:
            with open(filepath, 'w') as f:
                f.write(trace.to_json())
            logger.info(f"Exported trace {request_id} to {filepath}")
        else:
            logger.warning(f"No trace found for request_id={request_id}")


# Global orchestrator instance
_orchestrator = None

def get_orchestrator() -> Orchestrator:
    """Get global orchestrator instance (singleton)."""
    global _orchestrator
    if _orchestrator is None:
        _orchestrator = Orchestrator()
    return _orchestrator
