"""
Contract and policy validators.

Validates that:
1. Contracts match expected types and fields
2. Agent sequence follows handoff policy
3. Policy checks are satisfied (no unauthorized mutations, etc.)
"""

import re
from typing import Any, Dict, List, Optional
from pydantic import BaseModel, Field
from enum import Enum

# Import with fallback for standalone usage
try:
    from backend.agent import HandoffPolicy, AgentRole, PolicyViolationError
except ImportError:
    # Fallback implementations (LEGACY - see backend.config.agent_registry)
    class AgentRole(str, Enum):
        """Agent roles in the system."""
        INTENT_DETECTOR = "IntentDetector"
        GURU = "BioinformaticsGuru"
        PLANNER = "BioinformaticsExecutor"
        INFRA = "InfrastructureExpert"
        CODEGEN = "CodeGenerator"
        BROKER = "ExecutionBroker"
        VISUALIZER = "DataVisualizer"
    
    class PolicyViolationError(Exception):
        """Raised when a policy is violated."""
        pass
    
    class HandoffPolicy:
        """Minimal handoff policy for validation."""
        ALLOWED_HANDOFFS = {
            AgentRole.INTENT_DETECTOR: [AgentRole.GURU, AgentRole.PLANNER],
            AgentRole.GURU: [AgentRole.PLANNER],
            AgentRole.PLANNER: [AgentRole.INFRA],
            AgentRole.INFRA: [AgentRole.CODEGEN, AgentRole.BROKER],
            AgentRole.CODEGEN: [AgentRole.BROKER],
            AgentRole.BROKER: [AgentRole.VISUALIZER],
            AgentRole.VISUALIZER: [],
        }
        
        INTENT_ROUTING = {
            "ask": AgentRole.GURU,
            "qa": AgentRole.GURU,
            "execute": AgentRole.PLANNER,
        }
        
        def validate_handoff(self, from_agent, to_agent, user_consent=False):
            allowed = self.ALLOWED_HANDOFFS.get(from_agent, [])
            if to_agent not in allowed:
                raise PolicyViolationError(
                    f"Illegal handoff: {from_agent.value} → {to_agent.value}"
                )
        
        def validate_workflow_sequence(self, agent_sequence, intent):
            if not agent_sequence:
                raise PolicyViolationError("Empty agent sequence")
            if agent_sequence[0] != AgentRole.INTENT_DETECTOR:
                raise PolicyViolationError("First agent must be IntentDetector")
            for i in range(len(agent_sequence) - 1):
                self.validate_handoff(agent_sequence[i], agent_sequence[i + 1])

# Agent name mapping for validation (OLD → NEW canonical names)
# See docs/AGENT_NAMING_MIGRATION.md
CANONICAL_AGENT_NAMES = {
    "BioinformaticsGuru": "WorkflowPlannerAgent",
    "BioinformaticsExecutor": "ImplementationAgent",
    "InfrastructureExpert": "InfrastructureDecisionAgent",
    "ExecutionBroker": "ExecutionBroker",  # Not an agent, but kept for compatibility
    "IntentDetector": "IntentDetector",
    "DataVisualizer": "DataVisualizer",
    "CodeGenerator": "ToolGeneratorAgent",
}

from .scenario import (
    ExpectedBehavior,
    ContractExpectation,
    FieldValidation,
    ValidationOperator,
    PolicyCheck,
)
from .tracer import AgentCallTracer


class ValidationIssue(BaseModel):
    """A validation issue found during checking."""
    severity: str  # "error", "warning"
    category: str  # "contract", "policy", "sequence", "timing"
    message: str
    details: Optional[Dict[str, Any]] = None


class ValidationResult(BaseModel):
    """Result of validating a scenario execution."""
    passed: bool
    issues: List[ValidationIssue] = Field(default_factory=list)
    summary: str
    
    def add_error(self, category: str, message: str, details: Optional[Dict[str, Any]] = None) -> None:
        """Add an error issue."""
        self.issues.append(ValidationIssue(
            severity="error",
            category=category,
            message=message,
            details=details
        ))
        self.passed = False
    
    def add_warning(self, category: str, message: str, details: Optional[Dict[str, Any]] = None) -> None:
        """Add a warning issue."""
        self.issues.append(ValidationIssue(
            severity="warning",
            category=category,
            message=message,
            details=details
        ))
    
    def get_errors(self) -> List[ValidationIssue]:
        """Get only error-level issues."""
        return [i for i in self.issues if i.severity == "error"]
    
    def get_warnings(self) -> List[ValidationIssue]:
        """Get only warning-level issues."""
        return [i for i in self.issues if i.severity == "warning"]


class ContractValidator:
    """
    Validates contracts against expected behavior.
    
    This is the core validation engine that checks:
    - Agent sequence matches expectations
    - Contracts have correct types and fields
    - Policy constraints are satisfied
    """
    
    def __init__(self):
        self.handoff_policy = HandoffPolicy()
    
    def validate(
        self,
        tracer: AgentCallTracer,
        expected: ExpectedBehavior,
        intent: Optional[str] = None
    ) -> ValidationResult:
        """
        Validate a traced execution against expected behavior.
        
        Args:
            tracer: Tracer with recorded agent calls
            expected: Expected behavior from scenario
            intent: Intent classification (for policy validation)
        
        Returns:
            ValidationResult with pass/fail and detailed issues
        """
        result = ValidationResult(passed=True, summary="")
        
        # 1. Validate agent sequence
        self._validate_agent_sequence(tracer, expected, intent, result)
        
        # 2. Validate contracts
        self._validate_contracts(tracer, expected, result)
        
        # 3. Validate policy checks
        self._validate_policy_checks(tracer, expected, result)
        
        # 4. Validate semantic invariants (cross-contract checks)
        self._validate_semantic_invariants(tracer, result)
        
        # 5. Validate timing (if specified)
        self._validate_timing(tracer, expected, result)
        
        # Generate summary
        error_count = len(result.get_errors())
        warning_count = len(result.get_warnings())
        
        if error_count == 0 and warning_count == 0:
            result.summary = "All validations passed ✓"
        elif error_count == 0:
            result.summary = f"Passed with {warning_count} warning(s) ⚠"
        else:
            result.summary = f"Failed with {error_count} error(s) and {warning_count} warning(s) ✗"
        
        return result
    
    def _validate_agent_sequence(
        self,
        tracer: AgentCallTracer,
        expected: ExpectedBehavior,
        intent: Optional[str],
        result: ValidationResult
    ) -> None:
        """Validate the agent call sequence."""
        actual_sequence = tracer.get_agent_sequence()
        expected_sequence = expected.agent_sequence
        
        # Check sequence length
        if len(actual_sequence) != len(expected_sequence):
            result.add_error(
                "sequence",
                f"Agent sequence length mismatch: expected {len(expected_sequence)}, got {len(actual_sequence)}",
                {
                    "expected": expected_sequence,
                    "actual": actual_sequence
                }
            )
            return  # Don't check further if lengths don't match
        
        # Check each agent in sequence
        for i, (expected_agent, actual_agent) in enumerate(zip(expected_sequence, actual_sequence)):
            if expected_agent != actual_agent:
                result.add_error(
                    "sequence",
                    f"Agent mismatch at position {i}: expected {expected_agent}, got {actual_agent}",
                    {
                        "position": i,
                        "expected": expected_agent,
                        "actual": actual_agent,
                        "full_sequence": actual_sequence
                    }
                )
        
        # Validate handoff policy compliance (if we have intent)
        # Skip if policy_checks is explicitly empty (tool mapping mode)
        skip_policy = expected.policy_checks is not None and len(expected.policy_checks) == 0
        
        if intent and len(actual_sequence) > 0 and not skip_policy:
            try:
                # Convert string names to AgentRole enum
                agent_roles = []
                for agent_name in actual_sequence:
                    try:
                        role = AgentRole(agent_name)
                        agent_roles.append(role)
                    except ValueError:
                        result.add_error(
                            "policy",
                            f"Unknown agent role: {agent_name}",
                            {"agent_name": agent_name}
                        )
                        return
                
                # Validate the sequence
                self.handoff_policy.validate_workflow_sequence(agent_roles, intent)
            except PolicyViolationError as e:
                result.add_error(
                    "policy",
                    f"Handoff policy violation: {str(e)}",
                    {
                        "sequence": actual_sequence,
                        "intent": intent
                    }
                )
    
    def _validate_contracts(
        self,
        tracer: AgentCallTracer,
        expected: ExpectedBehavior,
        result: ValidationResult
    ) -> None:
        """Validate output contracts from agents."""
        actual_contracts = tracer.get_contracts()
        
        for expected_contract in expected.contracts:
            # Find the contract from this agent
            agent_contracts = [
                c for c in actual_contracts
                if c["agent"] == expected_contract.agent
            ]
            
            if not agent_contracts:
                result.add_error(
                    "contract",
                    f"No contract found from agent: {expected_contract.agent}",
                    {"expected_agent": expected_contract.agent}
                )
                continue
            
            # Take the most recent contract from this agent
            contract = agent_contracts[-1]
            
            # Check contract type
            if contract["contract_type"] != expected_contract.output_type:
                result.add_error(
                    "contract",
                    f"Contract type mismatch for {expected_contract.agent}: expected {expected_contract.output_type}, got {contract['contract_type']}",
                    {
                        "agent": expected_contract.agent,
                        "expected_type": expected_contract.output_type,
                        "actual_type": contract["contract_type"]
                    }
                )
                continue
            
            # Validate individual fields
            contract_data = contract["data"]
            for validation in expected_contract.validations:
                self._validate_field(
                    contract_data,
                    validation,
                    expected_contract.agent,
                    result
                )
    
    def _validate_field(
        self,
        contract_data: Dict[str, Any],
        validation: FieldValidation,
        agent: str,
        result: ValidationResult
    ) -> None:
        """Validate a single field in a contract."""
        field = validation.field
        
        # Handle nested fields (e.g., "metadata.scale")
        if "." in field:
            parts = field.split(".")
            value = contract_data
            for part in parts:
                if isinstance(value, dict) and part in value:
                    value = value[part]
                else:
                    result.add_error(
                        "contract",
                        f"Nested field '{field}' not found in contract from {agent}",
                        {
                            "agent": agent,
                            "field": field,
                            "available_fields": list(contract_data.keys())
                        }
                    )
                    return
        else:
            # Check if field exists
            if field not in contract_data:
                result.add_error(
                    "contract",
                    f"Field '{field}' not found in contract from {agent}",
                    {
                        "agent": agent,
                        "field": field,
                        "available_fields": list(contract_data.keys())
                    }
                )
                return
            
            value = contract_data[field]
        operator, expected_value = validation.get_operator_and_value()
        
        # Perform validation based on operator
        try:
            if operator == ValidationOperator.EQUALS:
                if value != expected_value:
                    result.add_error(
                        "contract",
                        f"Field '{field}' in {agent}: expected '{expected_value}', got '{value}'",
                        {"agent": agent, "field": field, "expected": expected_value, "actual": value}
                    )
            
            elif operator == ValidationOperator.CONTAINS:
                if isinstance(value, str):
                    missing = [item for item in expected_value if item.lower() not in value.lower()]
                    if missing:
                        result.add_error(
                            "contract",
                            f"Field '{field}' in {agent}: missing required substrings: {missing}",
                            {"agent": agent, "field": field, "missing": missing, "value": value}
                        )
                elif isinstance(value, list):
                    missing = [item for item in expected_value if item not in value]
                    if missing:
                        result.add_error(
                            "contract",
                            f"Field '{field}' in {agent}: missing required items: {missing}",
                            {"agent": agent, "field": field, "missing": missing, "value": value}
                        )
            
            elif operator == ValidationOperator.GREATER_THAN:
                if not (isinstance(value, (int, float)) and value > expected_value):
                    result.add_error(
                        "contract",
                        f"Field '{field}' in {agent}: expected > {expected_value}, got {value}",
                        {"agent": agent, "field": field, "threshold": expected_value, "actual": value}
                    )
            
            elif operator == ValidationOperator.LESS_THAN:
                if not (isinstance(value, (int, float)) and value < expected_value):
                    result.add_error(
                        "contract",
                        f"Field '{field}' in {agent}: expected < {expected_value}, got {value}",
                        {"agent": agent, "field": field, "threshold": expected_value, "actual": value}
                    )
            
            elif operator == ValidationOperator.MIN_LENGTH:
                if not (hasattr(value, '__len__') and len(value) >= expected_value):
                    result.add_error(
                        "contract",
                        f"Field '{field}' in {agent}: length {len(value) if hasattr(value, '__len__') else 'N/A'} is less than minimum {expected_value}",
                        {"agent": agent, "field": field, "min_length": expected_value, "actual_length": len(value) if hasattr(value, '__len__') else None}
                    )
            
            elif operator == ValidationOperator.MAX_LENGTH:
                if not (hasattr(value, '__len__') and len(value) <= expected_value):
                    result.add_error(
                        "contract",
                        f"Field '{field}' in {agent}: length {len(value)} exceeds maximum {expected_value}",
                        {"agent": agent, "field": field, "max_length": expected_value, "actual_length": len(value)}
                    )
            
            elif operator == ValidationOperator.REGEX_MATCH:
                if not isinstance(value, str) or not re.match(expected_value, value):
                    result.add_error(
                        "contract",
                        f"Field '{field}' in {agent}: does not match regex pattern '{expected_value}'",
                        {"agent": agent, "field": field, "pattern": expected_value, "value": value}
                    )
            
            elif operator == ValidationOperator.NOT_EMPTY:
                if not value or (hasattr(value, '__len__') and len(value) == 0):
                    result.add_error(
                        "contract",
                        f"Field '{field}' in {agent}: is empty",
                        {"agent": agent, "field": field}
                    )
        
        except Exception as e:
            result.add_error(
                "contract",
                f"Validation error for field '{field}' in {agent}: {str(e)}",
                {"agent": agent, "field": field, "error": str(e)}
            )
    
    def _validate_policy_checks(
        self,
        tracer: AgentCallTracer,
        expected: ExpectedBehavior,
        result: ValidationResult
    ) -> None:
        """Validate policy checks."""
        for check in expected.policy_checks:
            if check == PolicyCheck.NO_EXECUTION_SIDE_EFFECTS:
                # Check that ExecutionBroker was not called (unless explicitly expected)
                if "ExecutionBroker" in tracer.get_agent_sequence():
                    if "ExecutionBroker" not in expected.agent_sequence:
                        result.add_error(
                            "policy",
                            "Execution side effects detected: ExecutionBroker was called",
                            {"check": check.value}
                        )
            
            elif check == PolicyCheck.NO_INFRASTRUCTURE_DECISIONS:
                # Check that InfrastructureExpert was not called (unless explicitly expected)
                if "InfrastructureExpert" in tracer.get_agent_sequence():
                    if "InfrastructureExpert" not in expected.agent_sequence:
                        result.add_error(
                            "policy",
                            "Infrastructure decisions detected: InfrastructureExpert was called",
                            {"check": check.value}
                        )
            
            elif check == PolicyCheck.STRICT_ROLE_SEPARATION:
                # This is validated by the handoff policy check above
                pass
            
            elif check == PolicyCheck.NO_CONTRACT_MUTATION:
                # Would require comparing input/output contracts
                # Skip for now - can be added later
                pass
            
            elif check == PolicyCheck.CONFIDENCE_THRESHOLDS:
                # Validate that confidence levels meet expectations
                if expected.intent and expected.intent.confidence_min:
                    intent_contracts = [
                        c for c in tracer.get_contracts()
                        if c["agent"] == "IntentDetector"
                    ]
                    if intent_contracts:
                        contract = intent_contracts[0]
                        if "confidence" in contract["data"]:
                            confidence = contract["data"]["confidence"]
                            if confidence < expected.intent.confidence_min:
                                result.add_error(
                                    "policy",
                                    f"Intent confidence {confidence} below minimum {expected.intent.confidence_min}",
                                    {
                                        "check": check.value,
                                        "confidence": confidence,
                                        "min_required": expected.intent.confidence_min
                                    }
                                )
    
    def _validate_semantic_invariants(
        self,
        tracer: AgentCallTracer,
        result: ValidationResult
    ) -> None:
        """
        Validate semantic invariants across contracts.
        
        These checks enforce separation of concerns and consistency:
        1. InfraDecision must not include execution commands (infra agent decides WHERE, not HOW)
        2. ExecutionToolSpec.infrastructure must match InfraDecision.infrastructure (consistency)
        3. Low confidence or warnings should trigger appropriate handling
        
        These invariants prevent agents from overstepping their responsibilities.
        """
        contracts = tracer.get_contracts()
        
        # Find InfraDecision and ExecutionToolSpec contracts
        infra_decision = None
        execution_spec = None
        
        for contract in contracts:
            if contract["contract_type"] == "InfraDecision" or contract["agent"] in ["InfrastructureExpert", "InfrastructureDecisionAgent"]:
                infra_decision = contract
            elif contract["contract_type"] == "ExecutionToolSpec" or contract["agent"] in ["ImplementationAgent", "CodeGenerator"]:
                execution_spec = contract
        
        # Invariant 1: InfraDecision must not include execution commands
        if infra_decision and "data" in infra_decision:
            infra_data = infra_decision["data"]
            
            # Check for prohibited fields that belong in ExecutionToolSpec, not InfraDecision
            prohibited_fields = [
                "execution_command",
                "execution_commands",
                "commands",
                "command",
                "script",
                "container_spec",
                "docker_image",
                "retry_policy",
                "timeout",
            ]
            
            found_prohibited = []
            for field in prohibited_fields:
                if field in infra_data:
                    found_prohibited.append(field)
            
            if found_prohibited:
                result.add_error(
                    "semantic_invariant",
                    f"InfraDecision contains execution details (violations: {', '.join(found_prohibited)}). "
                    "Infrastructure agent must only decide WHERE to execute, not HOW.",
                    {
                        "agent": infra_decision["agent"],
                        "prohibited_fields": found_prohibited,
                        "principle": "Separation of Concerns: InfraDecision (WHERE) vs ExecutionToolSpec (HOW)"
                    }
                )
            
            # Check for execution-related strings in reasoning/summary
            if "reasoning" in infra_data:
                reasoning = str(infra_data["reasoning"]).lower()
                execution_keywords = ["docker run", "container", "script", "command", "execute", "retry"]
                found_keywords = [kw for kw in execution_keywords if kw in reasoning]
                if found_keywords:
                    result.add_warning(
                        "semantic_invariant",
                        f"InfraDecision reasoning mentions execution details: {', '.join(found_keywords)}. "
                        "Consider keeping reasoning focused on infrastructure choice.",
                        {
                            "agent": infra_decision["agent"],
                            "keywords": found_keywords
                        }
                    )
        
        # Invariant 2: ExecutionToolSpec.infrastructure must match InfraDecision.infrastructure
        if infra_decision and execution_spec:
            infra_data = infra_decision.get("data", {})
            exec_data = execution_spec.get("data", {})
            
            # Extract infrastructure choices
            infra_choice = None
            exec_choice = None
            
            # Try various field names (handle both OLD and NEW schemas)
            for field in ["infrastructure", "recommended_environment", "environment"]:
                if field in infra_data:
                    infra_choice = infra_data[field]
                    break
            
            for field in ["infrastructure", "environment"]:
                if field in exec_data:
                    exec_choice = exec_data[field]
                    break
            
            if infra_choice and exec_choice:
                # Normalize case for comparison
                infra_choice_norm = str(infra_choice).upper()
                exec_choice_norm = str(exec_choice).upper()
                
                if infra_choice_norm != exec_choice_norm:
                    result.add_error(
                        "semantic_invariant",
                        f"Infrastructure mismatch: InfraDecision recommends '{infra_choice}' but "
                        f"ExecutionToolSpec uses '{exec_choice}'. These must match.",
                        {
                            "infra_agent": infra_decision["agent"],
                            "exec_agent": execution_spec["agent"],
                            "infra_choice": infra_choice,
                            "exec_choice": exec_choice,
                            "principle": "Consistency: Implementation must respect infrastructure decision"
                        }
                    )
        
        # Invariant 3: Low confidence or high warning counts should trigger clarification
        for contract in contracts:
            if "data" in contract:
                data = contract["data"]
                agent = contract["agent"]
                
                # Check confidence scores
                confidence = None
                for field in ["confidence_score", "confidence"]:
                    if field in data:
                        confidence = data[field]
                        break
                
                if confidence is not None:
                    # Low confidence threshold: < 0.5 suggests uncertainty
                    if confidence < 0.5:
                        result.add_warning(
                            "semantic_invariant",
                            f"{agent} has low confidence ({confidence:.2f} < 0.5). "
                            "Consider asking clarifying questions or showing warnings to user.",
                            {
                                "agent": agent,
                                "confidence": confidence,
                                "threshold": 0.5,
                                "principle": "Uncertainty Handling: Low confidence should trigger user interaction"
                            }
                        )
                    
                    # Very low confidence: < 0.3 is critical
                    if confidence < 0.3:
                        result.add_error(
                            "semantic_invariant",
                            f"{agent} has critically low confidence ({confidence:.2f} < 0.3). "
                            "Execution should be blocked or require explicit user confirmation.",
                            {
                                "agent": agent,
                                "confidence": confidence,
                                "threshold": 0.3,
                                "principle": "Safety: Very low confidence should block automatic execution"
                            }
                        )
                
                # Check warnings count
                if "warnings" in data and isinstance(data["warnings"], list):
                    warning_count = len(data["warnings"])
                    
                    # High warning count suggests problems
                    if warning_count >= 3:
                        result.add_warning(
                            "semantic_invariant",
                            f"{agent} produced {warning_count} warnings. "
                            "Consider surfacing these to the user or requiring acknowledgment.",
                            {
                                "agent": agent,
                                "warning_count": warning_count,
                                "warnings": data["warnings"][:3],  # Show first 3
                                "principle": "Transparency: Multiple warnings should be communicated to user"
                            }
                        )
    
    def _validate_timing(
        self,
        tracer: AgentCallTracer,
        expected: ExpectedBehavior,
        result: ValidationResult
    ) -> None:
        """Validate timing constraints."""
        total_duration = tracer.get_total_duration_ms()
        
        if expected.min_duration_ms and total_duration < expected.min_duration_ms:
            result.add_warning(
                "timing",
                f"Execution too fast: {total_duration}ms < minimum {expected.min_duration_ms}ms",
                {"actual": total_duration, "min": expected.min_duration_ms}
            )
        
        if expected.max_duration_ms and total_duration > expected.max_duration_ms:
            result.add_warning(
                "timing",
                f"Execution too slow: {total_duration}ms > maximum {expected.max_duration_ms}ms",
                {"actual": total_duration, "max": expected.max_duration_ms}
            )
