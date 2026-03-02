"""
Tests for semantic invariant validation.

Validates that the ContractValidator correctly enforces:
1. Separation of concerns (InfraDecision: WHERE, ExecutionToolSpec: HOW)
2. Infrastructure consistency (ExecutionToolSpec matches InfraDecision)
3. Confidence/warning thresholds
"""

import pytest
import sys
from pathlib import Path

# Add demo_scenarios to path
sys.path.insert(0, str(Path(__file__).parent.parent / "demo_scenarios"))

from framework.validator import ContractValidator, ValidationResult
from framework.tracer import AgentCallTracer


class TestSemanticInvariant1_SeparationOfConcerns:
    """Test that InfraDecision cannot include execution details."""
    
    def test_infra_decision_with_execution_command_fails(self):
        """InfraDecision with execution_command should fail validation."""
        tracer = AgentCallTracer()
        
        # Record InfraDecision with prohibited execution_command
        tracer.start_agent_call("InfrastructureExpert")
        tracer.record_contract({
            "agent": "InfrastructureExpert",
            "contract_type": "InfraDecision",
            "data": {
                "infrastructure": "EMR",
                "execution_command": "fastqc -o output input.fq",  # ❌ PROHIBITED
                "reasoning": "Large files require EMR"
            }
        })
        tracer.end_agent_call()
        
        # Validate
        validator = ContractValidator()
        from framework.scenario import ExpectedBehavior
        expected = ExpectedBehavior(
            agent_sequence=["InfrastructureExpert"],
            contracts=[],
            policy_checks=[]
        )
        
        result = validator.validate(tracer, expected)
        
        # Should fail
        assert not result.passed, "Should fail when InfraDecision includes execution_command"
        
        # Should have semantic_invariant error
        errors = [e for e in result.issues if e.category == "semantic_invariant" and e.severity == "error"]
        assert len(errors) > 0, "Should have semantic_invariant error"
        
        # Error message should mention execution details
        error_messages = [e.message for e in errors]
        assert any("execution" in msg.lower() for msg in error_messages), \
            f"Error should mention execution details, got: {error_messages}"
    
    def test_infra_decision_with_container_spec_fails(self):
        """InfraDecision with container_spec should fail validation."""
        tracer = AgentCallTracer()
        
        # Record InfraDecision with prohibited container_spec
        tracer.start_agent_call("InfrastructureExpert")
        tracer.record_contract({
            "agent": "InfrastructureExpert",
            "contract_type": "InfraDecision",
            "data": {
                "infrastructure": "EMR",
                "container_spec": {  # ❌ PROHIBITED
                    "image": "biocontainers/fastqc:0.11.9"
                },
                "reasoning": "Large files require EMR"
            }
        })
        tracer.end_agent_call()
        
        # Validate
        validator = ContractValidator()
        from framework.scenario import ExpectedBehavior
        expected = ExpectedBehavior(
            agent_sequence=["InfrastructureExpert"],
            contracts=[],
            policy_checks=[]
        )
        
        result = validator.validate(tracer, expected)
        
        # Should fail
        assert not result.passed, "Should fail when InfraDecision includes container_spec"
        
        # Should have semantic_invariant error
        errors = [e for e in result.issues if e.category == "semantic_invariant"]
        assert len(errors) > 0, "Should have semantic_invariant error"
    
    def test_infra_decision_without_execution_details_passes(self):
        """InfraDecision without execution details should pass validation."""
        tracer = AgentCallTracer()
        
        # Record valid InfraDecision
        tracer.start_agent_call("InfrastructureExpert")
        tracer.record_contract({
            "agent": "InfrastructureExpert",
            "contract_type": "InfraDecision",
            "data": {
                "infrastructure": "EMR",  # ✅ Only infrastructure choice
                "confidence_score": 0.85,
                "reasoning": "Large files require EMR"
            }
        })
        tracer.end_agent_call()
        
        # Validate
        validator = ContractValidator()
        from framework.scenario import ExpectedBehavior
        expected = ExpectedBehavior(
            agent_sequence=["InfrastructureExpert"],
            contracts=[],
            policy_checks=[]
        )
        
        result = validator.validate(tracer, expected)
        
        # Should pass (or only have warnings, not errors)
        errors = [e for e in result.issues if e.severity == "error"]
        assert len(errors) == 0, f"Should pass without errors, got: {[e.message for e in errors]}"


class TestSemanticInvariant2_InfrastructureConsistency:
    """Test that ExecutionToolSpec.infrastructure matches InfraDecision.infrastructure."""
    
    def test_infrastructure_mismatch_fails(self):
        """ExecutionToolSpec with different infrastructure should fail."""
        tracer = AgentCallTracer()
        
        # Record InfraDecision: EMR
        tracer.start_agent_call("InfrastructureExpert")
        tracer.record_contract({
            "agent": "InfrastructureExpert",
            "contract_type": "InfraDecision",
            "data": {
                "infrastructure": "EMR",
                "confidence_score": 0.85
            }
        })
        tracer.end_agent_call()
        
        # Record ExecutionToolSpec: Local (mismatch!)
        tracer.start_agent_call("ImplementationAgent")
        tracer.record_contract({
            "agent": "ImplementationAgent",
            "contract_type": "ExecutionToolSpec",
            "data": {
                "infrastructure": "Local",  # ❌ MISMATCH!
                "tool_name": "fastqc"
            }
        })
        tracer.end_agent_call()
        
        # Validate
        validator = ContractValidator()
        from framework.scenario import ExpectedBehavior
        expected = ExpectedBehavior(
            agent_sequence=["InfrastructureExpert", "ImplementationAgent"],
            contracts=[],
            policy_checks=[]
        )
        
        result = validator.validate(tracer, expected)
        
        # Should fail
        assert not result.passed, "Should fail when infrastructure mismatches"
        
        # Should have semantic_invariant error about mismatch
        errors = [e for e in result.issues if e.category == "semantic_invariant" and "mismatch" in e.message.lower()]
        assert len(errors) > 0, f"Should have infrastructure mismatch error, got: {[e.message for e in result.issues]}"
    
    def test_infrastructure_match_passes(self):
        """ExecutionToolSpec with matching infrastructure should pass."""
        tracer = AgentCallTracer()
        
        # Record InfraDecision: EMR
        tracer.start_agent_call("InfrastructureExpert")
        tracer.record_contract({
            "agent": "InfrastructureExpert",
            "contract_type": "InfraDecision",
            "data": {
                "infrastructure": "EMR",
                "confidence_score": 0.85
            }
        })
        tracer.end_agent_call()
        
        # Record ExecutionToolSpec: EMR (match!)
        tracer.start_agent_call("ImplementationAgent")
        tracer.record_contract({
            "agent": "ImplementationAgent",
            "contract_type": "ExecutionToolSpec",
            "data": {
                "infrastructure": "EMR",  # ✅ MATCHES!
                "tool_name": "fastqc",
                "confidence_score": 0.82
            }
        })
        tracer.end_agent_call()
        
        # Validate
        validator = ContractValidator()
        from framework.scenario import ExpectedBehavior
        expected = ExpectedBehavior(
            agent_sequence=["InfrastructureExpert", "ImplementationAgent"],
            contracts=[],
            policy_checks=[]
        )
        
        result = validator.validate(tracer, expected)
        
        # Should pass (or only have warnings, not errors related to infrastructure)
        errors = [e for e in result.issues if e.severity == "error" and "infrastructure" in e.message.lower()]
        assert len(errors) == 0, f"Should pass without infrastructure errors, got: {[e.message for e in errors]}"


class TestSemanticInvariant3_ConfidenceThresholds:
    """Test confidence and warning threshold validation."""
    
    def test_critically_low_confidence_fails(self):
        """Confidence < 0.3 should trigger error."""
        tracer = AgentCallTracer()
        
        # Record contract with very low confidence
        tracer.start_agent_call("InfrastructureExpert")
        tracer.record_contract({
            "agent": "InfrastructureExpert",
            "contract_type": "InfraDecision",
            "data": {
                "infrastructure": "EMR",
                "confidence_score": 0.25  # ❌ Too low!
            }
        })
        tracer.end_agent_call()
        
        # Validate
        validator = ContractValidator()
        from framework.scenario import ExpectedBehavior
        expected = ExpectedBehavior(
            agent_sequence=["InfrastructureExpert"],
            contracts=[],
            policy_checks=[]
        )
        
        result = validator.validate(tracer, expected)
        
        # Should fail
        assert not result.passed, "Should fail when confidence < 0.3"
        
        # Should have semantic_invariant error about confidence
        errors = [e for e in result.issues if e.category == "semantic_invariant" and "critically low confidence" in e.message]
        assert len(errors) > 0, f"Should have critically low confidence error, got: {[e.message for e in result.issues]}"
    
    def test_low_confidence_warning(self):
        """Confidence < 0.5 should trigger warning."""
        tracer = AgentCallTracer()
        
        # Record contract with low confidence
        tracer.start_agent_call("InfrastructureExpert")
        tracer.record_contract({
            "agent": "InfrastructureExpert",
            "contract_type": "InfraDecision",
            "data": {
                "infrastructure": "EMR",
                "confidence_score": 0.45  # Low but not critical
            }
        })
        tracer.end_agent_call()
        
        # Validate
        validator = ContractValidator()
        from framework.scenario import ExpectedBehavior
        expected = ExpectedBehavior(
            agent_sequence=["InfrastructureExpert"],
            contracts=[],
            policy_checks=[]
        )
        
        result = validator.validate(tracer, expected)
        
        # Should have warning (but not error)
        warnings = [w for w in result.issues if w.severity == "warning" and "low confidence" in w.message]
        assert len(warnings) > 0, "Should have low confidence warning"
        
        errors = [e for e in result.issues if e.severity == "error" and "confidence" in e.message]
        assert len(errors) == 0, "Should not have confidence error for 0.45"
    
    def test_many_warnings_triggers_warning(self):
        """3+ warnings should trigger validation warning."""
        tracer = AgentCallTracer()
        
        # Record contract with many warnings
        tracer.start_agent_call("InfrastructureExpert")
        tracer.record_contract({
            "agent": "InfrastructureExpert",
            "contract_type": "InfraDecision",
            "data": {
                "infrastructure": "EMR",
                "confidence_score": 0.6,
                "warnings": [
                    "File size unknown",
                    "Location unknown",
                    "Cannot estimate costs",
                    "Confidence low"
                ]  # 4 warnings >= 3
            }
        })
        tracer.end_agent_call()
        
        # Validate
        validator = ContractValidator()
        from framework.scenario import ExpectedBehavior
        expected = ExpectedBehavior(
            agent_sequence=["InfrastructureExpert"],
            contracts=[],
            policy_checks=[]
        )
        
        result = validator.validate(tracer, expected)
        
        # Should have warning about many warnings
        warnings = [w for w in result.issues if w.severity == "warning" and "warnings" in w.message.lower()]
        assert len(warnings) > 0, f"Should have warning about multiple warnings, got: {[w.message for w in result.issues]}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
