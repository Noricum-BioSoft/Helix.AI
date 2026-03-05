"""
Unit tests for agent responsibilities policy enforcement.

Tests that the policy checks and handoff rules from agents/agent-responsibilities.md
are correctly enforced. These tests should fail if any overlap or policy bypass is introduced.

Test categories:
1. Handoff policy legality (who can call whom)
2. Contract "must not" checks (what each agent may/may not output)
3. Happy-path sequences (valid workflows)
"""

import pytest
from typing import Dict, List
from backend.agent import HandoffPolicy, AgentRole, PolicyViolationError
from backend.policy_checks import (
    check_planner_output,
    check_plan_not_mutated,
    check_codegen_output,
    check_visualizer_output,
    serialize_plan_steps,
    PolicyViolationError as PolicyCheckError,
)
from shared.contracts import (
    WorkflowPlan,
    PlanMetadata,
    WorkflowInput,
    ExecutionSpec,
    ExecutionTarget,
    VisualizationArtifacts,
    VisualizationArtifact,
)


# =============================================================================
# Test 1: Handoff Policy Legality
# =============================================================================

class TestHandoffPolicyLegality:
    """Test that handoff rules are enforced per agent-responsibilities.md."""
    
    def setup_method(self):
        """Create policy instance for each test."""
        self.policy = HandoffPolicy()
    
    def test_intent_detector_must_be_first(self):
        """Intent Detector must always be the first agent."""
        # Valid: IntentDetector is first
        self.policy.validate_workflow_sequence(
            [AgentRole.INTENT_DETECTOR, AgentRole.GURU],
            intent="ask"
        )
        
        # Invalid: Any other agent is first
        with pytest.raises(PolicyViolationError, match="First agent must be IntentDetector"):
            self.policy.validate_workflow_sequence(
                [AgentRole.GURU, AgentRole.PLANNER],
                intent="ask"
            )
    
    def test_ask_intent_routes_to_guru(self):
        """Intent 'ask' must route to BioinformaticsGuru."""
        assert self.policy.get_next_agent_for_intent("ask") == AgentRole.GURU
        assert self.policy.get_next_agent_for_intent("qa") == AgentRole.GURU
        
        # Valid sequence
        self.policy.validate_workflow_sequence(
            [AgentRole.INTENT_DETECTOR, AgentRole.GURU],
            intent="ask"
        )
        
        # Invalid: ask intent goes to Planner instead of Guru
        with pytest.raises(PolicyViolationError, match="second agent must be BioinformaticsGuru"):
            self.policy.validate_workflow_sequence(
                [AgentRole.INTENT_DETECTOR, AgentRole.PLANNER],
                intent="ask"
            )
    
    def test_execute_intent_routes_to_planner(self):
        """Intent 'execute' must route to BioinformaticsExecutor (Planner)."""
        assert self.policy.get_next_agent_for_intent("execute") == AgentRole.PLANNER
        
        # Valid sequence
        self.policy.validate_workflow_sequence(
            [AgentRole.INTENT_DETECTOR, AgentRole.PLANNER],
            intent="execute"
        )
        
        # Invalid: execute intent goes to Guru instead of Planner
        with pytest.raises(PolicyViolationError, match="second agent must be BioinformaticsExecutor"):
            self.policy.validate_workflow_sequence(
                [AgentRole.INTENT_DETECTOR, AgentRole.GURU],
                intent="execute"
            )
    
    def test_guru_cannot_handoff_to_broker(self):
        """Guru must not directly invoke Execution Broker."""
        # Guru can escalate to Planner (with user consent)
        self.policy.validate_handoff(AgentRole.GURU, AgentRole.PLANNER, user_consent=True)
        
        # Guru cannot call Broker
        with pytest.raises(PolicyViolationError, match="Illegal handoff"):
            self.policy.validate_handoff(AgentRole.GURU, AgentRole.BROKER)
    
    def test_infra_must_not_be_invoked_before_planner(self):
        """Infrastructure Expert must not be called before Planner produces a plan."""
        # Valid: Planner → Infra
        self.policy.validate_handoff(AgentRole.PLANNER, AgentRole.INFRA)
        
        # Invalid: Intent → Infra (skipping Planner)
        with pytest.raises(PolicyViolationError, match="Illegal handoff"):
            self.policy.validate_handoff(AgentRole.INTENT_DETECTOR, AgentRole.INFRA)
    
    def test_broker_must_be_after_infra_or_codegen(self):
        """Execution Broker must be invoked after infrastructure selection."""
        # Valid: Infra → Broker
        self.policy.validate_handoff(AgentRole.INFRA, AgentRole.BROKER)
        
        # Valid: CodeGen → Broker
        self.policy.validate_handoff(AgentRole.CODEGEN, AgentRole.BROKER)
        
        # Invalid: Planner → Broker (skipping Infra)
        with pytest.raises(PolicyViolationError, match="Illegal handoff"):
            self.policy.validate_handoff(AgentRole.PLANNER, AgentRole.BROKER)
    
    def test_visualizer_must_be_after_broker(self):
        """Visualizer must be invoked after Execution Broker produces results."""
        # Valid: Broker → Visualizer
        self.policy.validate_handoff(AgentRole.BROKER, AgentRole.VISUALIZER)
        
        # Invalid: Planner → Visualizer (skipping execution)
        with pytest.raises(PolicyViolationError, match="Illegal handoff"):
            self.policy.validate_handoff(AgentRole.PLANNER, AgentRole.VISUALIZER)
    
    def test_codegen_is_optional(self):
        """Code Generator can be skipped if not needed."""
        # Valid: Infra → Broker (skip CodeGen)
        self.policy.validate_handoff(AgentRole.INFRA, AgentRole.BROKER)
        
        # Valid: Infra → CodeGen → Broker
        self.policy.validate_handoff(AgentRole.INFRA, AgentRole.CODEGEN)
        self.policy.validate_handoff(AgentRole.CODEGEN, AgentRole.BROKER)
    
    def test_guru_escalation_requires_user_consent(self):
        """Guru → Planner escalation requires explicit user consent."""
        # Without consent: should fail
        with pytest.raises(PolicyViolationError, match="requires explicit user consent"):
            self.policy.validate_handoff(AgentRole.GURU, AgentRole.PLANNER, user_consent=False)
        
        # With consent: should pass
        self.policy.validate_handoff(AgentRole.GURU, AgentRole.PLANNER, user_consent=True)


# =============================================================================
# Test 2: Contract "Must Not" Checks
# =============================================================================

class TestPlannerContractChecks:
    """Test that Planner cannot output infrastructure decisions."""
    
    def test_planner_valid_output(self):
        """Valid planner output should pass validation."""
        valid_plan = WorkflowPlan(
            version="v1",
            metadata=PlanMetadata(scale="medium"),
            inputs=[WorkflowInput(uri="s3://bucket/input.fastq")],
            steps=[
                {"id": "step1", "tool_name": "fastqc", "arguments": {"input": "s3://bucket/input.fastq"}},
                {"id": "step2", "tool_name": "read_merging", "arguments": {"input_r1": "$ref:step1.result"}},
            ]
        )
        
        # Should not raise
        check_planner_output(valid_plan)
    
    def test_planner_cannot_choose_environment(self):
        """Planner output with 'environment' field should fail."""
        invalid_plan = {
            "version": "v1",
            "steps": [{"id": "step1", "tool_name": "fastqc", "arguments": {}}],
            "environment": "EMR",  # NOT ALLOWED
        }
        
        with pytest.raises(PolicyCheckError, match="choose infrastructure"):
            check_planner_output(invalid_plan)
    
    def test_planner_cannot_include_instance_type(self):
        """Planner output with instance types should fail."""
        invalid_plan = {
            "version": "v1",
            "steps": [
                {
                    "id": "step1",
                    "tool_name": "fastqc",
                    "arguments": {},
                    "instance_type": "m5.xlarge"  # NOT ALLOWED
                }
            ],
        }
        
        with pytest.raises(PolicyCheckError, match="Infrastructure keyword"):
            check_planner_output(invalid_plan)
    
    def test_planner_cannot_include_emr_config(self):
        """Planner output with EMR configuration should fail."""
        invalid_plan = {
            "version": "v1",
            "steps": [{"id": "step1", "tool_name": "fastqc", "arguments": {}}],
            "emr_steps": [  # NOT ALLOWED
                {"Name": "Run FastQC", "ActionOnFailure": "CONTINUE"}
            ],
        }
        
        with pytest.raises(PolicyCheckError, match="Infrastructure keyword"):
            check_planner_output(invalid_plan)
    
    def test_planner_cannot_set_execution_target(self):
        """Planner output with execution target should fail."""
        invalid_plan = {
            "version": "v1",
            "steps": [{"id": "step1", "tool_name": "fastqc", "arguments": {}}],
            "target": "EC2",  # NOT ALLOWED
        }
        
        with pytest.raises(PolicyCheckError, match="Execution target"):
            check_planner_output(invalid_plan)


class TestInfrastructureExpertContractChecks:
    """Test that Infrastructure Expert cannot mutate plan steps."""
    
    def test_infra_cannot_mutate_steps(self):
        """Infrastructure Expert must not modify workflow steps."""
        original_steps = [
            {"id": "step1", "tool_name": "fastqc", "arguments": {"input": "file.fastq"}},
            {"id": "step2", "tool_name": "read_merging", "arguments": {}},
        ]
        
        mutated_steps = [
            {"id": "step1", "tool_name": "fastqc", "arguments": {"input": "file.fastq"}},
            {"id": "step2", "tool_name": "DIFFERENT_TOOL", "arguments": {}},  # Changed!
        ]
        
        with pytest.raises(PolicyCheckError, match="must not mutate workflow plan steps"):
            check_plan_not_mutated(original_steps, mutated_steps, "InfrastructureExpert")
    
    def test_infra_cannot_add_steps(self):
        """Infrastructure Expert must not add steps to the plan."""
        original_steps = [
            {"id": "step1", "tool_name": "fastqc", "arguments": {}},
        ]
        
        steps_with_addition = [
            {"id": "step1", "tool_name": "fastqc", "arguments": {}},
            {"id": "step2", "tool_name": "new_step", "arguments": {}},  # Added!
        ]
        
        with pytest.raises(PolicyCheckError, match="must not mutate workflow plan steps"):
            check_plan_not_mutated(original_steps, steps_with_addition, "InfrastructureExpert")
    
    def test_infra_cannot_remove_steps(self):
        """Infrastructure Expert must not remove steps from the plan."""
        original_steps = [
            {"id": "step1", "tool_name": "fastqc", "arguments": {}},
            {"id": "step2", "tool_name": "read_merging", "arguments": {}},
        ]
        
        steps_with_removal = [
            {"id": "step1", "tool_name": "fastqc", "arguments": {}},
            # step2 removed!
        ]
        
        with pytest.raises(PolicyCheckError, match="must not mutate workflow plan steps"):
            check_plan_not_mutated(original_steps, steps_with_removal, "InfrastructureExpert")
    
    def test_infra_unchanged_steps_pass(self):
        """Infrastructure Expert with unchanged steps should pass."""
        steps = [
            {"id": "step1", "tool_name": "fastqc", "arguments": {"input": "file.fastq"}},
            {"id": "step2", "tool_name": "read_merging", "arguments": {}},
        ]
        
        # Should not raise
        check_plan_not_mutated(steps, steps, "InfrastructureExpert")
    
    def test_serialization_handles_pydantic_models(self):
        """serialize_plan_steps should handle Pydantic models."""
        from backend.plan_ir import PlanStep
        
        steps = [
            PlanStep(id="step1", tool_name="fastqc", arguments={"input": "file.fastq"}),
            PlanStep(id="step2", tool_name="read_merging", arguments={}),
        ]
        
        # Should not raise
        serialized = serialize_plan_steps(steps)
        assert "fastqc" in serialized
        assert "read_merging" in serialized


class TestCodeGeneratorContractChecks:
    """Test that Code Generator cannot change scientific intent."""
    
    def test_codegen_valid_output(self):
        """Valid CodeGen output should pass validation."""
        original_plan = WorkflowPlan(
            version="v1",
            steps=[
                {"id": "step1", "tool_name": "fastqc", "arguments": {}},
            ]
        )
        
        execution_spec = ExecutionSpec(
            target=ExecutionTarget.EMR,
            container_image="my-image:latest",
            command=["python", "run_workflow.py"],
        )
        
        # Should not raise
        check_codegen_output(original_plan, execution_spec)
    
    def test_codegen_with_none_output_passes(self):
        """CodeGen with None output (not needed) should pass."""
        original_plan = WorkflowPlan(
            version="v1",
            steps=[{"id": "step1", "tool_name": "fastqc", "arguments": {}}]
        )
        
        # Should not raise
        check_codegen_output(original_plan, None)
    
    def test_codegen_cannot_rewrite_steps(self):
        """Code Generator must not rewrite workflow steps."""
        original_plan = {
            "version": "v1",
            "steps": [
                {"id": "step1", "tool_name": "fastqc", "arguments": {}},
            ]
        }
        
        # ExecutionSpec should NOT have 'steps' field (that's a plan, not a spec)
        execution_spec_with_steps = {
            "target": "emr",
            "steps": [  # NOT ALLOWED in ExecutionSpec
                {"id": "step1", "tool_name": "DIFFERENT_TOOL", "arguments": {}},
            ]
        }
        
        with pytest.raises(PolicyCheckError, match="must not change scientific intent"):
            check_codegen_output(original_plan, execution_spec_with_steps)
    
    def test_codegen_cannot_choose_infrastructure(self):
        """Code Generator must not add infrastructure decision fields."""
        original_plan = WorkflowPlan(
            version="v1",
            steps=[{"id": "step1", "tool_name": "fastqc", "arguments": {}}]
        )
        
        # ExecutionSpec with infra fields
        execution_spec_with_infra = {
            "target": "emr",
            "command": ["python", "run.py"],
            "instance_type": "m5.xlarge",  # NOT ALLOWED (should come from InfraDecision)
        }
        
        with pytest.raises(PolicyCheckError, match="must not choose infrastructure"):
            check_codegen_output(original_plan, execution_spec_with_infra)


class TestVisualizerContractChecks:
    """Test that Visualizer cannot introduce workflow steps or infra changes."""
    
    def test_visualizer_valid_output(self):
        """Valid Visualizer output should pass validation."""
        valid_output = VisualizationArtifacts(
            summary="QC analysis complete",
            artifacts=[
                VisualizationArtifact(
                    artifact_type="plot",
                    uri="s3://bucket/qc_plot.png",
                    description="Quality scores"
                )
            ]
        )
        
        # Should not raise
        check_visualizer_output(valid_output)
    
    def test_visualizer_with_none_output_passes(self):
        """Visualizer with None output should pass."""
        # Should not raise
        check_visualizer_output(None)
    
    def test_visualizer_cannot_include_steps(self):
        """Visualizer must not include workflow steps in output."""
        invalid_output = {
            "summary": "Results",
            "artifacts": [],
            "steps": [  # NOT ALLOWED
                {"id": "step1", "tool_name": "fastqc", "arguments": {}}
            ]
        }
        
        with pytest.raises(PolicyCheckError, match="must not redesign workflow"):
            check_visualizer_output(invalid_output)
    
    def test_visualizer_cannot_choose_environment(self):
        """Visualizer must not include infrastructure decisions."""
        invalid_output = {
            "summary": "Results",
            "artifacts": [],
            "environment": "EMR",  # NOT ALLOWED
        }
        
        with pytest.raises(PolicyCheckError, match="must not choose infrastructure"):
            check_visualizer_output(invalid_output)
    
    def test_visualizer_cannot_include_commands(self):
        """Visualizer must not generate execution code."""
        invalid_output = {
            "summary": "Results",
            "artifacts": [],
            "command": ["python", "run.py"],  # NOT ALLOWED
        }
        
        with pytest.raises(PolicyCheckError, match="must not generate execution code"):
            check_visualizer_output(invalid_output)


# =============================================================================
# Test 3: Happy-Path Sequences
# =============================================================================

class TestHappyPathSequences:
    """Test that valid workflow sequences pass all checks."""
    
    def setup_method(self):
        """Create policy instance for each test."""
        self.policy = HandoffPolicy()
    
    def test_ask_workflow_happy_path(self):
        """Valid ask workflow: Intent → Guru."""
        sequence = [
            AgentRole.INTENT_DETECTOR,
            AgentRole.GURU,
        ]
        
        # Should not raise
        self.policy.validate_workflow_sequence(sequence, intent="ask")
    
    def test_execute_workflow_minimal_happy_path(self):
        """Valid minimal execute workflow: Intent → Planner → Infra → Broker → Visualizer."""
        sequence = [
            AgentRole.INTENT_DETECTOR,
            AgentRole.PLANNER,
            AgentRole.INFRA,
            AgentRole.BROKER,
            AgentRole.VISUALIZER,
        ]
        
        # Should not raise
        self.policy.validate_workflow_sequence(sequence, intent="execute")
    
    def test_execute_workflow_with_codegen_happy_path(self):
        """Valid execute workflow with CodeGen: Intent → Planner → Infra → CodeGen → Broker → Visualizer."""
        sequence = [
            AgentRole.INTENT_DETECTOR,
            AgentRole.PLANNER,
            AgentRole.INFRA,
            AgentRole.CODEGEN,
            AgentRole.BROKER,
            AgentRole.VISUALIZER,
        ]
        
        # Should not raise
        self.policy.validate_workflow_sequence(sequence, intent="execute")
    
    def test_full_execute_workflow_end_to_end(self):
        """End-to-end execute workflow with all contract checks."""
        # 1. Intent Detector
        sequence = [AgentRole.INTENT_DETECTOR]
        
        # 2. Planner produces a plan
        sequence.append(AgentRole.PLANNER)
        
        planner_output = WorkflowPlan(
            version="v1",
            metadata=PlanMetadata(scale="medium"),
            inputs=[WorkflowInput(uri="s3://bucket/input.fastq")],
            steps=[
                {"id": "step1", "tool_name": "fastqc", "arguments": {"input": "s3://bucket/input.fastq"}},
            ]
        )
        
        # Check planner output
        check_planner_output(planner_output)
        
        # 3. Infrastructure Expert selects infrastructure (without mutating plan)
        sequence.append(AgentRole.INFRA)
        
        steps_before_infra = planner_output.steps
        steps_after_infra = planner_output.steps  # Unchanged
        
        # Check infra didn't mutate plan
        check_plan_not_mutated(steps_before_infra, steps_after_infra, "InfrastructureExpert")
        
        # 4. Code Generator produces ExecutionSpec
        sequence.append(AgentRole.CODEGEN)
        
        execution_spec = ExecutionSpec(
            target=ExecutionTarget.EMR,
            container_image="fastqc:latest",
            command=["fastqc", "--input", "s3://bucket/input.fastq"],
        )
        
        # Check codegen output
        check_codegen_output(planner_output, execution_spec)
        
        # 5. Execution Broker runs the job
        sequence.append(AgentRole.BROKER)
        
        # 6. Visualizer produces artifacts
        sequence.append(AgentRole.VISUALIZER)
        
        visualizer_output = VisualizationArtifacts(
            summary="QC complete",
            artifacts=[
                VisualizationArtifact(
                    artifact_type="plot",
                    uri="s3://bucket/qc_plot.png"
                )
            ]
        )
        
        # Check visualizer output
        check_visualizer_output(visualizer_output)
        
        # Validate full sequence
        self.policy.validate_workflow_sequence(sequence, intent="execute")


# =============================================================================
# Test 4: Edge Cases and Error Messages
# =============================================================================

class TestEdgeCasesAndErrorMessages:
    """Test edge cases and verify error messages are actionable."""
    
    def test_error_messages_include_fix_instructions(self):
        """Policy violation errors should include 'How to fix' instructions."""
        invalid_plan = {
            "version": "v1",
            "steps": [],
            "environment": "EMR",
        }
        
        try:
            check_planner_output(invalid_plan)
            pytest.fail("Should have raised PolicyCheckError")
        except PolicyCheckError as e:
            error_msg = str(e)
            assert "How to fix:" in error_msg
            assert "Infrastructure Expert" in error_msg
    
    def test_empty_agent_sequence_fails(self):
        """Empty agent sequence should fail validation."""
        policy = HandoffPolicy()
        
        with pytest.raises(PolicyViolationError, match="Empty agent sequence"):
            policy.validate_workflow_sequence([], intent="ask")
    
    def test_unknown_intent_fails(self):
        """Unknown intent should fail with clear error."""
        policy = HandoffPolicy()
        
        with pytest.raises(PolicyViolationError, match="Unknown intent"):
            policy.get_next_agent_for_intent("unknown_intent")
    
    def test_plan_mutation_diff_in_error_message(self):
        """Plan mutation error should include diff for debugging."""
        original_steps = [{"id": "step1", "tool_name": "tool_a", "arguments": {}}]
        mutated_steps = [{"id": "step1", "tool_name": "tool_b", "arguments": {}}]
        
        try:
            check_plan_not_mutated(original_steps, mutated_steps, "InfrastructureExpert")
            pytest.fail("Should have raised PolicyCheckError")
        except PolicyCheckError as e:
            error_msg = str(e)
            # Should include some diff information
            assert "tool_a" in error_msg or "tool_b" in error_msg or "Line" in error_msg


# =============================================================================
# Test 5: Integration with Existing Code
# =============================================================================

class TestIntegrationWithExistingCode:
    """Test that policy checks integrate with existing agent code."""
    
    def test_handoff_policy_accessible_from_agent_module(self):
        """HandoffPolicy should be importable from backend.agent."""
        from backend.agent import HandoffPolicy, AgentRole, PolicyViolationError
        
        policy = HandoffPolicy()
        assert policy is not None
        assert hasattr(policy, 'validate_handoff')
        assert hasattr(policy, 'validate_workflow_sequence')
    
    def test_policy_checks_use_same_error_type(self):
        """Policy checks should use consistent error types."""
        from backend.policy_checks import PolicyViolationError as CheckError
        from backend.agent import PolicyViolationError as AgentError
        
        # Should be the same class or compatible
        # (We define it in both modules for now, but they're the same exception)
        assert CheckError.__name__ == AgentError.__name__
    
    def test_contracts_importable_from_shared(self):
        """All contracts should be importable from shared.contracts."""
        from shared.contracts import (
            IntentResult,
            WorkflowPlan,
            InfraDecisionContract,
            ExecutionSpec,
            ExecutionResult,
            VisualizationArtifacts,
        )
        
        # Should all be Pydantic models
        assert hasattr(IntentResult, 'parse_obj')
        assert hasattr(WorkflowPlan, 'parse_obj')
        assert hasattr(ExecutionSpec, 'parse_obj')
