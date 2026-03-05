"""
Tests for agent handoff policy enforcement.

Validates that the handoff policy correctly enforces routing rules from
agent-responsibilities.md.
"""

import pytest
from backend.agent import (
    HandoffPolicy,
    PolicyViolationError,
    AgentRole,
)


class TestHandoffPolicy:
    """Test agent handoff policy enforcement."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.policy = HandoffPolicy()
    
    def test_intent_detector_must_be_first(self):
        """Intent Detector must always be the first agent."""
        # Valid: Intent Detector is first
        sequence = [AgentRole.INTENT_DETECTOR, AgentRole.GURU]
        self.policy.validate_workflow_sequence(sequence, intent="ask")
        
        # Invalid: Starting with another agent
        with pytest.raises(PolicyViolationError, match="First agent must be IntentDetector"):
            sequence = [AgentRole.GURU, AgentRole.INTENT_DETECTOR]
            self.policy.validate_workflow_sequence(sequence, intent="ask")
    
    def test_ask_intent_routes_to_guru(self):
        """Ask intent must route to Guru."""
        next_agent = self.policy.get_next_agent_for_intent("ask")
        assert next_agent == AgentRole.GURU
        
        next_agent = self.policy.get_next_agent_for_intent("qa")
        assert next_agent == AgentRole.GURU
    
    def test_execute_intent_routes_to_planner(self):
        """Execute intent must route to Planner."""
        next_agent = self.policy.get_next_agent_for_intent("execute")
        assert next_agent == AgentRole.PLANNER
    
    def test_unknown_intent_raises_error(self):
        """Unknown intent should raise PolicyViolationError."""
        with pytest.raises(PolicyViolationError, match="Unknown intent"):
            self.policy.get_next_agent_for_intent("unknown")
    
    def test_ask_workflow_sequence(self):
        """Valid ask workflow: IntentDetector → Guru."""
        sequence = [AgentRole.INTENT_DETECTOR, AgentRole.GURU]
        # Should not raise
        self.policy.validate_workflow_sequence(sequence, intent="ask")
    
    def test_execute_workflow_sequence(self):
        """Valid execute workflow: IntentDetector → Planner → Infra → Broker → Visualizer."""
        sequence = [
            AgentRole.INTENT_DETECTOR,
            AgentRole.PLANNER,
            AgentRole.INFRA,
            AgentRole.BROKER,
            AgentRole.VISUALIZER,
        ]
        # Should not raise
        self.policy.validate_workflow_sequence(sequence, intent="execute")
    
    def test_execute_workflow_with_codegen(self):
        """Valid execute workflow with CodeGen: Planner → Infra → CodeGen → Broker → Visualizer."""
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
    
    def test_guru_cannot_call_broker(self):
        """Guru cannot directly call Broker (illegal handoff)."""
        with pytest.raises(PolicyViolationError, match="Illegal handoff.*BioinformaticsGuru.*ExecutionBroker"):
            self.policy.validate_handoff(
                from_agent=AgentRole.GURU,
                to_agent=AgentRole.BROKER
            )
    
    def test_infra_cannot_call_broker_directly(self):
        """
        Infra can call Broker directly (CodeGen is optional).
        
        This is actually allowed per the policy.
        """
        # Should NOT raise - Infra can skip CodeGen and go to Broker
        self.policy.validate_handoff(
            from_agent=AgentRole.INFRA,
            to_agent=AgentRole.BROKER
        )
    
    def test_infra_cannot_skip_to_visualizer(self):
        """Infra cannot skip Broker and go to Visualizer."""
        with pytest.raises(PolicyViolationError, match="Illegal handoff.*InfrastructureExpert.*DataVisualizer"):
            self.policy.validate_handoff(
                from_agent=AgentRole.INFRA,
                to_agent=AgentRole.VISUALIZER
            )
    
    def test_planner_must_go_to_infra(self):
        """Planner must hand off to Infra next."""
        # Valid: Planner → Infra
        self.policy.validate_handoff(
            from_agent=AgentRole.PLANNER,
            to_agent=AgentRole.INFRA
        )
        
        # Invalid: Planner → Broker (skipping Infra)
        with pytest.raises(PolicyViolationError, match="Illegal handoff.*BioinformaticsExecutor.*ExecutionBroker"):
            self.policy.validate_handoff(
                from_agent=AgentRole.PLANNER,
                to_agent=AgentRole.BROKER
            )
    
    def test_guru_escalation_requires_user_consent(self):
        """Guru → Planner escalation requires explicit user consent."""
        # Without consent: should raise
        with pytest.raises(PolicyViolationError, match="requires explicit user consent"):
            self.policy.validate_handoff(
                from_agent=AgentRole.GURU,
                to_agent=AgentRole.PLANNER,
                user_consent=False
            )
        
        # With consent: should pass
        self.policy.validate_handoff(
            from_agent=AgentRole.GURU,
            to_agent=AgentRole.PLANNER,
            user_consent=True
        )
    
    def test_get_allowed_next_agents(self):
        """Test getting allowed next agents."""
        # Intent Detector can go to Guru or Planner
        allowed = self.policy.get_allowed_next_agents(AgentRole.INTENT_DETECTOR)
        assert AgentRole.GURU in allowed
        assert AgentRole.PLANNER in allowed
        
        # Planner must go to Infra
        allowed = self.policy.get_allowed_next_agents(AgentRole.PLANNER)
        assert allowed == [AgentRole.INFRA]
        
        # Visualizer is terminal (no next agents)
        allowed = self.policy.get_allowed_next_agents(AgentRole.VISUALIZER)
        assert allowed == []
    
    def test_visualizer_is_terminal(self):
        """Visualizer should not hand off to any other agent."""
        allowed = self.policy.get_allowed_next_agents(AgentRole.VISUALIZER)
        assert len(allowed) == 0
        
        # Try to call any agent from Visualizer - should all fail
        with pytest.raises(PolicyViolationError):
            self.policy.validate_handoff(
                from_agent=AgentRole.VISUALIZER,
                to_agent=AgentRole.BROKER
            )
    
    def test_codegen_must_go_to_broker(self):
        """CodeGen must hand off to Broker."""
        # Valid: CodeGen → Broker
        self.policy.validate_handoff(
            from_agent=AgentRole.CODEGEN,
            to_agent=AgentRole.BROKER
        )
        
        # Invalid: CodeGen → Visualizer (skipping Broker)
        with pytest.raises(PolicyViolationError):
            self.policy.validate_handoff(
                from_agent=AgentRole.CODEGEN,
                to_agent=AgentRole.VISUALIZER
            )
    
    def test_empty_sequence_raises_error(self):
        """Empty agent sequence should raise error."""
        with pytest.raises(PolicyViolationError, match="Empty agent sequence"):
            self.policy.validate_workflow_sequence([], intent="ask")
    
    def test_invalid_second_agent_for_intent(self):
        """Second agent must match the intent."""
        # Execute intent should go to Planner, not Guru
        with pytest.raises(PolicyViolationError, match="second agent must be"):
            sequence = [AgentRole.INTENT_DETECTOR, AgentRole.GURU]
            self.policy.validate_workflow_sequence(sequence, intent="execute")
        
        # Ask intent should go to Guru, not Planner
        with pytest.raises(PolicyViolationError, match="second agent must be"):
            sequence = [AgentRole.INTENT_DETECTOR, AgentRole.PLANNER]
            self.policy.validate_workflow_sequence(sequence, intent="ask")


class TestCommandProcessorHandoffPolicy:
    """Test CommandProcessor integration with handoff policy."""
    
    @pytest.mark.asyncio
    async def test_process_validates_intent_detector_first(self):
        """CommandProcessor should validate IntentDetector is called first."""
        from backend.agent import CommandProcessor
        from unittest.mock import Mock
        
        # Mock intent classifier
        mock_intent = Mock()
        mock_intent.intent = "qa"
        mock_intent.reason = "test"
        
        def mock_classifier(text):
            return mock_intent
        
        processor = CommandProcessor(
            agent=None,
            intent_classifier=mock_classifier,
            is_mock_mode=True
        )
        
        # Process a command - should register IntentDetector first
        result = await processor.process("What is DNA?", "test_session")
        
        # Verify IntentDetector was registered first
        assert len(processor.agent_sequence) >= 1
        assert processor.agent_sequence[0] == AgentRole.INTENT_DETECTOR
    
    @pytest.mark.asyncio
    async def test_process_validates_intent_routing(self):
        """CommandProcessor should route to correct agent based on intent."""
        from backend.agent import CommandProcessor
        from unittest.mock import Mock
        
        # Test QA intent → Guru
        mock_intent = Mock()
        mock_intent.intent = "qa"
        mock_intent.reason = "test"
        
        def mock_classifier(text):
            return mock_intent
        
        processor = CommandProcessor(
            agent=None,
            intent_classifier=mock_classifier,
            is_mock_mode=True
        )
        
        result = await processor.process("What is DNA?", "test_session")
        
        # Verify IntentDetector → Guru
        assert len(processor.agent_sequence) >= 2
        assert processor.agent_sequence[0] == AgentRole.INTENT_DETECTOR
        assert processor.agent_sequence[1] == AgentRole.GURU
    
    @pytest.mark.asyncio
    async def test_process_validates_execute_routing(self):
        """CommandProcessor should route execute intent to Planner."""
        from backend.agent import CommandProcessor
        from unittest.mock import Mock
        
        # Test execute intent → Planner
        mock_intent = Mock()
        mock_intent.intent = "execute"
        mock_intent.reason = "test"
        
        def mock_classifier(text):
            return mock_intent
        
        processor = CommandProcessor(
            agent=None,
            intent_classifier=mock_classifier,
            is_mock_mode=True
        )
        
        result = await processor.process("Run FastQC", "test_session")
        
        # Verify IntentDetector → Planner
        assert len(processor.agent_sequence) >= 2
        assert processor.agent_sequence[0] == AgentRole.INTENT_DETECTOR
        assert processor.agent_sequence[1] == AgentRole.PLANNER
