"""
Tests for Phase 3 Multi-Agent Components:
- ExecutionToolSpec contract validation
- ImplementationAgent
- Orchestrator

**SYSTEM TESTED:** EXPERIMENTAL Orchestrator (backend/orchestrator.py) 🟡

This module validates the PROTOTYPE 2-agent pipeline with Pydantic contracts
and provenance hashing. This system is NOT used in production.

For tests of the PRIMARY orchestrator (backend/agent.py) used in production,
see tests/demo_scenarios/test_scenarios.py.

See docs/ORCHESTRATION_DUALITY.md for details on the two orchestration systems.

---

These tests verify multi-agent coordination and observability.
"""

import pytest
import json
from unittest.mock import AsyncMock, patch

from backend.contracts.workflow_plan import WorkflowPlan, DataInput, OperationSpec
from backend.contracts.infra_decision import InfraDecision, FileAnalysis, ComputationalRequirements, CostAnalysis
from backend.contracts.execution_spec import (
    ExecutionToolSpec,
    ContainerSpec,
    CommandSpec,
    RetryPolicy,
    ResourceRequirements,
    OutputSpec
)
from backend.implementation_agent import plan_implementation, _heuristic_implementation_plan
from backend.orchestrator import Orchestrator, AgentInvocation, OrchestratorTrace


# Helper functions for creating valid test objects
def create_test_workflow_plan(description="Test workflow for testing purposes") -> WorkflowPlan:
    """Create a valid WorkflowPlan for testing."""
    return WorkflowPlan(
        description=description,
        data_inputs=[DataInput(uri="s3://bucket/test.fq", format="fastq")],
        operations=[OperationSpec(operation_name="test", tool_name="test", parameters={})],
        expected_data_volume_mb=100.0,
        expected_compute_intensity="Medium"
    )


def create_test_infra_decision(infrastructure="EMR") -> InfraDecision:
    """Create a valid InfraDecision for testing."""
    return InfraDecision(
        infrastructure=infrastructure,
        confidence_score=0.8,
        decision_summary="Infrastructure decision for testing with sufficient detail to meet minimum length requirements.",
        reasoning="This is a test infrastructure decision with detailed reasoning that meets the minimum character requirement for validation.",
        file_analysis=FileAnalysis(
            total_size_bytes=100 * 1024 * 1024,  # 100MB in bytes
            total_size_mb=100.0,  # 100MB
            file_count=1,
            largest_file_bytes=100 * 1024 * 1024,
            largest_file_mb=100.0,
            all_in_s3=True
        ),
        computational_requirements=ComputationalRequirements(
            estimated_cpu_hours=1.0,
            estimated_memory_gb=4.0,
            parallelizable=True
        ),
        cost_analysis=CostAnalysis(
            estimated_cost_range_usd=(1.0, 5.0),
            cost_class="Low",
            cost_assumptions="Test cost assumptions for validation purposes",
            cost_confidence=0.7
        )
    )


class TestExecutionToolSpecContract:
    """Test ExecutionToolSpec Pydantic validation."""
    
    def test_valid_execution_spec(self):
        """Valid execution spec should validate."""
        spec = ExecutionToolSpec(
            tool_name="fastqc",
            infrastructure="EMR",
            container_spec=ContainerSpec(
                image="biocontainers/fastqc:0.11.9",
                image_type="docker",
                pull_policy="IfNotPresent"
            ),
            commands=[
                CommandSpec(
                    name="run_fastqc",
                    command="fastqc -o /output /input_R1.fq",
                    inputs=["s3://bucket/input_R1.fq"],
                    outputs=["s3://bucket/output/fastqc.html"],
                    timeout_minutes=30.0
                )
            ],
            retry_policy=RetryPolicy(max_retries=2),
            resource_requirements=ResourceRequirements(
                min_cpu_cores=4.0,
                min_memory_gb=8.0
            ),
            expected_outputs=[
                OutputSpec(
                    uri="s3://bucket/output/fastqc.html",
                    format="html",
                    required=True
                )
            ],
            confidence_score=0.85,
            reasoning="FastQC is well-tested on EMR with biocontainers support and distributed processing capabilities."
        )
        
        assert spec.tool_name == "fastqc"
        assert spec.infrastructure == "EMR"
        assert spec.confidence_score == 0.85
        assert len(spec.commands) == 1
        assert spec.commands[0].name == "run_fastqc"
    
    def test_empty_commands_fails(self):
        """ExecutionSpec with no commands should fail validation."""
        with pytest.raises(Exception):  # Pydantic ValidationError
            ExecutionToolSpec(
                tool_name="test",
                infrastructure="Local",
                commands=[],  # Invalid: empty
                confidence_score=0.5,
                reasoning="Test spec"
            )
    
    def test_confidence_out_of_range(self):
        """Confidence score must be 0-1."""
        with pytest.raises(Exception):  # Pydantic ValidationError
            ExecutionToolSpec(
                tool_name="test",
                infrastructure="Local",
                commands=[
                    CommandSpec(name="test", command="echo test")
                ],
                confidence_score=1.5,  # Invalid: >1.0
                reasoning="Test spec"
            )
    
    def test_gpu_consistency_warning(self):
        """GPU required without gpu_count should add warning."""
        spec = ExecutionToolSpec(
            tool_name="test",
            infrastructure="EC2",
            commands=[
                CommandSpec(name="test", command="echo test")
            ],
            resource_requirements=ResourceRequirements(
                gpu_required=True,
                gpu_count=None  # Missing
            ),
            confidence_score=0.5,
            reasoning="Test specification for validation purposes with minimal configuration and standard settings."
        )
        
        # Should auto-set gpu_count to 1 and add warning
        assert spec.resource_requirements.gpu_count == 1
        assert any("gpu_count" in w.lower() for w in spec.warnings)
    
    def test_batch_requires_container_warning(self):
        """Batch without container should add warning."""
        spec = ExecutionToolSpec(
            tool_name="test",
            infrastructure="Batch",
            container_spec=None,  # Missing
            commands=[
                CommandSpec(name="test", command="echo test")
            ],
            confidence_score=0.5,
            reasoning="Test specification for validation purposes with minimal configuration and standard settings."
        )
        
        assert any("container" in w.lower() for w in spec.warnings)
    
    def test_retry_policy_validation(self):
        """RetryPolicy should validate ranges."""
        # Valid
        policy = RetryPolicy(
            max_retries=3,
            retry_on=["exit_code", "timeout"],
            backoff_multiplier=2.0,
            initial_delay_seconds=10.0
        )
        assert policy.max_retries == 3
        
        # Invalid: max_retries > 5
        with pytest.raises(Exception):
            RetryPolicy(max_retries=10)  # Out of range
        
        # Invalid: backoff_multiplier < 1.0
        with pytest.raises(Exception):
            RetryPolicy(backoff_multiplier=0.5)  # Too low
    
    def test_resource_requirements_validation(self):
        """ResourceRequirements should validate non-negative values."""
        # Valid
        resources = ResourceRequirements(
            min_cpu_cores=4.0,
            min_memory_gb=8.0,
            min_disk_gb=50.0,
            gpu_required=True,
            gpu_count=2
        )
        assert resources.min_cpu_cores == 4.0
        
        # Invalid: negative CPU
        with pytest.raises(Exception):
            ResourceRequirements(min_cpu_cores=-1.0)
        
        # Invalid: negative memory
        with pytest.raises(Exception):
            ResourceRequirements(min_memory_gb=-4.0)
    
    def test_json_serialization(self):
        """ExecutionSpec should serialize to/from JSON."""
        spec = ExecutionToolSpec(
            tool_name="test",
            infrastructure="Local",
            commands=[
                CommandSpec(name="test", command="echo test")
            ],
            confidence_score=0.5,
            reasoning="Test specification for validation purposes with minimal configuration and standard settings."
        )
        
        # Serialize
        json_str = spec.model_dump_json()
        assert json_str
        
        # Deserialize
        spec_dict = json.loads(json_str)
        spec2 = ExecutionToolSpec(**spec_dict)
        assert spec2.tool_name == spec.tool_name
        assert spec2.infrastructure == spec.infrastructure


class TestImplementationAgent:
    """Test ImplementationAgent logic."""
    
    def test_heuristic_fallback(self):
        """Heuristic fallback should return valid ExecutionSpec."""
        workflow_plan = create_test_workflow_plan("Quality control analysis using FastQC")
        infra_decision = create_test_infra_decision("EMR")
        
        # Call heuristic fallback
        exec_spec = _heuristic_implementation_plan(workflow_plan, infra_decision, "test_req")
        
        # Should return valid ExecutionSpec
        assert isinstance(exec_spec, ExecutionToolSpec)
        assert exec_spec.infrastructure == "EMR"
        # Heuristic fallback should be "reasonable but not high" confidence.
        assert exec_spec.confidence_score == 0.6
        assert len(exec_spec.warnings) > 0  # Should have warnings
        assert "heuristic" in exec_spec.warnings[0].lower()
        assert exec_spec.request_id == "test_req"
    
    def test_heuristic_includes_container_for_batch(self):
        """Heuristic should include container for Batch/Lambda."""
        workflow_plan = create_test_workflow_plan("Test operation for containerization")
        infra_decision = create_test_infra_decision("Batch")
        
        exec_spec = _heuristic_implementation_plan(workflow_plan, infra_decision)
        
        # Should include container spec
        assert exec_spec.container_spec is not None
        assert exec_spec.container_spec.image_type == "docker"
    
    @pytest.mark.asyncio
    async def test_plan_implementation_with_mock_llm(self):
        """Test plan_implementation with mocked LLM."""
        workflow_plan = create_test_workflow_plan("Run FastQC quality control")
        infra_decision = create_test_infra_decision("EMR")
        
        # Mock LLM response with valid JSON
        mock_response_json = {
            "tool_name": "fastqc",
            "infrastructure": "EMR",
            "container_spec": {
                "image": "biocontainers/fastqc:0.11.9",
                "image_type": "docker",
                "pull_policy": "IfNotPresent",
                "env_vars": {},
                "mount_paths": [],
                "working_dir": None
            },
            "commands": [
                {
                    "name": "run_fastqc",
                    "command": "fastqc /data/file.fq",
                    "inputs": ["s3://bucket/file.fq"],
                    "outputs": ["s3://bucket/output/fastqc.html"],
                    "timeout_minutes": 30.0
                }
            ],
            "retry_policy": {
                "max_retries": 2,
                "retry_on": ["exit_code", "timeout"],
                "backoff_multiplier": 2.0,
                "initial_delay_seconds": 10.0
            },
            "resource_requirements": {
                "min_cpu_cores": 4.0,
                "min_memory_gb": 8.0,
                "gpu_required": False
            },
            "expected_outputs": [],
            "confidence_score": 0.85,
            "reasoning": "FastQC is well-supported in biocontainers with extensive testing and proven reliability on EMR clusters."
        }
        
        mock_llm_response = AsyncMock()
        mock_llm_response.content = json.dumps(mock_response_json)
        
        with patch('backend.implementation_agent._get_llm') as mock_get_llm:
            mock_llm = AsyncMock()
            # Make invoke return the mock response directly (not as a coroutine)
            mock_llm.invoke = lambda x: mock_llm_response
            mock_get_llm.return_value = mock_llm
            
            exec_spec = await plan_implementation(workflow_plan, infra_decision, "test_req")
            
            # Should use LLM response, not heuristic fallback
            assert exec_spec.confidence_score >= 0.8  # LLM returned 0.85
            assert exec_spec.infrastructure == "EMR"
            assert exec_spec.request_id == "test_req"


class TestOrchestrator:
    """Test Orchestrator coordination logic."""
    
    def test_generate_request_id(self):
        """Request IDs should be unique."""
        orchestrator = Orchestrator()
        
        req_id_1 = orchestrator.generate_request_id()
        req_id_2 = orchestrator.generate_request_id()
        
        assert req_id_1.startswith("req_")
        assert req_id_2.startswith("req_")
        assert req_id_1 != req_id_2
    
    def test_agent_invocation_to_dict(self):
        """AgentInvocation should serialize to dict."""
        invocation = AgentInvocation(
            agent_name="TestAgent",
            request_id="req_123",
            start_time=1000.0,
            end_time=1005.0,
            duration_ms=5000.0,
            success=True,
            confidence_score=0.8,
            warnings=["Warning 1"]
        )
        
        data = invocation.to_dict()
        assert data["agent_name"] == "TestAgent"
        assert data["request_id"] == "req_123"
        assert data["success"] is True
        assert data["confidence_score"] == 0.8
        assert len(data["warnings"]) == 1
    
    def test_orchestrator_trace_to_json(self):
        """OrchestratorTrace should serialize to JSON."""
        trace = OrchestratorTrace(
            request_id="req_123",
            user_command="Test command",
            start_time=1000.0,
            end_time=1010.0,
            duration_ms=10000.0,
            success=True
        )
        
        trace.invocations.append(AgentInvocation(
            agent_name="Agent1",
            request_id="req_123",
            start_time=1000.0,
            end_time=1005.0,
            duration_ms=5000.0,
            success=True
        ))
        
        json_str = trace.to_json()
        assert json_str
        
        # Parse back
        data = json.loads(json_str)
        assert data["request_id"] == "req_123"
        assert data["success"] is True
        assert len(data["invocations"]) == 1
    
    @pytest.mark.asyncio
    async def test_invoke_infrastructure_agent_success(self):
        """Test infrastructure agent invocation with mock."""
        orchestrator = Orchestrator()
        
        workflow_plan = create_test_workflow_plan("Test workflow for infrastructure decision")
        mock_infra_decision = create_test_infra_decision("Local")
        
        with patch('backend.orchestrator.decide_infrastructure', new_callable=AsyncMock) as mock_agent:
            mock_agent.return_value = mock_infra_decision
            
            infra_decision = await orchestrator.invoke_infrastructure_agent(
                workflow_plan=workflow_plan,
                request_id="test_req"
            )
            
            assert infra_decision.infrastructure == "Local"
            assert infra_decision.confidence_score == 0.8  # From helper function
    
    @pytest.mark.asyncio
    async def test_invoke_infrastructure_agent_failure(self):
        """Test infrastructure agent invocation failure handling."""
        orchestrator = Orchestrator()
        
        workflow_plan = create_test_workflow_plan("Test workflow for testing")
        
        trace = OrchestratorTrace(
            request_id="test_req",
            user_command="Test",
            start_time=1000.0
        )
        
        # Mock the agent to raise exception
        with patch('backend.orchestrator.decide_infrastructure', new_callable=AsyncMock) as mock_agent:
            mock_agent.side_effect = Exception("Agent failed")
            
            with pytest.raises(Exception, match="Agent failed"):
                await orchestrator.invoke_infrastructure_agent(
                    workflow_plan=workflow_plan,
                    request_id="test_req",
                    trace=trace
                )
            
            # Should have recorded failure in trace
            assert len(trace.invocations) == 1
            assert trace.invocations[0].success is False
            assert "Agent failed" in trace.invocations[0].error
    
    @pytest.mark.asyncio
    async def test_run_pipeline_success(self):
        """Test full pipeline execution."""
        orchestrator = Orchestrator()
        
        workflow_plan = create_test_workflow_plan("Run FastQC pipeline end-to-end")
        mock_infra_decision = create_test_infra_decision("EMR")
        
        mock_exec_spec = ExecutionToolSpec(
            tool_name="fastqc",
            infrastructure="EMR",
            commands=[CommandSpec(name="test", command="fastqc file.fq")],
            confidence_score=0.85,
            reasoning="FastQC is well-tested on EMR with biocontainers support and distributed processing."
        )
        
        with patch('backend.orchestrator.decide_infrastructure', new_callable=AsyncMock) as mock_infra:
            with patch('backend.orchestrator.plan_implementation', new_callable=AsyncMock) as mock_impl:
                mock_infra.return_value = mock_infra_decision
                mock_impl.return_value = mock_exec_spec
                
                trace = await orchestrator.run_pipeline(
                    command="Run FastQC",
                    workflow_plan=workflow_plan
                )
                
                # Should succeed
                assert trace.success is True
                assert len(trace.invocations) == 2  # 2 agents
                assert trace.infra_decision == mock_infra_decision
                assert trace.execution_spec == mock_exec_spec
                assert trace.request_id.startswith("req_")
    
    @pytest.mark.asyncio
    async def test_run_pipeline_failure(self):
        """Test pipeline failure handling."""
        orchestrator = Orchestrator()
        workflow_plan = create_test_workflow_plan("Test pipeline failure scenarios")
        
        # Mock infrastructure agent to fail
        with patch('backend.orchestrator.decide_infrastructure', new_callable=AsyncMock) as mock_infra:
            mock_infra.side_effect = Exception("Infrastructure agent failed")
            
            with pytest.raises(Exception, match="Infrastructure agent failed"):
                await orchestrator.run_pipeline(
                    command="Test",
                    workflow_plan=workflow_plan
                )
            
            # Trace should be stored even for failures
            traces = orchestrator.get_all_traces()
            assert len(traces) == 1
            assert traces[0].success is False


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
