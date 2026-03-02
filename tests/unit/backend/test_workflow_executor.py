"""
Unit tests for workflow_executor module.

Tests the multi-step workflow execution infrastructure including:
- Multi-step workflow detection
- Tool availability checking
- Step execution
- Data passing between steps
"""

import pytest
import asyncio
from unittest.mock import Mock, patch, AsyncMock
from pathlib import Path

from backend.workflow_executor import (
    WorkflowExecutor,
    WorkflowExecutionResult,
    get_workflow_executor
)
from backend.contracts.workflow_plan import WorkflowPlan, OperationSpec, DataInput


class TestWorkflowExecutionResult:
    """Test WorkflowExecutionResult class."""
    
    def test_initialization(self):
        """Test result initialization with default values."""
        result = WorkflowExecutionResult()
        
        assert result.status == "pending"
        assert result.steps_completed == 0
        assert result.total_steps == 0
        assert result.step_results == []
        assert result.errors == []
        assert result.outputs == {}
        assert result.start_time is not None
        assert result.end_time is None
    
    def test_to_dict(self):
        """Test conversion to dictionary."""
        result = WorkflowExecutionResult()
        result.workflow_id = "test_wf_123"
        result.status = "completed"
        result.steps_completed = 3
        result.total_steps = 3
        result.end_time = result.start_time + 10.5
        
        result_dict = result.to_dict()
        
        assert result_dict["workflow_id"] == "test_wf_123"
        assert result_dict["status"] == "completed"
        assert result_dict["steps_completed"] == 3
        assert result_dict["total_steps"] == 3
        assert result_dict["duration_seconds"] == 10.5
        assert "start_time" in result_dict
        assert "end_time" in result_dict


class TestWorkflowExecutor:
    """Test WorkflowExecutor class."""
    
    @pytest.fixture
    def executor(self):
        """Create a workflow executor instance."""
        return WorkflowExecutor()
    
    @pytest.fixture
    def sample_workflow_plan(self):
        """Create a sample workflow plan for testing."""
        return WorkflowPlan(
            description="Test RNA-seq workflow",
            data_inputs=[
                DataInput(
                    uri="s3://test-bucket/R1.fq",
                    size_bytes=1000000,
                    description="Read 1"
                ),
                DataInput(
                    uri="s3://test-bucket/R2.fq",
                    size_bytes=1000000,
                    description="Read 2"
                )
            ],
            operations=[
                OperationSpec(
                    operation_name="quality_control",
                    tool_name="fastqc",
                    parameters={},
                    expected_output_size_mb=10.0,
                    parallelizable=True
                ),
                OperationSpec(
                    operation_name="alignment",
                    tool_name="STAR",
                    parameters={"reference_genome": "hg38"},
                    expected_output_size_mb=2000.0,
                    parallelizable=False
                ),
                OperationSpec(
                    operation_name="quantification",
                    tool_name="featureCounts",
                    parameters={"feature_type": "exon"},
                    expected_output_size_mb=5.0,
                    parallelizable=False
                )
            ],
            expected_compute_intensity="High",
            session_id="test_session"
        )
    
    # =========================================================================
    # Multi-Step Detection Tests
    # =========================================================================
    
    def test_is_multi_step_command_with_workflow_keyword(self, executor):
        """Test detection with 'workflow' keyword."""
        command = "Run RNA-seq workflow on my samples"
        assert executor._is_multi_step_command(command) is True
    
    def test_is_multi_step_command_with_pipeline_keyword(self, executor):
        """Test detection with 'pipeline' keyword."""
        command = "Execute alignment pipeline"
        assert executor._is_multi_step_command(command) is True
    
    def test_is_multi_step_command_with_step_indicators(self, executor):
        """Test detection with step indicators."""
        command = "Step 1: run FastQC, Step 2: align with STAR"
        assert executor._is_multi_step_command(command) is True
    
    def test_is_multi_step_command_with_arrow_separator(self, executor):
        """Test detection with arrow separator."""
        command = "Run FastQC → STAR → featureCounts"
        assert executor._is_multi_step_command(command) is True
    
    def test_is_multi_step_command_with_multiple_tools(self, executor):
        """Test detection with multiple tools mentioned."""
        command = "Run FastQC and STAR on my data"
        assert executor._is_multi_step_command(command) is True
    
    def test_is_multi_step_command_single_tool(self, executor):
        """Test that single tool commands are not detected as multi-step."""
        command = "Run FastQC on my samples"
        assert executor._is_multi_step_command(command) is False
    
    def test_is_multi_step_command_no_indicators(self, executor):
        """Test that commands without indicators are not detected."""
        command = "Analyze my genomic data"
        assert executor._is_multi_step_command(command) is False
    
    # =========================================================================
    # Tool Availability Tests
    # =========================================================================
    
    @pytest.mark.asyncio
    async def test_check_tool_availability_existing_tool(self, executor):
        """Test checking availability of existing tool (FastQC)."""
        is_available = await executor._check_tool_availability("fastqc")
        assert is_available is True
    
    @pytest.mark.asyncio
    async def test_check_tool_availability_missing_tool(self, executor):
        """Test checking availability of missing tool."""
        is_available = await executor._check_tool_availability("nonexistent_tool_xyz")
        assert is_available is False
    
    @pytest.mark.asyncio
    async def test_check_tool_availability_caching(self, executor):
        """Test that tool availability is cached."""
        # First check
        is_available1 = await executor._check_tool_availability("fastqc")
        
        # Second check should use cache
        is_available2 = await executor._check_tool_availability("fastqc")
        
        assert is_available1 == is_available2
        assert "fastqc" in executor.tool_registry
    
    @pytest.mark.asyncio
    async def test_check_tool_availability_case_insensitive(self, executor):
        """Test that tool checking handles different naming conventions."""
        # FastQC is actually fastqc_quality_analysis in agent_tools
        is_available = await executor._check_tool_availability("fastqc")
        assert is_available is True
    
    # =========================================================================
    # Workflow Execution Tests
    # =========================================================================
    
    @pytest.mark.asyncio
    async def test_execute_workflow_initialization(self, executor, sample_workflow_plan):
        """Test workflow execution initialization."""
        with patch.object(executor, '_execute_step', new_callable=AsyncMock) as mock_execute:
            mock_execute.return_value = {
                "status": "completed",
                "tool_name": "fastqc",
                "outputs": {}
            }
            
            result = await executor.execute_workflow(
                workflow_plan=sample_workflow_plan,
                session_id="test_session"
            )
            
            assert result.workflow_id is not None
            assert result.workflow_id.startswith("wf_test_session_")
            assert result.total_steps == 3
            assert result.start_time is not None
    
    @pytest.mark.asyncio
    async def test_execute_workflow_success(self, executor, sample_workflow_plan):
        """Test successful workflow execution with all steps completing."""
        # Create a simple plan with just FastQC (which exists)
        simple_plan = WorkflowPlan(
            description="Simple test workflow",
            data_inputs=sample_workflow_plan.data_inputs,
            operations=[sample_workflow_plan.operations[0]],  # Just FastQC
            expected_compute_intensity="Low",
            session_id="test_session"
        )
        
        with patch.object(executor, '_execute_step', new_callable=AsyncMock) as mock_execute:
            mock_execute.return_value = {
                "status": "completed",
                "tool_name": "fastqc",
                "execution_mode": "local",
                "outputs": {"output_path": "s3://test-bucket/results/"}
            }
            
            result = await executor.execute_workflow(
                workflow_plan=simple_plan,
                session_id="test_session"
            )
            
            assert result.status == "completed"
            assert result.steps_completed == 1
            assert result.total_steps == 1
            assert len(result.step_results) == 1
            assert len(result.errors) == 0
            assert result.end_time is not None
    
    @pytest.mark.asyncio
    async def test_execute_workflow_step_failure(self, executor, sample_workflow_plan):
        """Test workflow execution when a step fails."""
        with patch.object(executor, '_execute_step', new_callable=AsyncMock) as mock_execute:
            # First step succeeds, second step fails
            mock_execute.side_effect = [
                {"status": "completed", "outputs": {}},
                {"status": "failed", "error": "Tool not found"}
            ]
            
            result = await executor.execute_workflow(
                workflow_plan=sample_workflow_plan,
                session_id="test_session"
            )
            
            assert result.status == "failed"
            assert result.steps_completed == 2  # Got to second step
            assert len(result.errors) > 0
            assert "Tool not found" in result.errors[0]
    
    @pytest.mark.asyncio
    async def test_execute_workflow_exception_handling(self, executor, sample_workflow_plan):
        """Test that workflow handles exceptions gracefully."""
        with patch.object(executor, '_execute_step', new_callable=AsyncMock) as mock_execute:
            mock_execute.side_effect = Exception("Unexpected error")
            
            result = await executor.execute_workflow(
                workflow_plan=sample_workflow_plan,
                session_id="test_session"
            )
            
            assert result.status == "failed"
            assert len(result.errors) > 0
            assert "Unexpected error" in result.errors[0]
    
    # =========================================================================
    # Step Execution Tests
    # =========================================================================
    
    @pytest.mark.asyncio
    async def test_execute_step_with_existing_tool(self, executor):
        """Test executing a step with an existing tool."""
        operation = OperationSpec(
            operation_name="quality_control",
            tool_name="fastqc",
            parameters={},
            expected_output_size_mb=10.0
        )
        
        with patch.object(executor, '_check_tool_availability', return_value=True):
            with patch.object(executor, '_execute_tool', new_callable=AsyncMock) as mock_exec:
                mock_exec.return_value = {
                    "status": "completed",
                    "execution_mode": "local",
                    "outputs": {"output_path": "s3://test/results/"}
                }
                
                result = await executor._execute_step(
                    operation=operation,
                    step_index=1,
                    previous_outputs={},
                    session_id="test",
                    session_context={}
                )
                
                assert result["status"] == "completed"
                assert result["tool_generated"] is False
                assert "outputs" in result
    
    @pytest.mark.asyncio
    async def test_execute_step_with_missing_tool(self, executor):
        """Test executing a step when tool needs generation."""
        operation = OperationSpec(
            operation_name="alignment",
            tool_name="STAR",
            parameters={},
            expected_output_size_mb=2000.0
        )
        
        with patch.object(executor, '_check_tool_availability', return_value=False):
            # Tool generation is done via _generate_and_execute_tool (legacy _generate_tool is deprecated).
            with patch.object(executor, '_generate_and_execute_tool', new_callable=AsyncMock) as mock_gen:
                mock_gen.return_value = {
                    "status": "failed",
                    "execution_mode": "generated_tool",
                    "message": "Tool generation failed: Tool generation not fully wired",
                    "outputs": {},
                }
                
                result = await executor._execute_step(
                    operation=operation,
                    step_index=1,
                    previous_outputs={},
                    session_id="test",
                    session_context={}
                )
                
                assert result["status"] == "failed"
                assert "Tool generation failed" in result["error"]
    
    # =========================================================================
    # Data Passing Tests
    # =========================================================================
    
    def test_prepare_execution_params_basic(self, executor):
        """Test basic parameter preparation."""
        operation = OperationSpec(
            operation_name="test_op",
            tool_name="test_tool",
            parameters={"param1": "value1"},
            expected_output_size_mb=10.0
        )
        
        params = executor._prepare_execution_params(
            operation=operation,
            previous_outputs={},
            session_context={}
        )
        
        assert params["param1"] == "value1"
    
    def test_prepare_execution_params_with_previous_outputs(self, executor):
        """Test parameter preparation with outputs from previous step."""
        operation = OperationSpec(
            operation_name="quantification",
            tool_name="featureCounts",
            parameters={},
            expected_output_size_mb=5.0
        )
        
        previous_outputs = {
            "bam_file": "s3://test-bucket/aligned.bam"
        }
        
        params = executor._prepare_execution_params(
            operation=operation,
            previous_outputs=previous_outputs,
            session_context={}
        )
        
        assert "input_bam" in params
        assert params["input_bam"] == "s3://test-bucket/aligned.bam"
    
    def test_prepare_execution_params_with_fastq_files(self, executor):
        """Test parameter preparation with FASTQ files from previous step."""
        operation = OperationSpec(
            operation_name="alignment",
            tool_name="STAR",
            parameters={},
            expected_output_size_mb=2000.0
        )
        
        previous_outputs = {
            "fastq_files": ["s3://test/R1.fq", "s3://test/R2.fq"]
        }
        
        params = executor._prepare_execution_params(
            operation=operation,
            previous_outputs=previous_outputs,
            session_context={}
        )
        
        assert "input_files" in params
        assert params["input_files"] == ["s3://test/R1.fq", "s3://test/R2.fq"]
    
    def test_extract_outputs_basic(self, executor):
        """Test extracting outputs from tool result."""
        tool_result = {
            "output": "s3://test-bucket/results/",
            "status": "completed"
        }
        
        outputs = executor._extract_outputs(tool_result)
        
        assert "output_path" in outputs
        assert outputs["output_path"] == "s3://test-bucket/results/"
    
    def test_extract_outputs_with_bam(self, executor):
        """Test extracting BAM file output."""
        tool_result = {
            "output_bam": "s3://test-bucket/aligned.bam",
            "status": "completed"
        }
        
        outputs = executor._extract_outputs(tool_result)
        
        assert "bam_file" in outputs
        assert outputs["bam_file"] == "s3://test-bucket/aligned.bam"
    
    # =========================================================================
    # Singleton Tests
    # =========================================================================
    
    def test_get_workflow_executor_singleton(self):
        """Test that get_workflow_executor returns singleton instance."""
        executor1 = get_workflow_executor()
        executor2 = get_workflow_executor()
        
        assert executor1 is executor2


class TestWorkflowExecutorIntegration:
    """Integration tests with actual backend modules."""
    
    @pytest.mark.asyncio
    async def test_fastqc_tool_availability(self):
        """Test that FastQC is detected as available."""
        executor = WorkflowExecutor()
        is_available = await executor._check_tool_availability("fastqc")
        
        # FastQC should be available in agent_tools
        assert is_available is True
    
    @pytest.mark.asyncio
    async def test_star_tool_unavailability(self):
        """Test that STAR is correctly detected as unavailable."""
        executor = WorkflowExecutor()
        is_available = await executor._check_tool_availability("STAR")
        
        # STAR is not implemented yet
        assert is_available is False
    
    @pytest.mark.asyncio
    async def test_workflow_plan_creation(self):
        """Test creating a workflow plan for testing."""
        workflow_plan = WorkflowPlan(
            description="Test workflow",
            data_inputs=[
                DataInput(uri="s3://test/file.fq", size_bytes=1000)
            ],
            operations=[
                OperationSpec(
                    operation_name="qc",
                    tool_name="fastqc",
                    parameters={},
                    expected_output_size_mb=10.0
                )
            ],
            expected_compute_intensity="Low",
            session_id="test"
        )
        
        assert workflow_plan is not None
        assert len(workflow_plan.operations) == 1
        assert workflow_plan.operations[0].tool_name == "fastqc"


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v", "-s"])
