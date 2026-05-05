"""
End-to-end tests for multi-step workflow execution.

Tests the complete workflow integration from command → detection → planning → execution.
"""

import pytest
import asyncio
from unittest.mock import patch, AsyncMock, Mock

from backend.agent import handle_command, CommandProcessor
from backend.workflow_planner_agent import plan_workflow
from backend.workflow_executor import get_workflow_executor


class TestMultiStepWorkflowE2E:
    """End-to-end tests for multi-step workflow execution."""
    
    @pytest.mark.asyncio
    async def test_multi_step_detection_in_command_processor(self):
        """Test that multi-step workflows are detected by command processor."""
        command = "Run complete RNA-seq workflow: FastQC → STAR → featureCounts"
        
        # The command processor should detect this as multi-step
        executor = get_workflow_executor()
        is_multi_step = executor._is_multi_step_command(command)
        
        assert is_multi_step is True
    
    @pytest.mark.asyncio
    async def test_workflow_planning_from_command(self):
        """Test that workflow planning works from natural language command."""
        command = """Run RNA-seq analysis:
        Step 1: Quality control with FastQC
        Step 2: Alignment with STAR
        Step 3: Quantification with featureCounts
        
        Input: s3://test-bucket/sample_R1.fq
        """
        
        session_context = {"session_id": "test_e2e"}
        
        # Plan workflow
        plan_result = await plan_workflow(command, session_context)
        
        assert plan_result["status"] == "success"
        workflow_plan = plan_result["workflow_plan"]
        
        # Should create multi-step plan
        assert len(workflow_plan.operations) >= 2
        
        # Check that operations are in the plan
        operation_names = [op.operation_name for op in workflow_plan.operations]
        assert any("quality" in name.lower() or "qc" in name.lower() for name in operation_names)
    
    @pytest.mark.asyncio
    async def test_handle_command_routes_to_workflow_handler(self):
        """Test that handle_command routes multi-step commands to workflow handler."""
        command = "Execute complete RNA-seq workflow: FastQC then STAR alignment"

        # In HELIX_MOCK_MODE the live LLM-based intent classifier is disabled
        # (raises by design — see backend/intent_classifier.py). This test only
        # verifies the multi-step routing branch, so stub classify_intent to
        # return a deterministic "execute" decision.
        from backend.intent_classifier import IntentDecision

        with patch('backend.agent.CommandProcessor._handle_multi_step_workflow', new_callable=AsyncMock) as mock_handler, \
             patch('backend.intent_classifier.classify_intent',
                   return_value=IntentDecision(intent="execute", reason="test_mock")), \
             patch('backend.agent._get_agent', return_value=Mock()):
            mock_handler.return_value = {
                "status": "workflow_executed",
                "success": True,
                "workflow_id": "test_wf_123"
            }

            result = await handle_command(command, session_id="test_e2e")

            # Should have called the multi-step workflow handler
            assert mock_handler.called or "workflow" in str(result)
    
    @pytest.mark.asyncio
    async def test_single_tool_not_routed_to_workflow(self):
        """Test that single tool commands are not routed to workflow handler."""
        command = "Run FastQC on my samples"
        
        executor = get_workflow_executor()
        is_multi_step = executor._is_multi_step_command(command)
        
        # Should NOT be detected as multi-step
        assert is_multi_step is False
    
    @pytest.mark.asyncio
    async def test_workflow_execution_with_fastqc_only(self):
        """Test executing a workflow with only FastQC (which exists)."""
        from backend.contracts.workflow_plan import WorkflowPlan, OperationSpec, DataInput
        
        # Create a simple workflow with just FastQC
        workflow_plan = WorkflowPlan(
            description="Simple FastQC workflow",
            data_inputs=[
                DataInput(
                    uri="s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq",
                    size_bytes=617923,
                    description="Test R1"
                )
            ],
            operations=[
                OperationSpec(
                    operation_name="quality_control",
                    tool_name="fastqc",
                    parameters={
                        "input_r1": "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq",
                        "output": "s3://noricum-ngs-data/results/test-workflow/fastqc/"
                    },
                    expected_output_size_mb=10.0
                )
            ],
            expected_compute_intensity="Low",
            session_id="test_e2e_fastqc"
        )
        
        # Execute workflow
        executor = get_workflow_executor()
        
        # Mock the tool execution to avoid actual FastQC run
        with patch.object(executor, '_execute_tool', new_callable=AsyncMock) as mock_exec:
            mock_exec.return_value = {
                "status": "completed",
                "execution_mode": "local",
                "message": "FastQC completed",
                "outputs": {
                    "output_path": "s3://noricum-ngs-data/results/test-workflow/fastqc/"
                }
            }
            
            result = await executor.execute_workflow(
                workflow_plan=workflow_plan,
                session_id="test_e2e_fastqc",
                session_context={}
            )
            
            assert result.status == "completed"
            assert result.steps_completed == 1
            assert result.total_steps == 1
            assert len(result.errors) == 0


class TestWorkflowDetectionPatterns:
    """Test various command patterns for workflow detection."""
    
    def test_detect_arrow_separator(self):
        """Test detection with → separator."""
        executor = get_workflow_executor()
        
        commands = [
            "Run FastQC → STAR → featureCounts",
            "Execute QC → alignment → quantification",
            "Process: trim → align → count"
        ]
        
        for cmd in commands:
            assert executor._is_multi_step_command(cmd) is True
    
    def test_detect_then_keyword(self):
        """Test detection with 'then' keyword."""
        executor = get_workflow_executor()
        
        commands = [
            "Run FastQC then align with STAR",
            "First trim, then align",
            "QC first, then quantify"
        ]
        
        for cmd in commands:
            assert executor._is_multi_step_command(cmd) is True
    
    def test_detect_step_numbering(self):
        """Test detection with step numbering."""
        executor = get_workflow_executor()
        
        commands = [
            "Step 1: FastQC, Step 2: STAR",
            "1. Quality control 2. Alignment",
            "First do QC, second align"
        ]
        
        for cmd in commands:
            assert executor._is_multi_step_command(cmd) is True
    
    def test_detect_workflow_keywords(self):
        """Test detection with workflow/pipeline keywords."""
        executor = get_workflow_executor()
        
        commands = [
            "Run RNA-seq workflow",
            "Execute alignment pipeline",
            "Complete analysis workflow",
            "Full processing pipeline"
        ]
        
        for cmd in commands:
            assert executor._is_multi_step_command(cmd) is True
    
    def test_detect_multiple_tools(self):
        """Test detection when multiple tools are mentioned."""
        executor = get_workflow_executor()
        
        commands = [
            "Run FastQC and STAR",
            "Use trimmomatic and bwa for processing",
            "Apply FastQC, STAR, and featureCounts"
        ]
        
        for cmd in commands:
            assert executor._is_multi_step_command(cmd) is True
    
    def test_not_detect_single_tool(self):
        """Test that single tool commands are not detected."""
        executor = get_workflow_executor()
        
        commands = [
            "Run FastQC",
            "Execute alignment",
            "Perform quality control",
            "Analyze with STAR"
        ]
        
        for cmd in commands:
            assert executor._is_multi_step_command(cmd) is False


class TestWorkflowPlannerIntegration:
    """Test integration with workflow planner."""
    
    @pytest.mark.asyncio
    async def test_rna_seq_playbook_selection(self):
        """Test that RNA-seq playbook is selected for RNA-seq commands."""
        command = """Analyze RNA-seq data:
        - Input: s3://bucket/sample_R1.fastq
        - Organism: human
        - Reference: hg38
        """
        
        session_context = {"session_id": "test_rna_seq"}
        plan_result = await plan_workflow(command, session_context)
        
        assert plan_result["status"] == "success"
        assert plan_result["workflow_type"] == "rna_seq_bulk"
    
    @pytest.mark.asyncio
    async def test_simple_operation_playbook_fallback(self):
        """Test that simple operations fall back to simple playbook."""
        command = "Run quality control on my data"
        
        session_context = {"session_id": "test_simple"}
        plan_result = await plan_workflow(command, session_context)
        
        assert plan_result["status"] == "success"
        # Should fall back to simple operation playbook
        workflow_plan = plan_result["workflow_plan"]
        assert len(workflow_plan.operations) >= 1
    
    @pytest.mark.asyncio
    async def test_workflow_plan_extraction_of_s3_paths(self):
        """Test that S3 paths are correctly extracted from commands."""
        command = """Run analysis on:
        - s3://bucket/file1.fastq
        - s3://bucket/file2.fastq
        """
        
        session_context = {"session_id": "test_s3"}
        plan_result = await plan_workflow(command, session_context)
        
        if plan_result["status"] == "success":
            workflow_plan = plan_result["workflow_plan"]
            # Should have extracted input files
            assert len(workflow_plan.data_inputs) >= 1


class TestWorkflowExecutorEdgeCases:
    """Test edge cases and error handling."""
    
    @pytest.mark.asyncio
    async def test_workflow_with_no_operations(self):
        """Test handling of empty workflow."""
        from backend.contracts.workflow_plan import WorkflowPlan
        
        workflow_plan = WorkflowPlan(
            description="Empty workflow",
            data_inputs=[],
            operations=[],
            expected_compute_intensity="Low",
            session_id="test_empty"
        )
        
        executor = get_workflow_executor()
        result = await executor.execute_workflow(
            workflow_plan=workflow_plan,
            session_id="test_empty",
            session_context={}
        )
        
        # Should complete immediately with no steps
        assert result.status == "completed"
        assert result.steps_completed == 0
        assert result.total_steps == 0
    
    @pytest.mark.asyncio
    async def test_workflow_with_step_failure_stops_execution(self):
        """Test that workflow stops after step failure."""
        from backend.contracts.workflow_plan import WorkflowPlan, OperationSpec
        
        workflow_plan = WorkflowPlan(
            description="Failing workflow",
            data_inputs=[],
            operations=[
                OperationSpec(
                    operation_name="step1",
                    tool_name="nonexistent_tool",
                    parameters={},
                    expected_output_size_mb=1.0
                ),
                OperationSpec(
                    operation_name="step2",
                    tool_name="fastqc",
                    parameters={},
                    expected_output_size_mb=1.0
                )
            ],
            expected_compute_intensity="Low",
            session_id="test_failure"
        )
        
        executor = get_workflow_executor()
        result = await executor.execute_workflow(
            workflow_plan=workflow_plan,
            session_id="test_failure",
            session_context={}
        )
        
        # Should fail on first step
        assert result.status == "failed"
        assert result.steps_completed == 1  # Attempted first step
        assert len(result.errors) > 0
        assert len(result.step_results) == 1  # Only first step result recorded
    
    @pytest.mark.asyncio
    async def test_tool_availability_caching_across_steps(self):
        """Test that tool availability checks are cached."""
        executor = get_workflow_executor()
        
        # Check same tool multiple times
        is_available1 = await executor._check_tool_availability("fastqc")
        is_available2 = await executor._check_tool_availability("fastqc")
        is_available3 = await executor._check_tool_availability("fastqc")
        
        # All should return same result
        assert is_available1 == is_available2 == is_available3
        
        # Should be cached
        assert "fastqc" in executor.tool_registry


class TestDataPassingBetweenSteps:
    """Test data passing between workflow steps."""
    
    def test_fastq_to_bam_passing(self):
        """Test passing FASTQ files from QC to alignment."""
        executor = get_workflow_executor()
        
        from backend.contracts.workflow_plan import OperationSpec
        
        operation = OperationSpec(
            operation_name="alignment",
            tool_name="STAR",
            parameters={"reference": "hg38"},
            expected_output_size_mb=2000.0
        )
        
        previous_outputs = {
            "fastq_files": ["s3://bucket/R1.fq", "s3://bucket/R2.fq"]
        }
        
        params = executor._prepare_execution_params(
            operation=operation,
            previous_outputs=previous_outputs,
            session_context={}
        )
        
        assert "input_files" in params
        assert params["input_files"] == ["s3://bucket/R1.fq", "s3://bucket/R2.fq"]
        assert params["reference"] == "hg38"
    
    def test_bam_to_counts_passing(self):
        """Test passing BAM file from alignment to quantification."""
        executor = get_workflow_executor()
        
        from backend.contracts.workflow_plan import OperationSpec
        
        operation = OperationSpec(
            operation_name="quantification",
            tool_name="featureCounts",
            parameters={},
            expected_output_size_mb=5.0
        )
        
        previous_outputs = {
            "bam_file": "s3://bucket/aligned.bam"
        }
        
        params = executor._prepare_execution_params(
            operation=operation,
            previous_outputs=previous_outputs,
            session_context={}
        )
        
        assert "input_bam" in params
        assert params["input_bam"] == "s3://bucket/aligned.bam"
    
    def test_output_extraction(self):
        """Test extracting outputs from tool results."""
        executor = get_workflow_executor()
        
        tool_result = {
            "output": "s3://bucket/results/",
            "output_bam": "s3://bucket/output.bam",
            "counts_file": "s3://bucket/counts.txt",
            "results_path": "s3://bucket/all_results/"
        }
        
        outputs = executor._extract_outputs(tool_result)
        
        assert outputs["output_path"] == "s3://bucket/results/"
        assert outputs["bam_file"] == "s3://bucket/output.bam"
        assert outputs["counts_file"] == "s3://bucket/counts.txt"
        assert outputs["results_path"] == "s3://bucket/all_results/"


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v", "-s"])
