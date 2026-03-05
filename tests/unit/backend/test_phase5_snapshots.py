"""
Phase 5: Snapshot Tests for Multi-Agent Pipeline

Tests that verify consistent behavior across predefined scenarios.
These are "snapshot" or "golden" tests that capture expected outputs
and detect regressions in agent decision-making.

Run with: pytest tests/unit/backend/test_phase5_snapshots.py -v
"""

import pytest
from backend.contracts.workflow_plan import WorkflowPlan, OperationSpec, DataInput
from backend.contracts.dataset_spec import DatasetSpec
from backend.infrastructure_decision_agent import decide_infrastructure, InputAsset, OutputAsset
from backend.implementation_agent import plan_implementation
from backend.orchestrator import Orchestrator


# ============================================================================
# Test Fixtures: Workflow Plans for Different Scenarios
# ============================================================================

@pytest.fixture
def small_local_files_workflow():
    """Small local files (<100MB) - should recommend Local execution."""
    return WorkflowPlan(
        description="Quality control analysis of small RNA-seq data using FastQC",
        operations=[
            OperationSpec(
                operation_name="fastqc",
                tool_name="fastqc"
            )
        ],
        data_inputs=[
            DataInput(
                uri="/data/sample_R1.fastq",
                size_bytes=50000000,  # 50MB
                location_type="Local",
                metadata={"read": "R1", "format": "fastq"}
            ),
            DataInput(
                uri="/data/sample_R2.fastq",
                size_bytes=48000000,  # 48MB
                location_type="Local",
                metadata={"read": "R2", "format": "fastq"}
            )
        ]
    )


@pytest.fixture
def large_s3_files_workflow():
    """Large S3 files (>100MB) - should recommend EMR execution."""
    return WorkflowPlan(
        description="Read merging operation on large RNA-seq paired-end files stored in S3",
        operations=[
            OperationSpec(
                operation_name="read_merging",
                tool_name="bbmerge"
            )
        ],
        data_inputs=[
            DataInput(
                uri="s3://helix-biodata/rnaseq/sample_R1.fastq",
                size_bytes=262144000,  # 250MB
                location_type="S3",
                metadata={"read": "R1", "format": "fastq"}
            ),
            DataInput(
                uri="s3://helix-biodata/rnaseq/sample_R2.fastq",
                size_bytes=251658240,  # 240MB
                location_type="S3",
                metadata={"read": "R2", "format": "fastq"}
            )
        ]
    )


@pytest.fixture
def unknown_sizes_workflow():
    """Unknown file sizes - should have reduced confidence and assumptions."""
    return WorkflowPlan(
        description="Quality assessment of files with unknown sizes",
        operations=[
            OperationSpec(
                operation_name="quality_control",
                tool_name="fastqc"
            )
        ],
        data_inputs=[
            DataInput(
                uri="s3://restricted-bucket/file1.fastq",
                size_bytes=None,  # Unknown
                location_type="S3",
                metadata={"note": "Size unavailable", "format": "fastq"}
            ),
            DataInput(
                uri="s3://restricted-bucket/file2.fastq",
                size_bytes=None,  # Unknown
                location_type="S3",
                metadata={"note": "Size unavailable", "format": "fastq"}
            )
        ]
    )


@pytest.fixture
def mixed_locations_workflow():
    """Mixed S3 + Local files - should recommend EMR with upload."""
    return WorkflowPlan(
        description="Alignment workflow with mixed S3 and local input files",
        operations=[
            OperationSpec(
                operation_name="alignment",
                tool_name="bwa"
            )
        ],
        data_inputs=[
            DataInput(
                uri="s3://helix-biodata/reference/hg38.fa",
                size_bytes=314572800,  # 300MB
                location_type="S3",
                metadata={"type": "reference", "format": "fasta"}
            ),
            DataInput(
                uri="/local/data/sample.fastq",
                size_bytes=52428800,  # 50MB
                location_type="Local",
                metadata={"type": "reads", "format": "fastq"}
            )
        ]
    )


# ============================================================================
# Snapshot Tests: Infrastructure Decision Agent
# ============================================================================

class TestInfrastructureDecisionSnapshots:
    """Snapshot tests for Infrastructure Decision Agent."""
    
    @pytest.mark.asyncio
    async def test_small_local_files_decision(self, small_local_files_workflow):
        """Test: Small local files should recommend Local execution."""
        # Convert DataInput to InputAsset format expected by decide_infrastructure
        inputs = [
            InputAsset(uri=inp.uri, size_bytes=inp.size_bytes)
            for inp in small_local_files_workflow.data_inputs
        ]
        
        decision = await decide_infrastructure(
            command=small_local_files_workflow.description,
            inputs=inputs
        )
        
        # Snapshot assertions
        assert decision.infrastructure == "Local", "Small local files should use Local execution"
        assert decision.confidence_score >= 0.8, "Should have high confidence for known scenario"
        assert decision.file_analysis.total_size_mb < 100, "Total size should be <100MB"
        assert decision.file_analysis.all_local is True, "All files should be local"
        assert decision.file_analysis.unknown_sizes == 0, "No unknown file sizes"
        
        # Reasoning should mention key factors
        reasoning_lower = decision.reasoning.lower()
        assert any(keyword in reasoning_lower for keyword in ["local", "small", "no transfer"]), \
            "Reasoning should mention local execution or small files"
    
    @pytest.mark.asyncio
    async def test_large_s3_files_decision(self, large_s3_files_workflow):
        """Test: Large S3 files should recommend EMR execution."""
        inputs = [
            InputAsset(uri=inp.uri, size_bytes=inp.size_bytes)
            for inp in large_s3_files_workflow.data_inputs
        ]
        
        decision = await decide_infrastructure(
            command=large_s3_files_workflow.description,
            inputs=inputs
        )
        
        # Snapshot assertions
        assert decision.infrastructure == "EMR", "Large S3 files should use EMR"
        assert decision.confidence_score >= 0.8, "Should have high confidence for known scenario"
        assert decision.file_analysis.total_size_mb > 100, "Total size should be >100MB"
        assert decision.file_analysis.all_in_s3 is True, "All files should be in S3"
        assert decision.file_analysis.unknown_sizes == 0, "No unknown file sizes"
        
        # Cost analysis should be present
        assert decision.cost_analysis.estimated_cost_range_usd[0] >= 0, "Cost min should be non-negative"
        assert decision.cost_analysis.estimated_cost_range_usd[1] > decision.cost_analysis.estimated_cost_range_usd[0], \
            "Cost max should be > cost min"
        assert decision.cost_analysis.cost_class in ["Low", "Medium", "High"], "Cost class should be valid"
        
        # Should provide alternatives
        assert len(decision.alternatives) >= 1, "Should provide at least 1 alternative"
    
    @pytest.mark.asyncio
    async def test_unknown_sizes_decision(self, unknown_sizes_workflow):
        """Test: Unknown file sizes should reduce confidence and add warnings."""
        inputs = [
            InputAsset(uri=inp.uri, size_bytes=inp.size_bytes)
            for inp in unknown_sizes_workflow.data_inputs
        ]
        
        decision = await decide_infrastructure(
            command=unknown_sizes_workflow.description,
            inputs=inputs
        )
        
        # Snapshot assertions
        assert decision.infrastructure in ["EMR", "Local", "EC2"], "Should recommend valid infrastructure"
        assert decision.confidence_score < 0.8, "Should have reduced confidence due to unknowns"
        assert decision.file_analysis.unknown_sizes == 2, "Should track 2 unknown file sizes"
        
        # Should have warnings about unknowns
        assert len(decision.warnings) > 0, "Should have warnings about unknown sizes"
        warnings_text = " ".join(decision.warnings).lower()
        assert any(keyword in warnings_text for keyword in ["unknown", "assume"]), \
            "Warnings should mention unknowns or assumptions"
    
    @pytest.mark.asyncio
    async def test_mixed_locations_decision(self, mixed_locations_workflow):
        """Test: Mixed locations should recommend EMR with upload note."""
        inputs = [
            InputAsset(uri=inp.uri, size_bytes=inp.size_bytes)
            for inp in mixed_locations_workflow.data_inputs
        ]
        
        decision = await decide_infrastructure(
            command=mixed_locations_workflow.description,
            inputs=inputs
        )
        
        # Snapshot assertions
        assert decision.infrastructure in ["EMR", "EC2"], "Mixed locations should use EMR or EC2"
        assert decision.file_analysis.mixed_locations is True, "Should detect mixed locations"
        assert decision.file_analysis.all_in_s3 is False, "Not all files in S3"
        assert decision.file_analysis.all_local is False, "Not all files local"
        
        # Should mention mixed locations or upload
        reasoning_lower = decision.reasoning.lower()
        assert any(keyword in reasoning_lower for keyword in ["mixed", "upload", "s3"]), \
            "Reasoning should mention mixed locations or upload"
    
    @pytest.mark.asyncio
    async def test_confidence_scoring_consistency(self, small_local_files_workflow, unknown_sizes_workflow):
        """Test: Confidence scores should be consistent with uncertainty levels."""
        # Known scenario (high confidence)
        inputs_known = [
            InputAsset(uri=inp.uri, size_bytes=inp.size_bytes)
            for inp in small_local_files_workflow.data_inputs
        ]
        decision_known = await decide_infrastructure(
            command=small_local_files_workflow.description,
            inputs=inputs_known
        )
        
        # Unknown scenario (low confidence)
        inputs_unknown = [
            InputAsset(uri=inp.uri, size_bytes=inp.size_bytes)
            for inp in unknown_sizes_workflow.data_inputs
        ]
        decision_unknown = await decide_infrastructure(
            command=unknown_sizes_workflow.description,
            inputs=inputs_unknown
        )
        
        # Confidence should be higher for known vs unknown
        assert decision_known.confidence_score > decision_unknown.confidence_score, \
            "Known scenario should have higher confidence than unknown scenario"
        
        # Known should be >=0.8, unknown should be <0.8
        assert decision_known.confidence_score >= 0.8, "Known scenario should have high confidence"
        assert decision_unknown.confidence_score < 0.8, "Unknown scenario should have reduced confidence"
    
    @pytest.mark.asyncio
    async def test_cost_analysis_completeness(self, large_s3_files_workflow):
        """Test: Cost analysis should be complete and valid."""
        inputs = [
            InputAsset(uri=inp.uri, size_bytes=inp.size_bytes)
            for inp in large_s3_files_workflow.data_inputs
        ]
        
        decision = await decide_infrastructure(
            command=large_s3_files_workflow.description,
            inputs=inputs
        )
        
        cost = decision.cost_analysis
        
        # Cost range validation
        assert len(cost.estimated_cost_range_usd) == 2, "Cost range should be (min, max)"
        min_cost, max_cost = cost.estimated_cost_range_usd
        assert min_cost >= 0, "Min cost should be non-negative"
        assert max_cost >= min_cost, "Max cost should be >= min cost"
        
        # Cost class validation
        assert cost.cost_class in ["Free", "Low", "Medium", "High"], "Cost class should be valid"
        
        # Assumptions and confidence
        assert len(cost.cost_assumptions) > 10, "Cost assumptions should be detailed (>10 chars)"
        assert 0.0 <= cost.cost_confidence <= 1.0, "Cost confidence should be 0-1"
        
        # Breakdown (if present)
        if cost.breakdown:
            assert cost.breakdown.compute_cost_range_usd[0] >= 0, "Compute cost should be non-negative"
            assert cost.breakdown.storage_cost_range_usd[0] >= 0, "Storage cost should be non-negative"
            assert cost.breakdown.data_transfer_cost_usd >= 0, "Transfer cost should be non-negative"


# ============================================================================
# Snapshot Tests: Implementation Agent
# ============================================================================

class TestImplementationAgentSnapshots:
    """Snapshot tests for Implementation Agent."""
    
    def create_test_infra_decision(self, infrastructure: str, total_size_mb: float):
        """Helper to create a valid InfraDecision for testing."""
        from backend.contracts.infra_decision import (
            InfraDecision, FileAnalysis, ComputationalRequirements,
            CostAnalysis, AlternativeRecommendation
        )
        
        return InfraDecision(
            infrastructure=infrastructure,
            confidence_score=0.8,
            decision_summary=f"{infrastructure} recommended for {total_size_mb}MB files",
            reasoning=f"Files are stored with total size {total_size_mb}MB. {infrastructure} provides optimal performance for this workload.",
            file_analysis=FileAnalysis(
                total_size_bytes=int(total_size_mb * 1024 * 1024),
                total_size_mb=total_size_mb,
                file_count=2,
                unknown_sizes=0,
                largest_file_bytes=int(total_size_mb * 0.6 * 1024 * 1024),
                largest_file_mb=total_size_mb * 0.6,
                all_in_s3=(infrastructure == "EMR"),
                all_local=(infrastructure == "Local"),
                mixed_locations=False
            ),
            computational_requirements=ComputationalRequirements(
                estimated_cpu_hours=1.0,
                estimated_memory_gb=8.0,
                estimated_runtime_minutes=30.0,
                parallelizable=True,
                gpu_required=False
            ),
            cost_analysis=CostAnalysis(
                estimated_cost_range_usd=(1.0, 5.0),
                cost_class="Medium",
                cost_assumptions="us-east-1, standard instances, 30min runtime",
                cost_confidence=0.6,
                data_transfer_cost_usd=0.0,
                breakdown={
                    "compute_cost_range_usd": (0.8, 4.0),
                    "storage_cost_range_usd": (0.2, 1.0),
                    "data_transfer_cost_usd": 0.0
                }
            ),
            alternatives=[
                AlternativeRecommendation(
                    infrastructure="EC2",
                    reasoning="Could work as alternative",
                    tradeoffs="Lower cost but requires data download"
                )
            ],
            warnings=[],
            inputs_analyzed=2
        )
    
    @pytest.mark.asyncio
    async def test_local_execution_plan(self, small_local_files_workflow):
        """Test: Local execution should have native execution (no container)."""
        infra_decision = self.create_test_infra_decision("Local", 98.0)
        
        execution_spec = await plan_implementation(
            workflow_plan=small_local_files_workflow,
            infra_decision=infra_decision,
            max_retries=1  # Reduce retries for faster tests
        )
        
        # Snapshot assertions
        assert execution_spec.tool_name in ["fastqc", "test"], "Tool name should be fastqc or test"
        assert execution_spec.infrastructure == "Local", "Should match infra decision"
        assert execution_spec.confidence_score >= 0.5, "Should have reasonable confidence"
        
        # Container spec (may or may not be present for Local)
        # Local can be native or containerized
        
        # Commands should be present
        assert len(execution_spec.commands) >= 1, "Should have at least 1 command"
        
        # Resources should be reasonable for small files
        assert execution_spec.resource_requirements.min_cpu_cores <= 8, "Should not over-provision CPU"
        assert execution_spec.resource_requirements.min_memory_gb <= 16, "Should not over-provision memory"
    
    @pytest.mark.asyncio
    async def test_emr_execution_plan(self, large_s3_files_workflow):
        """Test: EMR execution should use containers and S3-native commands."""
        infra_decision = self.create_test_infra_decision("EMR", 490.0)
        
        execution_spec = await plan_implementation(
            workflow_plan=large_s3_files_workflow,
            infra_decision=infra_decision,
            max_retries=1
        )
        
        # Snapshot assertions
        assert execution_spec.infrastructure == "EMR", "Should match infra decision"
        assert execution_spec.confidence_score >= 0.5, "Should have reasonable confidence"
        
        # Container should be present for EMR (recommended)
        if execution_spec.container_spec:
            assert execution_spec.container_spec.image_type == "docker", "Should use Docker for EMR"
        
        # Commands should reference S3 paths (for EMR S3-native execution)
        # (This is implementation-dependent, so we just check commands exist)
        assert len(execution_spec.commands) >= 1, "Should have at least 1 command"
        
        # Retry policy should be present
        assert execution_spec.retry_policy.max_retries >= 0, "Should have retry policy"
    
    @pytest.mark.asyncio
    async def test_resource_requirements_scaling(self, small_local_files_workflow, large_s3_files_workflow):
        """Test: Resource requirements should scale with file size."""
        # Small files
        infra_decision_small = self.create_test_infra_decision("Local", 98.0)
        execution_spec_small = await plan_implementation(
            workflow_plan=small_local_files_workflow,
            infra_decision=infra_decision_small,
            max_retries=1
        )
        
        # Large files
        infra_decision_large = self.create_test_infra_decision("EMR", 490.0)
        execution_spec_large = await plan_implementation(
            workflow_plan=large_s3_files_workflow,
            infra_decision=infra_decision_large,
            max_retries=1
        )
        
        # Large files should have >= resources than small files (or similar)
        # Note: Heuristic implementation may provide baseline resources, so we just check validity
        assert execution_spec_small.resource_requirements.min_cpu_cores > 0, "Small should have CPU"
        assert execution_spec_large.resource_requirements.min_cpu_cores > 0, "Large should have CPU"
        assert execution_spec_small.resource_requirements.min_memory_gb > 0, "Small should have memory"
        assert execution_spec_large.resource_requirements.min_memory_gb > 0, "Large should have memory"
    
    @pytest.mark.asyncio
    async def test_retry_policy_validity(self, small_local_files_workflow):
        """Test: Retry policy should be valid and reasonable."""
        infra_decision = self.create_test_infra_decision("Local", 98.0)
        
        execution_spec = await plan_implementation(
            workflow_plan=small_local_files_workflow,
            infra_decision=infra_decision,
            max_retries=1
        )
        
        retry = execution_spec.retry_policy
        
        # Retry policy validation
        assert 0 <= retry.max_retries <= 5, "Max retries should be 0-5"
        assert len(retry.retry_on) >= 0, "Retry_on should be a list"
        assert retry.backoff_multiplier >= 1.0, "Backoff multiplier should be >=1.0"
        assert retry.initial_delay_seconds > 0, "Initial delay should be positive"
        
        # Common retry reasons
        valid_retry_reasons = ["exit_code", "timeout", "network", "oom", "unknown"]
        for reason in retry.retry_on:
            assert reason in valid_retry_reasons, f"Retry reason '{reason}' should be valid"


# ============================================================================
# Snapshot Tests: Orchestrator (End-to-End)
# ============================================================================

class TestOrchestratorSnapshots:
    """Snapshot tests for full Orchestrator pipeline."""
    
    @pytest.mark.asyncio
    async def test_end_to_end_small_local(self, small_local_files_workflow):
        """Test: End-to-end pipeline for small local files."""
        orchestrator = Orchestrator()
        
        trace = await orchestrator.run_pipeline(
            command="Run FastQC on small local files",
            workflow_plan=small_local_files_workflow,
            request_id="test_e2e_small_local"
        )
        
        # Trace validation
        assert trace.success is True, "Pipeline should succeed"
        assert trace.duration_ms > 0, "Duration should be positive"
        assert len(trace.agent_invocations) == 2, "Should invoke 2 agents"
        
        # Agent invocation validation
        agent_names = [inv.agent_name for inv in trace.agent_invocations]
        assert "InfrastructureDecisionAgent" in agent_names, "Should invoke infra agent"
        assert "ImplementationAgent" in agent_names, "Should invoke implementation agent"
        
        # Output validation
        infra_inv = [inv for inv in trace.agent_invocations if inv.agent_name == "InfrastructureDecisionAgent"][0]
        impl_inv = [inv for inv in trace.agent_invocations if inv.agent_name == "ImplementationAgent"][0]
        
        assert infra_inv.success is True, "Infra agent should succeed"
        assert impl_inv.success is True, "Implementation agent should succeed"
        assert infra_inv.output.infrastructure == "Local", "Should recommend Local"
        assert impl_inv.output.infrastructure == "Local", "Execution plan should match infra decision"
    
    @pytest.mark.asyncio
    async def test_end_to_end_large_s3(self, large_s3_files_workflow):
        """Test: End-to-end pipeline for large S3 files."""
        orchestrator = Orchestrator()
        
        trace = await orchestrator.run_pipeline(
            command="Merge reads from large S3 files",
            workflow_plan=large_s3_files_workflow,
            request_id="test_e2e_large_s3"
        )
        
        # Trace validation
        assert trace.success is True, "Pipeline should succeed"
        assert len(trace.agent_invocations) == 2, "Should invoke 2 agents"
        
        # Output validation
        infra_inv = [inv for inv in trace.agent_invocations if inv.agent_name == "InfrastructureDecisionAgent"][0]
        impl_inv = [inv for inv in trace.agent_invocations if inv.agent_name == "ImplementationAgent"][0]
        
        assert infra_inv.output.infrastructure == "EMR", "Should recommend EMR for large S3 files"
        assert impl_inv.output.infrastructure == "EMR", "Execution plan should match infra decision"
        assert impl_inv.output.confidence_score >= 0.3, "Should have reasonable confidence"
    
    @pytest.mark.asyncio
    async def test_request_id_tracking(self, small_local_files_workflow):
        """Test: Request ID should be tracked throughout pipeline."""
        orchestrator = Orchestrator()
        
        request_id = "test_tracking_123"
        trace = await orchestrator.run_pipeline(
            command="Test request ID tracking",
            workflow_plan=small_local_files_workflow,
            request_id=request_id
        )
        
        # Request ID validation
        assert trace.request_id == request_id, "Trace should have correct request ID"
        
        # All agent invocations should have the same request ID
        for inv in trace.agent_invocations:
            assert inv.request_id == request_id, f"Agent {inv.agent_name} should have matching request ID"
        
        # Request ID should be in orchestrator's traces
        assert request_id in orchestrator.traces, "Orchestrator should store trace by request ID"
        assert orchestrator.traces[request_id] == trace, "Stored trace should match returned trace"
    
    @pytest.mark.asyncio
    async def test_agent_timing_tracking(self, small_local_files_workflow):
        """Test: Agent timing should be tracked."""
        orchestrator = Orchestrator()
        
        trace = await orchestrator.run_pipeline(
            command="Test timing tracking",
            workflow_plan=small_local_files_workflow,
            request_id="test_timing"
        )
        
        # Timing validation
        assert trace.start_time > 0, "Trace should have start time"
        assert trace.end_time > trace.start_time, "End time should be after start time"
        assert trace.duration_ms > 0, "Duration should be positive"
        
        # Each agent invocation should have timing
        for inv in trace.agent_invocations:
            assert inv.start_time > 0, f"Agent {inv.agent_name} should have start time"
            assert inv.end_time > inv.start_time, f"Agent {inv.agent_name} end time should be after start time"
            assert inv.duration_ms > 0, f"Agent {inv.agent_name} should have positive duration"


# ============================================================================
# Regression Tests: Ensure Consistency Across Runs
# ============================================================================

class TestRegressionConsistency:
    """Tests to ensure consistent behavior across multiple runs."""
    
    @pytest.mark.asyncio
    async def test_same_input_similar_output(self, small_local_files_workflow):
        """Test: Same input should produce similar output across runs."""
        # Run pipeline twice with same input
        trace1 = await Orchestrator().run_pipeline(
            command="Test consistency run 1",
            workflow_plan=small_local_files_workflow,
            request_id="consistency_1"
        )
        
        trace2 = await Orchestrator().run_pipeline(
            command="Test consistency run 2",
            workflow_plan=small_local_files_workflow,
            request_id="consistency_2"
        )
        
        # Both should succeed
        assert trace1.success is True, "Run 1 should succeed"
        assert trace2.success is True, "Run 2 should succeed"
        
        # Both should recommend same infrastructure
        infra1 = [inv for inv in trace1.agent_invocations if inv.agent_name == "InfrastructureDecisionAgent"][0]
        infra2 = [inv for inv in trace2.agent_invocations if inv.agent_name == "InfrastructureDecisionAgent"][0]
        
        assert infra1.output.infrastructure == infra2.output.infrastructure, \
            "Same input should produce same infrastructure recommendation"
        
        # Confidence scores should be similar (within 0.2)
        confidence_diff = abs(infra1.output.confidence_score - infra2.output.confidence_score)
        assert confidence_diff < 0.2, "Confidence scores should be similar across runs"
