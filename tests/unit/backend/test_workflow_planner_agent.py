"""
Unit tests for Workflow Planner Agent.

Tests workflow selection, task decomposition, clarification handling,
and feasibility assessment.
"""

import pytest
from backend.workflow_planner_agent import (
    plan_workflow,
    RNASeqPlaybook,
    ScRNASeqPlaybook,
    WGSWESPlaybook,
    PhylogeneticsPlaybook,
    SimpleOperationPlaybook,
    _extract_file_uris,
    _extract_parameter,
    _is_fastq,
    _is_fasta,
    get_playbook_for_workflow_type
)
from backend.contracts.workflow_plan import DataInput


# Fixtures

@pytest.fixture
def session_context():
    """Standard session context for tests."""
    return {
        "session_id": "test_session_123",
        "metadata": {
            "s3_bucket": "test-bucket",
            "uploaded_files": []
        }
    }


# Test file detection helpers

def test_is_fastq():
    """Test FASTQ file detection."""
    assert _is_fastq("sample.fastq")
    assert _is_fastq("sample.fq")
    assert _is_fastq("sample.fastq.gz")
    assert _is_fastq("sample.fq.gz")
    assert _is_fastq("s3://bucket/sample_R1.fastq")
    assert not _is_fastq("sample.fasta")
    assert not _is_fastq("sample.bam")


def test_is_fasta():
    """Test FASTA file detection."""
    assert _is_fasta("sample.fasta")
    assert _is_fasta("sample.fa")
    assert _is_fasta("sample.fna")
    assert _is_fasta("sample.faa")
    assert _is_fasta("sample.fasta.gz")
    assert not _is_fasta("sample.fastq")
    assert not _is_fasta("sample.bam")


# Test file URI extraction

def test_extract_file_uris_from_s3(session_context):
    """Test extracting S3 URIs from command."""
    command = "Analyze s3://bucket/sample_R1.fq and s3://bucket/sample_R2.fq"
    uris = _extract_file_uris(command, session_context)
    
    assert len(uris) == 2
    assert "s3://bucket/sample_R1.fq" in uris
    assert "s3://bucket/sample_R2.fq" in uris


def test_extract_file_uris_from_local_path(session_context):
    """Test extracting local file paths from command."""
    command = "Analyze /data/sample.fastq"
    uris = _extract_file_uris(command, session_context)
    
    assert len(uris) == 1
    assert "/data/sample.fastq" in uris


def test_extract_file_uris_from_session(session_context):
    """Test extracting file URIs from session context."""
    session_context["metadata"]["uploaded_files"] = [
        {"s3_key": "uploads/file1.fq", "s3_bucket": "test-bucket"},
        {"s3_key": "uploads/file2.fq", "s3_bucket": "test-bucket"}
    ]
    
    command = "Analyze uploaded files"
    uris = _extract_file_uris(command, session_context)
    
    assert len(uris) == 2
    assert "s3://test-bucket/uploads/file1.fq" in uris
    assert "s3://test-bucket/uploads/file2.fq" in uris


# Test parameter extraction

def test_extract_reference_genome():
    """Test extracting reference genome from command."""
    command = "RNA-seq analysis with reference genome hg38"
    genome = _extract_parameter("reference_genome", command, {})
    assert genome == "hg38"
    
    command = "Align to mm10"
    genome = _extract_parameter("reference_genome", command, {})
    assert genome == "mm10"


def test_extract_organism():
    """Test extracting organism from command."""
    command = "RNA-seq analysis for human samples"
    organism = _extract_parameter("organism", command, {})
    assert organism == "human"
    
    command = "Mouse genome analysis"
    organism = _extract_parameter("organism", command, {})
    assert organism == "mouse"


def test_extract_strandedness():
    """Test extracting strandedness from command."""
    command = "RNA-seq with unstranded library"
    strand = _extract_parameter("strandedness", command, {})
    assert strand == "unstranded"
    
    command = "Forward stranded RNA-seq"
    strand = _extract_parameter("strandedness", command, {})
    assert strand == "forward"
    
    command = "Reverse strand library"
    strand = _extract_parameter("strandedness", command, {})
    assert strand == "reverse"


# Test playbook matching

def test_rnaseq_playbook_matches():
    """Test RNA-seq playbook matching."""
    inputs = [
        DataInput(uri="s3://bucket/sample_R1.fastq", size_bytes=1000000),
        DataInput(uri="s3://bucket/sample_R2.fastq", size_bytes=1000000)
    ]
    
    # Should match with RNA-seq keywords
    assert RNASeqPlaybook.matches(inputs, "RNA-seq analysis", {})
    assert RNASeqPlaybook.matches(inputs, "Analyze RNA seq data", {})
    assert RNASeqPlaybook.matches(inputs, "Differential expression analysis", {})
    
    # Should not match without FASTQ files
    fasta_inputs = [DataInput(uri="s3://bucket/sample.fasta", size_bytes=1000000)]
    assert not RNASeqPlaybook.matches(fasta_inputs, "RNA-seq analysis", {})


def test_scrnaseq_playbook_matches():
    """Test scRNA-seq playbook matching."""
    inputs = [DataInput(uri="s3://bucket/sample.h5ad", size_bytes=1000000)]
    
    # Should match with scRNA-seq keywords
    assert ScRNASeqPlaybook.matches(inputs, "Single-cell RNA-seq analysis", {})
    assert ScRNASeqPlaybook.matches(inputs, "scRNA-seq clustering", {})
    
    # Should match with h5ad format even without keywords
    assert ScRNASeqPlaybook.matches(inputs, "Analyze this data", {})


def test_wgs_playbook_matches():
    """Test WGS/WES playbook matching."""
    inputs = [
        DataInput(uri="s3://bucket/sample_R1.fastq", size_bytes=1000000),
        DataInput(uri="s3://bucket/sample_R2.fastq", size_bytes=1000000)
    ]
    
    # Should match with WGS/variant keywords
    assert WGSWESPlaybook.matches(inputs, "WGS variant calling", {})
    assert WGSWESPlaybook.matches(inputs, "Whole genome sequencing", {})
    assert WGSWESPlaybook.matches(inputs, "Call variants from these samples", {})
    
    # Should not match RNA-seq commands
    assert not WGSWESPlaybook.matches(inputs, "RNA-seq analysis", {})


def test_phylogenetics_playbook_matches():
    """Test phylogenetics playbook matching."""
    inputs = [DataInput(uri="s3://bucket/sequences.fasta", size_bytes=1000000)]
    
    # Should match with phylo/tree keywords
    assert PhylogeneticsPlaybook.matches(inputs, "Build phylogenetic tree", {})
    assert PhylogeneticsPlaybook.matches(inputs, "Phylogenetic analysis", {})
    assert PhylogeneticsPlaybook.matches(inputs, "Create MSA and tree", {})
    
    # Should not match without FASTA
    fastq_inputs = [DataInput(uri="s3://bucket/sample.fastq", size_bytes=1000000)]
    assert not PhylogeneticsPlaybook.matches(fastq_inputs, "Build phylogenetic tree", {})


def test_simple_operation_playbook_matches():
    """Test SimpleOperationPlaybook always matches (fallback)."""
    inputs = [DataInput(uri="s3://bucket/file.txt", size_bytes=1000000)]
    assert SimpleOperationPlaybook.matches(inputs, "Do something", {})


# Test playbook workflow plan creation

@pytest.mark.asyncio
async def test_plan_rnaseq_workflow_success(session_context):
    """Test successful RNA-seq workflow planning."""
    command = "RNA-seq analysis for human samples with reference genome hg38 from s3://bucket/sample_R1.fastq and s3://bucket/sample_R2.fastq"
    
    result = await plan_workflow(command, session_context)
    
    assert result["status"] == "success"
    assert "workflow_plan" in result
    
    plan = result["workflow_plan"]
    assert plan.description
    assert len(plan.operations) > 0
    assert any(op.operation_name == "quality_control" for op in plan.operations)
    assert any(op.operation_name == "alignment" for op in plan.operations)
    assert any(op.operation_name == "quantification" for op in plan.operations)


@pytest.mark.asyncio
async def test_plan_rnaseq_workflow_missing_params(session_context):
    """Test RNA-seq workflow planning with missing parameters."""
    command = "RNA-seq analysis of s3://bucket/sample_R1.fastq and s3://bucket/sample_R2.fastq"
    
    result = await plan_workflow(command, session_context)
    
    # Should request clarification for missing organism/reference_genome
    assert result["status"] == "clarification_needed"
    assert "missing_parameters" in result
    
    missing = result["missing_parameters"]
    param_names = [p["parameter"] for p in missing]
    assert "organism" in param_names
    assert "reference_genome" in param_names


@pytest.mark.asyncio
async def test_plan_simple_operation(session_context):
    """Test simple single-operation planning."""
    command = "Merge s3://bucket/R1.fq and s3://bucket/R2.fq"
    
    result = await plan_workflow(command, session_context)
    
    assert result["status"] == "success"
    assert "workflow_plan" in result
    
    plan = result["workflow_plan"]
    assert len(plan.operations) == 1
    assert plan.operations[0].operation_name == "read_merging"


@pytest.mark.asyncio
async def test_plan_phylogenetics_workflow(session_context):
    """Test phylogenetics workflow planning."""
    command = "Build phylogenetic tree from s3://bucket/sequences.fasta"
    
    result = await plan_workflow(command, session_context)
    
    assert result["status"] == "success"
    assert "workflow_plan" in result
    
    plan = result["workflow_plan"]
    assert len(plan.operations) > 0
    assert any(op.operation_name == "multiple_sequence_alignment" for op in plan.operations)
    assert any(op.operation_name == "tree_inference" for op in plan.operations)


@pytest.mark.asyncio
async def test_plan_phylogenetics_with_existing_alignment(session_context):
    """Test phylogenetics planning when alignment exists in session."""
    # Add aligned sequences to session context
    session_context["aligned_sequences"] = "ATGC..."
    
    command = "Build phylogenetic tree from s3://bucket/sequences.fasta"
    
    result = await plan_workflow(command, session_context)
    
    assert result["status"] == "success"
    plan = result["workflow_plan"]
    
    # Should skip MSA step since alignment exists
    op_names = [op.operation_name for op in plan.operations]
    assert "multiple_sequence_alignment" not in op_names
    assert "tree_inference" in op_names


# Test playbook registry

def test_get_playbook_for_workflow_type():
    """Test getting playbook by workflow type."""
    assert get_playbook_for_workflow_type("rna_seq_bulk") == RNASeqPlaybook
    assert get_playbook_for_workflow_type("scrna_seq") == ScRNASeqPlaybook
    assert get_playbook_for_workflow_type("wgs_wes") == WGSWESPlaybook
    assert get_playbook_for_workflow_type("phylogenetics") == PhylogeneticsPlaybook
    assert get_playbook_for_workflow_type("simple_operation") == SimpleOperationPlaybook
    assert get_playbook_for_workflow_type("nonexistent") is None


# Test workflow plan structure

@pytest.mark.asyncio
async def test_workflow_plan_has_required_fields(session_context):
    """Test that workflow plans have all required fields."""
    command = "RNA-seq analysis for human hg38 samples from s3://bucket/sample.fastq"
    
    result = await plan_workflow(command, session_context)
    
    assert result["status"] == "success"
    plan = result["workflow_plan"]
    
    # Check required fields
    assert hasattr(plan, "description")
    assert hasattr(plan, "data_inputs")
    assert hasattr(plan, "operations")
    assert hasattr(plan, "constraints")
    assert hasattr(plan, "expected_compute_intensity")
    assert hasattr(plan, "session_id")
    
    # Check session_id is set
    assert plan.session_id == "test_session_123"


@pytest.mark.asyncio
async def test_workflow_plan_operations_have_required_fields(session_context):
    """Test that operations have all required fields."""
    command = "Merge s3://bucket/R1.fq and s3://bucket/R2.fq"
    
    result = await plan_workflow(command, session_context)
    
    assert result["status"] == "success"
    plan = result["workflow_plan"]
    
    for operation in plan.operations:
        assert hasattr(operation, "operation_name")
        assert hasattr(operation, "tool_name")
        assert hasattr(operation, "parameters")
        assert hasattr(operation, "expected_output_size_mb")
        assert hasattr(operation, "parallelizable")
        
        # Check operation_name is non-empty
        assert operation.operation_name


# Test edge cases

@pytest.mark.asyncio
async def test_plan_workflow_no_files(session_context):
    """Test workflow planning with no input files."""
    command = "RNA-seq analysis"
    
    result = await plan_workflow(command, session_context)
    
    # Should still work (fallback to simple operation)
    assert result["status"] in ["success", "clarification_needed"]


@pytest.mark.asyncio
async def test_plan_workflow_with_defaults(session_context):
    """Test that default parameters are applied correctly."""
    command = "RNA-seq analysis for human hg38 from s3://bucket/sample.fastq"
    
    result = await plan_workflow(command, session_context)
    
    assert result["status"] == "success"
    plan = result["workflow_plan"]
    
    # Check that strandedness default is applied
    alignment_op = next(
        (op for op in plan.operations if op.operation_name == "alignment"),
        None
    )
    if alignment_op:
        assert "strandedness" in alignment_op.parameters
        assert alignment_op.parameters["strandedness"] == "unstranded"


# Test operation detection

def test_simple_operation_detect_merge():
    """Test operation detection for merge."""
    inputs = [DataInput(uri="s3://bucket/R1.fq"), DataInput(uri="s3://bucket/R2.fq")]
    command = "Merge these reads"
    
    op_name, tool_name = SimpleOperationPlaybook.detect_operation(command, inputs)
    
    assert op_name == "read_merging"
    assert tool_name == "bbmerge"


def test_simple_operation_detect_trim():
    """Test operation detection for trim."""
    inputs = [DataInput(uri="s3://bucket/sample.fq")]
    command = "Trim low quality bases"
    
    op_name, tool_name = SimpleOperationPlaybook.detect_operation(command, inputs)
    
    assert op_name == "trimming"
    assert tool_name == "trimmomatic"


def test_simple_operation_detect_align():
    """Test operation detection for align."""
    inputs = [DataInput(uri="s3://bucket/sample.fq")]
    command = "Align to reference genome"
    
    op_name, tool_name = SimpleOperationPlaybook.detect_operation(command, inputs)
    
    assert op_name == "alignment"
    assert tool_name == "bwa_mem"


def test_simple_operation_detect_qc():
    """Test operation detection for QC."""
    inputs = [DataInput(uri="s3://bucket/sample.fq")]
    command = "Run quality control"
    
    op_name, tool_name = SimpleOperationPlaybook.detect_operation(command, inputs)
    
    assert op_name == "quality_control"
    assert tool_name == "fastqc"


# Run tests
if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
