"""
Workflow Planner Agent - Decides WHAT workflow to execute.

This agent analyzes user execution requests and produces a detailed workflow plan.
It selects appropriate bioinformatics workflows, decomposes them into atomic tasks,
and asks clarifying questions when needed.

Responsibilities:
- Workflow selection (choose appropriate playbook)
- Task decomposition (break into atomic operations)
- Clarification (ask minimum high-impact questions)
- Feasibility assessment (validate inputs and parameters)
- Workflow proposal (suggest workflows when user is uncertain)
"""

import logging
import re
from typing import Dict, Any, List, Optional, Tuple
from pathlib import Path

from backend.contracts.workflow_plan import (
    WorkflowPlan,
    DataInput,
    OperationSpec,
    ConstraintSpec
)
from backend.routing_keywords import SIMPLE_OPERATION_REGEX

logger = logging.getLogger(__name__)


class WorkflowPlaybook:
    """Base class for workflow playbooks."""
    
    workflow_type: str = "unknown"
    description: str = "Unknown workflow"
    
    @classmethod
    def matches(cls, inputs: List[DataInput], command: str, session_context: Dict[str, Any]) -> bool:
        """Check if this playbook matches the user request."""
        raise NotImplementedError
    
    @classmethod
    def required_parameters(cls) -> List[Dict[str, Any]]:
        """Return list of required parameters for this workflow."""
        return []
    
    @classmethod
    def create_workflow_plan(
        cls,
        inputs: List[DataInput],
        parameters: Dict[str, Any],
        session_context: Dict[str, Any]
    ) -> WorkflowPlan:
        """Create a workflow plan for this playbook."""
        raise NotImplementedError


class RNASeqPlaybook(WorkflowPlaybook):
    """RNA-seq (bulk) workflow playbook."""
    
    workflow_type = "rna_seq_bulk"
    description = "Bulk RNA-seq differential expression analysis"
    
    @classmethod
    def matches(cls, inputs: List[DataInput], command: str, session_context: Dict[str, Any]) -> bool:
        """Match if: FASTQ files + RNA-seq keywords, but NOT amplicon/microbiome data."""
        has_fastq = any(_is_fastq(inp.uri) for inp in inputs)
        command_lower = command.lower()

        # Explicitly exclude amplicon sequencing / microbiome pipelines.
        # "rRNA" contains "rna" but is 16S amplicon, not bulk RNA-seq.
        amplicon_exclusions = [
            "16s", "18s", "its ", "amplicon", "microbiome",
            "v3-v4", "v3–v4", "v4-v5", "hypervariable",
            "16s rrna", "rrna amplicon",
        ]
        if any(kw in command_lower for kw in amplicon_exclusions):
            return False

        rna_keywords = ["rna-seq", "rnaseq", "rna seq", "transcriptome", "expression", "differential"]
        has_rna_keyword = any(kw in command_lower for kw in rna_keywords)
        # Only match "rna" as standalone indicator when it is not preceded by "r" (i.e. avoid "rrna")
        import re as _re
        has_bare_rna = bool(_re.search(r'(?<![a-z])rna(?!-seq)', command_lower))
        return has_fastq and (has_rna_keyword or has_bare_rna)
    
    @classmethod
    def required_parameters(cls) -> List[Dict[str, Any]]:
        return [
            {
                "parameter": "organism",
                "required": True,
                "description": "Organism for the reference genome (e.g., 'human', 'mouse')"
            },
            {
                "parameter": "reference_genome",
                "required": True,
                "description": "Reference genome build (e.g., 'hg38', 'mm10')"
            },
            {
                "parameter": "strandedness",
                "required": False,
                "description": "Strandedness (unstranded, forward, reverse)",
                "default": "unstranded"
            }
        ]
    
    @classmethod
    def create_workflow_plan(
        cls,
        inputs: List[DataInput],
        parameters: Dict[str, Any],
        session_context: Dict[str, Any]
    ) -> WorkflowPlan:
        """Create RNA-seq workflow plan."""
        operations = [
            OperationSpec(
                operation_name="quality_control",
                tool_name="fastqc",
                parameters={},
                expected_output_size_mb=10.0,
                parallelizable=True
            ),
            OperationSpec(
                operation_name="trimming",
                tool_name="trimmomatic",
                parameters={"min_length": 36},
                expected_output_size_mb=sum(inp.size_bytes or 0 for inp in inputs) / (1024 * 1024) * 0.8,
                parallelizable=False
            ),
            OperationSpec(
                operation_name="alignment",
                tool_name="STAR",
                parameters={
                    "reference_genome": parameters.get("reference_genome"),
                    "strandedness": parameters.get("strandedness", "unstranded")
                },
                expected_output_size_mb=sum(inp.size_bytes or 0 for inp in inputs) / (1024 * 1024) * 2.0,
                parallelizable=False
            ),
            OperationSpec(
                operation_name="quantification",
                tool_name="featureCounts",
                parameters={
                    "feature_type": "exon",
                    "attribute": "gene_id"
                },
                expected_output_size_mb=5.0,
                parallelizable=False
            ),
            OperationSpec(
                operation_name="differential_expression",
                tool_name="DESeq2",
                parameters={},
                expected_output_size_mb=2.0,
                parallelizable=False
            ),
            OperationSpec(
                operation_name="visualization",
                tool_name="plotly",
                parameters={"plot_types": ["PCA", "volcano", "heatmap"]},
                expected_output_size_mb=1.0,
                parallelizable=False
            )
        ]
        
        return WorkflowPlan(
            description=f"RNA-seq analysis pipeline for {parameters.get('organism', 'unknown organism')}",
            data_inputs=inputs,
            operations=operations,
            constraints=ConstraintSpec(
                reproducibility_required=True,
                containerization_preferred=True
            ),
            expected_compute_intensity="High",
            session_id=session_context.get("session_id")
        )


class ScRNASeqPlaybook(WorkflowPlaybook):
    """Single-cell RNA-seq workflow playbook."""
    
    workflow_type = "scrna_seq"
    description = "Single-cell RNA-seq analysis"
    
    @classmethod
    def matches(cls, inputs: List[DataInput], command: str, session_context: Dict[str, Any]) -> bool:
        """Match if: h5ad/mtx files or scRNA-seq keywords."""
        command_lower = command.lower()
        sc_keywords = ["single-cell", "scrna", "sc-rna", "single cell", "scanpy", "seurat"]
        has_sc_keyword = any(kw in command_lower for kw in sc_keywords)
        
        has_sc_format = any(
            inp.uri.endswith(('.h5ad', '.mtx', '.loom'))
            for inp in inputs
        )
        
        return has_sc_keyword or has_sc_format
    
    @classmethod
    def required_parameters(cls) -> List[Dict[str, Any]]:
        return [
            {
                "parameter": "input_format",
                "required": False,
                "description": "Input format (h5ad, mtx, loom)",
                "default": "h5ad"
            }
        ]
    
    @classmethod
    def create_workflow_plan(
        cls,
        inputs: List[DataInput],
        parameters: Dict[str, Any],
        session_context: Dict[str, Any]
    ) -> WorkflowPlan:
        """Create scRNA-seq workflow plan."""
        operations = [
            OperationSpec(
                operation_name="quality_control",
                tool_name="scanpy_qc",
                parameters={"min_genes": 200, "min_cells": 3},
                expected_output_size_mb=10.0,
                parallelizable=False
            ),
            OperationSpec(
                operation_name="normalization",
                tool_name="scanpy_normalize",
                parameters={"method": "log"},
                expected_output_size_mb=20.0,
                parallelizable=False
            ),
            OperationSpec(
                operation_name="feature_selection",
                tool_name="scanpy_hvg",
                parameters={"n_top_genes": 2000},
                expected_output_size_mb=5.0,
                parallelizable=False
            ),
            OperationSpec(
                operation_name="dimensionality_reduction",
                tool_name="scanpy_pca",
                parameters={"n_comps": 50},
                expected_output_size_mb=15.0,
                parallelizable=False
            ),
            OperationSpec(
                operation_name="neighbor_graph",
                tool_name="scanpy_neighbors",
                parameters={"n_neighbors": 15},
                expected_output_size_mb=10.0,
                parallelizable=False
            ),
            OperationSpec(
                operation_name="embedding",
                tool_name="scanpy_umap",
                parameters={},
                expected_output_size_mb=5.0,
                parallelizable=False
            ),
            OperationSpec(
                operation_name="clustering",
                tool_name="scanpy_leiden",
                parameters={"resolution": 1.0},
                expected_output_size_mb=5.0,
                parallelizable=False
            ),
            OperationSpec(
                operation_name="marker_identification",
                tool_name="scanpy_rank_genes",
                parameters={},
                expected_output_size_mb=10.0,
                parallelizable=False
            )
        ]
        
        return WorkflowPlan(
            description="Single-cell RNA-seq analysis pipeline",
            data_inputs=inputs,
            operations=operations,
            constraints=ConstraintSpec(
                reproducibility_required=True,
                containerization_preferred=True
            ),
            expected_compute_intensity="High",
            session_id=session_context.get("session_id")
        )


class WGSWESPlaybook(WorkflowPlaybook):
    """Whole genome/exome sequencing workflow playbook."""
    
    workflow_type = "wgs_wes"
    description = "Whole genome or exome sequencing variant calling"
    
    @classmethod
    def matches(cls, inputs: List[DataInput], command: str, session_context: Dict[str, Any]) -> bool:
        """Match if: FASTQ files + WGS/WES/variant keywords."""
        has_fastq = any(_is_fastq(inp.uri) for inp in inputs)
        command_lower = command.lower()
        wgs_keywords = ["wgs", "wes", "whole genome", "whole exome", "variant", "snp", "indel", "mutation"]
        has_wgs_keyword = any(kw in command_lower for kw in wgs_keywords)
        return has_fastq and has_wgs_keyword
    
    @classmethod
    def required_parameters(cls) -> List[Dict[str, Any]]:
        return [
            {
                "parameter": "reference_genome",
                "required": True,
                "description": "Reference genome build (e.g., 'hg38', 'hg19')"
            }
        ]
    
    @classmethod
    def create_workflow_plan(
        cls,
        inputs: List[DataInput],
        parameters: Dict[str, Any],
        session_context: Dict[str, Any]
    ) -> WorkflowPlan:
        """Create WGS/WES workflow plan."""
        operations = [
            OperationSpec(
                operation_name="quality_control",
                tool_name="fastqc",
                parameters={},
                expected_output_size_mb=10.0,
                parallelizable=True
            ),
            OperationSpec(
                operation_name="alignment",
                tool_name="bwa_mem",
                parameters={"reference_genome": parameters.get("reference_genome")},
                expected_output_size_mb=sum(inp.size_bytes or 0 for inp in inputs) / (1024 * 1024) * 1.5,
                parallelizable=False
            ),
            OperationSpec(
                operation_name="duplicate_marking",
                tool_name="picard_markduplicates",
                parameters={},
                expected_output_size_mb=sum(inp.size_bytes or 0 for inp in inputs) / (1024 * 1024) * 1.5,
                parallelizable=False
            ),
            OperationSpec(
                operation_name="variant_calling",
                tool_name="gatk_haplotypecaller",
                parameters={},
                expected_output_size_mb=100.0,
                parallelizable=False
            ),
            OperationSpec(
                operation_name="variant_filtering",
                tool_name="gatk_variantfiltration",
                parameters={},
                expected_output_size_mb=50.0,
                parallelizable=False
            ),
            OperationSpec(
                operation_name="variant_annotation",
                tool_name="vep",
                parameters={},
                expected_output_size_mb=150.0,
                parallelizable=False
            )
        ]
        
        return WorkflowPlan(
            description=f"WGS/WES variant calling pipeline ({parameters.get('reference_genome', 'unknown reference')})",
            data_inputs=inputs,
            operations=operations,
            constraints=ConstraintSpec(
                reproducibility_required=True,
                containerization_preferred=True
            ),
            expected_compute_intensity="High",
            session_id=session_context.get("session_id")
        )


class PhylogeneticsPlaybook(WorkflowPlaybook):
    """Phylogenetic analysis workflow playbook."""
    
    workflow_type = "phylogenetics"
    description = "Phylogenetic tree inference from sequences"
    
    @classmethod
    def matches(cls, inputs: List[DataInput], command: str, session_context: Dict[str, Any]) -> bool:
        """Match if: FASTA files + phylo/tree keywords."""
        has_fasta = any(_is_fasta(inp.uri) for inp in inputs)
        command_lower = command.lower()
        phylo_keywords = ["phylogen", "tree", "msa", "alignment", "evolution", "ancestor"]
        has_phylo_keyword = any(kw in command_lower for kw in phylo_keywords)
        return has_fasta and has_phylo_keyword
    
    @classmethod
    def required_parameters(cls) -> List[Dict[str, Any]]:
        return []
    
    @classmethod
    def create_workflow_plan(
        cls,
        inputs: List[DataInput],
        parameters: Dict[str, Any],
        session_context: Dict[str, Any]
    ) -> WorkflowPlan:
        """Create phylogenetics workflow plan."""
        # Check if alignment already exists in session context
        has_alignment = "aligned_sequences" in session_context
        
        operations = []
        
        if not has_alignment:
            operations.append(
                OperationSpec(
                    operation_name="multiple_sequence_alignment",
                    tool_name="mafft",
                    parameters={},
                    expected_output_size_mb=sum(inp.size_bytes or 0 for inp in inputs) / (1024 * 1024) * 1.2,
                    parallelizable=False
                )
            )
        
        operations.extend([
            OperationSpec(
                operation_name="tree_inference",
                tool_name="iqtree",
                parameters={"bootstrap": 1000},
                expected_output_size_mb=5.0,
                parallelizable=False
            ),
            OperationSpec(
                operation_name="tree_visualization",
                tool_name="ete3",
                parameters={},
                expected_output_size_mb=1.0,
                parallelizable=False
            )
        ])
        
        return WorkflowPlan(
            description="Phylogenetic tree inference and visualization",
            data_inputs=inputs if not has_alignment else [],
            operations=operations,
            constraints=ConstraintSpec(
                reproducibility_required=True,
                containerization_preferred=False
            ),
            expected_compute_intensity="Medium",
            session_id=session_context.get("session_id")
        )


class ResultsBrowsingPlaybook(WorkflowPlaybook):
    """
    Browse existing results in S3 (list objects, fetch results.json, etc).

    This must take priority over bioinformatics execution when the user asks
    to *display/show/list* results at an S3 prefix.
    """

    workflow_type = "s3_results_browse"
    description = "Browse and display results from S3"

    @classmethod
    def matches(cls, inputs: List[DataInput], command: str, session_context: Dict[str, Any]) -> bool:
        command_lower = (command or "").lower()
        has_s3 = "s3://" in command_lower or any(
            isinstance(inp.uri, str) and inp.uri.startswith("s3://") for inp in (inputs or [])
        )
        if not has_s3:
            return False

        browse_verbs = [
            "display",
            "show",
            "list",
            "view",
            "browse",
            "fetch",
        ]
        browse_phrases = [
            "list objects",
            "list files",
            "fetch and show",
            "show results",
            "display the results",
        ]

        looks_like_browse = any(v in command_lower for v in browse_verbs) or any(
            p in command_lower for p in browse_phrases
        )
        looks_like_results = "results.json" in command_lower or "fastqc" in command_lower

        # If the user explicitly asks to run/perform analysis, don't treat as browsing.
        looks_like_execution = any(
            w in command_lower for w in ["run fastqc", "perform fastqc", "fastqc analysis", "analyze", "analysis on"]
        )

        return looks_like_browse and looks_like_results and not looks_like_execution

    @classmethod
    def required_parameters(cls) -> List[Dict[str, Any]]:
        return []

    @classmethod
    def create_workflow_plan(
        cls,
        inputs: List[DataInput],
        parameters: Dict[str, Any],
        session_context: Dict[str, Any],
    ) -> WorkflowPlan:
        uris = [inp.uri for inp in (inputs or []) if isinstance(inp.uri, str) and inp.uri.startswith("s3://")]
        prefix = None
        show = None

        for u in uris:
            if u.lower().endswith("results.json") or u.lower().endswith(".json"):
                show = u
            elif u.endswith("/"):
                prefix = u

        if not prefix and show:
            parts = show.rstrip("/").rsplit("/", 1)
            if len(parts) == 2:
                prefix = parts[0] + "/"

        if prefix and not prefix.endswith("/"):
            prefix = prefix + "/"

        op_params = {
            "prefix": prefix or (uris[0].rstrip("/") + "/" if uris else ""),
            "show": show,
            "recursive": True,
            "max_keys": 200,
        }

        operations = [
            OperationSpec(
                operation_name="browse_results",
                tool_name="s3_browse_results",
                parameters=op_params,
                expected_output_size_mb=1.0,
                parallelizable=False,
            )
        ]

        return WorkflowPlan(
            description="Browse results stored in S3",
            data_inputs=inputs,
            operations=operations,
            constraints=ConstraintSpec(),
            expected_compute_intensity="Low",
            session_id=session_context.get("session_id"),
        )


class SimpleOperationPlaybook(WorkflowPlaybook):
    """Simple single-operation playbook (fallback)."""
    
    workflow_type = "simple_operation"
    description = "Single bioinformatics operation"
    
    @classmethod
    def matches(cls, inputs: List[DataInput], command: str, session_context: Dict[str, Any]) -> bool:
        """Always matches as fallback."""
        return True
    
    @classmethod
    def detect_operation(cls, command: str, inputs: List[DataInput]) -> Tuple[str, Optional[str]]:
        """Detect operation type from command."""
        # IMPORTANT: avoid matching keywords inside file paths/URIs (e.g. "fastqc_new/")
        # which can contain substrings like "qc" that should NOT trigger QC/FastQC execution.
        command_lower = command.lower()
        scrubbed = re.sub(r"s3://[^\s]+", " ", command_lower)
        scrubbed = re.sub(r"(?:^|\s)(/[^\s]+)", " ", scrubbed)
        scrubbed = re.sub(r"\s+", " ", scrubbed).strip()
        
        # Operation keyword mapping
        # Use word-boundary matching for short keywords like "qc"
        for pattern, (op_name, tool_name) in SIMPLE_OPERATION_REGEX:
            if re.search(pattern, scrubbed):
                return op_name, tool_name
        
        # Default to generic operation
        return "bioinformatics_operation", None
    
    @classmethod
    def create_workflow_plan(
        cls,
        inputs: List[DataInput],
        parameters: Dict[str, Any],
        session_context: Dict[str, Any]
    ) -> WorkflowPlan:
        """Create simple single-operation workflow plan."""
        operation_name, tool_name = cls.detect_operation(
            parameters.get("command", ""),
            inputs
        )
        
        operations = [
            OperationSpec(
                operation_name=operation_name,
                tool_name=tool_name,
                parameters=parameters.get("operation_parameters", {}),
                expected_output_size_mb=sum(inp.size_bytes or 0 for inp in inputs) / (1024 * 1024),
                parallelizable=False
            )
        ]
        
        return WorkflowPlan(
            description=f"Simple {operation_name} operation",
            data_inputs=inputs,
            operations=operations,
            constraints=ConstraintSpec(),
            expected_compute_intensity="Low",
            session_id=session_context.get("session_id")
        )


# Playbook registry (ordered by specificity - most specific first)
PLAYBOOKS = [
    ResultsBrowsingPlaybook,
    RNASeqPlaybook,
    ScRNASeqPlaybook,
    WGSWESPlaybook,
    PhylogeneticsPlaybook,
    SimpleOperationPlaybook  # Fallback
]


def _is_fastq(uri: str) -> bool:
    """Check if URI is a FASTQ file."""
    return any(uri.lower().endswith(ext) for ext in ['.fastq', '.fq', '.fastq.gz', '.fq.gz'])


def _is_fasta(uri: str) -> bool:
    """Check if URI is a FASTA file."""
    return any(uri.lower().endswith(ext) for ext in ['.fasta', '.fa', '.fna', '.faa', '.fasta.gz', '.fa.gz'])


def _extract_file_uris(command: str, session_context: Dict[str, Any]) -> List[str]:
    """Extract file URIs from command and session context."""
    uris = []
    
    # Extract S3 URIs from command
    s3_pattern = r's3://[^\s]+'
    uris.extend(re.findall(s3_pattern, command))
    
    # Extract local paths from command (basic pattern)
    local_pattern = r'(?:^|\s)(/[^\s]+\.(?:fastq|fq|fasta|fa|h5ad|mtx|vcf|bam)(?:\.gz)?)'
    uris.extend(match.strip() for match in re.findall(local_pattern, command))
    
    # Check session context for uploaded files
    metadata = session_context.get("metadata", {})
    for uploaded_file in metadata.get("uploaded_files", []):
        if isinstance(uploaded_file, dict):
            file_uri = _get_file_uri_from_session(uploaded_file, session_context)
            if file_uri:
                uris.append(file_uri)
    
    return list(set(uris))  # Deduplicate


def _get_file_uri_from_session(file_info: Dict[str, Any], session_context: Dict[str, Any]) -> Optional[str]:
    """Get file URI from session file info."""
    if "s3_key" in file_info and "s3_bucket" in file_info:
        return f"s3://{file_info['s3_bucket']}/{file_info['s3_key']}"
    elif "s3_key" in file_info:
        # Try to get bucket from session metadata
        bucket = session_context.get("metadata", {}).get("s3_bucket")
        if bucket:
            return f"s3://{bucket}/{file_info['s3_key']}"
    return None


async def plan_workflow(
    command: str,
    session_context: Dict[str, Any],
    request_id: Optional[str] = None
) -> Dict[str, Any]:
    """
    Plan a workflow based on user command and session context.
    
    Args:
        command: User command (after intent classification)
        session_context: Session context with uploaded files, history, etc.
        request_id: Optional request ID for tracking
    
    Returns:
        Dict with one of:
        - {"status": "success", "workflow_plan": WorkflowPlan} - Plan created
        - {"status": "clarification_needed", "message": str, "missing_parameters": []} - Need more info
        - {"status": "infeasible", "error": str, "suggestions": []} - Cannot proceed
    """
    logger.info(f"[{request_id}] Planning workflow for command: {command[:100]}")
    
    # Extract file URIs from command and session context
    file_uris = _extract_file_uris(command, session_context)
    
    # Create DataInput objects
    inputs = []
    for uri in file_uris:
        inputs.append(DataInput(
            uri=uri,
            size_bytes=None,  # Size will be fetched by Infrastructure Agent
            description=f"Input file from {'session' if 's3://' in uri else 'command'}"
        ))
    
    # Select appropriate playbook
    selected_playbook = None
    for playbook_cls in PLAYBOOKS:
        if playbook_cls.matches(inputs, command, session_context):
            selected_playbook = playbook_cls
            logger.info(f"[{request_id}] Selected playbook: {playbook_cls.workflow_type}")
            break
    
    if not selected_playbook:
        # This should never happen since SimpleOperationPlaybook always matches
        return {
            "status": "infeasible",
            "error": "No suitable workflow found",
            "suggestions": ["Please provide more details about the analysis you want to perform"]
        }
    
    # Check for required parameters
    required_params = selected_playbook.required_parameters()
    missing_params = []
    extracted_params = {}
    
    for param_spec in required_params:
        param_name = param_spec["parameter"]
        
        # Try to extract parameter from command
        param_value = _extract_parameter(param_name, command, session_context)
        
        if param_value is None and param_spec.get("required", False):
            # Required parameter is missing
            missing_params.append(param_spec)
        elif param_value is None and "default" in param_spec:
            # Use default value
            extracted_params[param_name] = param_spec["default"]
        elif param_value is not None:
            extracted_params[param_name] = param_value
    
    # If missing required parameters, return clarification request
    if missing_params:
        # If the user provided an explicit multi-step plan (e.g. "Step 1:", "Step 2:"),
        # prefer producing a usable draft workflow with sensible defaults rather than
        # blocking on clarifications. This keeps E2E workflow planning deterministic.
        command_lower = (command or "").lower()
        has_explicit_steps = bool(re.search(r"\bstep\s+\d+\s*[:.)]", command_lower))
        if has_explicit_steps and selected_playbook is RNASeqPlaybook:
            extracted_params.setdefault("organism", "human")
            extracted_params.setdefault("reference_genome", "hg38")
            # Recompute missing params after defaults.
            missing_params = [
                p for p in missing_params
                if p.get("parameter") not in {"organism", "reference_genome"}
            ]

        if missing_params:
            return {
                "status": "clarification_needed",
                "message": _format_clarification_message(selected_playbook, missing_params, inputs),
                "missing_parameters": missing_params,
                "detected_workflow_type": selected_playbook.workflow_type,
                "detected_inputs": file_uris
            }
    
    # Create workflow plan
    try:
        extracted_params["command"] = command  # Pass original command for context
        workflow_plan = selected_playbook.create_workflow_plan(
            inputs=inputs,
            parameters=extracted_params,
            session_context=session_context
        )
        
        logger.info(
            f"[{request_id}] Created workflow plan: "
            f"{len(workflow_plan.operations)} operations, "
            f"{len(workflow_plan.data_inputs)} inputs"
        )
        
        return {
            "status": "success",
            "workflow_plan": workflow_plan,
            "workflow_type": selected_playbook.workflow_type
        }
        
    except Exception as e:
        logger.error(f"[{request_id}] Failed to create workflow plan: {e}", exc_info=True)
        return {
            "status": "infeasible",
            "error": f"Failed to create workflow plan: {str(e)}",
            "suggestions": ["Please check your input files and parameters"]
        }


def _extract_parameter(param_name: str, command: str, session_context: Dict[str, Any]) -> Optional[str]:
    """Extract parameter value from command or session context."""
    command_lower = command.lower()
    
    # Check for reference genome
    if param_name == "reference_genome":
        # Common patterns: "GRCh38", "GRCm38", "hg38", "mm10", "T2T-CHM13"
        # Try comprehensive pattern first (matches GRCh38, GRCm38, T2T-CHM13, etc.)
        genome_pattern = r'(?:genome|reference)[\s:]*([A-Za-z0-9][\w\-\.]+\d+)'
        match = re.search(genome_pattern, command, re.IGNORECASE)
        if match:
            return match.group(1)
        
        # Try direct mention patterns
        # Pattern 1: GRCh38, GRCm38, etc. (case-insensitive)
        grch_pattern = r'\b(GRC[hm]\d+(?:\.\w+)?)\b'
        match = re.search(grch_pattern, command, re.IGNORECASE)
        if match:
            return match.group(1)
        
        # Pattern 2: hg38, mm10, etc.
        direct_pattern = r'\b([a-z]{2}\d+)\b'
        match = re.search(direct_pattern, command_lower)
        if match:
            return match.group(1)
    
    # Check for organism
    if param_name == "organism":
        organisms = {
            "human": "human",
            "mouse": "mouse",
            "rat": "rat",
            "zebrafish": "zebrafish",
            "fly": "drosophila",
            "worm": "caenorhabditis",
            "yeast": "yeast"
        }
        for keyword, organism in organisms.items():
            if keyword in command_lower:
                return organism
    
    # Check for strandedness
    if param_name == "strandedness":
        if "unstranded" in command_lower or "non-stranded" in command_lower:
            return "unstranded"
        if "forward" in command_lower and "strand" in command_lower:
            return "forward"
        if "reverse" in command_lower and "strand" in command_lower:
            return "reverse"
    
    return None


def _format_clarification_message(
    playbook_cls: type,
    missing_params: List[Dict[str, Any]],
    inputs: List[DataInput]
) -> str:
    """Format clarification message for missing parameters."""
    workflow_name = playbook_cls.description
    
    # Build message
    msg_parts = [
        f"I detected {workflow_name} analysis with {len(inputs)} input file(s).",
        "To proceed, I need the following information:"
    ]
    
    for i, param in enumerate(missing_params, 1):
        msg_parts.append(f"{i}. **{param['parameter']}**: {param['description']}")
    
    msg_parts.append("")
    msg_parts.append("Please provide these parameters to continue.")
    
    return "\n".join(msg_parts)


# For backward compatibility and testing
def get_playbook_for_workflow_type(workflow_type: str) -> Optional[type]:
    """Get playbook class by workflow type."""
    for playbook_cls in PLAYBOOKS:
        if playbook_cls.workflow_type == workflow_type:
            return playbook_cls
    return None
