# /backend/agent_tools.py
# All @tool decorated functions for the Helix.AI agent

import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional

from langchain_core.tools import tool

logger = logging.getLogger(__name__)

# Ensure tools directory is in Python path for imports
# This allows importing modules from the tools/ directory
_TOOLS_DIR = Path(__file__).resolve().parent.parent / "tools"
if str(_TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(_TOOLS_DIR))


@tool
def toolbox_inventory() -> Dict:
    """List all tools Helix.AI has access to (registered MCP tools, discovered @tool functions, and local/EC2 CLI tools)."""
    from backend.tool_inventory import build_toolbox_inventory, format_toolbox_inventory_markdown
    inv = build_toolbox_inventory()
    return {
        "text": format_toolbox_inventory_markdown(inv),
        "input": {},
        "output": inv,
        "plot": {},
    }


@tool
def sequence_alignment(sequences: str) -> Dict:
    """Performs a sequence alignment on a given set of sequences.
    
    Use this tool ONLY when the user explicitly asks to align sequences without
    building a phylogenetic tree. If the user asks to visualize or create a
    phylogenetic tree, use phylogenetic_tree instead (which handles alignment internally).
    """
    
    # Import the alignment function
    from alignment import run_alignment_tool
    
    result = run_alignment_tool(sequences)
    
    return {
        "text": result.get("text", "Sequences aligned successfully."),
        "input": sequences,
        "output": result.get("alignment", []),
        "statistics": result.get("statistics", {}),
        "plot": {
            "data": [{"x": [1, 2, 3], "y": [3, 3, 3], "type": "bar"}],
            "layout": {"title": "Alignment Visualization"},
        },
    }


@tool
def mutate_sequence(sequence: str, num_variants: int = 96) -> Dict:
    """Mutates a given sequence and returns the specified number of variants. Use this when you need to CREATE new sequence variants from an input sequence."""
    
    # Import the mutation function
    from mutations import run_mutation_raw
    
    result = run_mutation_raw(sequence, num_variants)
    
    return {
        "text": result.get("text", "Sequence mutated successfully."),
        "input": {
            "sequence": sequence,
            "variants": num_variants,
        },
        "output": result.get("statistics", {}),
        "plot": result.get("plot", {
            "data": [{"x": [1, 2, 3], "y": [1, 1, 1], "type": "bar"}],
            "layout": {"title": "Mutation Visualization"},
        }),
    }


@tool
def dna_vendor_research(command: str, sequence_length: int = None, quantity: str = "standard") -> Dict:
    """Research DNA synthesis vendors and testing options for experimental validation. Use this when the user wants to ORDER sequences, find VENDORS, or research TESTING options. Keywords: order, vendor, synthesis, test, assay, expression, function, binding."""
    
    # Import the research function
    from dna_vendor_research import run_dna_vendor_research_raw
    
    result = run_dna_vendor_research_raw(command, sequence_length, quantity)
    
    return {
        "text": f"DNA vendor research completed: {result.get('message', 'Research done')}",
        "input": {
            "command": command,
            "sequence_length": sequence_length,
            "quantity": quantity
        },
        "output": result,
        "plot": {
            "data": [{"x": ["Vendors", "Testing Options"], "y": [result.get('total_vendors', 0), result.get('total_testing_options', 0)], "type": "bar"}],
            "layout": {"title": "DNA Vendor Research Results"},
        },
    }


@tool
def phylogenetic_tree(aligned_sequences: str) -> Dict:
    """Create phylogenetic tree visualization from sequences.
    
    This tool accepts either aligned or unaligned sequences in FASTA format.
    If sequences are unaligned, it will automatically align them first, then
    build the phylogenetic tree and create visualizations.
    
    Use this tool when the user asks to:
    - Visualize a phylogenetic tree
    - Create a phylogenetic tree
    - Build an evolutionary tree
    - Show phylogenetic relationships
    
    Do NOT use sequence_alignment for this - phylogenetic_tree handles alignment internally.
    """
    
    print(f"🔧 Agent phylogenetic_tree called with aligned_sequences: '{aligned_sequences}'")
    
    # Import the phylogenetic tree function
    from phylogenetic_tree import run_phylogenetic_tree
    
    result = run_phylogenetic_tree(aligned_sequences)
    
    # Log result summary without full data to avoid console spam
    if isinstance(result, dict):
        print(f"🔧 Agent phylogenetic_tree result: status={result.get('text', 'N/A')[:50]}...")
        print(f"🔧 Agent phylogenetic_tree has tree_newick: {bool(result.get('tree_newick'))}")
    else:
        print(f"🔧 Agent phylogenetic_tree result: {type(result).__name__}")
    
    # Return all fields from result, especially tree_newick and ete_visualization for frontend
    if isinstance(result, dict):
        return {
            "text": result.get("text", "Phylogenetic tree created successfully."),
            "input": aligned_sequences,
            "output": result.get("statistics", {}),
            "plot": result.get("plot", {
                "data": [{"x": [1, 2, 3], "y": [1, 1, 1], "type": "bar"}],
                "layout": {"title": "Phylogenetic Tree"},
            }),
            # Pass through tree data for frontend visualization
            "tree_newick": result.get("tree_newick"),
            "ete_visualization": result.get("ete_visualization"),
            "aligned_sequences": result.get("aligned_sequences"),
            "clustering_result": result.get("clustering_result"),
            "clustered_visualization": result.get("clustered_visualization"),
            "statistics": result.get("statistics", {}),
        }
    else:
        return {
            "text": "Phylogenetic tree created successfully.",
            "input": aligned_sequences,
            "output": {},
            "plot": {
                "data": [{"x": [1, 2, 3], "y": [1, 1, 1], "type": "bar"}],
                "layout": {"title": "Phylogenetic Tree"},
            },
        }


@tool
def sequence_selection(aligned_sequences: str, selection_type: str = "random", num_sequences: int = 1) -> Dict:
    """Select sequences from aligned sequences based on various criteria (random, best_conservation, lowest_gaps, highest_gc, longest, shortest)."""
    
    print(f"🔧 Agent sequence_selection called with:")
    print(f"🔧 aligned_sequences: '{aligned_sequences}'")
    print(f"🔧 selection_type: '{selection_type}'")
    print(f"🔧 num_sequences: {num_sequences}")
    
    # Import the sequence selection function
    from sequence_selection import run_sequence_selection_raw
    
    result = run_sequence_selection_raw(aligned_sequences, selection_type, num_sequences)
    
    # Log result summary without full data to avoid console spam
    if isinstance(result, dict):
        print(f"🔧 Agent sequence_selection result: status={result.get('text', 'N/A')[:50]}...")
    else:
        print(f"🔧 Agent sequence_selection result: {type(result).__name__}")
    
    return {
        "text": result.get("text", "Sequence selection completed successfully."),
        "input": {
            "aligned_sequences": aligned_sequences,
            "selection_type": selection_type,
            "num_sequences": num_sequences
        },
        "output": result.get("selected_sequences", []),
        "selection_criteria": result.get("selection_criteria", {}),
    }


@tool
def synthesis_submission(sequences: str, vendor_preference: Optional[str] = None, quantity: str = "standard", delivery_time: str = "standard") -> Dict:
    """Submit sequences for DNA synthesis and get pricing quote. Quantity options: standard, large, custom. Delivery options: rush, standard, economy."""
    
    # Import the synthesis submission function
    from synthesis_submission import run_synthesis_submission_raw
    
    result = run_synthesis_submission_raw(sequences, vendor_preference, quantity, delivery_time)
    
    return {
        "text": result.get("text", "Synthesis submission completed successfully."),
        "input": {
            "sequences": sequences,
            "vendor_preference": vendor_preference,
            "quantity": quantity,
            "delivery_time": delivery_time
        },
        "output": result.get("quote", {}),
        "validation_results": result.get("validation_results", {}),
    }


@tool
def create_session() -> Dict:
    """Create a new session for tracking user interactions. Use this when the user wants to start a new session or create a new experiment."""
    
    import uuid
    
    session_id = str(uuid.uuid4())
    
    return {
        "text": f"Session created successfully with ID: {session_id}",
        "input": {},
        "output": {"session_id": session_id},
        "plot": {}
    }


@tool
def plasmid_visualization(vector_name: str = None, cloning_sites: str = "", insert_sequence: str = "", 
                          full_plasmid_sequence: str = None, insert_position: int = None) -> Dict:
    """Generate plasmid visualization data.
    
    Supports two modes:
    1. Full plasmid: When full_plasmid_sequence is provided, visualize the complete plasmid
    2. Insert mode: When vector_name and insert_sequence are provided, create plasmid with insert
    """
    
    # Import the plasmid visualization function
    from plasmid_visualizer import run_plasmid_visualization_raw
    
    result = run_plasmid_visualization_raw(
        vector_name=vector_name,
        cloning_sites=cloning_sites,
        insert_sequence=insert_sequence,
        full_plasmid_sequence=full_plasmid_sequence,
        insert_position=insert_position
    )
    
    return {
        "text": result.get("text", "Plasmid visualization created successfully."),
        "input": {
            "vector_name": vector_name,
            "cloning_sites": cloning_sites,
            "insert_sequence": insert_sequence,
            "full_plasmid_sequence": full_plasmid_sequence,
            "insert_position": insert_position
        },
        "output": result.get("plasmid_data", {}),
        "visualization_type": result.get("visualization_type", "circular_plasmid"),
    }

@tool
def plasmid_for_representatives(representatives: List[str], aligned_sequences: str, vector_name: str = "pUC19", cloning_sites: str = "EcoRI, BamHI, HindIII") -> Dict:
    """Create plasmid visualizations for representative sequences from clustering analysis."""
    
    # Import the plasmid visualization function
    from plasmid_visualizer import create_plasmid_for_representatives
    
    result = create_plasmid_for_representatives(representatives, aligned_sequences, vector_name, cloning_sites)
    
    return {
        "text": result.get("text", "Plasmid visualizations for representatives created successfully."),
        "plasmid_results": result.get("plasmid_results", []),
        "total_representatives": result.get("total_representatives", 0),
        "vector_name": result.get("vector_name", vector_name),
        "cloning_sites": result.get("cloning_sites", cloning_sites),
    }


@tool
def single_cell_analysis(
    data_file: Optional[str] = None,
    data_format: str = "10x",
    steps: str = "all",
    question: Optional[str] = None
) -> Dict:
    """Answer questions about single-cell RNA-seq (scRNA-seq) sequencing and data analysis, OR perform actual analysis on data.
    
    USE THIS TOOL for ANY question about:
    - Single-cell sequencing methods and workflows
    - scRNA-seq data analysis
    - Seurat analysis pipelines
    - Cell type identification and annotation
    - Marker gene discovery
    - Differential expression in single-cell data
    - Pathway analysis for single-cell data
    - Batch correction methods
    - Single-cell data formats (10x, H5, Seurat objects)
    
    This tool provides comprehensive information about single-cell analysis capabilities and can also perform actual analysis when data is provided.
    
    Args:
        data_file: Path to input data file (10x format directory, H5, CSV, or Seurat RDS). Leave empty for informational questions.
        data_format: Format of input data ('10x', 'h5', 'csv', 'seurat'). Default: '10x'. Only needed if data_file is provided.
        steps: Comma-separated list of analysis steps: 'preprocessing', 'markers', 'differential', 'pathways', 'annotation', 'batch_correction', or 'all'. Default: 'all'. Only needed if data_file is provided.
        question: The user's question about single-cell analysis. Use this parameter to capture what the user is asking about.
    
    Returns:
        Dictionary containing detailed information about single-cell analysis capabilities, or analysis results if data_file is provided.
    """
    
    # Import the single-cell analysis function
    from single_cell_analysis import analyze_single_cell_data
    
    # If no data file is provided, return informational response about single-cell analysis
    if not data_file:
        # Build a comprehensive response about single-cell analysis
        intro = "I can help you with single-cell RNA-seq (scRNA-seq) analysis! "
        if question:
            intro = f"Regarding your question about single-cell sequencing and data analysis: "
        
        info_text = f"""{intro}Here's what I can help you with:

**Single-Cell Analysis Capabilities:**
- **Preprocessing**: Quality control, normalization, and filtering of single-cell data
- **Marker Gene Identification**: Find genes that distinguish cell clusters using tools like Seurat
- **Differential Expression Analysis**: Identify genes differentially expressed between cell types or conditions
- **Pathway Enrichment Analysis**: Discover biological pathways enriched in specific cell types
- **Cell-Type Annotation**: Automatically annotate cell types using reference datasets
- **Batch Correction**: Correct for batch effects in multi-sample datasets

**Supported Data Formats:**
- 10x Genomics format (standard output from Cell Ranger)
- H5/H5AD files (AnnData format)
- CSV files (count matrices)
- Seurat RDS objects

**Analysis Workflows:**
I use the scPipeline R package (built on Seurat) to perform comprehensive single-cell analysis. You can run individual steps or a complete end-to-end analysis.

**What I can do:**
- Answer questions about single-cell sequencing methods, workflows, and interpretation
- Explain Seurat analysis pipelines and best practices
- Help with cell type identification and marker gene discovery
- Perform actual analysis when you provide single-cell data
- Guide you through differential expression and pathway analysis

**To get started with analysis:**
- Upload your single-cell data (10x format directory, H5 file, or CSV)
- Ask me to analyze it, or specify which steps you want (e.g., "find marker genes", "perform differential expression analysis")

Feel free to ask me any specific questions about single-cell sequencing and data analysis!"""
        
        return {
            "text": info_text,
            "input": {"question": question or "general inquiry about single-cell analysis"},
            "output": {
                "capabilities": [
                    "Preprocessing and quality control",
                    "Marker gene identification",
                    "Differential expression analysis",
                    "Pathway enrichment analysis",
                    "Cell-type annotation",
                    "Batch correction"
                ],
                "supported_formats": ["10x", "h5", "h5ad", "csv", "seurat"],
                "analysis_package": "scPipeline (R) / Seurat",
                "can_perform_analysis": True,
                "can_answer_questions": True
            },
            "plot": {}
        }
    
    # If data file is provided, perform actual analysis
    steps_list = [s.strip() for s in steps.split(",")] if isinstance(steps, str) else steps
    
    result = analyze_single_cell_data(
        data_file=data_file,
        data_format=data_format,
        steps=steps_list
    )
    
    return {
        "text": result.get("text", f"Single-cell analysis completed. Steps: {', '.join(steps_list)}"),
        "input": {
            "data_file": data_file,
            "data_format": data_format,
            "steps": steps_list
        },
        "output": result,
        "plot": result.get("plot", {})
    }


@tool
def fetch_ncbi_sequence(
    accession: str,
    database: str = "nucleotide"
) -> Dict:
    """Fetch DNA/RNA/protein sequence from NCBI databases by accession number.
    
    USE THIS TOOL when users ask to:
    - Fetch sequences from NCBI
    - Get sequences by accession number
    - Retrieve reference sequences
    - Download GenBank/RefSeq sequences
    
    Examples:
    - "fetch sequence NC_000001.11"
    - "get protein sequence NP_000483.1"
    - "download BRCA1 gene sequence"
    
    Args:
        accession: NCBI accession number (e.g., "NC_000001.11", "NM_000492.3", "NP_000483.1")
        database: Database to search ('nucleotide' for DNA/RNA, 'protein' for proteins). Default: 'nucleotide'
    
    Returns:
        Dictionary containing sequence data and metadata.
    """
    
    # Import the NCBI tool function
    from ncbi_tools import fetch_sequence_from_ncbi
    
    result = fetch_sequence_from_ncbi(accession, database)
    
    if result.get("status") == "success":
        # For very large sequences (like full chromosomes), we don't need the full sequence
        # in the LLM call - just metadata and a small preview
        full_sequence = result.get("sequence", "")
        sequence_length = len(full_sequence)
        
        # For LLM: only include metadata and a tiny preview (20 bases) for validation
        # The full sequence is available in the response for frontend display
        if sequence_length > 20:
            # For large sequences, just show first 20 bases + metadata
            sequence_preview = full_sequence[:20] + f"... (full length: {sequence_length:,} bp)"
        else:
            sequence_preview = full_sequence
        
        return {
            "text": f"Fetched sequence {accession} ({result.get('length', 0):,} bp): {result.get('description', '')}",
            "input": {"accession": accession, "database": database},
            "output": {
                "sequence": sequence_preview,  # Minimal preview for LLM (20 bases + metadata)
                # Don't include full_sequence here - it causes LLM timeouts and huge responses
                # Full sequence is stored in tool execution result, not in LLM context
                "description": result.get("description"),
                "length": result.get("length"),
                "database": result.get("database"),
                "note": f"Full sequence ({sequence_length:,} bp) fetched and stored, preview shown above" if sequence_length > 20 else None
            },
            "plot": {}
        }
    else:
        return {
            "text": f"Error fetching sequence {accession}: {result.get('error', 'Unknown error')}",
            "input": {"accession": accession, "database": database},
            "output": result,
            "plot": {}
        }


@tool
def query_uniprot(
    query: str,
    format: str = "fasta",
    limit: int = 10
) -> Dict:
    """Query UniProt protein database for sequences and metadata."""
    
    from uniprot_tools import query_uniprot as query_uniprot_api
    
    result = query_uniprot_api(query, format=format, limit=limit)
    
    if result.get("status") == "success":
        summary = f"{result.get('count', 0)} result(s)" if "count" in result else "Query executed"
        return {
            "text": f"UniProt query for '{query}' completed: {summary}",
            "input": {"query": query, "format": format, "limit": limit},
            "output": result,
            "plot": {}
        }
    else:
        return {
            "text": f"Error querying UniProt for '{query}': {result.get('error', 'Unknown error')}",
            "input": {"query": query, "format": format, "limit": limit},
            "output": result,
            "plot": {}
        }


@tool
def lookup_go_term(go_id: str) -> Dict:
    """Lookup Gene Ontology (GO) term details by ID."""
    
    from go_tools import lookup_go_term as lookup_go
    
    result = lookup_go(go_id)
    
    if result.get("status") == "success":
        return {
            "text": f"GO term {go_id}: {result.get('name', '')}",
            "input": {"go_id": go_id},
            "output": result,
            "plot": {}
        }
    else:
        return {
            "text": f"Error looking up GO term {go_id}: {result.get('error', 'Unknown error')}",
            "input": {"go_id": go_id},
            "output": result,
            "plot": {}
        }


@tool
def bulk_rnaseq_analysis(
    count_matrix: str,
    sample_metadata: str,
    design_formula: str = "~condition",
    alpha: float = 0.05
) -> Dict:
    """Run bulk RNA-seq differential expression analysis using DESeq2."""
    
    from bulk_rnaseq import run_deseq2_analysis
    
    if not count_matrix or not sample_metadata:
        return {
            "text": "Bulk RNA-seq requires count matrix and sample metadata CSV paths.",
            "input": {
                "count_matrix": count_matrix,
                "sample_metadata": sample_metadata,
                "design_formula": design_formula,
                "alpha": alpha
            },
            "output": {"status": "error", "message": "Missing input files"},
            "plot": {}
        }
    
    result = run_deseq2_analysis(
        count_matrix=count_matrix,
        sample_metadata=sample_metadata,
        design_formula=design_formula,
        alpha=alpha
    )
    
    status = result.get("status", "success")
    summary = result.get("summary", {})
    text_summary = f"DESeq2 complete: {summary.get('significant_genes', 'n/a')} significant genes" if summary else result.get("message", "DESeq2 analysis completed")
    
    return {
        "text": text_summary if status == "success" else f"DESeq2 error: {result.get('message', 'Unknown error')}",
        "input": {
            "count_matrix": count_matrix,
            "sample_metadata": sample_metadata,
            "design_formula": design_formula,
            "alpha": alpha
        },
        "output": result,
        "plot": {}
    }


@tool
def fastqc_quality_analysis(
    input_r1: str,
    input_r2: str,
    output: Optional[str] = None
) -> Dict:
    """
    Run FastQC quality control analysis on paired-end FASTQ files stored in S3 using AWS EMR.
    
    This is an asynchronous long-running job that typically takes 10-30 minutes.
    The tool returns immediately with a job_id that can be used to track progress.
    
    Args:
        input_r1: S3 URI for the forward/read 1 FASTQ file (e.g., "s3://my-bucket/data/sample_R1.fastq")
        input_r2: S3 URI for the reverse/read 2 FASTQ file (e.g., "s3://my-bucket/data/sample_R2.fastq")
        output: Optional S3 URI for output directory (defaults to same directory as input files with /fastqc-results suffix)
    
    Returns:
        Dict with job_id and initial status. Format:
        {
            "type": "job",
            "job_id": "uuid-string",
            "status": "submitted",
            "message": "FastQC job submitted. Processing will take 10-30 minutes.",
            "input_r1": "s3://...",
            "input_r2": "s3://...",
            "output": "s3://..." (if provided)
        }
    """
    from backend.job_manager import get_job_manager
    
    job_manager = get_job_manager()
    
    try:
        job_id = job_manager.submit_fastqc_job(
            r1_path=input_r1,
            r2_path=input_r2,
            output_path=output
        )
        
        return {
            "type": "job",
            "job_id": job_id,
            "status": "submitted",
            "message": "FastQC job submitted. Processing will take 10-30 minutes.",
            "input_r1": input_r1,
            "input_r2": input_r2,
            "output": output
        }
    except Exception as e:
        logger.error(f"FastQC job submission failed: {e}")
        return {
            "type": "error",
            "status": "error",
            "message": f"Failed to submit FastQC job: {str(e)}",
            "error": str(e)
        }

