# /backend/agent.py

import asyncio
import os
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional

from langgraph.prebuilt import create_react_agent
from langgraph.checkpoint.memory import MemorySaver

from langchain_core.tools import tool

from pydantic import BaseModel, Field

# load the environment
from dotenv import load_dotenv

logger = logging.getLogger(__name__)

try:
    # Best-effort: in some sandbox/CI environments the .env file may be unreadable
    # (e.g., ignored files are blocked). This should never be fatal.
    load_dotenv()
except PermissionError as e:
    logger.warning(f"Could not read .env (permission denied); continuing without it: {e}")
except Exception as e:
    logger.warning(f"Could not load .env; continuing without it: {e}")

from langchain.globals import set_verbose, set_debug
from backend.prompts.templates import build_react_prompt
from backend.context_builder import build_context_snippet

# Disable verbose/debug logging by default to prevent large tool results (like sequences) from being printed
# Can be enabled via LANGCHAIN_VERBOSE=true or LANGCHAIN_DEBUG=true environment variables for debugging
verbose_enabled = os.getenv("LANGCHAIN_VERBOSE", "false").lower() == "true"
debug_enabled = os.getenv("LANGCHAIN_DEBUG", "false").lower() == "true"
set_verbose(verbose_enabled)
set_debug(debug_enabled)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
AGENT_PROMPT_PATH = PROJECT_ROOT / "agent.md"

try:
    BIOAGENT_SYSTEM_PROMPT = AGENT_PROMPT_PATH.read_text()
except Exception:
    BIOAGENT_SYSTEM_PROMPT = (
        "You are BioAgent, an autonomous bioinformatics assistant. "
        "Classify prompts, plan, execute with real tools, and return structured JSON "
        "for browser rendering. Decline non-bioinformatics or infeasible requests."
    )

class OutputFormatter(BaseModel):
    """Always use this tool to structure your response to the user."""

    input: str = Field(description="The input of the tool/function call")
    output: str = Field(description="The output of the tool/function call")
    plot: str = Field(description="An optional plot visualizing the output")


@tool
def toolbox_inventory() -> Dict:
    """List all tools Helix.AI has access to (registered MCP tools, discovered @tool functions, and local/EC2 CLI tools)."""
    from tool_inventory import build_toolbox_inventory, format_toolbox_inventory_markdown
    inv = build_toolbox_inventory()
    return {
        "text": format_toolbox_inventory_markdown(inv),
        "input": {},
        "output": inv,
        "plot": {},
    }


@tool
def sequence_alignment(sequences: str) -> Dict:
    """Performs a sequence alignment on a given set of sequences."""
    
    # Import the alignment function
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
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
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
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
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
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
    """Create phylogenetic tree visualization from aligned sequences."""
    
    print(f"🔧 Agent phylogenetic_tree called with aligned_sequences: '{aligned_sequences}'")
    
    # Import the phylogenetic tree function
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
    from phylogenetic_tree import run_phylogenetic_tree
    
    result = run_phylogenetic_tree(aligned_sequences)
    
    # Log result summary without full data to avoid console spam
    if isinstance(result, dict):
        print(f"🔧 Agent phylogenetic_tree result: status={result.get('text', 'N/A')[:50]}...")
    else:
        print(f"🔧 Agent phylogenetic_tree result: {type(result).__name__}")
    
    return {
        "text": result.get("text", "Phylogenetic tree created successfully."),
        "input": aligned_sequences,
        "output": result.get("statistics", {}),
        "plot": result.get("plot", {
            "data": [{"x": [1, 2, 3], "y": [1, 1, 1], "type": "bar"}],
            "layout": {"title": "Phylogenetic Tree"},
        }),
    }


@tool
def sequence_selection(aligned_sequences: str, selection_type: str = "random", num_sequences: int = 1) -> Dict:
    """Select sequences from aligned sequences based on various criteria (random, best_conservation, lowest_gaps, highest_gc, longest, shortest)."""
    
    print(f"🔧 Agent sequence_selection called with:")
    print(f"🔧 aligned_sequences: '{aligned_sequences}'")
    print(f"🔧 selection_type: '{selection_type}'")
    print(f"🔧 num_sequences: {num_sequences}")
    
    # Import the sequence selection function
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
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
def synthesis_submission(sequences: str, vendor_preference: str = None, quantity: str = "standard", delivery_time: str = "standard") -> Dict:
    """Submit sequences for DNA synthesis and get pricing quote. Quantity options: standard, large, custom. Delivery options: rush, standard, economy."""
    
    # Import the synthesis submission function
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
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
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
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
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
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
    data_file: str = None,
    data_format: str = "10x",
    steps: str = "all",
    question: str = None
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
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
    from single_cell_analysis import analyze_single_cell_data, process_single_cell_command
    
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
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
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
    
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
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
    
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
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
    
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
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


_llm = None
_agent = None
_memory = None


def _get_llm():
    """
    Lazily initialize and return the LLM for BioAgent.

    Important:
    - In tests/CI we often run with `HELIX_MOCK_MODE=1` (no LLM, no network).
    - Imports for OpenAI/DeepSeek can trigger SSL/cert loading that may be blocked
      in sandboxed environments; keep them local to this function.
    """
    if os.getenv("HELIX_MOCK_MODE") == "1":
        raise RuntimeError("LLM is disabled in HELIX_MOCK_MODE")

    openai_key = os.getenv("OPENAI_API_KEY", "").strip()
    deepseek_key = os.getenv("DEEPSEEK_API_KEY", "").strip()

    openai_enabled = openai_key and openai_key not in ["", "disabled", "your_openai_api_key_here", "none"]
    deepseek_enabled = deepseek_key and deepseek_key not in ["", "disabled", "your_deepseek_api_key_here", "none"]

    if openai_enabled:
        from langchain.chat_models import init_chat_model
        return init_chat_model("openai:gpt-4o", temperature=0)

    if deepseek_enabled:
        from langchain_deepseek import ChatDeepSeek
        return ChatDeepSeek(
            model="deepseek-chat",
            temperature=0,
            max_tokens=None,
            timeout=60.0,
            max_retries=1,
        )

    raise ValueError(
        "No API keys found. Please set either OPENAI_API_KEY or DEEPSEEK_API_KEY "
        "in your environment variables."
    )


def _get_agent():
    """Lazily create and return the LangGraph agent instance."""
    global _llm, _agent, _memory
    if _agent is not None:
        return _agent

    _llm = _get_llm()
    _memory = MemorySaver()
    _agent = create_react_agent(
        model=_llm,
        tools=[
            toolbox_inventory,
            sequence_alignment,
            mutate_sequence,
            dna_vendor_research,
            phylogenetic_tree,
            sequence_selection,
            synthesis_submission,
            plasmid_visualization,
            plasmid_for_representatives,
            single_cell_analysis,
            fetch_ncbi_sequence,
            query_uniprot,
            lookup_go_term,
            bulk_rnaseq_analysis,
            fastqc_quality_analysis,
        ],
        checkpointer=_memory,
    )
    return _agent

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
    from job_manager import get_job_manager
    
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


agent = None


async def handle_command(command: str, session_id: str = "default", session_context: Dict = None):
    """
    Use agent for tool mapping only - identify which tool to use and parameters.
    Tool execution is handled by the router/tool path in main_with_mcp.py.
    """
    global agent
    print(f"[handle_command] command: {command}")
    print(f"[handle_command] session_id: {session_id}")
    print(f"[handle_command] Agent mode: TOOL MAPPING ONLY (no execution)")

    # Deterministic short-circuit: capability/toolbox questions should never go to tool-generator-agent.
    # They should return the actual toolbox inventory.
    cmd_lower = (command or "").lower()
    if any(
        phrase in cmd_lower
        for phrase in [
            "what can you do",
            "what tools do you have",
            "what tools are available",
            "list tools",
            "show tools",
            "your toolbox",
            "toolbox",
            "capabilities",
        ]
    ):
        print("[handle_command] Matched capability/toolbox query -> toolbox_inventory")
        tool_mapping = {"tool_name": "toolbox_inventory", "parameters": {}}
        return {
            "tool_mapping": tool_mapping,
            "tool_name": "toolbox_inventory",
            "parameters": {},
            "status": "tool_mapped",
            "message": "Capability/toolbox query: returning toolbox inventory.",
        }

    # Note: Removed _maybe_handle_vendor_request pre-filter to let LLM handle all tool selection
    # The LLM agent is capable of correctly identifying tools based on context,
    # including distinguishing between vendor research and file paths containing keywords like "test"

    # Build session-aware context from mutated/aligned/selected sequences and uploaded files
    context_snippet = ""
    if session_context:
        context_snippet = build_context_snippet(session_context)
    
    # Prepend context to user command if available
    user_content = command
    if context_snippet:
        user_content = context_snippet + "\n\nUser request: " + command
        print(f"[handle_command] Prepend context snippet ({len(context_snippet)} chars)")
    
    # Log total message size to help diagnose timeout issues
    total_message_size = len(user_content) + len(BIOAGENT_SYSTEM_PROMPT)
    print(f"[handle_command] Total message size to LLM: {total_message_size:,} chars (~{total_message_size//4:,} tokens estimated)")

    # For other commands, use the agent with session-specific config
    input_message = {
        "role": "user",
        "content": user_content,
    }
    system_message = {
        "role": "system",
        "content": BIOAGENT_SYSTEM_PROMPT,
    }

    # Use a temporary session for tool mapping to avoid polluting the main session memory
    # This prevents issues where incomplete tool calls (AIMessage without ToolMessage) 
    # cause validation errors on subsequent commands
    import uuid
    temp_session_id = f"{session_id}_mapping_{uuid.uuid4().hex[:8]}"
    temp_session_config = {"configurable": {"thread_id": temp_session_id}}
    
    # In mock mode, do not run the real agent (it requires network + API keys).
    # Tests can monkeypatch `agent` with a dummy implementation to validate payloads.
    if os.getenv("HELIX_MOCK_MODE") == "1":
        if agent is None:
            return {
                "tool_mapping": {"tool_name": "toolbox_inventory", "parameters": {}},
                "tool_name": "toolbox_inventory",
                "parameters": {},
                "status": "tool_mapped",
                "message": "HELIX_MOCK_MODE=1: agent disabled; returning toolbox_inventory mapping.",
            }
        if hasattr(agent, "invoke"):
            # Allow tests to validate messages without requiring agent streaming.
            return agent.invoke({"messages": [system_message, input_message]}, config=None)

    # Ensure real agent exists for non-mock mode.
    if agent is None:
        agent = _get_agent()

    try:
        # Agent should only do tool mapping, not execution
        # Use streaming with "updates" mode to capture tool calls as they're made
        # This allows us to extract tool name and parameters before execution completes
        tool_mapping = None
        found_tool = False
        
        def capture_tool_call_sync():
            """Synchronous function to stream agent and capture first tool call."""
            nonlocal tool_mapping, found_tool
            try:
                # Stream with "updates" mode to see tool calls as they happen
                # Use temp session to avoid polluting main session memory
                for event in agent.stream({"messages": [system_message, input_message]}, temp_session_config, stream_mode="updates"):
                    # Event is a dict with node names as keys
                    if isinstance(event, dict):
                        for node_name, node_output in event.items():
                            # Look for "tools" node which executes tools
                            if node_name == "tools" or "tools" in node_name.lower():
                                # This node contains tool calls
                                if isinstance(node_output, dict) and "messages" in node_output:
                                    messages = node_output["messages"]
                                    for msg in messages:
                                        # Tool call messages have the tool name
                                        if hasattr(msg, 'name'):
                                            tool_name = msg.name
                                            # Get tool input/arguments
                                            tool_args = {}
                                            if hasattr(msg, 'content'):
                                                content = msg.content
                                                if isinstance(content, str):
                                                    import json
                                                    try:
                                                        tool_args = json.loads(content)
                                                    except:
                                                        # Try to extract from string
                                                        pass
                                            elif hasattr(msg, 'input'):
                                                tool_args = msg.input
                                            
                                            tool_mapping = {
                                                "tool_name": tool_name,
                                                "parameters": tool_args if isinstance(tool_args, dict) else {}
                                            }
                                            print(f"[handle_command] Agent mapped tool from stream: {tool_name} with params: {tool_args}")
                                            found_tool = True
                                            return tool_mapping
                            
                            # Also check for AIMessage with tool_calls in any node
                            if isinstance(node_output, dict) and "messages" in node_output:
                                messages = node_output["messages"]
                                for msg in messages:
                                    # Check for AIMessage with tool_calls
                                    if hasattr(msg, 'tool_calls') and msg.tool_calls:
                                        for tool_call in msg.tool_calls:
                                            # Handle different tool_call formats
                                            tool_name = None
                                            tool_args = {}
                                            
                                            if hasattr(tool_call, 'name'):
                                                tool_name = tool_call.name
                                            elif isinstance(tool_call, dict):
                                                tool_name = tool_call.get("name")
                                            
                                            if hasattr(tool_call, 'args'):
                                                tool_args = tool_call.args
                                            elif isinstance(tool_call, dict):
                                                tool_args = tool_call.get("args", {})
                                            
                                            if tool_name:
                                                tool_mapping = {
                                                    "tool_name": tool_name,
                                                    "parameters": tool_args if isinstance(tool_args, dict) else {}
                                                }
                                                print(f"[handle_command] Agent mapped tool from AIMessage: {tool_name} with params: {tool_args}")
                                                found_tool = True
                                                return tool_mapping
            except Exception as e:
                print(f"[handle_command] Error during tool call capture: {e}")
                import traceback
                traceback.print_exc()
            return tool_mapping
        
        # Try to get tool mapping with a shorter timeout (just for LLM decision + first tool call)
        try:
            tool_mapping = await asyncio.wait_for(
                asyncio.to_thread(capture_tool_call_sync),
                timeout=30.0  # Shorter timeout for tool mapping only
            )
        except asyncio.TimeoutError:
            print(f"⚠️  Agent tool mapping timed out after 30s")
            # Try to extract from temp session state as fallback
            extracted = _extract_tool_call_from_state(memory, temp_session_config)
            if extracted:
                tool_mapping = extracted
        
        if tool_mapping and tool_mapping.get("tool_name"):
            # Return tool mapping for router to execute
            print(f"[handle_command] Returning tool mapping (no execution): {tool_mapping['tool_name']}")
            return {
                "tool_mapping": tool_mapping,
                "tool_name": tool_mapping["tool_name"],
                "parameters": tool_mapping["parameters"],
                "status": "tool_mapped",
                "message": f"Agent identified tool: {tool_mapping['tool_name']}. Execution will be handled by router."
            }
        
        # Fallback: if streaming didn't capture tool call, use regular invoke
        # This means the agent might have completed without tool calls or streaming failed
        # Use temp session for fallback too to avoid polluting main session
        print("[handle_command] Streaming didn't capture tool call, using regular invoke as fallback...")
        result = await asyncio.wait_for(
            asyncio.to_thread(agent.invoke, {"messages": [system_message, input_message]}, temp_session_config),
            timeout=30.0
        )
        
        # Try to extract tool mapping from result messages
        if isinstance(result, dict) and "messages" in result:
            messages = result.get("messages", [])
            for msg in reversed(messages):
                if hasattr(msg, 'tool_calls') and msg.tool_calls:
                    for tool_call in msg.tool_calls:
                        tool_name = getattr(tool_call, 'name', None)
                        tool_args = getattr(tool_call, 'args', None) or {}
                        if tool_name:
                            print(f"[handle_command] Extracted tool mapping from result: {tool_name}")
                            return {
                                "tool_mapping": {"tool_name": tool_name, "parameters": tool_args},
                                "tool_name": tool_name,
                                "parameters": tool_args,
                                "status": "tool_mapped",
                                "message": f"Agent identified tool: {tool_name}. Execution will be handled by router."
                            }
        
        # If no tool mapping found, try tool-generator-agent
        msg_count = len(result.get("messages", [])) if isinstance(result, dict) else 0
        print(f"[handle_command] result: {msg_count} messages in response (no tool mapping extracted)")
        print(f"[handle_command] No tool mapping found, attempting tool-generator-agent...")
        
        try:
            from backend.intent_classifier import classify_intent

            intent = classify_intent(command)
            if intent.intent != "execute":
                print(f"[handle_command] Skipping tool-generator-agent (intent={intent.intent}, reason={intent.reason})")
                return result

            # Use relative import since we're in the backend directory
            from tool_generator_agent import generate_and_execute_tool
            
            tool_result = await generate_and_execute_tool(
                command=command,
                user_request=command,
                session_id=session_id
            )
            
            if tool_result.get("status") == "success":
                print(f"[handle_command] ✅ Tool-generator-agent successfully generated tool")
                return {
                    "status": "success",
                    "tool_generated": True,
                    "tool_name": "generated_tool",
                    "result": tool_result,
                    "message": tool_result.get("explanation", "Tool generated and executed successfully")
                }
            else:
                print(f"[handle_command] ⚠️  Tool-generator-agent failed: {tool_result.get('error', 'Unknown error')}")
                # Fall through to return original result
        except Exception as e:
            print(f"[handle_command] ❌ Tool-generator-agent exception: {e}")
            import traceback
            traceback.print_exc()
            # Fall through to return original result
        
        return result
    except asyncio.TimeoutError:
        print(f"⚠️  Agent tool mapping timed out after 30s")
        print("🔍 Attempting to extract tool call from agent state...")
        # Try to extract tool call from temp session state
        extracted = _extract_tool_call_from_state(memory, temp_session_config)
        if extracted:
            return {
                "tool_mapping": extracted,
                "tool_name": extracted["tool_name"],
                "parameters": extracted["parameters"],
                "status": "tool_mapped",
                "message": f"Agent identified tool: {extracted['tool_name']}. Execution will be handled by router."
            }
        # If extraction fails, raise to trigger fallback
        raise TimeoutError("Agent tool mapping timed out and could not extract tool call")
    except Exception as e:
        # If agent fails with connection error or other exception, try to extract tool call
        error_str = str(e)
        if "Connection error" in error_str or "APIConnectionError" in error_str or "timeout" in error_str.lower() or "INVALID_CHAT_HISTORY" in error_str:
            print(f"⚠️  Agent tool mapping failed: {e}")
            print("🔍 Attempting to extract tool call from agent state...")
            # Try to extract tool call from temp session state
            extracted = _extract_tool_call_from_state(memory, temp_session_config)
            if extracted:
                return {
                    "tool_mapping": extracted,
                    "tool_name": extracted["tool_name"],
                    "parameters": extracted["parameters"],
                    "status": "tool_mapped",
                    "message": f"Agent identified tool: {extracted['tool_name']}. Execution will be handled by router."
                }
        # For other errors or if extraction fails, re-raise to trigger fallback
        raise e


def _extract_tool_call_from_state(memory, session_config):
    """Extract tool call (name + params) from agent state before execution."""
    try:
        thread_id = session_config.get("configurable", {}).get("thread_id", "default")
        
        # Access memory storage
        if hasattr(memory, 'storage'):
            storage = memory.storage
        elif hasattr(memory, '_storage'):
            storage = memory._storage
        else:
            return None
        
        if thread_id in storage:
            checkpoints = storage.get(thread_id, [])
            if checkpoints:
                latest_checkpoint = checkpoints[-1] if isinstance(checkpoints, list) else checkpoints
                if isinstance(latest_checkpoint, tuple) and len(latest_checkpoint) >= 2:
                    checkpoint = latest_checkpoint[1]
                else:
                    checkpoint = latest_checkpoint
                
                if checkpoint and isinstance(checkpoint, dict):
                    if "channel_values" in checkpoint:
                        messages = checkpoint["channel_values"].get("messages", [])
                        # Look for AIMessage with tool_calls
                        for msg in reversed(messages):
                            if hasattr(msg, 'tool_calls') and msg.tool_calls:
                                for tool_call in msg.tool_calls:
                                    tool_name = getattr(tool_call, 'name', None)
                                    tool_args = getattr(tool_call, 'args', None) or {}
                                    if tool_name:
                                        return {
                                            "tool_name": tool_name,
                                            "parameters": tool_args
                                        }
    except Exception as e:
        print(f"[handle_command] Error extracting tool call from state: {e}")
    return None


def _extract_tool_results_from_state(memory, session_config, error_msg):
    """Extract tool results from agent state when LLM call fails."""
    try:
        thread_id = session_config.get("configurable", {}).get("thread_id", "default")
        
        # MemorySaver stores checkpoints in a dict keyed by thread_id
        # Access the internal storage (this is implementation-specific but works for MemorySaver)
        if hasattr(memory, 'storage'):
            storage = memory.storage
        elif hasattr(memory, '_storage'):
            storage = memory._storage
        else:
            # Try to get via the public API if available
            try:
                # MemorySaver might have a list() or get() method
                checkpoints = list(memory.list(session_config)) if hasattr(memory, 'list') else []
                if checkpoints:
                    # Get the latest checkpoint
                    latest = checkpoints[-1] if checkpoints else None
                    if latest:
                        state = memory.get(session_config, latest) if hasattr(memory, 'get') else None
                        if state and "channel_values" in state:
                            messages = state["channel_values"].get("messages", [])
                            return _extract_from_messages(messages, error_msg)
            except Exception:
                pass
            print("⚠️  Could not access memory storage")
            return None
        
        if thread_id in storage:
            # MemorySaver stores checkpoints as a list of (config, checkpoint) tuples
            # Get the latest checkpoint for this thread
            checkpoints = storage.get(thread_id, [])
            if checkpoints:
                latest_checkpoint = checkpoints[-1] if isinstance(checkpoints, list) else checkpoints
                if isinstance(latest_checkpoint, tuple) and len(latest_checkpoint) >= 2:
                    checkpoint = latest_checkpoint[1]
                else:
                    checkpoint = latest_checkpoint
                
                if checkpoint and isinstance(checkpoint, dict):
                    # Checkpoint structure: {"channel_values": {"messages": [...]}}
                    if "channel_values" in checkpoint:
                        messages = checkpoint["channel_values"].get("messages", [])
                        return _extract_from_messages(messages, error_msg)
        
        print("⚠️  Could not extract tool results from state")
        return None
    except Exception as extract_err:
        print(f"❌ Could not extract tool results: {extract_err}")
        import traceback
        traceback.print_exc()
        return None


def _extract_from_messages(messages, error_msg):
    """Extract tool results from message list."""
    # Look for the last tool message (which contains the actual results)
    for msg in reversed(messages):
        if hasattr(msg, 'type') and msg.type == 'tool':
            tool_name = getattr(msg, 'name', 'unknown')
            print(f"✅ Found tool result in agent state: {tool_name}")
            # Extract the tool result
            tool_result = msg.content if hasattr(msg, 'content') else str(msg)
            return {
                "messages": messages,
                "tool_result": tool_result,
                "tool_name": tool_name,
                "status": "partial_success",
                "error": error_msg
            }
    return None


if __name__ == "__main__":
    seq_align_results = asyncio.run(handle_command("align the given sequences"))
    print(seq_align_results)

    mut_results = asyncio.run(handle_command("mutate the given sequence."))
    print(mut_results)


def _maybe_handle_vendor_request(command: str):
    """Detect vendor-focused requests and handle them directly."""
    command_lower = command.lower()
    
    # Exclude "test" if it's part of a file path (e.g., s3://bucket/test/file.fq)
    # Check if "test" appears in an S3 path or file path context
    is_test_in_path = bool(re.search(r's3://[^/]+/[^/]*test[^/]*/', command_lower) or 
                           re.search(r'[/\\][^/\\]*test[^/\\]*[/\\]', command_lower))
    
    # FastQC detection (HIGH PRIORITY - before vendor check)
    if any(phrase in command_lower for phrase in ['fastqc', 'fastqc analysis', 'quality control', 'quality analysis', 'perform fastqc']):
        return None  # Let it be handled by the agent's tool selection
    
    vendor_keywords = ['order', 'vendor', 'synthesis', 'assay', 'expression', 'function', 'binding']
    # Only include "test" if it's not in a file path
    if not is_test_in_path:
        vendor_keywords.append('test')
    
    if not any(keyword in command_lower for keyword in vendor_keywords):
        return None

    print("[handle_command] Detected vendor research command, using dna_vendor_research tool")

    sequence_length = None
    if '96' in command or 'variants' in command:
        sequence_length = 1000

    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))
    from dna_vendor_research import run_dna_vendor_research_raw

    result = run_dna_vendor_research_raw(command, sequence_length, 'large')

    return {
        "text": f"DNA vendor research completed: {result.get('message', 'Research done')}",
        "input": command,
        "output": result,
        "plot": {
            "data": [
                {
                    "x": ["Vendors", "Testing Options"],
                    "y": [result.get('total_vendors', 0), result.get('total_testing_options', 0)],
                    "type": "bar",
                }
            ],
            "layout": {"title": "DNA Vendor Research Results"},
        },
    }
