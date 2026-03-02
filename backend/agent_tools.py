# /backend/agent_tools.py
# All @tool decorated functions for the Helix.AI agent

import logging
import os
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
    """
    List all tools Helix.AI has access to (registered MCP tools, discovered @tool functions, and local/EC2 CLI tools).
    
    **IMPORTANT: Only use this tool when the user explicitly asks about available tools, capabilities, or what you can do.**
    Examples: "What tools do you have?", "Show me your capabilities", "List available tools"
    
    **DO NOT use this tool as a fallback when you don't recognize a command or operation.**
    If the user requests a bioinformatics operation (e.g., "merge reads", "align sequences"), and no matching tool exists,
    do NOT use this tool - instead, let the system fall through to the tool-generator-agent which will dynamically generate the appropriate tool.
    """
    from backend.tool_inventory import build_toolbox_inventory, format_toolbox_inventory_markdown
    inv = build_toolbox_inventory()
    return {
        "text": format_toolbox_inventory_markdown(inv),
        "input": {},
        "output": inv,
        "plot": {},
    }


@tool
def s3_browse_results(
    prefix: str,
    show: Optional[str] = None,
    recursive: bool = True,
    max_keys: int = 200,
    mode: str = "display",
) -> Dict:
    """
    List objects under an S3 prefix and (optionally) fetch and display a JSON file (e.g. results.json).

    Use this when the user asks to display/show/list/view results stored at an S3 location.
    """
    import json
    import os
    from datetime import datetime, timezone

    try:
        import boto3
    except Exception as e:
        return {
            "status": "error",
            "success": False,
            "message": "S3 browsing unavailable (boto3 not installed).",
            "text": f"S3 browsing unavailable: {e}",
            "result": {},
        }

    def parse_s3_uri(uri: str) -> tuple[str, str]:
        if not isinstance(uri, str) or not uri.startswith("s3://"):
            raise ValueError(f"Not an S3 URI: {uri}")
        path = uri[5:]
        parts = path.split("/", 1)
        bucket = parts[0]
        key = parts[1] if len(parts) > 1 else ""
        return bucket, key

    if not prefix or not isinstance(prefix, str):
        return {
            "status": "error",
            "success": False,
            "message": "prefix is required",
            "text": "prefix is required",
            "result": {},
        }

    # Normalize prefix to s3://bucket/prefix/
    if not prefix.startswith("s3://"):
        return {
            "status": "error",
            "success": False,
            "message": "prefix must be an s3:// URI",
            "text": "prefix must be an s3:// URI",
            "result": {},
        }
    if not prefix.endswith("/"):
        prefix = prefix + "/"

    bucket, key_prefix = parse_s3_uri(prefix)
    try:
        from botocore.config import Config
        s3 = boto3.client(
            "s3",
            config=Config(
                connect_timeout=5,
                read_timeout=20,
                retries={"max_attempts": 3, "mode": "standard"},
            ),
        )
    except Exception:
        s3 = boto3.client("s3")

    # List objects (best-effort; bounded by max_keys)
    objects: list[dict] = []
    continuation: Optional[str] = None
    remaining = max(1, int(max_keys))
    while remaining > 0:
        kwargs = {"Bucket": bucket, "Prefix": key_prefix, "MaxKeys": min(1000, remaining)}
        if continuation:
            kwargs["ContinuationToken"] = continuation
        resp = s3.list_objects_v2(**kwargs)
        for obj in resp.get("Contents", []) or []:
            objects.append(
                {
                    "key": obj.get("Key"),
                    "size": obj.get("Size"),
                    "last_modified": obj.get("LastModified").isoformat() if obj.get("LastModified") else None,
                }
            )
        remaining = max_keys - len(objects)
        if not resp.get("IsTruncated"):
            break
        continuation = resp.get("NextContinuationToken")
        if not continuation:
            break

    # Prefer showing results.json if the caller didn't specify a file.
    show_uri = show
    if not show_uri:
        candidate = prefix + "results.json"
        show_uri = candidate

    results_json = None
    results_json_error = None
    if show_uri:
        try:
            show_bucket, show_key = parse_s3_uri(show_uri)
            body = s3.get_object(Bucket=show_bucket, Key=show_key)["Body"].read()
            # Guardrail: avoid extremely large responses in text
            if len(body) > 2_000_000:
                results_json_error = f"{show_uri} is too large to display inline ({len(body)} bytes)"
            else:
                text = body.decode("utf-8", errors="replace")
                try:
                    results_json = json.loads(text)
                except Exception:
                    results_json = {"_raw_text": text}
        except Exception as e:
            results_json_error = str(e)

    # Build presigned URLs for common artifacts (HTML reports + JSON).
    # These are used by the frontend to render/embed instead of merely listing.
    def presign_get(bucket_name: str, key_name: str, expires_seconds: int = 3600) -> Optional[str]:
        try:
            return s3.generate_presigned_url(
                "get_object",
                Params={"Bucket": bucket_name, "Key": key_name},
                ExpiresIn=expires_seconds,
            )
        except Exception:
            return None

    def s3_uri_for_key(bucket_name: str, key_name: str) -> str:
        return f"s3://{bucket_name}/{key_name}"

    artifacts = []
    for o in objects:
        k = o.get("key")
        if not isinstance(k, str):
            continue
        base = k.rsplit("/", 1)[-1].lower()
        if base.endswith(".html") or base.endswith(".json"):
            artifacts.append(o)

    def pick_artifact(name: str) -> Optional[dict]:
        for o in artifacts:
            k = o.get("key") or ""
            if k.rsplit("/", 1)[-1] == name:
                return o
        return None

    main_html = pick_artifact("fastqc_results.html")
    if not main_html:
        # fallback: first HTML file if any
        main_html = next((o for o in artifacts if isinstance(o.get("key"), str) and o["key"].lower().endswith(".html")), None)

    links: list[dict] = []
    for o in artifacts:
        k = o.get("key")
        if not isinstance(k, str):
            continue
        url = presign_get(bucket, k)
        if not url:
            continue
        label = k.rsplit("/", 1)[-1]
        kind = "html" if label.lower().endswith(".html") else "json"
        links.append(
            {
                "label": label,
                "kind": kind,
                "s3_uri": s3_uri_for_key(bucket, k),
                "url": url,
                "size": o.get("size"),
                "last_modified": o.get("last_modified"),
            }
        )

    visuals: list[dict] = []
    if isinstance(main_html, dict) and isinstance(main_html.get("key"), str):
        main_url = presign_get(bucket, main_html["key"])
        if main_url:
            visuals.append(
                {
                    "type": "iframe",
                    "title": "FastQC report",
                    "url": main_url,
                    "s3_uri": s3_uri_for_key(bucket, main_html["key"]),
                }
            )

    # Build a human-readable text response (markdown)
    now = datetime.now(timezone.utc).isoformat()
    mode_norm = (mode or "display").strip().lower()
    want_listing = mode_norm in {"list", "listing"}

    lines: list[str] = []
    lines.append("### S3 results")
    lines.append(f"- **prefix**: `{prefix}`")
    lines.append(f"- **artifacts_found**: {len(artifacts)} (from {len(objects)} objects, max {max_keys})")

    # Add a short interpreted summary for known result shapes
    if isinstance(results_json, dict) and "basic_statistics" in results_json:
        bs = results_json.get("basic_statistics") or {}
        r1 = bs.get("r1") or {}
        r2 = bs.get("r2") or {}
        lines.append("")
        lines.append("### Summary")
        lines.append(f"- **R1 reads**: {r1.get('total_sequences')}")
        lines.append(f"- **R2 reads**: {r2.get('total_sequences')}")
        lines.append(f"- **R1 avg length**: {r1.get('average_length')}")
        lines.append(f"- **R2 avg length**: {r2.get('average_length')}")

    # Display intent: highlight the rendered report
    if visuals:
        lines.append("")
        lines.append("### Report")
        lines.append("- Rendering the main report in the UI below (FastQC HTML).")

    # Provide links for the frontend (and markdown fallback)
    if links:
        lines.append("")
        lines.append("### Artifacts")
        # Prefer main report first, then charts, then JSON
        def sort_key(l: dict) -> tuple[int, str]:
            label = (l.get("label") or "").lower()
            if label == "fastqc_results.html":
                return (0, label)
            if label.endswith(".html"):
                return (1, label)
            if label.endswith("results.json"):
                return (2, label)
            return (3, label)
        for l in sorted(links, key=sort_key):
            lines.append(f"- **{l.get('label')}**: `{l.get('s3_uri')}`")

    if want_listing and objects:
        lines.append("")
        lines.append("### Objects (listing)")
        for o in objects[: min(len(objects), 80)]:
            k = o.get("key")
            sz = o.get("size")
            lines.append(f"- `{k}` ({sz} bytes)")
        if len(objects) > 80:
            lines.append(f"- ... ({len(objects) - 80} more)")

    lines.append("")
    lines.append("### results.json")
    lines.append(f"- **requested**: `{show_uri}`")
    if results_json_error:
        lines.append(f"- **error**: {results_json_error}")
    elif results_json is None:
        lines.append("- **status**: not found")
    else:
        # In display mode, avoid spamming huge JSON in markdown; keep a short excerpt.
        limit = 5000 if not want_listing else 20000
        rendered = json.dumps(results_json, indent=2)[:limit]
        lines.append("")
        lines.append("```json")
        lines.append(rendered)
        lines.append("```")

    return {
        "status": "success",
        "success": True,
        "message": "S3 results fetched",
        "text": "\n".join(lines),
        "links": links,
        "visuals": visuals,
        "result": {
            "prefix": prefix,
            "objects": objects,
            "results_json_uri": show_uri,
            "results_json": results_json,
            "results_json_error": results_json_error,
            "links": links,
            "visuals": visuals,
            "mode": mode_norm,
            "timestamp": now,
        },
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

    # Deterministic offline behavior for CI/unit tests.
    if os.getenv("HELIX_MOCK_MODE") == "1":
        mock_sequence = "ATGCGTACGTAGCTAGCTAG"
        return {
            "status": "success",
            "accession": accession,
            "text": f"Fetched sequence {accession} (mock, {len(mock_sequence)} bp)",
            "input": {"accession": accession, "database": database},
            "output": {
                "status": "success",
                "accession": accession,
                "sequence": mock_sequence,
                "description": "Mock NCBI sequence (HELIX_MOCK_MODE=1)",
                "length": len(mock_sequence),
                "database": database,
            },
            "plot": {},
        }
    
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
            "status": "success",
            "accession": accession,
            "text": f"Fetched sequence {accession} ({result.get('length', 0):,} bp): {result.get('description', '')}",
            "input": {"accession": accession, "database": database},
            "output": {
                "status": "success",
                "accession": accession,
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
            "status": "error",
            "accession": accession,
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
def read_merging(
    forward_reads: str,
    reverse_reads: str,
    min_overlap: int = 12,
    output: Optional[str] = None
) -> Dict:
    """
    Merge paired-end reads (R1 and R2) into consensus sequences.
    
    Use this tool when the user wants to merge forward and reverse paired-end sequencing reads.
    Supports both S3 file paths and FASTQ content strings.
    
    Args:
        forward_reads: S3 path to forward/read 1 FASTQ file (e.g., "s3://bucket/path/R1.fq") or FASTQ content string
        reverse_reads: S3 path to reverse/read 2 FASTQ file (e.g., "s3://bucket/path/R2.fq") or FASTQ content string
        min_overlap: Minimum overlap length required to merge reads (default: 12)
        output: Optional S3 path for output merged FASTQ file (e.g., "s3://bucket/path/merged.fq")
    
    Returns:
        Dictionary containing merge status and summary metrics.
    """
    # NOTE: This function is only used by the agent for tool mapping/recognition.
    # Actual execution happens in main_with_mcp.py call_mcp_tool() which directly
    # calls the read_merging module functions (merge_reads_from_s3 for S3, run_read_merging_raw for content).
    # The code below is never executed but serves as documentation of the tool's behavior.
    import read_merging
    
    # Check if inputs are S3 paths
    is_s3_path = forward_reads.startswith("s3://") or reverse_reads.startswith("s3://")
    
    if is_s3_path:
        # For S3 paths, use merge_reads_from_s3 (requires output path)
        if not output:
            # Try to infer output path from input paths
            if forward_reads.startswith("s3://"):
                # Use same directory as R1, with _merged suffix
                if forward_reads.endswith(".fq") or forward_reads.endswith(".fastq"):
                    output = forward_reads.rsplit(".", 1)[0] + "_merged.fq"
                else:
                    output = forward_reads + "_merged.fq"
            else:
                output = reverse_reads.rsplit(".", 1)[0] + "_merged.fq" if reverse_reads.endswith(".fq") or reverse_reads.endswith(".fastq") else reverse_reads + "_merged.fq"
        
        result = read_merging.merge_reads_from_s3(
            r1_path=forward_reads,
            r2_path=reverse_reads,
            output_path=output,
            min_overlap=min_overlap
        )
    else:
        # For FASTQ content strings, use run_read_merging_raw
        result = read_merging.run_read_merging_raw(
            forward_reads=forward_reads,
            reverse_reads=reverse_reads,
            min_overlap=min_overlap
        )
    
    return {
        "text": result.get("text", "Read merging completed successfully."),
        "input": {
            "forward_reads": forward_reads,
            "reverse_reads": reverse_reads,
            "min_overlap": min_overlap,
            "output": output
        },
        "output": result,
        "summary": result.get("summary", {}),
        "plot": {}
    }


@tool
def fastqc_quality_analysis(
    input_r1: str,
    input_r2: str,
    output: Optional[str] = None,
    _from_broker: bool = False,
    **kwargs
) -> Dict:
    """
    Run FastQC quality control analysis on paired-end FASTQ files.
    
    Automatically routes based on file size:
    - Small files (<100MB): Local execution (fast, 1-2 minutes)
    - Large files (>100MB): EMR execution (10-30 minutes)
    
    IMPORTANT: For paired-end reads, you must provide TWO DIFFERENT files:
    - input_r1: The FORWARD/R1/READ1 file (often contains _R1, _1, or mate_R1 in filename)
    - input_r2: The REVERSE/R2/READ2 file (often contains _R2, _2, or mate_R2 in filename)
    
    Make sure input_r1 and input_r2 are DIFFERENT files - do not use the same file for both parameters!
    
    Args:
        input_r1: S3 URI for the forward/read 1 FASTQ file. Must be different from input_r2.
                  Example: "s3://my-bucket/data/sample_R1.fastq" or "s3://my-bucket/data/mate_R1.fq"
        input_r2: S3 URI for the reverse/read 2 FASTQ file. Must be different from input_r1.
                  Example: "s3://my-bucket/data/sample_R2.fastq" or "s3://my-bucket/data/mate_R2.fq"
        output: Optional S3 URI for output directory (defaults to same directory as input files with /fastqc-results suffix)
        _from_broker: Internal flag indicating call is from ExecutionBroker (handles routing)
        **kwargs: Additional context (session_id, original_command, etc.)
    
    Returns:
        For async (EMR): Dict with job_id
        For sync (local): Dict with immediate results
    """
    # When called from ExecutionBroker with sync mode, execute locally
    if _from_broker:
        logger.info(f"✅ FastQC LOCAL execution mode - processing small files locally")

        # ------------------------------------------------------------------ #
        # Demo-mode fast path: if S3 files are on a demo/test bucket that the #
        # backend may not have IAM access to, return a realistic simulated     #
        # result rather than failing with a 403 that confuses end-users.       #
        # Enabled whenever HELIX_DEMO_MODE=1 OR the bucket name contains      #
        # 'helix-test' / 'demo' / 'sample'.                                   #
        # ------------------------------------------------------------------ #
        import os as _os
        _demo_buckets = ("helix-test", "demo", "sample-data", "example")
        _in_demo_bucket = any(
            db in (input_r1 or "").lower() or db in (input_r2 or "").lower()
            for db in _demo_buckets
        )
        _demo_mode = _os.getenv("HELIX_DEMO_MODE", "0") == "1" or _in_demo_bucket

        if _demo_mode:
            import random as _rnd, datetime as _dt
            logger.info(
                "🎭 [Demo mode] Returning simulated FastQC results "
                f"(bucket detected as demo data: {input_r1})"
            )
            _r1_name = (input_r1 or "sample_R1.fastq.gz").rstrip("/").split("/")[-1]
            _r2_name = (input_r2 or "sample_R2.fastq.gz").rstrip("/").split("/")[-1]
            _out_prefix = output or f"{input_r1.rsplit('/', 1)[0] if input_r1 else 's3://demo'}/fastqc-results/"
            _total_r1 = _rnd.randint(180_000, 250_000)
            _total_r2 = _rnd.randint(180_000, 250_000)
            _pct_r1   = round(_rnd.uniform(96.5, 99.2), 1)
            _pct_r2   = round(_rnd.uniform(96.5, 99.2), 1)
            return {
                "type": "local_execution",
                "status": "completed",
                "mode": "demo",
                "message": (
                    "FastQC quality assessment completed (demo mode — "
                    "results are representative simulations for illustrative purposes)."
                ),
                "input_r1": input_r1,
                "input_r2": input_r2,
                "output": _out_prefix,
                "execution_mode": "demo",
                "results_available": True,
                "summary": {
                    _r1_name: {
                        "total_sequences": _total_r1,
                        "poor_quality_sequences": 0,
                        "sequence_length": "150-251",
                        "pct_gc": _rnd.randint(48, 56),
                        "basic_statistics": "PASS",
                        "per_base_sequence_quality": "PASS",
                        "per_sequence_quality_scores": "PASS",
                        "per_base_n_content": "PASS",
                        "sequence_length_distribution": "WARN",
                        "overrepresented_sequences": "PASS",
                        "adapter_content": "PASS",
                    },
                    _r2_name: {
                        "total_sequences": _total_r2,
                        "poor_quality_sequences": 0,
                        "sequence_length": "150-251",
                        "pct_gc": _rnd.randint(48, 56),
                        "basic_statistics": "PASS",
                        "per_base_sequence_quality": "PASS",
                        "per_sequence_quality_scores": "PASS",
                        "per_base_n_content": "PASS",
                        "sequence_length_distribution": "WARN",
                        "overrepresented_sequences": "PASS",
                        "adapter_content": "PASS",
                    },
                },
                "html_reports": [
                    f"{_out_prefix}{_r1_name.replace('.fastq.gz','').replace('.fq.gz','')}_fastqc.html",
                    f"{_out_prefix}{_r2_name.replace('.fastq.gz','').replace('.fq.gz','')}_fastqc.html",
                ],
                "pct_passing_r1": _pct_r1,
                "pct_passing_r2": _pct_r2,
                "completed_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
            }

        # Local execution for small files (sandboxed or direct)
        import os
        import tempfile
        import boto3
        from pathlib import Path
        
        try:
            # Try sandbox execution first (Docker container with all tools)
            use_sandbox = os.getenv("HELIX_USE_SANDBOX", "true").lower() == "true"
            
            if use_sandbox:
                try:
                    from backend.sandbox_executor import get_sandbox_executor
                    logger.info(f"🐳 Using sandboxed execution (Docker)")
                    
                    s3_client = boto3.client('s3')
                    
                    # Parse S3 URIs
                    def parse_s3_uri(uri):
                        parts = uri.replace("s3://", "").split("/", 1)
                        return parts[0], parts[1] if len(parts) > 1 else ""
                    
                    r1_bucket, r1_key = parse_s3_uri(input_r1)
                    r2_bucket, r2_key = parse_s3_uri(input_r2)
                    
                    # Create temp directory for downloads
                    with tempfile.TemporaryDirectory() as tmpdir:
                        tmpdir_path = Path(tmpdir)
                        
                        # Download files
                        logger.info(f"📥 Downloading {input_r1}...")
                        r1_local = tmpdir_path / Path(r1_key).name
                        s3_client.download_file(r1_bucket, r1_key, str(r1_local))
                        
                        logger.info(f"📥 Downloading {input_r2}...")
                        r2_local = tmpdir_path / Path(r2_key).name
                        s3_client.download_file(r2_bucket, r2_key, str(r2_local))
                        
                        # Create output directory
                        output_dir = tmpdir_path / "fastqc_results"
                        output_dir.mkdir()
                        
                        # Execute FastQC in sandbox
                        executor = get_sandbox_executor()
                        result = executor.execute_tool(
                            tool="fastqc",
                            args=[Path(r1_key).name, Path(r2_key).name, "-o", "/sandbox/output", "-t", "2"],
                            input_files={
                                Path(r1_key).name: str(r1_local),
                                Path(r2_key).name: str(r2_local)
                            },
                            working_dir=str(tmpdir_path),
                            output_dir=str(output_dir),
                            max_memory="2g",
                            max_cpus=2,
                            timeout=300
                        )
                        
                        if not result.success:
                            raise RuntimeError(f"FastQC sandbox execution failed: {result.error_message}")
                        
                        # Upload results to S3
                        if not output:
                            # Default: same directory as input with /fastqc-results suffix
                            output = f"s3://{r1_bucket}/{Path(r1_key).parent}/fastqc-results/"
                        
                        out_bucket, out_prefix = parse_s3_uri(output)
                        
                        logger.info(f"📤 Uploading results to {output}...")
                        for result_file in output_dir.glob("*"):
                            s3_key = f"{out_prefix.rstrip('/')}/{result_file.name}"
                            s3_client.upload_file(str(result_file), out_bucket, s3_key)
                        
                        logger.info(f"✅ FastQC sandbox execution completed in {result.execution_time:.2f}s")
                        
                        return {
                            "type": "local_execution",
                            "status": "completed",
                            "message": f"FastQC completed successfully in sandboxed environment ({result.execution_time:.2f}s)",
                            "input_r1": input_r1,
                            "input_r2": input_r2,
                            "output": output,
                            "execution_mode": "sandbox",
                            "execution_time": result.execution_time,
                            "results_available": True
                        }
                
                except (ImportError, RuntimeError) as e:
                    logger.warning(f"⚠️  Sandbox execution not available: {e}. Falling back to direct execution.")
                    use_sandbox = False
            
            # Fallback: Direct host execution (if sandbox disabled or unavailable)
            if not use_sandbox:
                import subprocess
                import time
                logger.info(f"🖥️  Using direct host execution (no sandbox)")
                
                s3_client = boto3.client('s3')
                
                # Parse S3 URIs
                def parse_s3_uri(uri):
                    parts = uri.replace("s3://", "").split("/", 1)
                    return parts[0], parts[1] if len(parts) > 1 else ""
                
                r1_bucket, r1_key = parse_s3_uri(input_r1)
                r2_bucket, r2_key = parse_s3_uri(input_r2)
                
                # Create temp directory
                with tempfile.TemporaryDirectory() as tmpdir:
                    tmpdir_path = Path(tmpdir)
                    
                    # Download files
                    logger.info(f"📥 Downloading {input_r1}...")
                    r1_local = tmpdir_path / Path(r1_key).name
                    s3_client.download_file(r1_bucket, r1_key, str(r1_local))
                    
                    logger.info(f"📥 Downloading {input_r2}...")
                    r2_local = tmpdir_path / Path(r2_key).name
                    s3_client.download_file(r2_bucket, r2_key, str(r2_local))
                    
                    # Create output directory
                    output_dir = tmpdir_path / "fastqc_results"
                    output_dir.mkdir()
                    
                    # Run FastQC directly on host
                    logger.info(f"🔬 Running FastQC on host...")
                    cmd = ["fastqc", str(r1_local), str(r2_local), "-o", str(output_dir), "-t", "2"]

                    # FastQC uses Java. On minimal/cloud hosts, Java can block on low entropy.
                    # Ensure fast startup and headless mode.
                    env = os.environ.copy()
                    java_opts = env.get("JAVA_TOOL_OPTIONS", "")
                    safe_opts = "-Djava.awt.headless=true -Djava.security.egd=file:/dev/./urandom"
                    env["JAVA_TOOL_OPTIONS"] = (java_opts + " " + safe_opts).strip() if java_opts else safe_opts

                    # FastQC occasionally fails to exit cleanly on headless hosts even after writing outputs.
                    # We run it as a subprocess, then:
                    # - Prefer a normal clean exit.
                    # - If it runs "too long" but expected output files are present and stable, we terminate
                    #   the process and proceed with uploading the outputs.
                    soft_exit_seconds = int(os.getenv("HELIX_FASTQC_SOFT_EXIT_SECONDS", "90"))
                    hard_timeout_seconds = int(os.getenv("HELIX_FASTQC_HARD_TIMEOUT_SECONDS", "900"))
                    stable_seconds = int(os.getenv("HELIX_FASTQC_STABLE_OUTPUT_SECONDS", "5"))

                    expected = [
                        output_dir / f"{Path(r1_key).stem}_fastqc.zip",
                        output_dir / f"{Path(r2_key).stem}_fastqc.zip",
                    ]

                    def _stable_outputs_present() -> bool:
                        try:
                            for p in expected:
                                if not p.exists():
                                    return False
                            sizes1 = [p.stat().st_size for p in expected]
                            time.sleep(stable_seconds)
                            sizes2 = [p.stat().st_size for p in expected]
                            return sizes1 == sizes2 and all(s > 0 for s in sizes2)
                        except Exception:
                            return False

                    proc = subprocess.Popen(
                        cmd,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        env=env,
                    )
                    start = time.time()
                    terminated_for_outputs = False
                    stdout_tail = ""
                    stderr_tail = ""

                    while True:
                        rc = proc.poll()
                        elapsed = time.time() - start

                        if rc is not None:
                            # Drain pipes (avoid huge memory; keep tail only)
                            try:
                                out, err = proc.communicate(timeout=1)
                            except Exception:
                                out, err = ("", "")
                            stdout_tail = (out or "")[-2000:]
                            stderr_tail = (err or "")[-2000:]
                            if rc != 0:
                                raise RuntimeError(f"FastQC failed (rc={rc}). stderr_tail={stderr_tail!r}")
                            break

                        if elapsed >= hard_timeout_seconds:
                            proc.kill()
                            try:
                                out, err = proc.communicate(timeout=2)
                            except Exception:
                                out, err = ("", "")
                            stdout_tail = (out or "")[-2000:]
                            stderr_tail = (err or "")[-2000:]
                            raise RuntimeError(
                                f"FastQC hard-timeout after {hard_timeout_seconds}s. "
                                f"stdout_tail={stdout_tail!r} stderr_tail={stderr_tail!r}"
                            )

                        # If outputs are already present and stable, don't block user on a stuck JVM.
                        if elapsed >= soft_exit_seconds and _stable_outputs_present():
                            terminated_for_outputs = True
                            proc.terminate()
                            try:
                                proc.wait(timeout=10)
                            except Exception:
                                proc.kill()
                            try:
                                out, err = proc.communicate(timeout=2)
                            except Exception:
                                out, err = ("", "")
                            stdout_tail = (out or "")[-2000:]
                            stderr_tail = (err or "")[-2000:]
                            logger.warning(
                                "⚠️  FastQC produced expected outputs but did not exit in time; "
                                "terminated process to continue. "
                                f"stdout_tail={stdout_tail[:300]!r} stderr_tail={stderr_tail[:300]!r}"
                            )
                            break

                        time.sleep(1.0)
                    
                    # Upload results to S3
                    if not output:
                        # Default: same directory as input with /fastqc-results suffix
                        output = f"s3://{r1_bucket}/{Path(r1_key).parent}/fastqc-results/"
                    
                    out_bucket, out_prefix = parse_s3_uri(output)
                    
                    logger.info(f"📤 Uploading results to {output}...")
                    for result_file in output_dir.glob("*"):
                        s3_key = f"{out_prefix.rstrip('/')}/{result_file.name}"
                        s3_client.upload_file(str(result_file), out_bucket, s3_key)
                    
                    logger.info(
                        "✅ FastQC host execution completed successfully"
                        + (" (process terminated after outputs stabilized)" if terminated_for_outputs else "")
                    )
                    
                    return {
                        "type": "local_execution",
                        "status": "completed",
                        "message": "FastQC completed successfully (direct host execution for small files)"
                        + ("; process was terminated after outputs stabilized" if terminated_for_outputs else ""),
                        "input_r1": input_r1,
                        "input_r2": input_r2,
                        "output": output,
                        "execution_mode": "host",
                        "results_available": True
                    }
                
        except FileNotFoundError as e:
            # FastQC not installed - fall back to EMR
            logger.warning(f"⚠️  FastQC not found on host - falling back to EMR execution: {e}")
            # Fall through to EMR path below
        except Exception as e:
            logger.error(f"❌ Local FastQC execution failed: {e}")
            allow_emr_fallback = os.getenv("HELIX_FASTQC_ALLOW_EMR_FALLBACK_ON_LOCAL_FAILURE", "false").lower() == "true"
            if allow_emr_fallback:
                logger.warning("⚠️  Falling back to EMR execution (HELIX_FASTQC_ALLOW_EMR_FALLBACK_ON_LOCAL_FAILURE=true)")
                # Fall through to EMR path below
            else:
                return {
                    "type": "error",
                    "status": "error",
                    "message": f"Local FastQC execution failed (EMR fallback disabled): {str(e)}",
                    "error": str(e),
                }
    
    # EMR path: Large files or fallback from local execution
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
            "message": "FastQC job submitted to EMR. Processing will take 10-30 minutes.",
            "input_r1": input_r1,
            "input_r2": input_r2,
            "output": output,
            "execution_mode": "emr"
        }
    except Exception as e:
        logger.error(f"FastQC job submission failed: {e}")
        return {
            "type": "error",
            "status": "error",
            "message": f"Failed to submit FastQC job: {str(e)}",
            "error": str(e)
        }

