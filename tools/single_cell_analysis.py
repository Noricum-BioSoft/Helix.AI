"""
Single-cell RNA-seq analysis tool using scPipeline (R package).

This tool interfaces with the scPipeline R package to perform end-to-end
single-cell analysis including:
- Marker gene identification
- Differential expression analysis
- Pathway analysis
- Cell-type annotation
- Batch correction
"""

import subprocess
import json
import tempfile
import os
from pathlib import Path
from typing import Dict, Any, Optional, List
import pandas as pd
import base64
import io

# Get the tools directory for R scripts
TOOLS_DIR = Path(__file__).resolve().parent
R_SCRIPTS_DIR = TOOLS_DIR / "r_scripts"


def _ensure_r_scripts_dir():
    """Ensure the R scripts directory exists."""
    R_SCRIPTS_DIR.mkdir(exist_ok=True)


def _run_r_script(script_path: Path, args: List[str] = None) -> Dict[str, Any]:
    """
    Execute an R script and return the JSON result.
    
    Args:
        script_path: Path to the R script
        args: Additional arguments to pass to the script
    
    Returns:
        Dictionary containing the script output
    """
    if not script_path.exists():
        return {
            "status": "error",
            "message": f"R script not found: {script_path}"
        }
    
    cmd = ["Rscript", str(script_path)]
    if args:
        cmd.extend(args)
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=300  # 5 minute timeout
        )
        
        # Try to parse JSON output
        output = result.stdout.strip()
        if output:
            try:
                return json.loads(output)
            except json.JSONDecodeError:
                # If not JSON, return as text
                return {
                    "status": "success",
                    "output": output,
                    "stderr": result.stderr
                }
        else:
            return {
                "status": "success",
                "message": "Script executed successfully",
                "stderr": result.stderr
            }
    except subprocess.TimeoutExpired:
        return {
            "status": "error",
            "message": "R script execution timed out"
        }
    except subprocess.CalledProcessError as e:
        return {
            "status": "error",
            "message": f"R script execution failed: {e.stderr}",
            "stderr": e.stderr,
            "stdout": e.stdout
        }


def analyze_single_cell_data(
    data_file: Optional[str] = None,
    data_format: str = "10x",
    output_dir: Optional[str] = None,
    steps: Optional[List[str]] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Perform single-cell RNA-seq analysis using scPipeline.
    
    Args:
        data_file: Path to input data file (10x format directory, H5, or CSV)
        data_format: Format of input data ('10x', 'h5', 'csv', 'seurat')
        output_dir: Directory for output files
        steps: List of analysis steps to perform:
            - 'preprocessing': Quality control and normalization
            - 'markers': Find marker genes
            - 'differential': Differential expression analysis
            - 'pathways': Pathway enrichment analysis
            - 'annotation': Cell-type annotation
            - 'batch_correction': Batch effect correction
            - 'all': All steps (default)
        **kwargs: Additional parameters for specific steps
    
    Returns:
        Dictionary containing analysis results
    """
    _ensure_r_scripts_dir()
    
    # Default to all steps if not specified
    if steps is None:
        steps = ['all']
    
    # Create temporary output directory if not provided
    if output_dir is None:
        output_dir = tempfile.mkdtemp(prefix="scpipeline_")
    
    # Create main analysis script
    script_path = R_SCRIPTS_DIR / "scpipeline_analysis.R"
    
    # Prepare arguments
    args = [
        "--data-format", data_format,
        "--output-dir", output_dir,
        "--steps", ",".join(steps)
    ]
    
    if data_file:
        args.extend(["--data-file", data_file])
    
    # Add additional parameters
    for key, value in kwargs.items():
        args.extend([f"--{key.replace('_', '-')}", str(value)])
    
    # Run the analysis
    result = _run_r_script(script_path, args)
    
    # Add output directory to result
    result["output_dir"] = output_dir
    
    return result


def find_marker_genes(
    seurat_object_path: str,
    ident_column: str = "seurat_clusters",
    min_pct: float = 0.25,
    logfc_threshold: float = 0.25,
    only_pos: bool = True
) -> Dict[str, Any]:
    """
    Find marker genes for cell clusters using scPipeline.
    
    Args:
        seurat_object_path: Path to Seurat object RDS file
        ident_column: Column name for cluster identity
        min_pct: Minimum percentage of cells expressing the gene
        logfc_threshold: Minimum log fold change threshold
        only_pos: Only return positive markers
    
    Returns:
        Dictionary containing marker genes
    """
    _ensure_r_scripts_dir()
    
    script_path = R_SCRIPTS_DIR / "find_markers.R"
    
    args = [
        "--seurat-object", seurat_object_path,
        "--ident-column", ident_column,
        "--min-pct", str(min_pct),
        "--logfc-threshold", str(logfc_threshold),
        "--only-pos" if only_pos else "--all-markers"
    ]
    
    return _run_r_script(script_path, args)


def perform_differential_expression(
    seurat_object_path: str,
    group1: str,
    group2: str,
    ident_column: str = "seurat_clusters",
    min_pct: float = 0.1,
    logfc_threshold: float = 0.25
) -> Dict[str, Any]:
    """
    Perform differential expression analysis between two groups.
    
    Args:
        seurat_object_path: Path to Seurat object RDS file
        group1: First group identifier
        group2: Second group identifier
        ident_column: Column name for group identity
        min_pct: Minimum percentage of cells expressing the gene
        logfc_threshold: Minimum log fold change threshold
    
    Returns:
        Dictionary containing differential expression results
    """
    _ensure_r_scripts_dir()
    
    script_path = R_SCRIPTS_DIR / "differential_expression.R"
    
    args = [
        "--seurat-object", seurat_object_path,
        "--group1", group1,
        "--group2", group2,
        "--ident-column", ident_column,
        "--min-pct", str(min_pct),
        "--logfc-threshold", str(logfc_threshold)
    ]
    
    return _run_r_script(script_path, args)


def annotate_cell_types(
    seurat_object_path: str,
    reference: str = "HumanPrimaryCellAtlasData",
    method: str = "SingleR"
) -> Dict[str, Any]:
    """
    Annotate cell types using reference datasets.
    
    Args:
        seurat_object_path: Path to Seurat object RDS file
        reference: Reference dataset name (from celldex)
        method: Annotation method ('SingleR' or 'scType')
    
    Returns:
        Dictionary containing cell type annotations
    """
    _ensure_r_scripts_dir()
    
    script_path = R_SCRIPTS_DIR / "annotate_cell_types.R"
    
    args = [
        "--seurat-object", seurat_object_path,
        "--reference", reference,
        "--method", method
    ]
    
    return _run_r_script(script_path, args)


def perform_batch_correction(
    seurat_object_path: str,
    batch_column: str = "batch",
    method: str = "harmony"
) -> Dict[str, Any]:
    """
    Perform batch correction on single-cell data.
    
    Args:
        seurat_object_path: Path to Seurat object RDS file
        batch_column: Column name for batch information
        method: Batch correction method ('harmony', 'seurat', or 'batchelor')
    
    Returns:
        Dictionary containing batch-corrected results
    """
    _ensure_r_scripts_dir()
    
    script_path = R_SCRIPTS_DIR / "batch_correction.R"
    
    args = [
        "--seurat-object", seurat_object_path,
        "--batch-column", batch_column,
        "--method", method
    ]
    
    return _run_r_script(script_path, args)


def analyze_pathways(
    seurat_object_path: str,
    ident_column: str = "seurat_clusters",
    database: str = "Reactome"
) -> Dict[str, Any]:
    """
    Perform pathway enrichment analysis.
    
    Args:
        seurat_object_path: Path to Seurat object RDS file
        ident_column: Column name for cluster identity
        database: Pathway database ('Reactome', 'GO', or 'KEGG')
    
    Returns:
        Dictionary containing pathway analysis results
    """
    _ensure_r_scripts_dir()
    
    script_path = R_SCRIPTS_DIR / "pathway_analysis.R"
    
    args = [
        "--seurat-object", seurat_object_path,
        "--ident-column", ident_column,
        "--database", database
    ]
    
    return _run_r_script(script_path, args)


def process_single_cell_command(
    command: str,
    session_context: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Process a natural language command for single-cell analysis.
    
    This function parses commands and routes them to appropriate scPipeline functions.
    
    Args:
        command: Natural language command
        session_context: Session context with previous results
    
    Returns:
        Dictionary containing analysis results
    """
    command_lower = command.lower()
    
    # Extract data file from command or session context
    data_file = None
    if session_context:
        # Check for uploaded files
        uploaded_files = session_context.get("uploaded_files", [])
        for file_info in uploaded_files:
            if any(ext in file_info.get("name", "").lower() for ext in [".h5", ".h5ad", ".rds", "matrix.mtx"]):
                data_file = file_info.get("name")
                break
    
    # Determine analysis steps from command
    steps = []
    if "marker" in command_lower or "markers" in command_lower:
        steps.append("markers")
    if "differential" in command_lower or "deg" in command_lower:
        steps.append("differential")
    if "pathway" in command_lower or "enrichment" in command_lower:
        steps.append("pathways")
    if "annotate" in command_lower or "cell type" in command_lower:
        steps.append("annotation")
    if "batch" in command_lower or "correct" in command_lower:
        steps.append("batch_correction")
    if not steps or "all" in command_lower or "complete" in command_lower:
        steps = ["all"]
    
    # Run analysis
    result = analyze_single_cell_data(
        data_file=data_file,
        steps=steps,
        session_context=session_context
    )
    
    return {
        "status": "success" if result.get("status") == "success" else "error",
        "result": result,
        "text": f"Single-cell analysis completed. Steps: {', '.join(steps)}"
    }

