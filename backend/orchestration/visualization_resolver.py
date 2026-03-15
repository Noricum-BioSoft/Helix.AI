from __future__ import annotations

from typing import Any, Dict


def determine_visualization_type(tool: str, result: Any, prompt: str) -> str:
    """
    Artifact-first visualization resolver.
    Uses structured result metadata first; falls back to tool and prompt heuristics.
    """
    tool_lower = (tool or "").lower()
    prompt_lower = (prompt or "").lower()

    if isinstance(result, dict):
        artifact_kind = str(result.get("artifact_kind") or result.get("result_type") or "").lower()
        plot_family = str(result.get("plot_family") or "").lower()
        if artifact_kind in {"deg_table", "enrichment", "results_viewer"}:
            return "results_viewer"
        if plot_family in {"pca", "volcano", "heatmap"}:
            return "results_viewer"
        if artifact_kind in {"sequence", "fasta"}:
            return "sequence_viewer"
        if artifact_kind in {"alignment"}:
            return "alignment_viewer"
        if artifact_kind in {"phylogenetic_tree", "newick"}:
            return "phylogenetic_tree"
        if result.get("visualization_type"):
            return str(result["visualization_type"])

    if tool_lower == "agent":
        return "markdown"
    if tool_lower in ("s3_browse_results", "bulk_rnaseq_analysis"):
        return "results_viewer"
    if tool_lower in ("single_cell_analysis", "quality_assessment"):
        if isinstance(result, dict) and (result.get("visuals") or result.get("links")):
            return "results_viewer"

    if tool_lower in {"fetch_ncbi_sequence", "query_uniprot"}:
        return "sequence_viewer"
    if tool_lower == "sequence_alignment":
        return "alignment_viewer"
    if tool_lower == "lookup_go_term":
        return "go_term_viewer"
    if tool_lower == "plasmid_visualization":
        return "plasmid_viewer"
    if tool_lower == "phylogenetic_tree":
        return "phylogenetic_tree"
    if tool_lower in {"quality_assessment", "read_trimming"}:
        return "quality_plot"

    if isinstance(result, dict):
        if result.get("sequence") or result.get("sequences"):
            return "sequence_viewer"
        if result.get("alignment") or result.get("tree_data") or result.get("newick"):
            return "alignment_viewer" if result.get("alignment") else "phylogenetic_tree"
        if result.get("plot_data") or result.get("metrics"):
            return "quality_plot"

    if any(q in prompt_lower for q in ["what is", "what are", "explain", "tell me about", "describe", "how does"]):
        return "markdown"
    return "default"

