from __future__ import annotations

from typing import Any, Dict, Optional


def infer_action_type(command: str, tool_name: Optional[str] = None) -> str:
    c = (command or "").lower()
    t = (tool_name or "").lower()
    if any(k in c for k in ["what is", "explain", "describe"]):
        return "clarify_intent"
    if any(k in c for k in ["correct", "mislabeled", "should be"]):
        return "correct_metadata"
    if any(k in c for k in ["exclude", "subset", "focus only", "only female", "only male"]):
        return "subset_data"
    if any(k in c for k in ["compare", "difference", "changed significance"]):
        return "compare_versions"
    if any(k in c for k in ["plot", "volcano", "heatmap", "pca", "highlight", "color"]):
        return "generate_plot"
    if any(k in c for k in ["enrichment", "pathway", "go term"]):
        return "run_enrichment"
    if t in {"bulk_rnaseq_analysis", "single_cell_analysis"}:
        return "perform_group_comparison"
    if t in {"bio_rerun", "patch_and_rerun"}:
        return "rerun_downstream_steps"
    if t in {"bio_diff_runs"}:
        return "compare_versions"
    return "execute_tool"


def map_action_to_tool(action_type: str, fallback_tool: Optional[str], arguments: Dict[str, Any]) -> Optional[str]:
    """
    Resolve a generic action type to a concrete executor tool.
    """
    action = (action_type or "").lower()
    if action in {"perform_group_comparison"}:
        return fallback_tool or "bulk_rnaseq_analysis"
    if action in {"run_enrichment"}:
        return fallback_tool or "lookup_go_term"
    if action in {"generate_plot"}:
        return fallback_tool or "patch_and_rerun"
    if action in {"compare_versions"}:
        return fallback_tool or "bio_diff_runs"
    if action in {"rerun_downstream_steps", "correct_metadata", "subset_data"}:
        return fallback_tool or "bio_rerun"
    if action in {"clarify_intent"}:
        return fallback_tool or "handle_natural_command"
    return fallback_tool

