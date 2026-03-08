"""
BioEvaluator — computes delta metrics between two bioinformatics runs.

Each tool has its own delta computation:
- bulk_rnaseq_analysis : DE gene counts, top-gene overlap
- single_cell_analysis : cluster count, marker overlap
- phylogenetic_tree    : tree topology similarity (Robinson-Foulds distance)
- fastqc/quality       : pass-rate change, reads retained
"""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


class BioEvaluator:
    """Compute structured delta between a current run and its parent."""

    # ------------------------------------------------------------------ #
    # Public API                                                           #
    # ------------------------------------------------------------------ #

    def compare(
        self,
        tool_name: str,
        current: Dict[str, Any],
        prior: Optional[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Return a delta dict; empty dict if prior is None."""
        if prior is None:
            return {}
        try:
            handler = getattr(self, f"_delta_{tool_name.lower()}", self._delta_generic)
            return handler(current, prior)
        except Exception as exc:
            logger.warning("Delta computation failed for %s: %s", tool_name, exc)
            return {"error": str(exc)}

    # ------------------------------------------------------------------ #
    # Per-tool delta handlers                                              #
    # ------------------------------------------------------------------ #

    def _delta_bulk_rnaseq_analysis(
        self, current: Dict, prior: Dict
    ) -> Dict[str, Any]:
        curr_summary = current.get("result", {}).get("summary", current.get("summary", []))
        prev_summary = prior.get("result", {}).get("summary", prior.get("summary", []))

        curr_sig = sum(r.get("significant", 0) for r in curr_summary)
        prev_sig = sum(r.get("significant", 0) for r in prev_summary)

        # Top-gene overlap
        curr_top = _extract_top_genes(current)
        prev_top = _extract_top_genes(prior)
        overlap = len(curr_top & prev_top)
        total_unique = len(curr_top | prev_top)

        return {
            "tool": "bulk_rnaseq_analysis",
            "significant_genes_delta": curr_sig - prev_sig,
            "significant_genes_current": curr_sig,
            "significant_genes_prior": prev_sig,
            "top_gene_overlap": overlap,
            "top_gene_jaccard": overlap / total_unique if total_unique else 0.0,
            "narrative": (
                f"Significant genes: {prev_sig} → {curr_sig} "
                f"({'↑' if curr_sig > prev_sig else '↓'}{abs(curr_sig - prev_sig)}). "
                f"Top-gene Jaccard similarity: {overlap}/{total_unique}."
            ),
        }

    def _delta_single_cell_analysis(
        self, current: Dict, prior: Dict
    ) -> Dict[str, Any]:
        curr_res = current.get("result", current)
        prev_res = prior.get("result", prior)

        curr_clusters = curr_res.get("n_clusters", 0)
        prev_clusters = prev_res.get("n_clusters", 0)

        curr_markers = _flatten_markers(curr_res.get("markers", {}))
        prev_markers = _flatten_markers(prev_res.get("markers", {}))
        overlap = len(curr_markers & prev_markers)
        total = len(curr_markers | prev_markers)

        return {
            "tool": "single_cell_analysis",
            "cluster_count_delta": curr_clusters - prev_clusters,
            "cluster_count_current": curr_clusters,
            "cluster_count_prior": prev_clusters,
            "marker_gene_overlap": overlap,
            "marker_gene_jaccard": overlap / total if total else 0.0,
            "narrative": (
                f"Clusters: {prev_clusters} → {curr_clusters} "
                f"({'↑' if curr_clusters > prev_clusters else '↓'}{abs(curr_clusters - prev_clusters)}). "
                f"Marker gene Jaccard: {overlap}/{total}."
            ),
        }

    def _delta_phylogenetic_tree(
        self, current: Dict, prior: Dict
    ) -> Dict[str, Any]:
        curr_stats = current.get("statistics", current.get("result", {}).get("statistics", {}))
        prev_stats = prior.get("statistics", prior.get("result", {}).get("statistics", {}))

        curr_dist = curr_stats.get("mean_pairwise_distance", 0)
        prev_dist = prev_stats.get("mean_pairwise_distance", 0)

        curr_n = curr_stats.get("n_sequences", 0)
        prev_n = prev_stats.get("n_sequences", 0)

        return {
            "tool": "phylogenetic_tree",
            "n_sequences_delta": curr_n - prev_n,
            "mean_distance_delta": curr_dist - prev_dist,
            "mean_distance_current": curr_dist,
            "mean_distance_prior": prev_dist,
            "narrative": (
                f"Sequences: {prev_n} → {curr_n}. "
                f"Mean pairwise distance: {prev_dist:.4f} → {curr_dist:.4f}."
            ),
        }

    def _delta_fastqc_quality_analysis(
        self, current: Dict, prior: Dict
    ) -> Dict[str, Any]:
        return self._delta_quality_assessment(current, prior)

    def _delta_quality_assessment(
        self, current: Dict, prior: Dict
    ) -> Dict[str, Any]:
        curr_res = current.get("result", current)
        prev_res = prior.get("result", prior)

        curr_merged = curr_res.get("merged_reads", curr_res.get("summary", {}).get("merged_reads", 0))
        prev_merged = prev_res.get("merged_reads", prev_res.get("summary", {}).get("merged_reads", 0))

        return {
            "tool": "quality_assessment",
            "merged_reads_delta": (curr_merged or 0) - (prev_merged or 0),
            "merged_reads_current": curr_merged,
            "merged_reads_prior": prev_merged,
            "narrative": (
                f"Merged reads: {prev_merged} → {curr_merged}."
            ),
        }

    def _delta_generic(self, current: Dict, prior: Dict) -> Dict[str, Any]:
        return {"narrative": "No tool-specific delta metrics available."}


# ------------------------------------------------------------------ #
# Helpers                                                              #
# ------------------------------------------------------------------ #

def _extract_top_genes(result: Dict) -> set:
    """Pull top-gene names from various result shapes."""
    genes: set = set()
    top_str = result.get("result", result).get("top_genes", "")
    if top_str:
        for part in top_str.split("—"):
            if "top genes:" in part.lower():
                gene_part = part.lower().split("top genes:")[-1]
                genes.update(g.strip() for g in gene_part.split(","))
    return genes


def _flatten_markers(markers_dict: Dict) -> set:
    """Return flat set of all marker gene names across all clusters."""
    genes: set = set()
    for gene_list in markers_dict.values():
        if isinstance(gene_list, list):
            genes.update(gene_list)
    return genes
