"""
BioReviewer — generates human-readable narrative summaries from bio run results.

Does NOT call an LLM; uses templated text to stay deterministic and fast.
"""
from __future__ import annotations

import logging
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)


class BioReviewer:
    """Template-based narrative summariser for bioinformatics runs."""

    def summarize(
        self,
        tool_name: str,
        result: Dict[str, Any],
        delta: Optional[Dict[str, Any]] = None,
        objective: str = "",
    ) -> str:
        """Return a markdown narrative for the run result (and delta if provided)."""
        try:
            handler = getattr(self, f"_review_{tool_name.lower()}", self._review_generic)
            text = handler(result, objective)
        except Exception as exc:
            logger.warning("Review generation failed for %s: %s", tool_name, exc)
            text = f"Analysis completed (review generation failed: {exc})."

        if delta and delta.get("narrative"):
            text += f"\n\n### Comparison to Prior Run\n\n{delta['narrative']}"

        return text

    # ------------------------------------------------------------------ #
    # Per-tool review templates                                            #
    # ------------------------------------------------------------------ #

    def _review_bulk_rnaseq_analysis(
        self, result: Dict, objective: str
    ) -> str:
        res = result.get("result", result)
        n_samples = res.get("n_samples", "?")
        n_genes = res.get("n_genes_total", "?")
        summary = res.get("summary", [])
        mode = res.get("mode", "real")

        sig_total = sum(r.get("significant", 0) for r in summary)
        contrasts_text = "\n".join(
            f"- **{r['contrast']}**: {r['significant']} significant genes "
            f"({r['upregulated']} ↑, {r['downregulated']} ↓)"
            for r in summary
        )

        top = res.get("top_genes", "")
        top_section = f"\n### Top Differentially Expressed Genes\n{top}" if top else ""

        mode_note = (
            "\n> ⚠️ Synthetic demo data — representative results only."
            if mode == "synthetic" else ""
        )

        obj_note = f"\n**Objective:** {objective}" if objective else ""

        return (
            f"### Bulk RNA-seq DE Analysis Review{obj_note}\n\n"
            f"**Dataset:** {n_samples} samples, {n_genes} genes  \n"
            f"**Total significant genes (across all contrasts):** {sig_total}\n\n"
            f"#### Contrasts\n{contrasts_text}"
            f"{top_section}"
            f"{mode_note}"
        )

    def _review_single_cell_analysis(
        self, result: Dict, objective: str
    ) -> str:
        res = result.get("result", result)
        n_cells = res.get("n_cells", "?")
        n_genes = res.get("n_genes", "?")
        n_clusters = res.get("n_clusters", "?")
        resolution = res.get("resolution", "?")
        mode = res.get("mode", "real")

        ct_comp = res.get("cell_type_composition", {})
        ct_text = ""
        if ct_comp:
            total = sum(ct_comp.values())
            ct_text = "\n".join(
                f"- **{ct}**: {n} ({n / total * 100:.0f}%)"
                for ct, n in sorted(ct_comp.items(), key=lambda x: -x[1])[:5]
            )

        mode_note = (
            "\n> ⚠️ Synthetic PBMC data — cell-type composition is representative."
            if mode == "synthetic" else ""
        )

        return (
            f"### Single-Cell RNA-seq Analysis Review\n\n"
            f"**Cells:** {n_cells}  |  **Genes:** {n_genes}  |  "
            f"**Clusters:** {n_clusters}  |  **Resolution:** {resolution}\n\n"
            + (f"#### Top Cell Types\n{ct_text}\n\n" if ct_text else "")
            + mode_note
        )

    def _review_phylogenetic_tree(
        self, result: Dict, objective: str
    ) -> str:
        stats = result.get("statistics", result.get("result", {}).get("statistics", {}))
        n_seqs = stats.get("n_sequences", "?")
        mean_dist = stats.get("mean_pairwise_distance")
        mean_id = stats.get("mean_pairwise_identity_pct")

        dist_text = (
            f"  \n**Mean pairwise distance:** {mean_dist:.4f}"
            if mean_dist is not None else ""
        )
        id_text = (
            f"  \n**Mean pairwise identity:** {mean_id:.1f}%"
            if mean_id is not None else ""
        )

        return (
            f"### Phylogenetic Analysis Review\n\n"
            f"**Sequences:** {n_seqs}"
            f"{dist_text}{id_text}\n\n"
            f"Neighbor-joining tree constructed from pairwise Hamming distances. "
            f"Clusters visualised as hierarchical dendrogram."
        )

    def _review_quality_assessment(
        self, result: Dict, objective: str
    ) -> str:
        return self._review_fastqc_quality_analysis(result, objective)

    def _review_fastqc_quality_analysis(
        self, result: Dict, objective: str
    ) -> str:
        # Result can be nested: result["result"] or result directly
        res = result.get("result", result)
        if not isinstance(res, dict):
            res = result

        status = res.get("status") or result.get("status", "completed")
        mode = res.get("mode", "")

        # Extract QC summary — dict of {filename: {total_sequences, pct_passing, ...}}
        qc_summary = res.get("summary", {})
        if isinstance(qc_summary, dict):
            total_r1 = total_r2 = None
            pass_r1 = pass_r2 = None
            files = list(qc_summary.keys())
            if files:
                first = qc_summary[files[0]]
                total_r1 = first.get("total_sequences")
                pass_r1 = res.get("pct_passing_r1")
            if len(files) >= 2:
                second = qc_summary[files[1]]
                total_r2 = second.get("total_sequences")
                pass_r2 = res.get("pct_passing_r2")

            reads_line = ""
            if total_r1 is not None:
                reads_line = f"**R1:** {total_r1:,} reads"
                if pass_r1 is not None:
                    reads_line += f" ({pass_r1:.0f}% passing)"
            if total_r2 is not None:
                reads_line += f"  |  **R2:** {total_r2:,} reads"
                if pass_r2 is not None:
                    reads_line += f" ({pass_r2:.0f}% passing)"
        else:
            # Legacy merged_reads path (read_merging tool output)
            merged = res.get("merged_reads", "?")
            rate = res.get("merge_rate")
            rate_text = f" ({rate * 100:.1f}% merge rate)" if rate else ""
            reads_line = f"**Merged reads:** {merged}{rate_text}"

        mode_note = (
            "\n> ⚠️ Simulated FastQC results — real FASTQ files unavailable."
            if mode in ("demo", "simulated") else ""
        )

        return (
            f"### Amplicon QC Review\n\n"
            f"{reads_line}\n\n"
            f"FastQC quality assessment completed."
            f"{mode_note}"
        )

    def _review_generic(self, result: Dict, objective: str) -> str:
        status = result.get("status", "unknown")
        return (
            f"Analysis completed with status: **{status}**."
        )
