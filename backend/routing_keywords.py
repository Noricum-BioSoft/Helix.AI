"""
Shared routing keyword configuration.

Keeps keyword mappings in one place so command router, workflow planning, and
inline pipeline extraction do not drift.
"""

from __future__ import annotations

# Used by backend.agent._extract_inline_pipeline_plan
#
# Order matters: more specific patterns must come before broader fallbacks.
# Each entry maps a list of keyword fragments to a tool name that has a real
# handler in dispatch_tool (main.py).
INLINE_PIPELINE_TOOL_MAP = [
    # ── NGS QC / preprocessing ───────────────────────────────────────────────
    (["quality report", "qc report", "quality summary", "pipeline qc"], "quality_report"),
    (["fastqc", "quality assessment", "quality control"], "fastqc_quality_analysis"),
    (["trim", "adapter", "cutadapt", "trimmomatic", "fastp"], "read_trimming"),
    (["merge overlap", "merge paired", "overlapping paired", "flash", "pear", "pandaseq"], "read_merging"),
    (["merge", "overlap"], "read_merging"),

    # ── Phylogenetics ─────────────────────────────────────────────────────────
    (["phylogenetic tree", "phylogeny", "maximum-likelihood", "bootstrap", "iq-tree", "fasttree"], "phylogenetic_tree"),
    (["pairwise", "identity matrix", "pairwise amino", "pairwise identity"], "pairwise_identity"),
    # "annotate" alone is too broad — it matches "annotated by" in heatmap descriptions.
    # Use more specific phrases that only appear in true tree-annotation steps.
    (["annotate tree", "annotate the tree", "annotate.*voc", "annotate.*variant",
      "mutation site", "rbd mutation", "key mutation", "color-code by variant",
      "rectangular cladogram", "color.*lineage", "label.*voc", "label.*variant"], "annotate_tree"),
    (["multiple sequence alignment", "mafft", "muscle", "clustal"], "sequence_alignment"),
    (["retrieve", "fetch", "ncbi", "refseq", "genbank", "accession"], "fetch_ncbi_sequence"),
    (["align", " map ", "star ", "hisat", "bowtie"], "sequence_alignment"),

    # ── Bulk RNA-seq / DESeq2 ─────────────────────────────────────────────────
    # Map DESeq2-related steps to bulk_rnaseq_analysis which has a real handler.
    (["deseq2", "deseq", "edger", "limma", "lrt", "wald test", "likelihood ratio",
      "differential expression", "de genes", "de analysis", "pca plot", "pca ",
      "heatmap", "strip plot", "volcano", "count matrix", "load and validate",
      "sample metadata", "independent filtering", "benjamini", "padj", "log2fc",
      "temporal", "k-means", "kmeans", "cluster.*gene", "gene set", "gsea", "ora",
      "hallmark", "msigdb", "time.course", "timecourse", "time point"], "bulk_rnaseq_analysis"),

    # ── Single-cell RNA-seq ───────────────────────────────────────────────────
    (["single.cell", "scrna", "scanpy", "seurat", "cell ranger", "umap", "t-sne", "tsne",
      "pseudobulk", "cell type", "pbmc", "normalize.*count", "log.transform",
      "dot plot", "dotplot"], "single_cell_analysis"),

    # ── Quantification ────────────────────────────────────────────────────────
    (["quantif", "featurecount", "htseq", "salmon", "kallisto"], "bulk_rnaseq_analysis"),

    # ── Microbiome / 16S ──────────────────────────────────────────────────────
    (["diversity", "qiime", "kraken", "metaphlan", "microbiome", "16s", "otu", "asv"], "microbiome_analysis"),

    # ── Generic fallback: visualisation / summary reporting ──────────────────
    # "quality_report" requires FASTQ read-count metrics and should only be used
    # at the end of NGS preprocessing pipelines (fastqc → trim → merge).
    # Biological narrative / summary steps in RNA-seq or phylogenetics workflows
    # are gracefully handled as custom_step (simulated success with text output).
    (["quality report", "qc table", "pipeline qc summary", "pipeline summary table"], "quality_report"),
    # Everything else that looks like a text summary is a custom_step.
    (["biological interpretation", "biological narrative", "biological summary",
      "written summary", "written interpretation", "write a biological",
      "write a summary", "produce a written"], "custom_step"),
]

# Used by backend.workflow_planner_agent.SimpleOperationPlaybook.detect_operation
SIMPLE_OPERATION_REGEX = [
    (r"\bmerge\b", ("read_merging", "bbmerge")),
    (r"\btrim\b|\btrimming\b", ("trimming", "trimmomatic")),
    (r"\balign\b|\balignment\b", ("alignment", "bwa_mem")),
    (r"\bqc\b|\bquality\s*control\b", ("quality_control", "fastqc")),
    (r"\bquality\b", ("quality_control", "fastqc")),
    (r"\bblast\b", ("blast_search", "blastn")),
    (r"\bannotate\b|\bannotation\b", ("annotation", "prokka")),
]
