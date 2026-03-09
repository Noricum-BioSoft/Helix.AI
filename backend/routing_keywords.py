"""
Shared routing keyword configuration.

Keeps keyword mappings in one place so command router, workflow planning, and
inline pipeline extraction do not drift.
"""

from __future__ import annotations

# Used by backend.agent._extract_inline_pipeline_plan
INLINE_PIPELINE_TOOL_MAP = [
    (["quality report", "qc report", "quality summary"], "quality_report"),
    (["fastqc", "quality assessment", "quality control"], "fastqc_quality_analysis"),
    (["trim", "adapter", "cutadapt", "trimmomatic"], "read_trimming"),
    (["merge overlap", "merge paired", "overlapping paired", "flash", "pear"], "read_merging"),
    (["merge", "overlap"], "read_merging"),
    (["phylogenetic tree", "phylogeny", "maximum-likelihood", "bootstrap"], "phylogenetic_tree"),
    (["pairwise", "identity matrix", "pairwise amino", "pairwise identity"], "pairwise_identity"),
    (["annotate", "mutation site", "rbd mutation", "key mutation"], "annotate_tree"),
    (["multiple sequence alignment", "mafft", "muscle", "clustal"], "sequence_alignment"),
    (["retrieve", "fetch", "ncbi", "refseq", "genbank"], "fetch_ncbi_sequence"),
    (["align", "map ", "star ", "hisat", "bowtie"], "sequence_alignment"),
    (["quantif", "featurecount", "htseq"], "read_quantification"),
    (["differential", "deseq", "edger", "limma"], "differential_expression"),
    (["single.cell", "scanpy", "seurat", "cell ranger"], "single_cell_analysis"),
    (["diversity", "qiime", "kraken", "metaphlan"], "microbiome_analysis"),
    (["report", "summary", "csv", "visualiz", "plot"], "quality_report"),
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
