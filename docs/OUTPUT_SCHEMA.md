# Output Format Schema

## Overview

This document defines the strict JSON output format that all agents must return for frontend rendering. The schema ensures consistent, predictable responses across all bioinformatics operations.

---

## Top-Level Response Schema

Every final response must be a single JSON object with these fields:

```json
{
  "domain": "bioinformatics|non_bioinformatics",
  "task_type": "qa|analysis|mixed",
  "status": "success|partial_success|failed|declined",
  "user_friendly_summary": "",
  "details_markdown": "",
  "data_inventory": {
    "inputs_detected": [],
    "assumptions": [],
    "qc_warnings": []
  },
  "artifacts": [],
  "errors": [],
  "provenance": {
    "tools_used": [],
    "references": [],
    "commands_run": [],
    "environment": {
      "os": "",
      "container_or_env": ""
    }
  }
}
```

### Field Descriptions

| Field | Type | Required | Description |
|---|---|---|---|
| `domain` | string | Yes | Classification: "bioinformatics" or "non_bioinformatics" |
| `task_type` | string | Yes | Type: "qa" (question/answer), "analysis" (computation), or "mixed" |
| `status` | string | Yes | Outcome: "success", "partial_success", "failed", or "declined" |
| `user_friendly_summary` | string | Yes | Brief plain-text summary (1-3 sentences) for user display |
| `details_markdown` | string | No | Detailed markdown explanation with formatting |
| `data_inventory` | object | No | Summary of detected inputs, assumptions, and warnings |
| `artifacts` | array | No | List of output artifacts (plots, tables, sequences, etc.) |
| `errors` | array | No | List of errors encountered during execution |
| `provenance` | object | No | Complete execution provenance for reproducibility |

---

## Status Codes

### success
Operation completed successfully with no errors.

```json
{
  "status": "success",
  "user_friendly_summary": "RNA-seq analysis completed. Found 1,234 differentially expressed genes."
}
```

### partial_success
Operation completed but with warnings or incomplete results.

```json
{
  "status": "partial_success",
  "user_friendly_summary": "QC analysis completed for 8 out of 10 samples. 2 samples failed quality thresholds.",
  "errors": [
    {
      "code": "LOW_QUALITY",
      "message": "Sample3.fastq has mean quality score below threshold (Q20)",
      "severity": "warning"
    }
  ]
}
```

### failed
Operation failed due to an error.

```json
{
  "status": "failed",
  "user_friendly_summary": "Alignment failed: reference genome not found.",
  "errors": [
    {
      "code": "REFERENCE_NOT_FOUND",
      "message": "Reference genome 'hg39' does not exist. Did you mean 'hg38' or 'hg19'?",
      "severity": "error"
    }
  ]
}
```

### declined
Operation was declined (out of scope, unsafe, or inappropriate).

```json
{
  "status": "declined",
  "user_friendly_summary": "Request declined: not a bioinformatics query.",
  "details_markdown": "I specialize in bioinformatics analysis. Your request about web development is outside my scope."
}
```

---

## Artifact Schemas

All artifacts must include: `type`, `title`, `description`, and `data`.

### Supported Artifact Types

#### sequence
Biological sequences (DNA, RNA, protein).

```json
{
  "type": "sequence",
  "title": "Aligned Sequences",
  "description": "Multiple sequence alignment of 10 protein sequences",
  "data": {
    "format": "fasta",
    "content": ">seq1\nATGCGTAC...\n>seq2\nATGCGTAC..."
  }
}
```

**Supported formats:**
- `fasta`: Standard FASTA format
- `fastq`: FASTQ with quality scores
- `clustal`: Clustal alignment format
- `stockholm`: Stockholm alignment format

#### table
Tabular data (CSV/TSV-like).

```json
{
  "type": "table",
  "title": "Differential Expression Results",
  "description": "Top 100 differentially expressed genes (FDR < 0.05)",
  "data": {
    "columns": ["Gene", "Log2FC", "P-value", "FDR"],
    "rows": [
      ["GENE1", 2.5, 0.0001, 0.001],
      ["GENE2", -1.8, 0.0005, 0.003]
    ]
  }
}
```

**Optional fields:**
- `column_types`: Array of types ("string", "number", "boolean")
- `row_count`: Total number of rows (for pagination)
- `uri`: URI if table is too large for inline display

#### chart
Interactive visualizations.

```json
{
  "type": "chart",
  "title": "PCA Plot",
  "description": "Principal component analysis of samples colored by condition",
  "data": {
    "chart_kind": "scatter",
    "axes": {
      "x": {"label": "PC1 (45.2%)", "range": [-10, 10]},
      "y": {"label": "PC2 (23.1%)", "range": [-8, 8]}
    },
    "series": [
      {
        "name": "Control",
        "points": [[1.2, 3.4], [2.1, 4.5]],
        "color": "#1f77b4",
        "marker": "circle"
      },
      {
        "name": "Treatment",
        "points": [[-2.3, -1.5], [-3.1, -2.2]],
        "color": "#ff7f0e",
        "marker": "circle"
      }
    ]
  }
}
```

**Supported chart types:**
- `scatter`: Scatter plot
- `line`: Line plot
- `bar`: Bar chart
- `heatmap`: Heatmap
- `volcano`: Volcano plot
- `histogram`: Histogram

#### tree
Phylogenetic or hierarchical trees.

```json
{
  "type": "tree",
  "title": "Phylogenetic Tree",
  "description": "Maximum likelihood tree of 20 sequences (1000 bootstraps)",
  "data": {
    "format": "newick",
    "newick": "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);",
    "metadata": {
      "method": "IQ-TREE",
      "model": "GTR+G",
      "bootstrap": 1000,
      "log_likelihood": -1234.56
    }
  }
}
```

#### genome_track
Genomic annotations and coverage tracks.

```json
{
  "type": "genome_track",
  "title": "ChIP-seq Peaks",
  "description": "ChIP-seq peaks for transcription factor X",
  "data": {
    "format": "bed",
    "content_or_uri": "chr1\t1000\t2000\tpeak1\t100\n...",
    "reference_genome": "hg38"
  }
}
```

**Supported formats:**
- `bed`: BED format
- `bigwig`: BigWig coverage
- `bedgraph`: BedGraph coverage

#### plasmid_map
Plasmid maps with features.

```json
{
  "type": "plasmid_map",
  "title": "pUC19 Plasmid Map",
  "description": "Standard cloning vector with multiple cloning site",
  "data": {
    "length_bp": 2686,
    "features": [
      {
        "name": "lacZ",
        "type": "CDS",
        "start": 203,
        "end": 1276,
        "strand": "+",
        "color": "#0000FF"
      },
      {
        "name": "ampR",
        "type": "resistance",
        "start": 1615,
        "end": 2475,
        "strand": "-",
        "color": "#FF0000"
      }
    ],
    "sequence_uri_or_inline": "ATGCGTAC..."
  }
}
```

#### image
Static images (plots, diagrams, etc.).

```json
{
  "type": "image",
  "title": "Quality Distribution",
  "description": "Per-base quality score distribution across all reads",
  "data": {
    "uri": "data:image/png;base64,iVBORw0KGgoAAAANSU...",
    "caption": "FastQC quality distribution plot",
    "format": "png"
  }
}
```

**Supported formats:**
- Base64-encoded: `data:image/png;base64,...`
- External URI: `https://...` or `s3://...`

#### log
Execution logs and debugging information.

```json
{
  "type": "log",
  "title": "Execution Log",
  "description": "Detailed execution log for alignment step",
  "data": {
    "steps": [
      {
        "timestamp": "2026-01-18T10:00:00Z",
        "level": "INFO",
        "message": "Starting alignment with STAR v2.7.10a"
      },
      {
        "timestamp": "2026-01-18T10:05:00Z",
        "level": "WARNING",
        "message": "Low memory warning: using 90% of available RAM"
      }
    ],
    "warnings": ["Low memory", "Long runtime"],
    "debug": {
      "command": "STAR --genomeDir hg38 ...",
      "exit_code": 0,
      "runtime_seconds": 300
    }
  }
}
```

---

## Provenance Schema

Complete execution history for reproducibility.

```json
{
  "provenance": {
    "tools_used": [
      {
        "name": "STAR",
        "version": "2.7.10a",
        "purpose": "RNA-seq alignment",
        "inputs_summary": "sample_R1.fastq, sample_R2.fastq",
        "outputs_summary": "aligned.bam",
        "session_reference": "Used uploaded files from session"
      },
      {
        "name": "featureCounts",
        "version": "2.0.1",
        "purpose": "Gene quantification",
        "inputs_summary": "aligned.bam",
        "outputs_summary": "counts.tsv",
        "session_reference": "Used alignment from STAR step"
      }
    ],
    "references": [
      "Reference genome: hg38 (GRCh38.p13)",
      "GTF annotation: Ensembl v104"
    ],
    "commands_run": [
      "STAR --genomeDir /ref/hg38 --readFilesIn sample_R1.fastq sample_R2.fastq --outFileNamePrefix aligned",
      "featureCounts -a genes.gtf -o counts.tsv aligned.bam"
    ],
    "environment": {
      "os": "Linux 5.10.0",
      "container_or_env": "helix/bioinformatics:v1.2.3"
    }
  }
}
```

### Session Context References

When using data from session context, explicitly reference it in provenance:

```json
{
  "tools_used": [
    {
      "name": "IQ-TREE",
      "purpose": "Phylogenetic tree inference",
      "session_reference": "Used aligned_sequences from alignment_1 step (MAFFT, completed 2026-01-18T10:05:00Z)"
    }
  ]
}
```

---

## Data Inventory

Summary of inputs, assumptions, and quality warnings.

```json
{
  "data_inventory": {
    "inputs_detected": [
      "2 paired-end FASTQ files (sample_R1.fq, sample_R2.fq)",
      "Total size: 950 MB",
      "Estimated 50M reads per file"
    ],
    "assumptions": [
      "Assumed unstranded library (not specified by user)",
      "Used default alignment parameters",
      "Filtered genes with <10 total counts"
    ],
    "qc_warnings": [
      "Sample R2 has slightly lower quality scores in last 10 bases",
      "3 samples have <30M reads (recommended minimum is 30M for RNA-seq)"
    ]
  }
}
```

---

## Error Handling

Structured error reporting.

```json
{
  "errors": [
    {
      "code": "MISSING_PARAMETER",
      "message": "Reference genome not specified. Please provide 'hg38', 'hg19', or 'mm10'.",
      "severity": "error",
      "context": {
        "parameter": "reference_genome",
        "suggestions": ["hg38", "hg19", "mm10"]
      }
    },
    {
      "code": "LOW_QUALITY_WARNING",
      "message": "Sample3 has mean Q-score of 28 (below recommended Q30)",
      "severity": "warning",
      "context": {
        "sample": "Sample3",
        "metric": "mean_q_score",
        "value": 28,
        "threshold": 30
      }
    }
  ]
}
```

### Error Severity Levels

- `error`: Operation failed, cannot proceed
- `warning`: Operation succeeded but with caveats
- `info`: Informational message, no action needed

---

## Complete Example: RNA-seq Analysis

```json
{
  "domain": "bioinformatics",
  "task_type": "analysis",
  "status": "success",
  "user_friendly_summary": "RNA-seq differential expression analysis completed. Found 1,234 significantly differentially expressed genes (FDR < 0.05) between control and treatment groups.",
  "details_markdown": "# RNA-seq Analysis Results\n\n## Summary\n- **Total genes analyzed**: 20,000\n- **Differentially expressed**: 1,234 (6.2%)\n- **Upregulated**: 645 genes\n- **Downregulated**: 589 genes\n\n## Quality Control\nAll samples passed QC thresholds (>30M reads, Q-score >30).\n\n## Methods\n- Alignment: STAR v2.7.10a\n- Quantification: featureCounts v2.0.1\n- Differential expression: DESeq2 v1.34.0\n- Multiple testing correction: Benjamini-Hochberg (FDR < 0.05)",
  "data_inventory": {
    "inputs_detected": [
      "6 paired-end FASTQ files (3 control, 3 treatment)",
      "Total size: 5.8 GB",
      "Average 45M reads per sample"
    ],
    "assumptions": [
      "Used hg38 reference genome (user-specified)",
      "Unstranded library (user-specified)",
      "Filtered genes with <10 total counts across all samples"
    ],
    "qc_warnings": []
  },
  "artifacts": [
    {
      "type": "table",
      "title": "Top 50 Differentially Expressed Genes",
      "description": "Genes ranked by adjusted p-value",
      "data": {
        "columns": ["Gene", "BaseMean", "Log2FC", "P-value", "FDR"],
        "rows": [
          ["GENE1", 1234.5, 3.2, 1.2e-50, 2.4e-46],
          ["GENE2", 2345.6, -2.8, 3.4e-45, 3.4e-41]
        ]
      }
    },
    {
      "type": "chart",
      "title": "PCA Plot",
      "description": "Principal component analysis showing clear separation between control and treatment groups",
      "data": {
        "chart_kind": "scatter",
        "axes": {
          "x": {"label": "PC1 (52.3%)", "range": [-30, 30]},
          "y": {"label": "PC2 (21.4%)", "range": [-20, 20]}
        },
        "series": [
          {
            "name": "Control",
            "points": [[10.2, 5.3], [12.1, 6.5], [11.5, 4.8]],
            "color": "#1f77b4",
            "marker": "circle"
          },
          {
            "name": "Treatment",
            "points": [[-15.3, -8.2], [-14.1, -7.5], [-16.2, -9.1]],
            "color": "#ff7f0e",
            "marker": "circle"
          }
        ]
      }
    },
    {
      "type": "chart",
      "title": "Volcano Plot",
      "description": "Log2 fold change vs. -log10(FDR) for all genes",
      "data": {
        "chart_kind": "volcano",
        "axes": {
          "x": {"label": "Log2 Fold Change", "range": [-6, 6]},
          "y": {"label": "-log10(FDR)", "range": [0, 50]}
        },
        "series": [
          {
            "name": "Significant",
            "points": [[3.2, 46.6], [-2.8, 41.5]],
            "color": "#d62728",
            "marker": "circle"
          },
          {
            "name": "Not Significant",
            "points": [[0.5, 1.2], [-0.3, 0.8]],
            "color": "#7f7f7f",
            "marker": "circle"
          }
        ]
      }
    }
  ],
  "errors": [],
  "provenance": {
    "tools_used": [
      {
        "name": "STAR",
        "version": "2.7.10a",
        "purpose": "RNA-seq alignment to hg38 reference genome",
        "inputs_summary": "6 paired-end FASTQ files",
        "outputs_summary": "6 BAM files (aligned reads)"
      },
      {
        "name": "featureCounts",
        "version": "2.0.1",
        "purpose": "Gene-level quantification",
        "inputs_summary": "6 BAM files from STAR alignment",
        "outputs_summary": "Gene counts matrix (20,000 genes × 6 samples)"
      },
      {
        "name": "DESeq2",
        "version": "1.34.0",
        "purpose": "Differential expression analysis",
        "inputs_summary": "Gene counts matrix from featureCounts",
        "outputs_summary": "DE results table (1,234 significant genes)"
      }
    ],
    "references": [
      "Reference genome: Homo sapiens hg38 (GRCh38.p13)",
      "Gene annotation: Ensembl v104 GTF",
      "Statistical method: DESeq2 with Benjamini-Hochberg FDR correction"
    ],
    "commands_run": [
      "STAR --genomeDir /ref/hg38 --readFilesIn sample1_R1.fq sample1_R2.fq --outFileNamePrefix sample1_",
      "featureCounts -a /ref/ensembl_v104.gtf -o counts.tsv sample1.bam sample2.bam sample3.bam sample4.bam sample5.bam sample6.bam",
      "DESeq2::DESeq(dds, contrast=c('condition','treatment','control'))"
    ],
    "environment": {
      "os": "Linux 5.10.0-23-cloud-amd64",
      "container_or_env": "helix/rnaseq:v2.1.0"
    }
  }
}
```

---

## Best Practices

### 1. Always Include Provenance
- Record all tools, versions, and parameters
- Enable reproducibility and debugging

### 2. Use Appropriate Artifact Types
- Choose the right artifact type for your data
- Don't force data into the wrong schema

### 3. Keep Responses Lean
- Use URIs for large data (>1MB)
- Downsample plots if necessary (max 20k points)
- Limit inline tables to 500-2,000 rows

### 4. Provide Context in Descriptions
- Explain what each artifact shows
- Note important findings or warnings

### 5. Reference Session Data
- When using previous results, note it in provenance
- Help users understand workflow continuity

---

## Related Documents

- `docs/SESSION_MANAGEMENT.md` - Session context and data retrieval
- `docs/TASK_ARCHITECTURE.md` - Micro-/Macroflow pattern
- `backend/execution_broker.py` - Execution result formatting
- `agents/workflow-planner-agent.md` - Workflow planning and output
