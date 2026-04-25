/**
 * CapabilityGrid — shown in the empty-state conversation pane.
 *
 * Replaces the "Run a command to see your prompts here" placeholder with a
 * visual overview of what Helix.AI can do, grouped by capability category.
 * Each row is clickable and pre-fills the input with an example prompt.
 */
import React from "react";

interface Capability {
  icon: string;
  title: string;
  description: string;
  examplePrompt: string;
}

interface CapabilityCategory {
  label: string;
  color: string;
  items: Capability[];
}

const CATEGORIES: CapabilityCategory[] = [
  {
    label: "Genomics & Sequencing",
    color: "#EFF6FF",
    items: [
      {
        icon: "🧬",
        title: "Multiple Sequence Alignment",
        description: "Align DNA, RNA, or protein sequences",
        examplePrompt: "Align these sequences and show me the result",
      },
      {
        icon: "🌳",
        title: "Phylogenetic Trees",
        description: "Build evolutionary trees from FASTA sequences",
        examplePrompt: "Build a phylogenetic tree from my sequences",
      },
      {
        icon: "📥",
        title: "Fetch from NCBI",
        description: "Download sequences by accession number",
        examplePrompt: "Fetch the sequence for accession NM_001301717 from NCBI",
      },
      {
        icon: "✂️",
        title: "Read QC & Trimming",
        description: "FastQC quality control and adapter trimming",
        examplePrompt: "Run FastQC quality control on my FASTQ files",
      },
    ],
  },
  {
    label: "RNA-seq & Single-Cell",
    color: "#F0FDF4",
    items: [
      {
        icon: "📊",
        title: "Bulk RNA-seq (DESeq2)",
        description: "Differential expression analysis on count matrices",
        examplePrompt: "Run differential expression analysis on my RNA-seq data",
      },
      {
        icon: "🔬",
        title: "Single-Cell RNA-seq",
        description: "Clustering, UMAP, and marker gene analysis",
        examplePrompt: "Run single-cell RNA-seq analysis on my .h5ad file",
      },
      {
        icon: "🗺️",
        title: "Update Visualizations",
        description: "Edit plots, change colors, or adjust axes",
        examplePrompt: "Update the volcano plot — highlight the top 10 genes in red",
      },
      {
        icon: "🔄",
        title: "Iterative Refinement",
        description: "Re-run analysis with different parameters",
        examplePrompt: "Re-run with a stricter threshold: alpha=0.01",
      },
    ],
  },
  {
    label: "Tabular & Data Science",
    color: "#FFFBEB",
    items: [
      {
        icon: "📋",
        title: "Tabular Analysis",
        description: "Sort, filter, rank, or aggregate CSV/Excel files",
        examplePrompt: "Sort my CSV file by the expression column and show the top 20 rows",
      },
      {
        icon: "❓",
        title: "Ask About Your Data",
        description: "Answer questions about any uploaded file",
        examplePrompt: "What is the average logFC value in my CSV?",
      },
      {
        icon: "🧪",
        title: "Exploratory Analysis",
        description: "Automatic EDA pipeline on any dataset",
        examplePrompt: "Run an exploratory data analysis on my uploaded dataset",
      },
      {
        icon: "📈",
        title: "GO / Pathway Enrichment",
        description: "Enrich gene lists against GO or KEGG",
        examplePrompt: "Run GO enrichment analysis on my list of significant genes",
      },
    ],
  },
];

interface Props {
  onSelectPrompt: (prompt: string) => void;
}

export const CapabilityGrid: React.FC<Props> = ({ onSelectPrompt }) => {
  return (
    <div className="mt-2">
      <p className="text-muted mb-3" style={{ fontSize: "0.9rem" }}>
        Here's what Helix.AI can do — click any item to try it:
      </p>
      <div className="d-flex flex-column gap-3">
        {CATEGORIES.map((cat) => (
          <div key={cat.label}>
            <div
              className="fw-semibold mb-2"
              style={{ fontSize: "0.78rem", color: "#6B7280", textTransform: "uppercase", letterSpacing: "0.05em" }}
            >
              {cat.label}
            </div>
            <div className="row g-2">
              {cat.items.map((item) => (
                <div key={item.title} className="col-6">
                  <button
                    className="w-100 text-start border-0 rounded-3 p-2 d-flex align-items-start gap-2"
                    style={{
                      background: cat.color,
                      cursor: "pointer",
                      transition: "filter 0.15s",
                    }}
                    onMouseEnter={(e) => (e.currentTarget.style.filter = "brightness(0.95)")}
                    onMouseLeave={(e) => (e.currentTarget.style.filter = "none")}
                    onClick={() => onSelectPrompt(item.examplePrompt)}
                    title={`Try: ${item.examplePrompt}`}
                  >
                    <span style={{ fontSize: "1.1rem", lineHeight: 1.2, flexShrink: 0 }}>{item.icon}</span>
                    <div>
                      <div style={{ fontSize: "0.82rem", fontWeight: 600, color: "#1E293B", lineHeight: 1.2 }}>
                        {item.title}
                      </div>
                      <div style={{ fontSize: "0.74rem", color: "#64748B", marginTop: "2px", lineHeight: 1.3 }}>
                        {item.description}
                      </div>
                    </div>
                  </button>
                </div>
              ))}
            </div>
          </div>
        ))}
      </div>
    </div>
  );
};
