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
        examplePrompt: `Align these 4 CCR7 protein sequences and show the alignment:

>CCR7_human
MAAASNNTSSGSGTESNYYTTRESMLGGHSSVSTESEVAQNLSSMLLSNLTEQLALLVNQNRSIASAQEDQAFQIFQM
>CCR7_mouse
MAAASNNTSSGSGTENNYYTTRESMLGGHSSVSTESEVAQNLSSMLLSNLTEQLALLVNQNRSIASAQEDQAFQIFQM
>CCR7_rat
MAAASNNTSSGSGTENNYYTTRESMLGGHSSVSTESEVAQNLSSMLLSNLTEQLALLVNQNRSIASAQEDQAFQIFQM
>CCR7_chimp
MAAASNNTSSGSGTESNYYTTRESMLGGHSSVSTESEVAQNLSSMLLSNLTEQLALLVNQNRSIASAQEDQAFQIFQM`,
      },
      {
        icon: "🌳",
        title: "Phylogenetic Trees",
        description: "Build evolutionary trees from FASTA sequences",
        examplePrompt: `Build a phylogenetic tree from these CCR7 sequences across species:

>CCR7_human
MAAASNNTSSGSGTESNYYTTRESMLGGHSSVSTESEVAQNLSSMLLSNLTEQLALLVNQNRSIASAQEDQAFQIFQM
>CCR7_mouse
MAAASNNTSSGSGTENNYYTTRESMLGGHSSVSTESEVAQNLSSMLLSNLTEQLALLVNQNRSIASAQEDQAFQIFQM
>CCR7_rat
MAAASNNTSSGSGTENNYYTTRESMLGGHSSVSTESEVAQNLSSMLLSNLTEQLALLVNQNRSIASAQEDQAFQIFQM
>CCR7_chimp
MAAASNNTSSGSGTESNYYTTRESMLGGHSSVSTESEVAQNLSSMLLSNLTEQLALLVNQNRSIASAQEDQAFQIFQM
>CCR7_zebrafish
MPAASNNTSSGSGTESNYYTNRESMLGGHSSVSTESEVAQNLSSMLLSNLTEQLALLVNQNRSIASAQEDQAFQIFQM`,
      },
      {
        icon: "📥",
        title: "Fetch from NCBI",
        description: "Download sequences by accession number",
        examplePrompt:
          "Fetch the sequence for accession NM_001301717 from NCBI",
      },
      {
        icon: "✂️",
        title: "Read QC & Trimming",
        description: "FastQC quality control and adapter trimming",
        examplePrompt:
          "I uploaded a FASTQ file called sample_R1.fastq.gz. Please run FastQC quality control on it and show me the per-base quality scores, sequence duplication levels, and flag any adapter contamination.",
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
        examplePrompt: `Run differential expression analysis (DESeq2, tumor vs normal) on this counts matrix:

gene,tumor_1,tumor_2,tumor_3,tumor_4,normal_1,normal_2,normal_3,normal_4
TP53,120,145,98,132,890,920,845,912
MYC,2340,2890,2101,2567,456,389,421,445
BRCA1,89,76,92,85,456,489,412,478
EGFR,678,712,645,698,234,198,212,245
VEGFA,1230,1456,1098,1345,234,267,198,245
CDK4,567,623,489,534,123,145,112,134
PTEN,234,198,256,221,789,823,756,812
PIK3CA,890,923,845,878,456,489,412,478
AKT1,456,512,389,478,234,256,212,245
KRAS,1234,1389,1098,1267,345,389,312,356`,
      },
      {
        icon: "🔬",
        title: "Single-Cell RNA-seq",
        description: "Clustering, UMAP, and marker gene analysis",
        examplePrompt:
          "I uploaded a PBMC single-cell dataset (pbmc3k.h5ad, 2,700 cells). Please cluster the cells with resolution 0.5, compute a UMAP, and show me the top 5 marker genes per cluster with their p-values.",
      },
      {
        icon: "🗺️",
        title: "Update Visualizations",
        description: "Edit plots, change colors, or adjust axes",
        examplePrompt:
          "Update the last volcano plot: highlight all genes with |log2FC| > 2 in red, add dashed threshold lines at padj = 0.05, and label the top 5 significant genes",
      },
      {
        icon: "🔄",
        title: "Iterative Refinement",
        description: "Re-run analysis with different parameters",
        examplePrompt:
          "Re-run the last differential expression analysis with stricter thresholds: alpha = 0.01 and minimum log2FC = 1.5",
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
        examplePrompt: `What are the top 5 genes by absolute log2FC value in this differential expression table, and which of them are statistically significant (padj < 0.05)?

gene,log2FC,padj,baseMean
TP53,-2.14,0.00012,450.3
BRCA1,3.42,0.000045,892.1
MYC,1.83,0.0082,2103.5
EGFR,-0.91,0.048,310.7
VEGFA,2.71,0.00031,678.4
CDK4,1.21,0.019,543.2
PTEN,-1.54,0.0075,189.6
KRAS,0.88,0.062,789.3
CDKN2A,-2.33,0.00089,234.5
PIK3CA,1.67,0.0041,445.8
ERBB2,2.95,0.00018,567.9
AKT1,-1.12,0.027,334.6`,
      },
      {
        icon: "❓",
        title: "Ask About Your Data",
        description: "Answer questions about any uploaded file",
        examplePrompt: `How many genes in this table have an adjusted p-value below 0.01, and what is the average log2FC for those significant genes?

gene,log2FC,padj,baseMean
TP53,-2.14,0.00012,450.3
BRCA1,3.42,0.000045,892.1
MYC,1.83,0.0082,2103.5
EGFR,-0.91,0.048,310.7
VEGFA,2.71,0.00031,678.4
CDK4,1.21,0.019,543.2
PTEN,-1.54,0.0075,189.6
CDKN2A,-2.33,0.00089,234.5
PIK3CA,1.67,0.0041,445.8
ERBB2,2.95,0.00018,567.9`,
      },
      {
        icon: "🧪",
        title: "Exploratory Analysis",
        description: "Automatic EDA pipeline on any dataset",
        examplePrompt: `Run an exploratory data analysis on this clinical + expression dataset — describe the distributions, correlations, and any interesting patterns:

sample,age,BMI,tumor_stage,TP53_expr,MYC_expr,BRCA1_expr,outcome
P001,45,24.3,II,120,2340,89,survived
P002,62,28.7,III,89,2890,76,deceased
P003,38,22.1,I,145,2101,92,survived
P004,71,31.2,IV,98,2567,85,deceased
P005,55,26.8,II,132,1890,98,survived
P006,48,25.4,III,112,3201,71,deceased
P007,43,23.9,I,167,1567,112,survived
P008,66,29.1,IV,76,3456,68,deceased`,
      },
      {
        icon: "📈",
        title: "GO / Pathway Enrichment",
        description: "Enrich gene lists against GO or KEGG",
        examplePrompt:
          "Run GO enrichment analysis on these significantly upregulated genes from a cancer RNA-seq study: TP53, MYC, BRCA1, EGFR, VEGFA, CDK4, PIK3CA, AKT1, KRAS, ERBB2, CDKN2A, PTEN, RB1, MDM2, CCND1, BCL2, NOTCH1, WNT5A, STAT3, NF1",
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
