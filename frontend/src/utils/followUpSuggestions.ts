/**
 * Contextual follow-up suggestions.
 *
 * Maps the tool name returned by the backend to a short list of
 * natural-language next-step prompts the user can click to continue.
 * Each suggestion becomes a clickable chip below the response card.
 */

export interface FollowUpSuggestion {
  label: string;
  prompt: string;
  /** Emoji icon shown on the chip */
  icon: string;
}

const SUGGESTIONS: Record<string, FollowUpSuggestion[]> = {
  sequence_alignment: [
    { icon: "🌳", label: "Build phylogenetic tree", prompt: "Build a phylogenetic tree from this alignment" },
    { icon: "🔬", label: "Cluster sequences", prompt: "Cluster these sequences and get representative members" },
    { icon: "📥", label: "Fetch more from NCBI", prompt: "Fetch the reference sequence from NCBI for comparison" },
  ],
  fetch_ncbi_sequence: [
    { icon: "🧬", label: "Align sequences", prompt: "Align these sequences" },
    { icon: "🌳", label: "Build phylogenetic tree", prompt: "Build a phylogenetic tree from these sequences" },
    { icon: "🔬", label: "Find variants", prompt: "Show me variants relative to the reference" },
  ],
  phylogenetic_tree: [
    { icon: "🧬", label: "Align sequences", prompt: "Align the sequences used for this tree" },
    { icon: "🔬", label: "Cluster sequences", prompt: "Cluster these sequences and get representative members" },
    { icon: "📊", label: "Export tree data", prompt: "What is the structure of the phylogenetic tree?" },
  ],
  clustering_analysis: [
    { icon: "🌳", label: "Build phylogenetic tree", prompt: "Build a phylogenetic tree from these representatives" },
    { icon: "🧬", label: "Align representatives", prompt: "Align the representative sequences" },
    { icon: "🔢", label: "Select variants", prompt: "Select the top 5 diverse variants from the clusters" },
  ],
  bulk_rnaseq_analysis: [
    { icon: "🗺️", label: "Update volcano plot", prompt: "Update the volcano plot — highlight the top 10 significant genes in red" },
    { icon: "🔄", label: "Adjust parameters", prompt: "Re-run with a stricter threshold: alpha=0.01" },
    { icon: "📋", label: "GO enrichment", prompt: "Run GO enrichment analysis on the significant genes" },
    { icon: "📊", label: "Compare runs", prompt: "What changed between my runs?" },
  ],
  single_cell_analysis: [
    { icon: "🗺️", label: "Update UMAP", prompt: "Update the UMAP — color cells by a different gene" },
    { icon: "🔄", label: "Change resolution", prompt: "Re-run clustering with a higher resolution" },
    { icon: "📊", label: "Marker genes", prompt: "What are the top marker genes for each cluster?" },
    { icon: "📈", label: "Differential expression", prompt: "Run differential expression between cluster 0 and cluster 1" },
  ],
  fastqc_quality_analysis: [
    { icon: "✂️", label: "Trim reads", prompt: "Trim low-quality bases and adapters from my reads" },
    { icon: "📊", label: "Merge reads", prompt: "Merge the paired-end reads" },
    { icon: "🧬", label: "Align reads", prompt: "Align trimmed reads to the reference genome" },
  ],
  read_trimming: [
    { icon: "🔍", label: "Check quality again", prompt: "Run FastQC on the trimmed reads to confirm quality" },
    { icon: "📊", label: "Merge reads", prompt: "Merge the trimmed paired-end reads" },
    { icon: "🧬", label: "Align reads", prompt: "Align trimmed reads to the reference genome" },
  ],
  tabular_analysis: [
    { icon: "❓", label: "Ask a question", prompt: "What is the average value in the top-ranked rows?" },
    { icon: "📊", label: "Plot results", prompt: "Create a bar chart of the top 10 rows" },
    { icon: "🔄", label: "Different operation", prompt: "Filter rows where the p-value is below 0.05" },
  ],
  tabular_qa: [
    { icon: "🔢", label: "Sort data", prompt: "Sort the table by the most important column" },
    { icon: "📊", label: "Plot distribution", prompt: "Plot the distribution of values in the key column" },
    { icon: "🔎", label: "Filter rows", prompt: "Filter to show only rows where the value is above average" },
  ],
  ds_run_analysis: [
    { icon: "❓", label: "Ask about results", prompt: "What are the most important features in my dataset?" },
    { icon: "📊", label: "Plot feature importance", prompt: "Plot the feature importance from the analysis" },
    { icon: "🔄", label: "Different target", prompt: "Re-run the analysis targeting a different column" },
  ],
  patch_and_rerun: [
    { icon: "🔄", label: "Another change", prompt: "Make another change to the plot" },
    { icon: "📊", label: "Compare runs", prompt: "What changed between my runs?" },
    { icon: "💾", label: "Summarize results", prompt: "Summarize what the current plot shows" },
  ],
  bio_rerun: [
    { icon: "📊", label: "Compare runs", prompt: "What changed between my runs?" },
    { icon: "🗺️", label: "Update visualization", prompt: "Update the volcano plot with the new results" },
    { icon: "📋", label: "Summarize changes", prompt: "Summarize how the results changed after the parameter update" },
  ],
  bio_diff_runs: [
    { icon: "🔄", label: "Adjust parameters", prompt: "Re-run with different parameters to investigate further" },
    { icon: "📊", label: "Plot differences", prompt: "Plot the genes that changed most between the runs" },
  ],
  go_enrichment_analysis: [
    { icon: "📊", label: "Plot top pathways", prompt: "Plot the top 15 enriched pathways as a bar chart" },
    { icon: "🔍", label: "Filter by category", prompt: "Show only biological process terms" },
  ],
};

/** Fallback suggestions shown when the tool has no specific mapping. */
const DEFAULT_SUGGESTIONS: FollowUpSuggestion[] = [
  { icon: "🔄", label: "Refine results", prompt: "Refine the previous results" },
  { icon: "❓", label: "Explain results", prompt: "Explain what these results mean" },
  { icon: "📋", label: "What's next?", prompt: "What should I do next with these results?" },
];

/**
 * Return contextual follow-up suggestions for a given tool.
 * Returns an empty array for tools where no follow-up makes sense
 * (e.g. toolbox_inventory, session_run_io_summary).
 */
export function getFollowUpSuggestions(tool: string | undefined): FollowUpSuggestion[] {
  if (!tool) return [];
  const noFollowUp = new Set([
    "toolbox_inventory",
    "session_run_io_summary",
    "handle_natural_command",
    "__plan__",
    "local_demo_plot_script",
    "local_demo_scatter_plot",
  ]);
  if (noFollowUp.has(tool)) return [];
  return SUGGESTIONS[tool] ?? DEFAULT_SUGGESTIONS;
}

/**
 * Return a context-aware input placeholder string based on the last tool used.
 */
export function getContextualPlaceholder(lastTool: string | undefined): string {
  const placeholders: Record<string, string> = {
    sequence_alignment: "Try: build a phylogenetic tree from this alignment…",
    fetch_ncbi_sequence: "Try: align these sequences…",
    phylogenetic_tree: "Try: cluster the sequences and get representatives…",
    clustering_analysis: "Try: build a phylogenetic tree from the representatives…",
    bulk_rnaseq_analysis: "Try: update the volcano plot or re-run with alpha=0.01…",
    single_cell_analysis: "Try: update the UMAP or change the clustering resolution…",
    fastqc_quality_analysis: "Try: trim low-quality reads and adapters…",
    read_trimming: "Try: run FastQC again to verify quality…",
    tabular_analysis: "Try: ask a question about the results…",
    tabular_qa: "Try: sort or filter the data…",
    bio_rerun: "Try: compare runs to see what changed…",
    patch_and_rerun: "Try: make another change to the visualization…",
    go_enrichment_analysis: "Try: plot the top 15 enriched pathways…",
  };
  return placeholders[lastTool ?? ""] ?? "Ask anything, describe your data, or upload a file…";
}
