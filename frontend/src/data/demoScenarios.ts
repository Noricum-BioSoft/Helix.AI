/**
 * Helix.AI Demo Scenarios
 *
 * Five realistic bioinformatics prompts used for investor and user demonstrations.
 * Each scenario exercises a different tool or workflow, covering both
 * "needs_inputs" responses (no data files provided) and multi-step pipeline
 * execution (S3 paths supplied inline).
 */

export type ExpectedBehavior = 'needs_inputs' | 'executes_pipeline' | 'multi_step_plan';

export interface DemoOutput {
  label: string;
  type: 'plot' | 'csv' | 'fasta' | 'report' | 'newick';
}

/** Required or optional inputs the user must supply (e.g. S3 paths, parameters). */
export interface DemoInput {
  label: string;
  description?: string;
}

export interface DataPreviewTable {
  title: string;
  headers: string[];
  rows: string[][];
  note?: string;
}

export interface DemoScenario {
  id: string;
  title: string;
  subtitle: string;
  domain: string;
  domainColor: string;
  icon: string;
  tags: string[];
  expectedBehavior: ExpectedBehavior;
  behaviorLabel: string;
  tool: string;
  estimatedRuntime: string;
  /** Inputs the scenario needs (e.g. count matrix, sample metadata, S3 paths). */
  inputs: DemoInput[];
  outputs: DemoOutput[];
  prompt: string;
  /**
   * For `needs_inputs` scenarios: a ready-to-send follow-up message that
   * supplies the required data paths from the demo S3 bucket.
   */
  followUpPrompt?: string;
  /** Optional inline preview tables shown in the "Preview data" modal. */
  dataPreview?: DataPreviewTable[];
}

export const demoScenarios: DemoScenario[] = [
  // ── 1. Bulk RNA-seq · 2×2 factorial ─────────────────────────────────────────
  {
    id: 'bulk-rnaseq-factorial',
    title: 'Bulk RNA-seq: Infection × Time Factorial',
    subtitle: 'T. gondii infection in mouse brain — 2×2 design, DESeq2',
    domain: 'Bulk RNA-seq',
    domainColor: '#3A60A8',
    icon: '🧫',
    tags: ['Factorial Design', 'DESeq2', 'Interaction Term', 'Mouse'],
    expectedBehavior: 'needs_inputs',
    behaviorLabel: 'Requests Data',
    tool: 'bulk_rnaseq_analysis',
    estimatedRuntime: '5–10 min (with data)',
    inputs: [
      { label: 'Count matrix (CSV)', description: 'Gene × sample raw counts, e.g. S3 path' },
      { label: 'Sample metadata (CSV)', description: 'Sample IDs, factors (infection_status, time_point)' },
      { label: 'Design formula', description: 'e.g. ~infection_status + time_point + interaction' },
    ],
    outputs: [
      { label: 'Volcano plots (per contrast)', type: 'plot' },
      { label: 'PCA + sample-distance heatmap', type: 'plot' },
      { label: 'DE results tables', type: 'csv' },
    ],
    prompt: `You are analyzing an RNA-seq transcriptome dataset from a mouse study investigating the effects of Toxoplasma gondii infection on brain gene expression.

Study Design
This is a 2 × 2 factorial experimental design with the following factors:

Factor 1: Infection Status
  Infected
  Uninfected

Factor 2: Time Point
  11 days post infection (11 dpi)
  33 days post infection (33 dpi)

Each of the four experimental groups contains 3 biological replicates, for a total of 12 samples.

Objectives
  Model gene expression changes using an appropriate statistical framework for a two-factor design.
  Assess:
    The main effect of infection status
    The main effect of time point
    The interaction effect between infection status and time point
  Perform differential expression analysis.
  Apply appropriate normalization, statistical testing, and multiple testing correction.
  Generate clear, reproducible code (e.g., in R with DESeq2/edgeR or Python with appropriate libraries).
  Include exploratory data analysis (PCA, clustering, sample distance heatmap).
  Clearly structure outputs (results tables, contrasts, visualizations).

Assume raw count data and sample metadata are available.
Implement best practices for RNA-seq analysis.`,
    followUpPrompt: `Run bulk RNA-seq differential expression analysis with these inputs:
count_matrix: s3://noricum-ngs-data/demo/rnaseq/tgondii_counts.csv
sample_metadata: s3://noricum-ngs-data/demo/rnaseq/tgondii_metadata.csv
design_formula: ~infection_status + time_point + infection_status:time_point`,
    dataPreview: [
      {
        title: 'Count Matrix (tgondii_counts.csv)',
        headers: ['gene', 'inf_11dpi_r1', 'inf_11dpi_r2', 'inf_11dpi_r3', 'inf_33dpi_r1', 'inf_33dpi_r2', 'uninf_11dpi_r1', 'uninf_11dpi_r2', 'uninf_33dpi_r1', '…'],
        rows: [
          ['Gene0702', '1748', '1409', '1467', '1603', '1686', '1703', '1540', '1263', '…'],
          ['Gene0106', '1419', '1478', '1886', '2028', '1455', '1656', '1444', '1395', '…'],
          ['Gene0019', '1261', '1240', '1278', '1216', '1383', '1061', '1189', '1096', '…'],
          ['Gene0291', '1532', '1612', '1489', '1601', '1558', '1498', '1521', '1477', '…'],
          ['Gene0413', '1387', '1301', '1424', '1356', '1412', '1342', '1398', '1319', '…'],
        ],
        note: '250 genes × 12 samples  ·  s3://noricum-ngs-data/demo/rnaseq/tgondii_counts.csv',
      },
      {
        title: 'Sample Metadata (tgondii_metadata.csv)',
        headers: ['sample', 'infection_status', 'time_point', 'group'],
        rows: [
          ['inf_11dpi_rep1', 'infected', '11dpi', 'infected_11dpi'],
          ['inf_11dpi_rep2', 'infected', '11dpi', 'infected_11dpi'],
          ['inf_11dpi_rep3', 'infected', '11dpi', 'infected_11dpi'],
          ['inf_33dpi_rep1', 'infected', '33dpi', 'infected_33dpi'],
          ['inf_33dpi_rep2', 'infected', '33dpi', 'infected_33dpi'],
          ['uninf_11dpi_rep1', 'uninfected', '11dpi', 'uninfected_11dpi'],
          ['uninf_11dpi_rep2', 'uninfected', '11dpi', 'uninfected_11dpi'],
          ['uninf_33dpi_rep1', 'uninfected', '33dpi', 'uninfected_33dpi'],
          ['…', '…', '…', '…'],
        ],
        note: '12 samples  ·  2×2 factorial: 2 infection states × 2 time points × 3 replicates',
      },
    ],
  },

  // ── 2. Single-cell RNA-seq · SLE PBMC ────────────────────────────────────────
  {
    id: 'scrna-sle-pbmc',
    title: 'Single-Cell RNA-seq: SLE Immune Profiling',
    subtitle: 'Human PBMC study — lupus vs. healthy, cell-type DEG',
    domain: 'Single-Cell',
    domainColor: '#7B3FA8',
    icon: '🔬',
    tags: ['10x Genomics', 'UMAP', 'Pseudobulk', 'Immune Cells', 'Human'],
    expectedBehavior: 'needs_inputs',
    behaviorLabel: 'Requests Data',
    tool: 'single_cell_analysis',
    estimatedRuntime: '15–30 min (with data)',
    inputs: [
      { label: 'Gene–cell matrix (H5)', description: '10x HDF5 or similar; S3 path' },
      { label: 'Data format', description: 'e.g. 10x' },
      { label: 'Resolution (optional)', description: 'Clustering resolution, e.g. 0.5' },
      { label: 'Steps (optional)', description: 'e.g. all, or subset' },
    ],
    outputs: [
      { label: 'UMAP: cluster, cell type, disease', type: 'plot' },
      { label: 'Marker gene dot plot', type: 'plot' },
      { label: 'DEGs per cell type', type: 'csv' },
      { label: 'Cell-type proportion chart', type: 'plot' },
    ],
    prompt: `You are analyzing a single-cell RNA-seq dataset from a human peripheral blood mononuclear cell (PBMC) study investigating immune cell dysregulation in patients with systemic lupus erythematosus (SLE) compared to healthy controls.

Dataset
  Single-cell RNA-seq (10x Genomics Chromium v3) gene-expression matrix.
  Estimated 8,000–12,000 cells per sample; 5 SLE patients, 5 healthy donors.

Study Design
  Two-group comparison: SLE (n=5) vs. Healthy (n=5).
  Cells span major PBMC lineages: T cells, B cells, NK cells, monocytes, pDCs.

Objectives
  1. Perform quality control (mitochondrial gene %, doublet detection, UMI distribution)
     and filter low-quality cells.
  2. Normalize counts (SCTransform), select highly variable genes, and reduce
     dimensionality (PCA → UMAP).
  3. Cluster cells (Leiden algorithm, resolution 0.5) and annotate major cell types
     using canonical marker genes.
  4. Identify cell-type-specific differentially expressed genes between SLE and
     healthy donors using a pseudobulk approach (DESeq2 per cell type).
  5. Quantify changes in cell-type composition (cell-type frequency per donor).

Desired Outputs
  - UMAP plots colored by cluster, cell type, and disease status (PNG).
  - Dot plot of top marker genes per cluster.
  - CSV tables of DEGs per cell type (log2FC, adjusted p-value, mean expression).
  - Stacked bar chart of cell-type proportions per sample.`,
    followUpPrompt: `data_file: s3://noricum-ngs-data/demo/scrna/sle_pbmc_filtered_feature_bc_matrix.h5
data_format: 10x
resolution: 0.5
steps: all`,
    dataPreview: [
      {
        title: 'Gene–Cell Matrix (excerpt, UMI counts)',
        headers: ['barcode', 'CD3D', 'CD3E', 'MS4A1', 'CD14', 'FCGR3A', 'GNLY', 'IL7R', '…'],
        rows: [
          ['AAACCTGAGAAACCAT-1', '12', '9', '0', '0', '1', '0', '14', '…'],
          ['AAACCTGAGAAACCGC-1', '0', '1', '18', '0', '0', '0', '2', '…'],
          ['AAACCTGAGAAACCTA-1', '1', '0', '0', '34', '2', '0', '0', '…'],
          ['AAACCTGAGAAACGAG-1', '0', '0', '0', '0', '0', '41', '0', '…'],
          ['AAACCTGAGAAACGCC-1', '8', '11', '0', '1', '0', '0', '22', '…'],
        ],
        note: '~10,000 cells × 33,538 genes  ·  10x Genomics Chromium v3  ·  s3://noricum-ngs-data/demo/scrna/sle_pbmc_filtered_feature_bc_matrix.h5',
      },
      {
        title: 'Cell Metadata (donor & disease status)',
        headers: ['barcode', 'donor_id', 'disease', 'n_genes', 'n_umis', 'pct_mito'],
        rows: [
          ['AAACCTGAGAAACCAT-1', 'SLE_01', 'SLE', '842', '2314', '2.1%'],
          ['AAACCTGAGAAACCGC-1', 'SLE_01', 'SLE', '631', '1189', '1.8%'],
          ['AAACCTGAGAAACCTA-1', 'HC_01',  'Healthy', '714', '1843', '1.4%'],
          ['AAACCTGAGAAACGAG-1', 'HC_02',  'Healthy', '508', '982', '3.2%'],
          ['…', '…', '…', '…', '…', '…'],
        ],
        note: '5 SLE donors + 5 healthy controls  ·  2-group comparison: SLE vs. Healthy',
      },
    ],
  },

  // ── 3. Amplicon QC pipeline · with S3 paths ──────────────────────────────────
  {
    id: 'amplicon-qc-pipeline',
    title: 'Amplicon QC Pipeline: 16S Gut Microbiome',
    subtitle: 'IBD study — FastQC → trim → merge → quality report',
    domain: 'Sequencing QC',
    domainColor: '#2A8A5A',
    icon: '⚙️',
    tags: ['16S rRNA', 'FastQC', 'Trimming', 'Read Merging', 'Pipeline'],
    expectedBehavior: 'multi_step_plan',
    behaviorLabel: 'Multi-Step Plan',
    tool: 'quality_assessment',
    estimatedRuntime: '2–5 min',
    inputs: [
      { label: 'Forward reads (FASTQ)', description: 'S3 path to R1, e.g. s3://bucket/sample_R1.fastq.gz' },
      { label: 'Reverse reads (FASTQ)', description: 'S3 path to R2' },
      { label: 'Output prefix (optional)', description: 'S3 prefix for trimmed/merged outputs' },
    ],
    outputs: [
      { label: 'FastQC HTML reports', type: 'report' },
      { label: 'Merged amplicon FASTA', type: 'fasta' },
      { label: 'Quality summary (per-step stats)', type: 'csv' },
    ],
    prompt: `You are processing a 16S rRNA amplicon sequencing dataset from a gut microbiome study. Raw paired-end FASTQ files are on S3 and need a full preprocessing pipeline before downstream diversity analysis.

Dataset
  Illumina MiSeq 2×250 bp paired-end reads; V3–V4 hypervariable region.
  Forward reads: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq
  Reverse reads: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq
  Output prefix:  s3://noricum-ngs-data/test-output/amplicon-demo/

Study Design
  Case-control: 20 IBD patients vs. 20 healthy controls.
  Single sample shown here (sample01); pipeline will be applied to all 40.

Pipeline Steps
  1. Run FastQC quality assessment on both raw R1 and R2 files.
  2. Trim adapter sequences (CTGTCTCTTATACACATCT) and low-quality bases
     (Phred < 20) from both ends; minimum read length 150 bp.
  3. Merge overlapping paired-end reads with minimum overlap of 20 bp.
  4. Generate a quality report summarizing read counts before and after each step.

Desired Outputs
  - FastQC HTML reports for raw R1 and R2.
  - Trimmed FASTQ files saved to the output S3 prefix.
  - Merged FASTA file of consensus amplicon sequences.
  - Quality summary CSV (sample, raw_reads, post_trim_reads, merged_reads, merge_rate).`,
    dataPreview: [
      {
        title: 'Sample Manifest',
        headers: ['sample_id', 'R1_path', 'R2_path', 'group'],
        rows: [
          ['sample01', 's3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq', 's3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq', 'IBD'],
          ['sample02', 's3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq', 's3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq', 'Healthy'],
          ['…', '…', '…', '…'],
        ],
        note: '40 samples · 20 IBD patients + 20 healthy controls · 16S V3–V4 · MiSeq 2×250 bp',
      },
    ],
  },

  // ── 4. Longitudinal bulk RNA-seq · APAP liver injury ─────────────────────────
  {
    id: 'bulk-rnaseq-timecourse',
    title: 'Bulk RNA-seq Workflow Design: Liver Injury Recovery',
    subtitle: 'Ask Helix to propose the time-course analysis plan + outputs',
    domain: 'Workflow Design',
    domainColor: '#3A60A8',
    icon: '🧩',
    tags: ['Bulk RNA-seq', 'Time-Course', 'Workflow Planning', 'LRT', 'Figures'],
    expectedBehavior: 'multi_step_plan',
    behaviorLabel: 'Multi-Step Plan',
    tool: 'handle_natural_command',
    estimatedRuntime: '1–3 min (planning)',
    inputs: [
      { label: 'Study objective', description: 'What biological progression/recovery signal to quantify' },
      { label: 'Available data', description: 'Count matrix + sample metadata summary' },
      { label: 'Design constraints (optional)', description: 'Statistical method, runtime, reproducibility needs' },
    ],
    outputs: [
      { label: 'Proposed end-to-end workflow', type: 'report' },
      { label: 'Expected calculations and result tables', type: 'csv' },
      { label: 'Recommended plots/figures', type: 'plot' },
    ],
    prompt: `I have a bulk RNA-seq time-course study from a murine acetaminophen liver-injury model and I want you to design the analysis workflow before execution.

Study Context
  - 5 post-injury time points: 0 h, 6 h, 24 h, 72 h, 168 h
  - 4 biological replicates per time point (20 samples total)
  - Raw count matrix and sample metadata are available

Task
  1. Propose a complete, publication-ready workflow from count-level QC to biological interpretation.
  2. For each stage, specify:
     - key calculations/statistical tests,
     - expected tables/artifacts,
     - recommended plots/figures.
  3. Include explicit checkpoints for data quality and model validity.
  4. Recommend how to identify early-response, late-recovery, and persistent signatures.

Desired deliverable format
  - Ordered workflow steps
  - Calculations/methods per step
  - Output artifacts (tables + plots) per step
  - Brief rationale and risk notes`,
    dataPreview: [
      {
        title: 'Count Matrix (apap_timecourse_counts.csv)',
        headers: ['gene', '0h_rep1', '0h_rep2', '6h_rep1', '6h_rep2', '24h_rep1', '72h_rep1', '168h_rep1', '…'],
        rows: [
          ['Cyp2e1', '3241', '3187', '412', '389', '201', '1876', '3102', '…'],
          ['Hmox1',  '412',  '398',  '3842', '3967', '4201', '2103', '520', '…'],
          ['Lcn2',   '89',   '94',   '5621', '5489', '6102', '1234', '156', '…'],
          ['Alb',    '8921', '9103', '7234', '7056', '4210', '7891', '8742', '…'],
          ['Il6',    '45',   '52',   '2341', '2567', '3102', '891',  '67',  '…'],
        ],
        note: '1000 genes × 20 samples  ·  s3://noricum-ngs-data/demo/rnaseq/apap_timecourse_counts.csv',
      },
      {
        title: 'Sample Metadata (apap_timecourse_metadata.csv)',
        headers: ['sample', 'time_point', 'replicate', 'hours_post_APAP'],
        rows: [
          ['0h_rep1',  '0h',   '1', '0'],
          ['0h_rep2',  '0h',   '2', '0'],
          ['6h_rep1',  '6h',   '1', '6'],
          ['6h_rep2',  '6h',   '2', '6'],
          ['24h_rep1', '24h',  '1', '24'],
          ['72h_rep1', '72h',  '1', '72'],
          ['168h_rep1','168h', '1', '168'],
          ['…', '…', '…', '…'],
        ],
        note: '20 samples  ·  5 time points × 4 replicates  ·  APAP overdose model',
      },
    ],
  },

  // ── 5. Phylogenetic analysis · SARS-CoV-2 spike ──────────────────────────────
  {
    id: 'phylogenetics-sarscov2',
    title: 'Phylogenetics: SARS-CoV-2 Variant Divergence',
    subtitle: 'Spike protein across 8 VOCs — alignment + ML tree',
    domain: 'Phylogenetics',
    domainColor: '#A85A3A',
    icon: '🌳',
    tags: ['SARS-CoV-2', 'Spike Protein', 'ML Tree', 'NCBI', 'Variants'],
    expectedBehavior: 'multi_step_plan',
    behaviorLabel: 'Multi-Step Plan',
    tool: 'phylogenetic_tree',
    estimatedRuntime: '3–8 min',
    inputs: [
      { label: 'Sequences (FASTA or names)', description: 'User can paste FASTA or name variants; Helix can fetch from NCBI' },
      { label: 'Algorithm (optional)', description: 'Alignment (e.g. MAFFT), tree (e.g. ML with bootstrap)' },
    ],
    outputs: [
      { label: 'Multiple sequence alignment', type: 'fasta' },
      { label: 'Annotated phylogenetic tree', type: 'plot' },
      { label: 'Tree in Newick format', type: 'newick' },
      { label: 'Pairwise identity matrix', type: 'csv' },
    ],
    prompt: `You are conducting a comparative evolutionary analysis of the SARS-CoV-2 spike protein across major variants of concern (VOCs) to characterize mutational divergence from the ancestral Wuhan-Hu-1 reference sequence.

Dataset
  Full-length spike protein amino acid sequences for 8 SARS-CoV-2 variants:
  Wuhan-Hu-1 (reference), Alpha (B.1.1.7), Beta (B.1.351), Gamma (P.1),
  Delta (B.1.617.2), Omicron BA.1, Omicron BA.4/5, and XBB.1.5.
  Use this curated GenBank-allowed accession set for reproducible retrieval:
    Wuhan-Hu-1: MN908947.3
    Alpha (B.1.1.7): OQ898928.1
    Beta (B.1.351): OR353131.1
    Gamma (P.1): MW642250.1
    Delta (B.1.617.2): OR323381.1
    Omicron BA.1: PP847536.1
    Omicron BA.4/5: PP848071.1
    XBB.1.5: PP405604.1

Study Design
  Comparative sequence analysis — no experimental replicates.
  Unit of analysis: one representative consensus sequence per variant.

Objectives
  1. Retrieve spike protein sequences from NCBI for each variant (or accept
     user-provided FASTA).
  2. Perform multiple sequence alignment (MAFFT L-INS-i algorithm).
  3. Reconstruct a maximum-likelihood phylogenetic tree with 1,000 bootstrap
     replicates.
  4. Annotate key mutation sites on the tree (RBD mutations: K417N, E484K/A,
     N501Y, L452R; Furin cleavage site mutations).
  5. Calculate pairwise amino acid identity matrix across all variant pairs.

Desired Outputs
  - Multiple sequence alignment FASTA file.
  - Phylogenetic tree in Newick format with bootstrap support values.
  - Visualization of the annotated tree (PNG, rectangular cladogram style).
  - Pairwise identity matrix as CSV (variant × variant, % amino acid identity).
  - Table of key mutation presence/absence per variant.`,
    followUpPrompt: `Build a phylogenetic tree and pairwise amino acid identity matrix for the SARS-CoV-2 spike protein across 8 variants of concern using these pre-aligned sequences:
sequences: s3://noricum-ngs-data/demo/phylo/sars_cov2_spike.fasta
variants: Wuhan-Hu-1, Alpha B.1.1.7, Beta B.1.351, Gamma P.1, Delta B.1.617.2, Omicron BA.1, Omicron BA.4/5, XBB.1.5`,
    dataPreview: [
      {
        title: 'Variants to Analyse',
        headers: ['Variant', 'Lineage', 'NCBI Nucleotide Accession (Curated Demo Set)', 'Key RBD Mutations', 'First Detected'],
        rows: [
          ['Wuhan-Hu-1', 'B', 'MN908947.3', '— (reference)', 'Dec 2019'],
          ['Alpha', 'B.1.1.7', 'OQ898928.1', 'N501Y, P681H, Δ69-70', 'Sep 2020 (UK)'],
          ['Beta', 'B.1.351', 'OR353131.1', 'K417N, E484K, N501Y', 'May 2020 (SA)'],
          ['Gamma', 'P.1', 'MW642250.1', 'K417T, E484K, N501Y', 'Nov 2020 (Brazil)'],
          ['Delta', 'B.1.617.2', 'OR323381.1', 'L452R, T478K, P681R', 'Oct 2020 (India)'],
          ['Omicron BA.1', 'BA.1', 'PP847536.1', 'K417N, E484A, N501Y (+30)', 'Nov 2021 (SA)'],
          ['Omicron BA.4/5', 'BA.4/5', 'PP848071.1', 'L452R, F486V, R493Q', 'Jan 2022 (SA)'],
          ['XBB.1.5', 'XBB.1.5', 'PP405604.1', 'F486P (recombinant)', 'Oct 2022 (USA)'],
        ],
        note: 'Pre-aligned spike protein sequences  ·  s3://noricum-ngs-data/demo/phylo/sars_cov2_spike.fasta',
      },
    ],
  },
];

export const getDemoScenarioById = (id: string): DemoScenario | undefined =>
  demoScenarios.find((s) => s.id === id);

export const getDemoScenarioByTool = (toolName: string): DemoScenario | undefined =>
  demoScenarios.find((s) => {
    if (!s.followUpPrompt) return false;
    return s.tool === toolName || s.tool.includes(toolName);
  });

/**
 * Match a demo scenario by command text and tool name.
 * Used when scenarioId is not available — checks command keywords against
 * each scenario's distinctive terms to avoid returning the wrong scenario
 * when multiple demos share the same tool (e.g. Demo 1 and Demo 4 both use
 * bulk_rnaseq_analysis).
 *
 * Each demo gets a small list of unique keywords extracted from its prompt.
 */
export const getDemoScenarioByCommandAndTool = (
  command: string,
  toolName: string,
): DemoScenario | undefined => {
  const cmdLower = command.toLowerCase();

  // Keywords unique to each demo scenario (used for disambiguation)
  const KEYWORDS: Record<string, string[]> = {
    'bulk-rnaseq-factorial':     ['toxoplasma', 'tgondii', 't. gondii', 'gondii', 'infection'],
    'bulk-rnaseq-timecourse':    ['apap', 'acetaminophen', 'liver injury', 'time-course', 'timecourse', 'time course'],
    'scrna-sle-pbmc':            ['lupus', 'sle', 'pbmc', 'single-cell', 'single cell', '10x genomics', 'umap', 'leiden'],
    'amplicon-qc-pipeline':      ['16s', 'amplicon', 'microbiome', 'fastqc', 'gut'],
    'phylogenetics-sarscov2':    ['sars', 'covid', 'spike protein', 'phylogen', 'wuhan', 'omicron', 'variant'],
  };

  // Check each scenario whose tool matches
  for (const scenario of demoScenarios) {
    if (!scenario.followUpPrompt) continue;
    if (scenario.tool !== toolName && !scenario.tool.includes(toolName)) continue;

    const kws = KEYWORDS[scenario.id];
    if (kws && kws.some((kw) => cmdLower.includes(kw))) {
      return scenario;
    }
  }

  // Fall back to the first matching tool scenario
  return getDemoScenarioByTool(toolName);
};

/**
 * Match a demo scenario by command text only (no tool).
 * Use when the backend returns tool "__plan__" so we can still show "Load & run"
 * and data preview for known demo prompts.
 */
export const getDemoScenarioByCommand = (command: string): DemoScenario | undefined => {
  const cmdLower = command.toLowerCase();
  const KEYWORDS: Record<string, string[]> = {
    'bulk-rnaseq-factorial':     ['toxoplasma', 'tgondii', 't. gondii', 'gondii', 'infection'],
    'bulk-rnaseq-timecourse':    ['apap', 'acetaminophen', 'liver injury', 'time-course', 'timecourse', 'time course'],
    'scrna-sle-pbmc':            ['lupus', 'sle', 'pbmc', 'single-cell', 'single cell', '10x genomics', 'umap', 'leiden'],
    'amplicon-qc-pipeline':      ['16s', 'amplicon', 'microbiome', 'fastqc', 'gut'],
    'phylogenetics-sarscov2':    ['sars', 'covid', 'spike protein', 'phylogen', 'wuhan', 'omicron', 'variant'],
  };
  for (const scenario of demoScenarios) {
    if (!scenario.followUpPrompt) continue;
    const kws = KEYWORDS[scenario.id];
    if (kws && kws.some((kw) => cmdLower.includes(kw))) return scenario;
  }
  return undefined;
};
