# Helix.AI — User Scenarios for Scientists and Bioinformaticians

This document defines **representative scenarios** that scientists and bioinformaticians may run when using Helix.AI. The focus is on **reproducibility**, **moving code into production quickly**, and **correctness of results**. Each scenario is described so it can drive design decisions, acceptance criteria, and test cases.

---

## Domain context: Bioinformatics and data analysis

Bioinformatics and data analysis involve using computational tools, statistics, and algorithms to interpret complex biological data—such as DNA sequences, RNA-seq, and proteomics—to drive insights in genomics, drug discovery, and diagnostics. The field bridges biology and computer science, transforming raw, high-throughput sequencing data into actionable insights.

**Key aspects**

- **Data types:** Genomic, transcriptomic (RNA-seq), proteomic, and metabolomic data.
- **Workflow:** Analysis generally follows a pipeline: *primary* (base calling), *secondary* (alignment/assembly), and *tertiary* (interpretation and functional analysis).
- **Techniques:** Sequence alignment (e.g. BLAST), machine learning for predictive modeling, and statistical analysis for identifying patterns.

**Applications**

- **Clinical diagnostics:** Identifying genetic variants (SNPs, INDELs) for disease risk assessment (e.g. Alzheimer’s, cancer).
- **Drug discovery:** Identifying novel drug targets and biomarker discovery.
- **Evolutionary biology:** Tracking mutations and species evolution.

**Tools and skills**

- **Languages:** Python and R for data manipulation and visualization.
- **Platforms/software:** Linux environments, BLAST, Ensembl, DNAnexus, and similar.
- **Methods:** Handling large-scale datasets, NGS workflow management, and building reproducible, automated pipelines.

Bioinformatics is essential in modern medicine for personalized, patient-centered care and precision medicine, as it provides the analytical infrastructure for understanding complex biological systems. The scenarios in this document are chosen so that Helix.AI can support this breadth of tasks—from exploration and first-run analysis to iterative refinement, validation, and production handoff—while preserving correctness and reproducibility.

---

## 1. User Personas and Goals

| Persona | Primary goal | Secondary goals |
|--------|---------------|-----------------|
| **Wet-lab scientist** | Answer a biological question from data without writing code or configuring pipelines. | Get interpretable figures and tables; share a reproducible bundle with a collaborator or core facility. |
| **Bioinformatician** | Run a standard analysis fast, then iterate (parameters, design, visualizations) without re-running from scratch. | Export a clean script and bundle for version control, papers, or production runs on new cohorts. |
| **Core facility / platform team** | Offer a single interface for multiple assay types (bulk RNA-seq, scRNA-seq, QC) with consistent outputs. | Ensure every run is auditable and reproducible; hand off `analysis.py` or bundle to the lab for re-runs. |

**Cross-cutting requirements**
- **Correctness**: Results must come from real tool execution (ExecutionBroker), not LLM-generated prose. Parameters and design must be applied consistently.
- **Reproducibility**: Every run produces `analysis.py` and (optionally) `bundle.zip` so the same result can be re-produced without Helix.AI.
- **Production handoff**: The generated script has a clear parameter block; changing one value and re-running must be straightforward (e.g., new cohort, new alpha).

---

## 2. Scenario Categories

Scenarios are grouped by **user intent** and **session pattern**:

- **A. Exploration & tool discovery** — User asks what is possible or what inputs are needed.
- **B. First run (supply data)** — User provides data (paths or example data); system runs the analysis once.
- **C. Iterative refinement** — User asks to change a parameter or visualization; system patches the script and re-runs (no full re-prompt).
- **D. Validation & interpretation** — User asks about the results (Q&A) or requests a different view (e.g., axis scale, threshold).
- **E. Reproducibility & handoff** — User downloads script/bundle, re-runs locally or in another environment, or shares with a collaborator.

---

## 3. Scenario Definitions

### A. Exploration & tool discovery

| ID | Persona | Trigger | System behavior | Success criteria |
|----|---------|--------|------------------|------------------|
| **A1** | Scientist | “I have bulk RNA-seq count data and sample metadata. What can you do with it?” | Intent: execute. Router selects `bulk_rnaseq_analysis`. Response: **needs_inputs** with Required/Optional tables (GFM-rendered), “Provide the required inputs above…”. | User sees a clear table of parameters (count_matrix, sample_metadata, design_formula, alpha) and example values; no execution yet. |
| **A2** | Scientist | “What do I need to run single-cell RNA-seq analysis?” | **Intent can be ask or execute.** If *ask*: general answer about requirements (no tool). If *execute*: same as A1 for `single_cell_analysis` — needs_inputs with required (e.g. data_file, data_format), optional (resolution, steps). | User either gets a conceptual answer or a clear parameter table; in the execute path, understands required H5 path and optional parameters. |
| **A3** | Bioinformatician | “How do I run FastQC on my FASTQ files?” | Router → `fastqc_quality_analysis`; **needs_inputs** with input_r1, input_r2, optional output prefix. | User sees S3 path examples and knows they can paste paths and re-send. |
| **A4** | Scientist | “What’s the difference between bulk and single-cell RNA-seq?” | Intent: ask. BioinformaticsGuru answers in markdown (no tool execution). | Answer is accurate, references concepts (normalization, cell-level vs sample-level); no tables/plots. |

**Correctness**  
- Needs_inputs text must be generated from the tool’s **registered** required/optional inputs (no hallucinated parameters).  
- Q&A (A4, and A2 when intent is *ask*) must not trigger tool execution or alter session state.

---

### B. First run (supply data)

| ID | Persona | Trigger | System behavior | Success criteria |
|----|---------|--------|------------------|------------------|
| **B1** | Scientist | Uses Demo “Bulk RNA-seq: Infection × Time”, clicks “Load & Run” with pre-filled S3 paths. | Run `bulk_rnaseq_analysis` with count_matrix, sample_metadata, design_formula. ExecutionBroker runs tool; result: volcano plots, DE table, PCA. | Real plots and tables; run_id and session stored; “analysis.py” and “bundle.zip” links appear. |
| **B2** | Scientist | Same for “Single-Cell RNA-seq: SLE PBMC” with follow-up prompt (H5 path, 10x, resolution). | Run `single_cell_analysis`; result: UMAP, marker dot plot, DEG tables, cell-type proportions. | Real artifacts; script saved under session/run; download links present. |
| **B3** | Bioinformatician | “Run FastQC on s3://bucket/R1.fastq.gz and s3://bucket/R2.fastq.gz” (or pastes from manifest). | Run `fastqc_quality_analysis`; async job or inline; result: FastQC HTML report (or link). | Report is from real FastQC run; no hallucinated metrics. |
| **B4** | Scientist | Pastes FASTA sequences; “Build a phylogenetic tree from these sequences.” | Router → `phylogenetic_tree` (alignment + tree). Result: interactive tree (e.g. Newick + visualization). | Tree is computed from the provided sequences; Newick and viz are real. |
| **B5** | Scientist | “Trim my reads: R1 and R2 on S3, adapter AGATCGGAAGAGC, quality 20.” | Run `read_trimming` with forward_reads, reverse_reads, adapter, quality_threshold. | Trimmed output reflects real trimming; summary stats are correct. |

**Correctness**  
- All outputs (plots, tables, HTML, Newick) must be produced by the actual tool/script, not synthesized by the LLM.  
- If S3 is unreachable, bulk/scRNA tools should fall back to synthetic data **and** surface a clear notice (e.g. “Synthetic demo data used — real S3 data unavailable”).

---

### C. Iterative refinement (same analysis, change parameters)

| ID | Persona | Trigger | System behavior | Success criteria |
|----|---------|--------|------------------|------------------|
| **C1** | Bioinformatician | After B1: “Change the significance threshold to 0.01.” | Command routed to `patch_and_rerun`. System finds last scriptable run’s `analysis.py`, patches `ALPHA = 0.01`, re-runs script, returns new plots + parameter_diff + output_diff. | New run is a **child** of the previous (parent_run_id); DE table and volcano reflect alpha=0.01; user sees what changed. |
| **C2** | Bioinformatician | After B1: “Switch the volcano plot x-axis to linear fold change.” | `patch_and_rerun` with X_SCALE = "linear"; script patched and re-executed. | New volcano uses linear x-axis; run_id updated; analysis.py in new run dir has X_SCALE="linear". |
| **C3** | Scientist | After B1: “Use design formula ~infection_status + time_point only, no interaction.” | Either re-run with new design (new run) or patch DESIGN_FORMULA. | Result matches the requested design; contrast table and plots are consistent with the formula. |
| **C4** | Bioinformatician | After B2: “Increase clustering resolution to 0.8.” | If single_cell tool supports resolution in script: patch and re-run; else may need new run with new parameter. | Either patched re-run or clear “re-run with resolution=0.8” path; no silent ignore. |

**Correctness**  
- Parameter patches must apply to the **actual** script (e.g. ALPHA, X_SCALE, DESIGN_FORMULA) and re-execution must use the patched file.  
- Run ledger must record parent_run_id so the chain is traceable (reproducibility bundle can include iteration history).

---

### D. Validation & interpretation

| ID | Persona | Trigger | System behavior | Success criteria |
|----|---------|--------|------------------|------------------|
| **D1** | Scientist | After B1: “Which genes are most significant for the infection effect?” | Intent: ask. Answer derived from session/run context (e.g. DE table) or general knowledge; no re-execution. | Answer is consistent with the DE results (e.g. top genes by padj); no fabricated gene names. |
| **D2** | Scientist | “What does padj < 0.05 mean?” | Intent: ask. Markdown explanation (FDR, multiple testing). | Correct statistical interpretation; no tool run. |
| **D3** | Bioinformatician | After B1: “Show me the same volcano but with log2 scale on the x-axis again.” | Can be `patch_and_rerun` (X_SCALE = "log2") to get back to previous view. | Same as C2 but reverting; plot matches. |

**Correctness**  
- Q&A must not mutate run state or overwrite artifacts.  
- When user asks to “show again” or “change axis”, routing to `patch_and_rerun` must use the **last scriptable run** (or specified target_run) so the correct script is patched.

---

### E. Reproducibility & production handoff

| ID | Persona | Trigger | System behavior | Success criteria |
|----|---------|--------|------------------|------------------|
| **E1** | Scientist | After B1 or C1: Clicks “analysis.py” download. | GET /download/script?path=… returns the exact script that produced the current run (with parameter block). | File runs standalone (e.g. `python analysis.py`) and reproduces the same outputs when data is available. |
| **E2** | Bioinformatician | Clicks “bundle.zip” after an iterative chain (e.g. B1 → C1 → C2). | GET /download/bundle?session_id=…&run_id=… returns ZIP: analysis.py, plots, tables, README, run_manifest.json, iteration_history (ancestor runs). | Unzip → README describes how to run; manifest lists inputs/outputs; iteration_history allows replaying the chain. |
| **E3** | Core facility | “We need to run this exact analysis on the next cohort (new count matrix and metadata).” | User downloads analysis.py, edits COUNT_MATRIX and SAMPLE_METADATA (and RUN_DIR), runs locally or in their pipeline. | Script is self-contained; only parameter block and paths need editing; no Helix-specific runtime required. |
| **E4** | Bioinformatician | Shares bundle.zip with a collaborator who doesn’t use Helix.AI. | Collaborator extracts ZIP, reads README, runs analysis.py (or follows manifest). | Reproducibility is achievable without Helix.AI; README and manifest are accurate. |

**Correctness**  
- Script in the bundle must be the one that was **actually executed** for that run_id (or the final patched version in an iteration chain).  
- Manifest and README must list real artifact paths and parameters; no placeholder or wrong run_id.

---

## 4. Session and State Assumptions

- **One session** can contain multiple runs (e.g. B1, then C1, then C2). Runs are linked by `parent_run_id` where applicable.
- **New Session** (UI) creates a new backend session and clears in-app history; no cross-session state leakage.
- **needs_inputs** does not create a run; it only returns structured text and optional “Use example data” / “Load & Run” for demos.
- **Correctness** is enforced by: (1) ExecutionBroker as sole executor, (2) typed contracts and HandoffPolicy, (3) real artifact schemas, (4) script-first execution for scriptable tools (bulk RNA-seq, etc.) so the same code path that produced results is what gets patched and re-run.

---

## 5. Summary Table (by category)

| Category | Scenario IDs | Main guarantee |
|----------|--------------|----------------|
| **A. Exploration** | A1–A4 | Right tool and parameters surfaced; no spurious execution. |
| **B. First run** | B1–B5 | Real tool execution; run_id, script, and download links. |
| **C. Iteration** | C1–C4 | Parameter patch → re-run same script; parent_run_id and diffs. |
| **D. Validation** | D1–D3 | Q&A does not alter state; axis/param changes go to patch_and_rerun. |
| **E. Handoff** | E1–E4 | analysis.py and bundle.zip are complete and reproducible. |

---

## 6. Additional Scenarios (Existing Tools)

These scenarios use **tools that already exist** in Helix but were not listed in sections 3–5. **To support them:** (1) ensure the command router and agent map the relevant natural-language triggers to the right tool, (2) add the tool to `_TOOL_INPUT_REQUIREMENTS` in `main_with_mcp.py` if you want a structured “needs_inputs” response when required inputs are missing, (3) for multi-step flows (F3→F4, I1, I2), ensure session state is passed so downstream tools receive prior outputs. Reproducibility (script/bundle) applies where a scriptable run is produced (e.g. bulk RNA-seq, single-cell); sequence/plasmid lookups may return artifacts without a full analysis.py.

### F. Sequence, alignment & plasmid workflows

| ID | Persona | Trigger | System behavior | Success criteria |
|----|---------|--------|------------------|------------------|
| **F1** | Scientist | “Align these FASTA sequences” (pastes sequences). | Router → `sequence_alignment`. Result: aligned FASTA (or report). | Alignment is from real MSA tool; output is usable as input to phylogenetic_tree. |
| **F2** | Scientist | “Build a tree from these sequences” (unaligned FASTA). | Router → `phylogenetic_tree` (alignment done internally) or sequence_alignment → phylogenetic_tree. | Tree is real; Newick + visualization returned. |
| **F3** | Bioinformatician | “Generate 50 random variants of this sequence.” | Router → `mutate_sequence` (sequence, num_variants). Result: list of variant sequences. | Variants are generated by the tool; stored in session for downstream select_variants. |
| **F4** | Bioinformatician | After F3: “Select 10 variants with best conservation” (or lowest gaps, etc.). | Router → `select_variants` with session_id + selection_criteria. Uses session’s mutation_results. | Selection is from real variant_selection logic; 10 sequences returned. |
| **F5** | Scientist | “Show my gene in a plasmid map — vector pUC19, EcoRI/BamHI sites.” | Router → `plasmid_visualization` (vector_name, cloning_sites, optional insert_sequence). Result: plasmid map image. | Image is from real plasmid visualizer; correct vector and sites. |
| **F6** | Scientist | “Put these representative sequences into pUC19 and show the plasmid maps.” | Router → `plasmid_for_representatives` (representatives, aligned_sequences, vector_name, cloning_sites). Result: multiple plasmid visuals. | One map per representative; insert positions correct. |

**Support status:** Tools exist; ensure they appear in toolbox_inventory, are routable from natural language, and (where applicable) session state is passed correctly for multi-step (F3→F4).

---

### G. Literature & knowledge lookup

| ID | Persona | Trigger | System behavior | Success criteria |
|----|---------|--------|------------------|------------------|
| **G1** | Scientist | “Fetch the BRCA1 human sequence from NCBI.” | Router → `fetch_ncbi_sequence` (or equivalent). Result: FASTA or accession summary. | Sequence/content comes from NCBI API (or cached); no hallucinated sequence. |
| **G2** | Scientist | “Look up protein P53 in UniProt” or “What does UniProt say about TP53?” | Router → `query_uniprot`. Result: structured summary (function, domains, etc.). | Data from UniProt API; no fabricated annotations. |
| **G3** | Scientist | “What is GO:0006915?” | Router → `lookup_go_term`. Result: term name, definition, ontology. | Correct GO term; answer from go_tools. |

**Support status:** Tools exist; intent should be “execute” for lookup (not just Q&A) so the tool runs and returns real data. Q&A can still answer “what is GO?” without calling the tool.

---

### H. DNA synthesis & vendor

| ID | Persona | Trigger | System behavior | Success criteria |
|----|---------|--------|------------------|------------------|
| **H1** | Scientist | “Research vendors for synthesizing a 500 bp sequence” or “Who can make this oligo?” | Router → `dna_vendor_research`. Result: vendor options, lead times, cost hints. | Response is from real vendor-research logic (or structured template); no hallucinated vendors. |
| **H2** | Scientist | “Submit this sequence for synthesis to Twist” (or default vendor). | Router → `synthesis_submission` (sequences, vendor_preference, quantity). Result: submission confirmation or order ref. | Submission path is clearly indicated; if mock, user is not misled that a real order was placed. |

**Support status:** Tools exist; ensure routing and parameter extraction (sequence length, quantity, vendor) work from natural language; clarify mock vs real submission in UI.

---

### I. Multi-step pipelines (same session)

| ID | Persona | Trigger | System behavior | Success criteria |
|----|---------|--------|------------------|------------------|
| **I1** | Bioinformatician | “Align these sequences, then build a phylogenetic tree.” | Router/agent chains sequence_alignment → phylogenetic_tree (alignment output → tree input). Or single `phylogenetic_tree` if it accepts unaligned. | Tree is built from the user’s sequences; alignment is preserved in session if needed. |
| **I2** | Bioinformatician | “Generate 96 variants of this sequence, select the 12 with lowest gaps, and show them in plasmid maps.” | Chain: mutate_sequence → select_variants → plasmid_for_representatives. Session carries mutation_results between steps. | All three steps run; final output is 12 plasmid maps; reproducibility bundle can describe the chain. |

**Support status:** Tools exist; session state (mutated_sequences, mutation_results, alignment result) must be passed correctly. Agent or router must resolve “then” / “and then” to ordered tool calls.

---

## 7. Biotech & Biopharma — Extended Coverage (Roadmap)

To make Helix work for a **wide range of bioinformatics and data analysis requests** in biotech and biopharma, the following scenario groups are recommended. Each is marked as **Supported today** (existing tool + scenario in 3–6), **Partially supported** (some pieces exist), or **Planned** (future tools or integrations).

### J. Variant interpretation & VCF (Planned)

| ID | Persona | Trigger | Goal | Notes |
|----|---------|--------|------|------|
| **J1** | Clinical bioinformatician | “Annotate this VCF with gene and effect.” | VCF in → annotated VCF or table (gene, transcript, consequence). | Requires VCF parser + annotation (e.g. VEP, SnpEff, or API). |
| **J2** | Scientist | “Filter variants by ACMG classification” or “Which are pathogenic?” | Filter/flag by ACMG (pathogenic, likely pathogenic, VUS, etc.). | Often needed for clinical reports; depends on annotation + ACMG rules. |
| **J3** | Scientist | “Compare variant calls between this tumor and normal pair.” | Somatic variant calling or comparison of two VCFs. | Large scope; can start with “compare two VCFs” summary. |

**Support:** *Planned.* Today Helix has `variant_selection` in a **mutation workflow** (mutate_sequence → select_variants), not VCF-based variant annotation. J1–J3 require dedicated VCF/annotation tools or integrations.

---

### K. Microbiome & metagenomics (Partially supported → Planned)

| ID | Persona | Trigger | Goal | Notes |
|----|---------|--------|------|------|
| **K1** | Scientist | “Run FastQC and trim/merge my 16S amplicon reads.” | QC + trim + merge per sample. | **Supported today:** read_trimming, read_merging, fastqc_quality_analysis. |
| **K2** | Scientist | “Assign taxonomy to these ASV sequences” or “What taxa are in this 16S run?” | Taxonomy table (e.g. SILVA/GTDB), abundance matrix. | **Planned:** no built-in taxonomy assignment yet. |
| **K3** | Scientist | “Differential abundance of taxa between IBD and healthy.” | Statistical comparison of taxon (or gene) counts between groups. | **Planned:** conceptually similar to bulk RNA-seq DE but for count tables from microbiome; could be new tool or extension. |
| **K4** | Scientist | “Alpha and beta diversity for these samples.” | Diversity indices, PCoA, PERMANOVA. | **Planned:** requires abundance table + metadata; standard in microbiome pipelines. |

**Support:** K1 is supported. K2–K4 are common biotech/pharma needs (e.g. gut microbiome, biomanufacturing); add as roadmap items or tool-generation targets.

---

### L. Gene set enrichment & pathways (Partially supported)

| ID | Persona | Trigger | Goal | Notes |
|----|---------|--------|------|------|
| **L1** | Scientist | “What does GO term GO:0006915 mean?” | Term definition, hierarchy. | **Supported today:** lookup_go_term. |
| **L2** | Scientist | “Run GO enrichment on my DE gene list.” | Input: gene list (e.g. from B1). Output: enriched GO terms (or KEGG/Reactome). | **Partially supported:** single_cell_analysis can include “pathways” step; dedicated GO/KEGG enrichment on a **bulk DE list** may need a dedicated tool or patch to bulk_rnaseq output. |
| **L3** | Scientist | “Which KEGG pathways are overrepresented in these genes?” | Pathway enrichment table. | **Planned:** KEGG/Reactome enrichment tool or API integration. |

**Support:** L1 supported. L2 can be partially met via scRNA “pathways” or by adding an enrichment step to the bulk RNA-seq bundle; L3 is a natural extension.

---

### M. Biomarker discovery & validation (Planned)

| ID | Persona | Trigger | Goal | Notes |
|----|---------|--------|------|------|
| **M1** | Scientist | “Which genes best separate responders from non-responders?” | Feature selection or ranking (e.g. from DE, ROC, or ML). | **Planned:** could build on top of bulk or scRNA results; needs clear input (e.g. expression matrix + response labels). |
| **M2** | Scientist | “Validate this signature on the holdout cohort.” | Apply existing gene list/coefficients to new data; report performance. | **Planned:** requires “signature” as first-class object and a validation tool. |
| **M3** | Scientist | “ROC curve and AUC for my classifier.” | Binary (or multi-class) outcome + predictions → ROC, AUC, confusion matrix. | **Planned:** generic small ML/eval tool or integration. |

**Support:** *Planned.* Valuable for biopharma (e.g. companion diagnostics, stratification); no built-in tool yet.

---

### N. Compliance, audit & reporting (Partially supported)

| ID | Persona | Trigger | Goal | Notes |
|----|---------|--------|------|------|
| **N1** | Core facility | “Export an audit trail for this run.” | Timestamped log of inputs, parameters, tool, outputs, run_id. | **Partially supported:** run.json, session history, and bundle manifest provide much of this; a dedicated “export audit log” format (e.g. PDF or CSV) could be added. |
| **N2** | Scientist | “Regenerate the exact report from run X.” | Re-run report generation from stored artifacts (no re-execution of analysis). | **Partially supported:** bundle + analysis.py allow full re-run; “report only” from run_id is a thin wrapper (e.g. ds_report or similar). |
| **N3** | Regulator / QA | “Prove that this result came from this code and this data.” | Immutable link: run_id → script hash, data references, output hashes. | **Partially supported:** manifest and iteration_history in bundle; optional hashing of inputs/outputs would strengthen “prove” story. |

**Support:** N1–N3 are largely about **exposing and formatting** existing provenance (run.json, manifest, bundle); no new execution tools strictly required.

---

## 8. Domain → Scenario Map (Quick Reference)

| Domain | Scenario IDs | Support |
|--------|--------------|--------|
| **Bulk RNA-seq** | A1, B1, C1–C3, D1–D3, E1–E4 | Supported |
| **Single-cell RNA-seq** | A2, B2, C4, E* | Supported |
| **Sequencing QC (FastQC, trim, merge)** | A3, B3, B5, K1 | Supported |
| **Phylogenetics** | B4, F1–F2, I1 | Supported |
| **Sequence alignment** | F1, I1 | Supported |
| **Mutation & variant selection** | F3–F4, I2 | Supported |
| **Plasmid visualization** | F5–F6, I2 | Supported |
| **Literature / knowledge** | A4, G1–G3, L1 | Supported |
| **DNA synthesis** | H1–H2 | Supported |
| **Multi-step pipelines** | I1–I2 | Supported (session state) |
| **Variant/VCF annotation** | J1–J3 | Planned |
| **Microbiome (taxonomy, diversity, DA)** | K2–K4 | Planned (K1 supported) |
| **Enrichment (GO/KEGG on gene lists)** | L2–L3 | Partially / Planned |
| **Biomarker discovery & validation** | M1–M3 | Planned |
| **Compliance & audit** | N1–N3 | Partially (provenance exists) |

---

## 9. How to Use This Document

- **Product / design**: Use scenarios to define “must work” flows for launch and to prioritize UI (e.g. download links, New Session, needs_inputs table rendering). Sections 6–7 extend coverage to sequence/plasmid, literature, DNA synthesis, and multi-step flows; section 8 is the roadmap for biotech/biopharma breadth.
- **QA / testing**: A full **testbed** maps every scenario to concrete tests. See **[TESTBED.md](TESTBED.md)** for test IDs, layers (backend unit/integration, frontend, E2E), and assertions. Automated backend tests live in **`tests/testbed/`**; run with `pytest tests/testbed/ -v`. You can also run manual test cases (e.g. run B1 with demo data, then C1, then download bundle; run F3→F4→F6 for variant-to-plasmid chain).
- **Correctness**: Any change to ExecutionBroker, command routing, or script generation should be checked against B* and C* (real execution, correct script) and E* (bundle integrity). Same standard applies to F–I when those tools are invoked.
- **Reproducibility**: E1–E4 are the contract for “production fast” — the script and bundle must be sufficient to re-run outside Helix.AI.
- **Roadmap**: Use section 7 (J–N) and section 8 to prioritize new tools or integrations (VCF/ACMG, microbiome taxonomy/diversity, GO/KEGG enrichment, biomarker validation, audit export) so Helix supports a wide range of bioinformatics and data analysis requests in biotech and biopharma.
