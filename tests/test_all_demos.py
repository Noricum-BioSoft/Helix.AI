#!/usr/bin/env python3
"""
Test all 5 demo scenarios against the local (or remote) Helix.AI backend.

Usage:
  python tests/test_all_demos.py [--url http://localhost:8001]
"""

import sys, json, time, textwrap, argparse
import requests

BASE_URL = "http://localhost:8001"

DEMOS = [
    {
        "id": "1-bulk-rnaseq-factorial",
        "title": "Bulk RNA-seq: Infection × Time Factorial",
        "expected_behavior": "needs_inputs",
        "tool": "bulk_rnaseq_analysis",
        "prompt": """You are analyzing an RNA-seq transcriptome dataset from a mouse study investigating the effects of Toxoplasma gondii infection on brain gene expression.

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
Implement best practices for RNA-seq analysis.""",
    },
    {
        "id": "2-scrna-sle-pbmc",
        "title": "Single-Cell RNA-seq: SLE Immune Profiling",
        "expected_behavior": "needs_inputs",
        "tool": "single_cell_analysis",
        "prompt": """You are analyzing a single-cell RNA-seq dataset from a human peripheral blood mononuclear cell (PBMC) study investigating immune cell dysregulation in patients with systemic lupus erythematosus (SLE) compared to healthy controls.

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
  - Stacked bar chart of cell-type proportions per sample.""",
    },
    {
        "id": "3-amplicon-qc-pipeline",
        "title": "Amplicon QC Pipeline: 16S Gut Microbiome",
        "expected_behavior": "executes_pipeline",
        "tool": "fastqc_quality_analysis / read_trimming / read_merging",
        "prompt": """You are processing a 16S rRNA amplicon sequencing dataset from a gut microbiome study. Raw paired-end FASTQ files are on S3 and need a full preprocessing pipeline before downstream diversity analysis.

Dataset
  Illumina MiSeq 2×250 bp paired-end reads; V3–V4 hypervariable region.
  Forward reads: s3://helix-test-data/microbiome/run1/sample01_R1.fastq.gz
  Reverse reads: s3://helix-test-data/microbiome/run1/sample01_R2.fastq.gz
  Output prefix:  s3://helix-test-data/microbiome/run1/processed/

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
  - Quality summary CSV (sample, raw_reads, post_trim_reads, merged_reads, merge_rate).""",
    },
    {
        "id": "4-bulk-rnaseq-timecourse",
        "title": "Bulk RNA-seq: Liver Injury Time-Course",
        "expected_behavior": "needs_inputs",
        "tool": "bulk_rnaseq_analysis",
        "prompt": """You are analyzing a bulk RNA-seq time-course dataset from a murine model of acute liver injury (acetaminophen overdose), tracking transcriptional recovery from peak injury back to baseline.

Dataset
  Raw count matrix from strand-specific paired-end RNA-seq (75 bp).
  5 time points post-APAP administration × 4 biological replicates = 20 samples.

Study Design
  One-factor time-course (independent cross-sectional, not paired longitudinal):
    Time points: 0 h (baseline), 6 h, 24 h, 72 h, 168 h (7 days)
    n = 4 mice per time point

Objectives
  1. Model time as an ordered factor in DESeq2 using likelihood ratio tests
     (full model: ~time; reduced: ~1) to identify any time-regulated genes.
  2. Identify early-response genes (peak at 6–24 h), late-recovery genes
     (returning to baseline by 72–168 h), and persistently dysregulated genes.
  3. Fit smooth expression trajectories over time using spline regression for
     the top 500 most variable genes.
  4. Perform GO biological process enrichment for each temporal cluster.

Desired Outputs
  - Line plots of mean normalized expression ± SEM over time for the top 20
    injury-response and top 20 recovery genes (PNG).
  - Heatmap of all temporally regulated genes grouped by expression pattern.
  - CSV of LRT results (gene, LRT statistic, padj, peak time point).
  - GO enrichment table per temporal cluster (term, p.adjust, gene ratio).""",
    },
    {
        "id": "5-phylogenetics-sarscov2",
        "title": "Phylogenetics: SARS-CoV-2 Variant Divergence",
        "expected_behavior": "multi_step_plan",
        "tool": "fetch_ncbi_sequence → sequence_alignment → phylogenetic_tree",
        "prompt": """You are conducting a comparative evolutionary analysis of the SARS-CoV-2 spike protein across major variants of concern (VOCs) to characterize mutational divergence from the ancestral Wuhan-Hu-1 reference sequence.

Dataset
  Full-length spike protein amino acid sequences for 8 SARS-CoV-2 variants:
  Wuhan-Hu-1 (reference), Alpha (B.1.1.7), Beta (B.1.351), Gamma (P.1),
  Delta (B.1.617.2), Omicron BA.1, Omicron BA.4/5, and XBB.1.5.
  Sequences should be fetched from NCBI RefSeq where available.

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
  - Table of key mutation presence/absence per variant.""",
    },
]

SEP = "=" * 70

def send_prompt(url: str, prompt: str, timeout: int = 90) -> dict:
    payload = {"command": prompt, "session_id": f"demo-test-{int(time.time())}"}
    t0 = time.time()
    try:
        resp = requests.post(f"{url}/execute", json=payload, timeout=timeout)
        elapsed = time.time() - t0
        return {
            "status_code": resp.status_code,
            "elapsed": round(elapsed, 2),
            "data": resp.json() if resp.headers.get("content-type", "").startswith("application/json") else resp.text,
        }
    except requests.exceptions.Timeout:
        elapsed = time.time() - t0
        return {"status_code": 504, "elapsed": round(elapsed, 2), "error": "TIMEOUT"}
    except Exception as e:
        elapsed = time.time() - t0
        return {"status_code": 0, "elapsed": round(elapsed, 2), "error": str(e)}


def extract_message(data: dict) -> str:
    """Pull message text from various response shapes."""
    if not isinstance(data, dict):
        return str(data)
    for key in ("response", "message", "content", "result"):
        val = data.get(key)
        if val and isinstance(val, str) and len(val) > 10:
            return val
        if val and isinstance(val, dict):
            sub = extract_message(val)
            if sub:
                return sub
    return str(data)


def analyze_result(demo: dict, result: dict) -> tuple[bool, str]:
    """Returns (passed, note)."""
    if result.get("error") == "TIMEOUT":
        return False, f"TIMEOUT after {result['elapsed']}s"

    sc = result.get("status_code")
    if sc != 200:
        return False, f"HTTP {sc}"

    data = result.get("data", {})
    if isinstance(data, str):
        return False, f"Non-JSON response: {data[:120]}"

    # Any backend-reported failure is a fail (unless it's a valid needs_inputs response)
    if isinstance(data, dict) and data.get("success") is False:
        err = data.get("error", "unknown error")
        # needs_inputs responses might have success=False in some paths; check content first
        msg = extract_message(data)
        if "needs_inputs" not in msg.lower() and "input" not in msg.lower():
            return False, f"Backend error: {err}"

    msg = extract_message(data)
    msg_lower = msg.lower()

    expected = demo["expected_behavior"]

    if expected == "needs_inputs":
        bad_phrases = ["completed successfully", "analysis complete", "results saved"]
        good_phrases = ["input", "data", "file", "provide", "missing", "count matrix",
                        "fastq", "path", "require", "need", "please", "workflow"]
        has_bad = any(p in msg_lower for p in bad_phrases)
        has_good = any(p in msg_lower for p in good_phrases)
        if has_bad:
            return False, f"Incorrectly claimed execution success: '{msg[:120]}'"
        if not has_good:
            return False, f"Did not ask for inputs: '{msg[:120]}'"
        return True, "Correctly returned needs_inputs / plan response"

    elif expected == "executes_pipeline":
        # Should return HTTP 200 with some meaningful response (job ID, plan, or results)
        if data.get("success") is False:
            err = data.get("error", "unknown")
            return False, f"Pipeline failed: {err[:120]}"
        good_phrases = ["job", "fastqc", "trim", "merg", "pipeline", "step", "processing",
                        "workflow", "output", "needs_inputs", "input", "file", "report"]
        has_good = any(p in msg_lower for p in good_phrases)
        if not has_good:
            return False, f"No meaningful pipeline response: '{msg[:120]}'"
        return True, f"Got pipeline response in {result['elapsed']}s"

    elif expected == "multi_step_plan":
        bad_phrases = ["clustering completed successfully", "no module named"]
        good_phrases = ["step", "align", "tree", "ncbi", "sequence", "pipeline",
                        "fetch", "phylogen", "variant", "blast", "mafft", "workflow"]
        has_bad = any(p in msg_lower for p in bad_phrases)
        has_good = any(p in msg_lower for p in good_phrases)
        if has_bad:
            return False, f"Wrong routing or error: '{msg[:120]}'"
        if data.get("success") is False:
            return False, f"Backend error: {data.get('error', 'unknown')[:120]}"
        if has_good:
            return True, "Got multi-step plan response"
        return False, f"Missing plan keywords in response: '{msg[:120]}'"

    return True, "No specific check"


def run_tests(url: str):
    print(f"\n{SEP}")
    print(f"  Helix.AI Demo Scenario Tests  —  {url}")
    print(SEP)

    results = []
    for i, demo in enumerate(DEMOS, 1):
        print(f"\n[{i}/5] {demo['title']}")
        print(f"      Expected: {demo['expected_behavior']} | Tool: {demo['tool']}")
        print(f"      Sending request...", end="", flush=True)

        result = send_prompt(url, demo["prompt"], timeout=90)
        passed, note = analyze_result(demo, result)

        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"\r      {status}  [{result['elapsed']}s]  {note}")

        data = result.get("data", {})
        if isinstance(data, dict):
            msg = data.get("response") or data.get("message") or data.get("content") or ""
            if msg:
                preview = textwrap.shorten(str(msg), width=200, placeholder="…")
                print(f"      Response preview: {preview}")

        results.append({"demo": demo["id"], "passed": passed, "note": note, "elapsed": result["elapsed"]})

    print(f"\n{SEP}")
    n_pass = sum(1 for r in results if r["passed"])
    print(f"  Results: {n_pass}/{len(results)} passed")
    print(SEP)
    for r in results:
        icon = "✅" if r["passed"] else "❌"
        print(f"  {icon}  {r['demo']:40s}  {r['note']}")
    print(SEP + "\n")
    return all(r["passed"] for r in results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", default="http://localhost:8001")
    args = parser.parse_args()
    ok = run_tests(args.url)
    sys.exit(0 if ok else 1)
