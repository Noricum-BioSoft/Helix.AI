#!/usr/bin/env python3
"""50-workflow benchmark for Helix.AI advisory-rendering quality check."""
import requests, time, json, sys, os

BASE    = os.environ.get("HELIX_BASE_URL", "http://localhost:8001")
TIMEOUT = int(os.environ.get("HELIX_BENCH_TIMEOUT", "120"))

PROBES = [
    ("WF01","bulk-rnaseq",         "Which genes are differentially expressed between treated and control samples? I have FASTQ files ready."),
    ("WF02","scrna-seq",           "I have a 10x Genomics scRNA-seq dataset. Can you identify cell types and create a UMAP?"),
    ("WF03","scrna-condition",     "How do T-cell subtypes respond to anti-PD1 treatment in my scRNA-seq data?"),
    ("WF04","atac-seq",            "I have ATAC-seq data. Which regulatory regions are more open in tumor vs normal?"),
    ("WF05","chip-seq",            "I have H3K27ac ChIP-seq data. Find binding sites and motifs for this histone mark."),
    ("WF06","cut-and-run",         "Analyze my CUT&RUN experiment for CTCF binding. Show signal tracks and differential binding."),
    ("WF07","wgs-germline",        "I have whole genome sequencing. Find SNPs, indels, and structural variants in the germline."),
    ("WF08","wes-rare",            "My patient has rare disease. Which variants in the WES data might explain the phenotype?"),
    ("WF09","tumor-normal",        "Identify somatic mutations and driver genes in my tumor-normal WES pair."),
    ("WF10","cnv",                 "Find copy number amplifications and deletions in my WGS data."),
    ("WF11","sv",                  "Detect structural variants including inversions and translocations in my long-read sequencing."),
    ("WF12","assembly",            "I have Oxford Nanopore reads. Please assemble the genome and assess completeness."),
    ("WF13","metagenomics-tax",    "Which organisms are present in my gut microbiome shotgun metagenomics sample?"),
    ("WF14","metagenomics-fn",     "What metabolic pathways are enriched in my microbiome compared to controls?"),
    ("WF15","16s",                 "Analyze my 16S rRNA sequencing data for microbial community composition and alpha/beta diversity."),
    ("WF16","viral-seq",           "What SARS-CoV-2 variants and lineages are in my sequencing data?"),
    ("WF17","phylogenetics",       "Build a phylogenetic tree from my pathogen genomes to detect transmission clusters."),
    ("WF18","proteomics",          "Which proteins are differentially abundant between my two proteomics conditions?"),
    ("WF19","phosphoproteomics",   "Which kinases are active based on my phosphoproteomics data?"),
    ("WF20","metabolomics",        "Find biomarker metabolites that distinguish disease from control in my metabolomics dataset."),
    ("WF21","spatial-tx",          "I have Visium spatial transcriptomics data. Map gene expression across the tissue section."),
    ("WF22","imaging",             "Segment cells and extract phenotypic features from my confocal microscopy images."),
    ("WF23","crispr-ko",           "Which gene knockouts cause fitness defects in my MAGeCK CRISPR screen?"),
    ("WF24","crispra",             "Find hits from my CRISPRa screen that activate resistance to drug treatment."),
    ("WF25","perturb-seq",         "Analyze transcriptional responses to CRISPR perturbations in my Perturb-seq experiment."),
    ("WF26","timecourse",          "Identify dynamically expressed gene clusters across my 6-timepoint RNA-seq experiment."),
    ("WF27","multi-omics",         "Integrate my RNA-seq and ATAC-seq data to find correlated expression and accessibility changes."),
    ("WF28","ewas",                "Which CpG sites are differentially methylated between cases and controls in my EWAS?"),
    ("WF29","wgbs",                "Profile DNA methylation from my WGBS data and find differentially methylated regions."),
    ("WF30","eqtl",                "Map expression quantitative trait loci using my genotype and RNA-seq data."),
    ("WF31","gwas",                "Run a genome-wide association study on my SNP array data for blood pressure phenotype."),
    ("WF32","prs",                 "Calculate polygenic risk scores from my genotype data using published GWAS weights."),
    ("WF33","classifier",          "Build a machine learning classifier to predict cancer subtype from gene expression."),
    ("WF34","drug-response",       "Predict drug sensitivity from omics features in my cancer cell line data."),
    ("WF35","cell-screen",         "Analyze dose-response curves from my compound screening viability assay."),
    ("WF36","grn",                 "Infer a gene regulatory network from my scRNA-seq data to identify master regulators."),
    ("WF37","splicing",            "Find differentially spliced exons between tumor and normal in my RNA-seq data."),
    ("WF38","isoform",             "Discover novel transcript isoforms from my PacBio long-read RNA-seq data."),
    ("WF39","ribo-seq",            "Calculate translation efficiency by combining my Ribo-seq and RNA-seq data."),
    ("WF40","mirna",               "Find differentially expressed miRNAs and predict their mRNA targets."),
    ("WF41","tcr-bcr",             "Profile T-cell receptor repertoire diversity and clonal expansion in my TCR-seq data."),
    ("WF42","hla",                 "Type HLA alleles from my WGS data."),
    ("WF43","neoantigen",          "Predict neoantigens from somatic mutations and HLA type for my tumor sample."),
    ("WF44","de-novo",             "Assemble a de novo transcriptome from my RNA-seq data from a non-model organism."),
    ("WF45","comparative",         "Compare gene families and synteny across my 5 fungal genomes."),
    ("WF46","annotation",          "Annotate genes and functions in my newly assembled bacterial genome."),
    ("WF47","edna",                "Which species are present in my eDNA water sample using amplicon sequencing?"),
    ("WF48","fastq-qc",            "Run quality control on my FASTQ files and tell me if the data is usable."),
    ("WF49","public-data",         "Download and reanalyze the RNA-seq dataset GSE123456 from GEO."),
    ("WF50","integrative",         "Integrate evidence from RNA-seq, ATAC-seq, and published GWAS to prioritize candidate genes for my disease."),
]


def classify_rendering(text: str) -> str:
    if not text:
        return "empty"
    try:
        parsed = json.loads(text)
        ht = parsed.get("helix_type")
        return f"advisory:{ht}" if ht else "json"
    except Exception:
        return "text"


def fresh_session() -> str:
    return requests.post(f"{BASE}/create_session", timeout=30).json()["session_id"]


def main() -> int:
    print(f"Backend: {BASE}  timeout={TIMEOUT}s  (fresh session per probe)\n")
    sys.stdout.flush()

    results, passed, failed, errors = [], 0, 0, 0

    for wf_id, wf_type, prompt in PROBES:
        try:
            sid = fresh_session()
            t0 = time.time()
            r  = requests.post(
                f"{BASE}/execute",
                json={"command": prompt, "session_id": sid},
                timeout=TIMEOUT,
            )
            elapsed = time.time() - t0
            d       = r.json()
            tool    = d.get("tool", "?")
            status  = d.get("status", "?")
            text    = d.get("text") or ""
            rendering = classify_rendering(text)

            ok = status not in ("error", "failed") and bool(text)
            if ok:
                passed += 1
            else:
                failed += 1

            results.append(
                dict(
                    id=wf_id, type=wf_type, tool=tool, status=status,
                    elapsed=round(elapsed, 1), rendering=rendering, ok=ok,
                )
            )
            mark = "✓" if ok else "✗"
            print(f"{mark} {wf_id} ({elapsed:.1f}s)  tool={tool}  status={status}  rendering={rendering}")
        except Exception as exc:
            errors += 1
            results.append(
                dict(
                    id=wf_id, type=wf_type, tool="ERROR", status="error",
                    elapsed=0, rendering="error", ok=False,
                )
            )
            print(f"✗ {wf_id} ERROR: {exc}")
        sys.stdout.flush()

    print(f"\n{'='*60}")
    print(f"PASSED: {passed}/50  |  FAILED: {failed}  |  ERRORS: {errors}")

    advisory_ok = sum(1 for r in results if r["rendering"].startswith("advisory:advisory"))
    text_ok     = sum(1 for r in results if r["rendering"] == "text")
    print(f"Rendering breakdown → canonical advisory: {advisory_ok}  plain text: {text_ok}  other: {50 - advisory_ok - text_ok}")

    out = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "passed": passed, "failed": failed, "errors": errors,
        "advisory_canonical": advisory_ok,
        "results": results,
    }
    os.makedirs("artifacts/benchmark_results", exist_ok=True)
    path = f"artifacts/benchmark_results/benchmark_{time.strftime('%Y%m%d_%H%M%S')}.json"
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"Saved → {path}")

    return 0 if passed >= 45 else 1


if __name__ == "__main__":
    sys.exit(main())
