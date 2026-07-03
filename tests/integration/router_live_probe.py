#!/usr/bin/env python3
"""
Live probe of the deployed router via /debug/route.

Hits every category of input the router must handle:
  1. The 5 demo initial prompts
  2. The 3 demo followUpPrompts (the "Load & run" path)
  3. Multi-step composite requests
  4. Unknown intents
  5. Read-only / Q&A
  6. Edge cases (empty, very long, weird unicode)

Each case prints PASS/FAIL and the exact routing decision.
Intended to be run end-to-end as a single script — no pytest dependencies.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import List, Tuple

import base64
import urllib.request

DEFAULT_BASE = os.getenv("HELIX_BASE_URL", "https://helix-beta.noricum-biosoft.com")


def _auth_headers() -> dict:
    """Return an Authorization header if HELIX_BASIC_AUTH ('user:pass') is set."""
    creds = os.getenv("HELIX_BASIC_AUTH", "")
    if not creds:
        return {}
    token = base64.b64encode(creds.encode()).decode()
    return {"Authorization": f"Basic {token}"}


# ─────────────────────────────────────────────────────────────────────────────
# Test cases — (label, command, expected_pattern, allowed_fallbacks)
# expected_pattern: substring that must appear in result["final"]["tool"]
# allowed_fallbacks: alternate acceptable tools (oneof)
# ─────────────────────────────────────────────────────────────────────────────

TestCase = Tuple[str, str, List[str], str]


CASES: List[TestCase] = [
    # ── Demo initial prompts (multi-paragraph study descriptions) ─────────
    ("demo1.bulk_rnaseq.initial",
     "You are analyzing an RNA-seq transcriptome dataset from a mouse study investigating Toxoplasma gondii infection. 2x2 factorial: Infection (Inf/Uninf) × Time (11dpi/33dpi). Perform DE analysis. Include PCA, sample distance heatmap.",
     ["bulk_rnaseq_analysis", "multi_step_workflow"],
     "bulk RNA-seq specialist tool or multi-step plan"),

    ("demo2.scrna.initial",
     "You are analyzing a single-cell RNA-seq dataset from PBMCs comparing SLE patients to healthy controls. 10x v3, 5 vs 5. Objectives: QC, normalize, PCA→UMAP, Leiden cluster, annotate cell types, pseudobulk DEGs.",
     ["single_cell_analysis", "multi_step_workflow"],
     "scRNA specialist tool or multi-step plan"),

    ("demo3.amplicon.initial",
     "You are processing a 16S rRNA amplicon sequencing dataset. Illumina MiSeq 2x250 paired-end. Run FastQC, trim adapters, merge paired-end reads, generate QC report.",
     ["multi_step_workflow", "metagenomics_16s", "fastqc_quality_analysis", "read_trimming"],
     "16S specialist or multi-step plan"),

    ("demo4.apap.initial",
     "I have a bulk RNA-seq time-course study from a murine acetaminophen liver-injury model. 5 time points, 4 reps each. Design the analysis workflow before execution. Include checkpoints for QC and biological interpretation.",
     ["multi_step_workflow", "bulk_rnaseq_analysis", "handle_natural_command"],
     "workflow design — agent or multi-step"),

    ("demo5.phylo.initial",
     "You are conducting a comparative evolutionary analysis of the SARS-CoV-2 spike protein across major variants of concern. Retrieve sequences, MSA, ML tree with bootstrap, RBD mutation annotation, pairwise identity matrix.",
     ["multi_step_workflow", "phylogenetic_tree"],
     "multi-step or phylo direct"),

    # ── Demo followUpPrompts (the "Load & run" path) ─────────────────────
    ("demo1.bulk_rnaseq.followup",
     "Run bulk RNA-seq differential expression analysis with these inputs:\ncount_matrix: s3://noricum-ngs-data/demo/rnaseq/tgondii_counts.csv\nsample_metadata: s3://noricum-ngs-data/demo/rnaseq/tgondii_metadata.csv\ndesign_formula: ~infection_status + time_point",
     ["bulk_rnaseq_analysis"],
     "MUST match deterministic — bulk RNA-seq inputs block"),

    ("demo2.scrna.followup",
     "data_file: s3://noricum-ngs-data/demo/scrna/sle_pbmc_filtered_feature_bc_matrix.h5\ndata_format: 10x\nresolution: 0.5\nsteps: all",
     ["single_cell_analysis"],
     "MUST match deterministic — scRNA 10x inputs block"),

    ("demo5.phylo.followup",
     "Build a phylogenetic tree from sequences:\nsequences: s3://noricum-ngs-data/demo/phylo/sars_cov2_spike.fasta",
     ["phylogenetic_tree"],
     "MUST match deterministic — phylo with sequences:"),

    # ── Multi-step composite requests ────────────────────────────────────
    ("multistep.fetch_align_tree",
     "Fetch SARS-CoV-2 spike sequences from NCBI then align with MAFFT then build a phylogenetic tree",
     ["multi_step_workflow"],
     "multi-step composite (LLM)"),

    ("multistep.qc_trim_merge",
     "Run FastQC then trim adapters then merge paired-end reads",
     ["multi_step_workflow"],
     "multi-step composite (LLM)"),

    # ── Unknown intent / nonsense ────────────────────────────────────────
    ("unknown.gibberish",
     "qwerty foobar baz quux not a real bioinformatics request",
     ["unknown_intent", "handle_natural_command"],
     "unknown / no-tool match"),

    # ── Q&A / informational ──────────────────────────────────────────────
    ("qa.what_is_deseq",
     "What is DESeq2 and when should I use it for RNA-seq?",
     ["handle_natural_command", "unknown_intent"],
     "informational query — handle_natural_command"),

    ("qa.list_tools",
     "What tools do you have available?",
     ["toolbox_inventory", "handle_natural_command"],
     "tool list — toolbox_inventory or NLP"),

    # ── Specific tool requests ───────────────────────────────────────────
    ("specific.fastqc",
     "Run FastQC on s3://my-bucket/r1.fastq.gz and s3://my-bucket/r2.fastq.gz",
     ["fastqc_quality_analysis"],
     "FastQC specialist"),

    ("specific.read_trimming",
     "Trim adapter AGATCGGAAGAGC with quality threshold 25 from /data/sample.fq",
     ["read_trimming"],
     "read trimming specialist"),

    ("specific.fetch_ncbi",
     "Fetch sequence MN908947.3 from NCBI",
     ["fetch_ncbi_sequence"],
     "NCBI fetch specialist"),

    # ── Explicit tool directive (deterministic shortcut) ─────────────────
    ("explicit.tool_directive",
     "tool: phylogenetic_tree",
     ["phylogenetic_tree"],
     "MUST match deterministic — explicit `tool:` directive"),

    # ── Tabular Q&A ──────────────────────────────────────────────────────
    ("tabular_qa.de_question",
     "What are the top 10 genes by absolute log2FC that are statistically significant (padj < 0.05) in the uploaded DE results table?",
     ["tabular_qa", "tabular_analysis", "handle_natural_command"],
     "DE table Q&A — tabular_qa/tabular_analysis or NLP fallback"),

    ("tabular_qa.rank_genes",
     "Rank these genes by absolute fold-change and flag the ones with padj < 0.05:\n"
     "gene,log2FC,padj\nTP53,-2.14,0.00012\nBRCA1,3.42,0.000045\nMYC,1.83,0.0082\nKRAS,0.88,0.062",
     ["tabular_qa", "tabular_analysis", "handle_natural_command"],
     "inline DE table ranking — tabular_analysis (explicit rank op) preferred"),

    # ── Edge cases ───────────────────────────────────────────────────────
    ("edge.short",
     "align",
     ["sequence_alignment", "handle_natural_command", "unknown_intent"],
     "very short prompt"),

    ("edge.weird_unicode",
     "Build a phylogenetic tree 🧬 from these sequences: >s1 ATCG >s2 ATGC",
     ["phylogenetic_tree", "multi_step_workflow"],
     "unicode + inline FASTA"),

    ("edge.long_pasted_block",
     "I have an RNA-seq dataset. " + ("Some background context. " * 60) + " Please run differential expression with count_matrix: s3://b/c.csv sample_metadata: s3://b/m.csv design_formula: ~x",
     ["bulk_rnaseq_analysis", "multi_step_workflow"],
     "long pasted block with hidden inputs"),
]


# ─────────────────────────────────────────────────────────────────────────────
# HTTP helper (urllib only — no extra deps)
# ─────────────────────────────────────────────────────────────────────────────

def post_json(url: str, payload: dict, timeout: int = 45) -> dict:
    body = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url,
        data=body,
        headers={"Content-Type": "application/json", **_auth_headers()},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", default=DEFAULT_BASE)
    parser.add_argument("--out", default=None, help="Optional JSON output path")
    args = parser.parse_args()

    base = args.base.rstrip("/")
    print(f"# Live router probe — base={base}")
    print(f"# {len(CASES)} cases")
    print()

    results = []
    passed = 0
    failed = 0
    for label, command, allowed, note in CASES:
        try:
            t0 = time.time()
            r = post_json(f"{base}/debug/route", {"command": command})
            elapsed_ms = int((time.time() - t0) * 1000)
        except Exception as exc:
            print(f"FAIL {label}: HTTP error: {exc}")
            failed += 1
            results.append({"label": label, "ok": False, "error": str(exc)})
            continue

        final_tool = (r.get("final") or {}).get("tool")
        deterministic = bool(r.get("deterministic_match"))
        is_fallback = bool(r.get("is_router_fallback"))

        ok = final_tool in allowed
        status = "PASS" if ok else "FAIL"
        det_marker = " [deterministic]" if deterministic else " [llm]"
        print(f"{status} {label:<40} -> {final_tool}{det_marker}  ({elapsed_ms} ms)")
        if not ok:
            print(f"     expected one of: {allowed}")
            print(f"     note: {note}")

        results.append({
            "label": label,
            "ok": ok,
            "final_tool": final_tool,
            "expected_oneof": allowed,
            "deterministic": deterministic,
            "is_fallback": is_fallback,
            "elapsed_ms": elapsed_ms,
            "command_excerpt": command[:90],
        })
        if ok:
            passed += 1
        else:
            failed += 1

    print()
    print(f"Results: {passed} passed, {failed} failed (out of {len(CASES)})")

    # Summary breakdown
    deterministic_count = sum(1 for r in results if r.get("deterministic"))
    llm_count = sum(1 for r in results if not r.get("deterministic"))
    avg_det = (
        sum(r["elapsed_ms"] for r in results if r.get("deterministic")) / max(deterministic_count, 1)
    )
    avg_llm = (
        sum(r["elapsed_ms"] for r in results if not r.get("deterministic")) / max(llm_count, 1)
    )
    print(f"  Deterministic matches: {deterministic_count}  avg latency {avg_det:.0f} ms")
    print(f"  LLM matches:           {llm_count}  avg latency {avg_llm:.0f} ms")

    if args.out:
        with open(args.out, "w") as f:
            json.dump({
                "base": base,
                "passed": passed,
                "failed": failed,
                "deterministic_count": deterministic_count,
                "llm_count": llm_count,
                "avg_deterministic_ms": avg_det,
                "avg_llm_ms": avg_llm,
                "results": results,
            }, f, indent=2)
        print(f"  Wrote {args.out}")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
