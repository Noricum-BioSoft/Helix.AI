#!/usr/bin/env python3
"""
End-to-end probe of the deployed beta for every demo flow.

For each demo:
  1. Submit the initial study-description prompt
  2. (If applicable) submit the followUpPrompt
  3. Capture the response and assert no known error strings + non-empty content

This exercises the full stack: nginx → SSE → router → approval gate →
agent / specialist tool → response renderer.
"""
from __future__ import annotations

import base64
import json
import os
import sys
import time
import urllib.request
from typing import Any, Dict, List, Optional

DEFAULT_BASE = os.getenv("HELIX_BASE_URL", "https://helix-beta.noricum-biosoft.com")
TIMEOUT_S = int(os.getenv("HELIX_PROBE_TIMEOUT_S", "90"))


def _auth_headers() -> dict:
    """Return an Authorization header if HELIX_BASIC_AUTH ('user:pass') is set."""
    creds = os.getenv("HELIX_BASIC_AUTH", "")
    if not creds:
        return {}
    token = base64.b64encode(creds.encode()).decode()
    return {"Authorization": f"Basic {token}"}


_FORBIDDEN = [
    "Error: No sequences provided",
    "Error: At least 2 sequences",
    "at least 2 sequences are required",
    "internal server error",
    "HTTP 500",
    "HTTP 404",
    "Pipeline Execution Complete\nThe approved plan has been executed.\n\nSteps:\n\nstep1 (handle_natural_command",
]


def _post_json(url: str, payload: dict, timeout: int = TIMEOUT_S) -> dict:
    body = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url, data=body, headers={"Content-Type": "application/json", **_auth_headers()}, method="POST"
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def _post_sse(url: str, payload: dict, timeout: int = TIMEOUT_S) -> Optional[Dict[str, Any]]:
    """POST to /execute/stream and return the final 'result' event data."""
    body = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url, data=body, headers={"Content-Type": "application/json", **_auth_headers()}, method="POST"
    )
    result = None
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        buffer = b""
        for chunk in iter(lambda: resp.read(4096), b""):
            buffer += chunk
            while b"\n\n" in buffer:
                event_bytes, buffer = buffer.split(b"\n\n", 1)
                for line in event_bytes.splitlines():
                    if not line.startswith(b"data: "):
                        continue
                    try:
                        event = json.loads(line[6:].decode("utf-8"))
                    except Exception:
                        continue
                    if event.get("type") == "result":
                        return event.get("data", {})
                    if event.get("type") == "error":
                        return {"error": event.get("detail", "sse_error")}
    return result


def _flatten(result: dict) -> str:
    parts: List[str] = []
    for k in ("text", "output", "result", "advisory", "message", "explanation"):
        v = result.get(k)
        if isinstance(v, str):
            parts.append(v)
        elif isinstance(v, dict):
            parts.append(json.dumps(v))
    inner = result.get("result") or result.get("data") or {}
    if isinstance(inner, dict):
        for k in ("text", "output", "advisory"):
            if isinstance(inner.get(k), str):
                parts.append(inner[k])
    return " ".join(parts)


def _is_acceptable_plan_card(result: dict, text_lower: str) -> bool:
    """A multi-step request that produces a plan card waiting for approval is a PASS.

    The system is correctly staging — the user will then click 'Approve & Run'.
    """
    if "pipeline plan" in text_lower or "approve" in text_lower or "execute pipeline" in text_lower:
        return True
    tool = (result or {}).get("tool", "")
    if tool == "__plan__":
        return True
    return False


def _create_session(base: str) -> str:
    r = _post_json(f"{base}/create_session", {})
    return r["session_id"]


# Minimal DE table matching the tabular_qa demo scenario
_DEMO_DE_CSV = """\
gene,log2FC,padj,baseMean
PDCD1,-3.21,0.0000081,445.3
CD274,2.87,0.000023,678.4
HAVCR2,-2.64,0.000041,312.7
TIGIT,-2.45,0.000098,289.1
LAG3,-2.33,0.00015,334.6
CD226,1.98,0.0012,567.9
CTLA4,-1.87,0.0021,234.5
BTLA,-1.72,0.0038,198.2
CD80,2.41,0.00067,412.8
CD86,2.15,0.0015,389.3
ICOS,1.63,0.0088,276.4
KRAS,0.88,0.062,789.3
AKT1,-0.54,0.21,334.6
""".encode("utf-8")


def _upload_csv_to_session(base: str, session_id: str, filename: str, content: bytes, timeout: int = TIMEOUT_S) -> dict:
    """Upload a CSV file to a session via multipart/form-data POST."""
    boundary = "----HelixProbeBoundary7MA4YWxkTrZu0gW"
    body = (
        f"--{boundary}\r\n"
        f'Content-Disposition: form-data; name="files"; filename="{filename}"\r\n'
        f"Content-Type: text/csv\r\n\r\n"
    ).encode("utf-8") + content + f"\r\n--{boundary}--\r\n".encode("utf-8")

    req = urllib.request.Request(
        f"{base}/session/{session_id}/uploads",
        data=body,
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}", **_auth_headers()},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


# ─────────────────────────────────────────────────────────────────────────────
# Demo definitions — one per scenario card
# ─────────────────────────────────────────────────────────────────────────────

DEMOS = [
    # ── Demo 6: tabular_qa — upload CSV then ask (no inline data) ────────────
    {
        "id": "demo6.tabular_qa",
        "upload_csv": ("pan_cancer_tcell_de_results.csv", _DEMO_DE_CSV),
        "initial": (
            "What are the top 10 therapeutic T cell targets ranked by absolute log2FC "
            "that are statistically significant (padj < 0.05)? Use the uploaded DE table."
        ),
        "expect_initial_keywords": ["brca1", "pdcd1", "cd274", "havcr2", "tigit", "lag3",
                                     "gene", "log2fc", "padj", "significant", "rank"],
        "followup": None,
    },
    {
        "id": "demo1.bulk_rnaseq",
        "initial": (
            "You are analyzing an RNA-seq transcriptome dataset from a mouse study investigating "
            "Toxoplasma gondii infection. 2x2 factorial: Infection (Inf/Uninf) × Time (11dpi/33dpi). "
            "Perform DE analysis. Include PCA, sample distance heatmap."
        ),
        "expect_initial_keywords": ["rna", "differential", "gene", "expression", "factorial"],
        "followup": (
            "Run bulk RNA-seq differential expression analysis with these inputs:\n"
            "count_matrix: s3://noricum-ngs-data/demo/rnaseq/tgondii_counts.csv\n"
            "sample_metadata: s3://noricum-ngs-data/demo/rnaseq/tgondii_metadata.csv\n"
            "design_formula: ~infection_status + time_point + infection_status:time_point"
        ),
        "expect_followup_keywords": ["rnaseq", "deseq", "differential", "expression", "gene"],
    },
    {
        "id": "demo2.scrna",
        "initial": (
            "You are analyzing a single-cell RNA-seq dataset from PBMCs comparing SLE patients "
            "to healthy controls. 10x v3, 5 vs 5. Objectives: QC, normalize, PCA→UMAP, Leiden cluster, "
            "annotate cell types, pseudobulk DEGs."
        ),
        "expect_initial_keywords": ["single", "cell", "umap", "cluster", "scrna"],
        "followup": (
            "data_file: s3://noricum-ngs-data/demo/scrna/sle_pbmc_filtered_feature_bc_matrix.h5\n"
            "data_format: 10x\nresolution: 0.5\nsteps: all"
        ),
        "expect_followup_keywords": ["cell", "cluster", "umap", "scrna", "single"],
    },
    {
        "id": "demo3.amplicon",
        "initial": (
            "You are processing a 16S rRNA amplicon sequencing dataset. Illumina MiSeq 2x250 paired-end. "
            "Run FastQC, trim adapters, merge paired-end reads, generate QC report."
        ),
        "expect_initial_keywords": ["16s", "amplicon", "quality", "fastqc", "trim", "merge"],
        "followup": None,
    },
    {
        "id": "demo4.apap",
        "initial": (
            "I have a bulk RNA-seq time-course study from a murine acetaminophen liver-injury model. "
            "5 time points, 4 reps each. Design the analysis workflow before execution. Include checkpoints "
            "for QC and biological interpretation."
        ),
        "expect_initial_keywords": ["workflow", "time", "course", "design", "analysis"],
        "followup": None,
    },
    {
        "id": "demo5.phylo",
        "initial": (
            "You are conducting a comparative evolutionary analysis of the SARS-CoV-2 spike protein across "
            "major variants of concern. Retrieve sequences, MSA, ML tree with bootstrap, RBD mutation "
            "annotation, pairwise identity matrix."
        ),
        "expect_initial_keywords": ["phylogen", "tree", "spike", "variant", "alignment"],
        "followup": (
            "Build a phylogenetic tree and pairwise amino acid identity matrix for the SARS-CoV-2 spike "
            "protein across 8 variants of concern using these pre-aligned sequences:\n"
            "sequences: s3://noricum-ngs-data/demo/phylo/sars_cov2_spike.fasta"
        ),
        "expect_followup_keywords": ["tree", "phylogen", "spike", "variant", "identity", "newick"],
    },
]


def main() -> int:
    base = os.getenv("HELIX_BASE_URL", DEFAULT_BASE).rstrip("/")
    out_path = sys.argv[1] if len(sys.argv) > 1 else None

    print(f"# Demo end-to-end probe — base={base}")
    print()

    all_results = []
    failures = 0

    for demo in DEMOS:
        sess = _create_session(base)
        print(f"## {demo['id']}  (session={sess[:8]})")

        # ── 0. optional file upload (before the initial prompt) ───────
        if demo.get("upload_csv"):
            filename, content = demo["upload_csv"]
            try:
                up = _upload_csv_to_session(base, sess, filename, content)
                n_files = len(up.get("files", []))
                print(f"  upload:   OK  ({n_files} file(s) registered)")
            except Exception as exc:
                print(f"  upload:   FAIL  {exc}")
                failures += 1
                all_results.append({"demo": demo["id"], "stage": "upload", "ok": False, "error": str(exc)})
                continue

        # ── 1. initial prompt ─────────────────────────────────────────
        t0 = time.time()
        try:
            r = _post_sse(f"{base}/execute/stream",
                          {"command": demo["initial"], "session_id": sess})
        except Exception as exc:
            print(f"  initial: HTTP error {exc}")
            failures += 1
            all_results.append({"demo": demo["id"], "stage": "initial", "ok": False, "error": str(exc)})
            continue
        dt = int((time.time() - t0) * 1000)

        text = _flatten(r or {}).lower()
        forbidden_hit = next((p for p in _FORBIDDEN if p.lower() in text), None)
        kw_match = any(kw in text for kw in demo["expect_initial_keywords"])
        plan_card_ok = _is_acceptable_plan_card(r or {}, text)
        ok = forbidden_hit is None and (kw_match or plan_card_ok)
        status = "PASS" if ok else "FAIL"

        tool = (r or {}).get("tool", "?")
        print(f"  initial:  {status}  tool={tool}  ({dt} ms)")
        if not ok:
            print(f"    forbidden_hit: {forbidden_hit!r}  kw_match: {kw_match}")
            print(f"    text(0..200): {text[:200]!r}")
        all_results.append({
            "demo": demo["id"], "stage": "initial", "ok": ok, "tool": tool,
            "elapsed_ms": dt, "forbidden_hit": forbidden_hit, "kw_match": kw_match,
        })
        if not ok:
            failures += 1

        # ── 2. followup ──────────────────────────────────────────────
        if demo["followup"]:
            t0 = time.time()
            try:
                r2 = _post_sse(f"{base}/execute/stream",
                               {"command": demo["followup"], "session_id": sess})
            except Exception as exc:
                print(f"  followup: HTTP error {exc}")
                failures += 1
                all_results.append({"demo": demo["id"], "stage": "followup", "ok": False, "error": str(exc)})
                continue
            dt2 = int((time.time() - t0) * 1000)
            text2 = _flatten(r2 or {}).lower()
            forbidden2 = next((p for p in _FORBIDDEN if p.lower() in text2), None)
            kw2 = any(kw in text2 for kw in demo["expect_followup_keywords"])
            plan_card2 = _is_acceptable_plan_card(r2 or {}, text2)
            ok2 = forbidden2 is None and (kw2 or plan_card2)
            status2 = "PASS" if ok2 else "FAIL"
            tool2 = (r2 or {}).get("tool", "?")
            print(f"  followup: {status2}  tool={tool2}  ({dt2} ms)")
            if not ok2:
                print(f"    forbidden_hit: {forbidden2!r}  kw_match: {kw2}")
                print(f"    text(0..200): {text2[:200]!r}")
            all_results.append({
                "demo": demo["id"], "stage": "followup", "ok": ok2, "tool": tool2,
                "elapsed_ms": dt2, "forbidden_hit": forbidden2, "kw_match": kw2,
            })
            if not ok2:
                failures += 1
        print()

    print(f"Results: {len(all_results) - failures} passed, {failures} failed (of {len(all_results)})")

    if out_path:
        with open(out_path, "w") as f:
            json.dump({"base": base, "results": all_results, "failures": failures}, f, indent=2)
        print(f"Wrote {out_path}")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
