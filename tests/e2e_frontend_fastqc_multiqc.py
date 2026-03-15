#!/usr/bin/env python3
"""
End-to-end test: 16S demo + FastQC/MultiQC follow-ups via localhost frontend.

Prerequisites:
  1. Backend: HELIX_MOCK_MODE=1 uvicorn backend.main_with_mcp:app --host 0.0.0.0 --port 8001
  2. Frontend: cd frontend && npm run dev  (runs on http://localhost:5173)

Run: python tests/e2e_frontend_fastqc_multiqc.py
"""

import asyncio
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from playwright.sync_api import sync_playwright, expect


DEMO_PROMPT = """You are processing a 16S rRNA amplicon sequencing dataset from a gut microbiome study. Raw paired-end FASTQ files are on S3 and need a full preprocessing pipeline before downstream diversity analysis.

Input files:
- R1: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq
- R2: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq

Pipeline steps (run in order):
1. FastQC quality assessment on both R1 and R2 files.
2. Read trimming to remove low-quality bases and adapter sequences.
3. Read merging of trimmed paired-end reads.
4. Final quality report.

Output prefix: s3://noricum-ngs-data/test-output/amplicon-demo/"""

# (prompt, must_contain, must_not_contain)
FOLLOWUPS = [
    ("Run FastQC quality assessment on the merged reads.", "FASTA", None),
    ("perform FastQC analysis on s3://noricum-ngs-data/test-output/amplicon-demo/merged.fasta", "FASTA", None),
    ("Run FastQC on the raw reads from the pipeline.", None, "FastQC cannot run on merged"),  # Should run, not merged error
    ("run MultiQC on s3://noricum-ngs-data/test-output/amplicon-demo/merged.fasta", "not supported", None),
]


def run():
    frontend_url = os.environ.get("HELIX_FRONTEND_URL", "http://localhost:5173")
    print(f"Frontend: {frontend_url}")
    print("=" * 60)

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=os.environ.get("HEADLESS", "1") == "1")
        page = browser.new_page()
        page.set_default_timeout(60000)

        try:
            page.goto(frontend_url)
            page.wait_for_load_state("networkidle", timeout=15000)
        except Exception as e:
            print(f"Failed to load frontend: {e}")
            print("Ensure frontend is running: cd frontend && npm run dev")
            browser.close()
            return 1

        # Find prompt input (main command textarea)
        textarea = page.locator("textarea").first
        textarea.wait_for(state="visible", timeout=5000)

        def submit_and_wait():
            # Submit: try Ctrl+Enter (or Cmd+Enter on Mac), fallback to button click
            textarea.press("Control+Enter")
            # Wait for loading to finish
            page.wait_for_function(
                "() => !document.body.innerText.includes('Processing...') && !document.body.innerText.includes('Submitting…')",
                timeout=120000,
            )
            page.wait_for_timeout(2000)

        # 1. Run demo
        print("\n1. Submitting 16S amplicon demo...")
        textarea.fill(DEMO_PROMPT)
        submit_and_wait()
        print("   Demo completed.")

        # 2-5. Follow-up prompts
        # Each: (prompt, must_contain, must_not_contain) - use None to skip that check
        results = []
        for i, (prompt, must_contain, must_not_contain) in enumerate(FOLLOWUPS, 2):
            print(f"\n{i}. Submitting: {prompt[:60]}...")
            # Count history items before submit
            history_count_before = page.locator(".mt-4").count()
            textarea.fill(prompt)
            submit_and_wait()

            # Get only the latest response (last history item)
            history_items = page.locator(".mt-4")
            if history_items.count() > history_count_before:
                last_item = history_items.nth(history_items.count() - 1)
                response_text = last_item.inner_text()
            else:
                response_text = page.locator("body").inner_text()

            ok = True
            if must_contain:
                ok = must_contain in response_text
            if must_not_contain and ok:
                ok = must_not_contain not in response_text
            results.append((prompt[:50], ok, response_text[-500:] if len(response_text) > 500 else response_text))
            print(f"   ok={ok} | preview={response_text[-150:][:120]}...")

        browser.close()

    # Summary
    print("\n" + "=" * 60)
    all_ok = all(r[1] for r in results)
    print(f"Result: {'PASS' if all_ok else 'FAIL'}")
    for prompt, ok, _ in results:
        print(f"  {prompt}... : {'✓' if ok else '✗'}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(run())
