#!/usr/bin/env python3
"""
Test FastQC and MultiQC follow-up requests after the 16S amplicon demo.

Scenario:
1. Run 16S amplicon demo (creates session with pipeline history)
2. "Run FastQC quality assessment on the merged reads" -> expect clear explanation (merged FASTA incompatible)
3. "perform FastQC analysis on s3://noricum-ngs-data/test-output/amplicon-demo/merged.fasta" -> same
4. "Run FastQC on the raw reads" -> expect session-aware R1/R2 from pipeline
5. "run MultiQC on s3://noricum-ngs-data/test-output/amplicon-demo/merged.fasta" -> expect unsupported + alternatives
"""

import asyncio
import json
import os
import sys
from pathlib import Path

# Add project root
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# Use local sessions
os.environ.setdefault("HELIX_MOCK_MODE", "1")
os.environ.setdefault("HELIX_ANALYSIS_USE_SANDBOX", "false")
os.environ.setdefault("HELIX_SANDBOX_HOST_FALLBACK", "1")


async def run_test():
    from backend.command_router import CommandRouter
    from backend.history_manager import history_manager

    # Create session and seed with 16S demo history (simulating completed pipeline)
    session_id = "test-fastqc-multiqc-session"
    history_manager.ensure_session_exists(session_id)

    demo_prompt = """You are processing a 16S rRNA amplicon sequencing dataset from a gut microbiome study. Raw paired-end FASTQ files are on S3 and need a full preprocessing pipeline before downstream diversity analysis.

Input files:
- R1: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq
- R2: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq

Pipeline steps (run in order):
1. FastQC quality assessment on both R1 and R2 files.
2. Read trimming to remove low-quality bases and adapter sequences.
3. Read merging of trimmed paired-end reads.
4. Final quality report.

Output prefix: s3://noricum-ngs-data/test-output/amplicon-demo/"""

    # Seed session with demo as if it was just run (agent result with plan)
    history_manager.add_history_entry(
        session_id,
        demo_prompt,
        "agent",
        {
            "status": "success",
            "type": "plan_result",
            "steps": [
                {
                    "tool_name": "fastqc_quality_analysis",
                    "arguments": {
                        "input_r1": "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq",
                        "input_r2": "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq",
                    },
                },
                {
                    "tool_name": "read_merging",
                    "arguments": {
                        "forward_reads": "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq",
                        "reverse_reads": "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq",
                        "output": "s3://noricum-ngs-data/test-output/amplicon-demo/merged.fasta",
                    },
                },
            ],
        },
    )

    session_context = history_manager.sessions.get(session_id, {})
    router = CommandRouter()

    results = []

    # Test 1: "Run FastQC on merged reads" -> expect session_resolution_error (merged FASTA incompatible)
    cmd1 = "Run FastQC quality assessment on the merged reads."
    tool1, params1 = router.route_command(cmd1, session_context)
    results.append({
        "cmd": cmd1,
        "tool": tool1,
        "params": {k: v for k, v in (params1 or {}).items() if k != "session_context"},
        "expect": "session_resolution_error or needs_inputs with clear merged-FASTA message",
    })
    assert tool1 == "fastqc_quality_analysis", f"Expected fastqc_quality_analysis, got {tool1}"
    err = params1.get("session_resolution_error")
    assert err and "FastQC cannot run on merged" in err and "FASTA" in err, (
        f"Expected session_resolution_error about merged FASTA, got: {err}"
    )
    print("✓ Test 1: 'Run FastQC on merged reads' -> clear merged-FASTA explanation")

    # Test 2: "perform FastQC analysis on merged.fasta" -> same
    cmd2 = "perform FastQC analysis on s3://noricum-ngs-data/test-output/amplicon-demo/merged.fasta"
    tool2, params2 = router.route_command(cmd2, session_context)
    results.append({"cmd": cmd2, "tool": tool2, "params": params2})
    assert tool2 == "fastqc_quality_analysis"
    err2 = params2.get("session_resolution_error")
    assert err2 and "FASTA" in err2, f"Expected FASTA incompatibility, got: {err2}"
    print("✓ Test 2: 'FastQC on merged.fasta' -> FASTA incompatibility message")

    # Test 3: "Run FastQC on the raw reads" -> session-aware R1/R2
    cmd3 = "Run FastQC on the raw reads from the pipeline."
    tool3, params3 = router.route_command(cmd3, session_context)
    results.append({"cmd": cmd3, "tool": tool3, "params": params3})
    assert tool3 == "fastqc_quality_analysis"
    r1 = params3.get("input_r1", "")
    r2 = params3.get("input_r2", "")
    assert r1 and r2 and "test_mate_R1" in r1 and "test_mate_R2" in r2, (
        f"Expected session-aware R1/R2, got input_r1={r1!r} input_r2={r2!r}"
    )
    assert not params3.get("needs_inputs"), "Should have resolved R1/R2 from session"
    print("✓ Test 3: 'FastQC on raw reads' -> session-aware R1/R2 from pipeline")

    # Test 4: "run MultiQC" -> unsupported_tool with alternatives
    cmd4 = "run MultiQC on s3://noricum-ngs-data/test-output/amplicon-demo/merged.fasta"
    tool4, params4 = router.route_command(cmd4, session_context)
    results.append({"cmd": cmd4, "tool": tool4, "params": params4})
    assert tool4 == "unsupported_tool", f"Expected unsupported_tool, got {tool4}"
    assert params4.get("requested_tool") == "multiqc"
    print("✓ Test 4: 'run MultiQC' -> unsupported_tool with alternatives")

    # Test 5: Dispatch FastQC with invalid inputs -> validation error with clear message
    from backend.main import dispatch_tool
    bad_result = await dispatch_tool("fastqc_quality_analysis", {
        "input_r1": "s3://bucket/merged.fasta",
        "input_r2": "s3://bucket/merged.fasta",  # Same file + FASTA format
        "session_id": session_id,
    })
    text = bad_result.get("text") or ""
    assert ("FASTA" in text or "both" in text.lower()) and "What to do" in text, (
        f"Expected validation with clear guidance, got: {bad_result}"
    )
    print("✓ Test 5: FastQC dispatch with invalid inputs -> validation error with guidance")

    # Test 6: Dispatch unsupported MultiQC
    multi_result = await dispatch_tool("multiqc", {"command": "run MultiQC"})
    assert "not supported" in (multi_result.get("text") or "").lower()
    assert "Alternatives" in (multi_result.get("text") or "")
    print("✓ Test 6: MultiQC dispatch -> unsupported message with alternatives")

    return results


def test_fastqc_multiqc_followups():
    """Pytest entry point."""
    asyncio.run(run_test())


if __name__ == "__main__":
    out = asyncio.run(run_test())
    print("\n" + "=" * 60)
    print("All 6 tests passed.")
    print("=" * 60)
