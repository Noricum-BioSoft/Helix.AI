"""Local Nextflow pipeline executor for Helix.AI.

Runs nf-core / custom Nextflow workflows as subprocesses on the local machine
using Docker containers.  No Seqera Platform dependency.

Execution tiers
---------------
Tier 3a  local_docker   — Nextflow + Docker on this server (small datasets)
Tier 3b  aws_batch      — Nextflow + AWS Batch (deferred; not yet implemented)

Decision heuristic (choose_tier):
- sample_count <= TIER3A_MAX_SAMPLES AND total_size_gb < TIER3A_MAX_GB → local_docker
- otherwise → aws_batch (returns a friendly error until Tier 3b is implemented)

Weblog integration
------------------
Every launched pipeline is started with::

    nextflow run ... -with-weblog http://localhost:{PORT}/internal/nextflow/events

This causes Nextflow to POST JSON events to Helix at each workflow lifecycle
change (started, process_submitted, process_completed, completed, failed).
The weblog receiver in main.py forwards events to NextflowEventBus, which
fans them out to SSE subscribers.

Job IDs are injected as a Nextflow param (``--helix_job_id``) and as the
``runName`` so the weblog receiver can correlate events back to job_id.
"""

from __future__ import annotations

import asyncio
import json
import logging
import os
import shutil
import subprocess
import tempfile
import uuid
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

NEXTFLOW_BIN = os.getenv(
    "HELIX_NEXTFLOW_BIN",
    str(Path.home() / ".local" / "bin" / "nextflow"),
)
HELIX_PORT = int(os.getenv("PORT", "8001"))
WORK_BASE = Path(os.getenv("HELIX_NEXTFLOW_WORK_DIR", "/tmp/helix_nextflow"))
RESULTS_BASE = Path(os.getenv("HELIX_NEXTFLOW_RESULTS_DIR", "/tmp/helix_results"))
NF_CONFIG = Path(__file__).resolve().parent.parent / "workflows" / "nextflow.config"

# Tier-3a limits: run locally if dataset is small enough
TIER3A_MAX_SAMPLES = int(os.getenv("HELIX_NF_LOCAL_MAX_SAMPLES", "4"))
TIER3A_MAX_GB = float(os.getenv("HELIX_NF_LOCAL_MAX_GB", "10.0"))

# ---------------------------------------------------------------------------
# Pipeline registry — maps Helix tool names to Nextflow pipeline identifiers
# ---------------------------------------------------------------------------

# Format: (pipeline_path_or_nfcore_name, default_profile)
HELIX_TO_PIPELINE: Dict[str, Tuple[str, str]] = {
    "chip_seq_analysis":      ("workflows/chip_seq.nf",   "docker"),
    "atac_seq_analysis":      ("nf-core/atacseq",         "docker"),
    "genome_assembly":        ("nf-core/mag",             "docker"),
    "variant_calling":        ("nf-core/sarek",           "docker"),
    "metagenomics_16s":       ("nf-core/ampliseq",        "docker"),
    "metagenomics_shotgun":   ("nf-core/taxprofiler",     "docker"),
    "rna_splicing_isoform":   ("nf-core/rnasplice",       "docker"),
    "crispr_screen_analysis": ("nf-core/crisprseq",       "docker"),
}

# ---------------------------------------------------------------------------
# Prerequisites check
# ---------------------------------------------------------------------------

_prereq_cache: Optional[Dict[str, bool]] = None


def check_prerequisites() -> Dict[str, bool]:
    """Return availability of each required binary.

    Results are cached after first call.
    """
    global _prereq_cache
    if _prereq_cache is not None:
        return _prereq_cache

    results: Dict[str, bool] = {}
    checks = {
        "java":     ["java", "-version"],
        "nextflow": [NEXTFLOW_BIN, "-version"],
        "docker":   ["docker", "--version"],
    }
    for name, cmd in checks.items():
        try:
            result = subprocess.run(cmd, capture_output=True, timeout=10)
            results[name] = result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            results[name] = False

    _prereq_cache = results
    missing = [k for k, v in results.items() if not v]
    if missing:
        logger.warning("Nextflow prerequisites missing: %s", missing)
    else:
        logger.info("Nextflow prerequisites OK: %s", list(results.keys()))
    return results


def _prereqs_ok() -> Tuple[bool, str]:
    """Return (ok, error_message)."""
    status = check_prerequisites()
    missing = [k for k, v in status.items() if not v]
    if not missing:
        return True, ""
    return False, (
        f"Nextflow pipeline execution requires {', '.join(missing)} to be installed. "
        f"Run `scripts/setup-nextflow.sh` to install missing dependencies."
    )


# ---------------------------------------------------------------------------
# Tier decision
# ---------------------------------------------------------------------------

def choose_tier(arguments: Dict[str, Any]) -> str:
    """Choose execution tier based on dataset size hints in *arguments*.

    Returns 'local_docker' or 'aws_batch'.
    The caller is responsible for handling 'aws_batch' (not yet implemented).
    """
    sample_count = int(arguments.get("sample_count") or arguments.get("n_samples") or 1)
    size_gb = float(arguments.get("estimated_size_gb") or 0.0)

    if sample_count <= TIER3A_MAX_SAMPLES and size_gb < TIER3A_MAX_GB:
        return "local_docker"
    return "aws_batch"


# ---------------------------------------------------------------------------
# Parameter builders — translate Helix tool arguments to nf-core params
# ---------------------------------------------------------------------------

def _build_chip_seq_params(arguments: Dict[str, Any], result_dir: Path) -> Dict[str, Any]:
    return {
        "treatment_bam": arguments.get("treatment_bam") or arguments.get("input_bam", ""),
        "control_bam":   arguments.get("control_bam", ""),
        "genome_size":   arguments.get("genome") or arguments.get("genome_size", "hs"),
        "peak_type":     arguments.get("peak_type", "narrow"),
        "outdir":        str(result_dir),
    }


def _build_atac_seq_params(arguments: Dict[str, Any], result_dir: Path) -> Dict[str, Any]:
    """nf-core/atacseq parameter mapping."""
    return {
        "input":   arguments.get("bam_files") or arguments.get("input", ""),
        "genome":  arguments.get("genome", "GRCh38"),
        "outdir":  str(result_dir),
    }


def _build_genome_assembly_params(arguments: Dict[str, Any], result_dir: Path) -> Dict[str, Any]:
    """nf-core/mag parameter mapping."""
    return {
        "input":       arguments.get("reads") or arguments.get("input", ""),
        "outdir":      str(result_dir),
    }


def _build_variant_calling_params(arguments: Dict[str, Any], result_dir: Path) -> Dict[str, Any]:
    """nf-core/sarek parameter mapping."""
    return {
        "input":   arguments.get("bam_files") or arguments.get("input", ""),
        "genome":  arguments.get("reference_genome") or arguments.get("genome", "GRCh38"),
        "tools":   "haplotypecaller",
        "outdir":  str(result_dir),
    }


def _build_metagenomics_16s_params(arguments: Dict[str, Any], result_dir: Path) -> Dict[str, Any]:
    """nf-core/ampliseq parameter mapping."""
    return {
        "input":        arguments.get("fastq_dir") or arguments.get("input", ""),
        "metadata":     arguments.get("sample_metadata", ""),
        "FW_primer":    arguments.get("primer_f", "GTGYCAGCMGCCGCGGTAA"),
        "RV_primer":    arguments.get("primer_r", "GGACTACNVGGGTWTCTAAT"),
        "outdir":       str(result_dir),
    }


def _build_metagenomics_shotgun_params(arguments: Dict[str, Any], result_dir: Path) -> Dict[str, Any]:
    """nf-core/taxprofiler parameter mapping."""
    return {
        "input":   arguments.get("fastq_files") or arguments.get("input", ""),
        "outdir":  str(result_dir),
    }


def _build_rna_splicing_params(arguments: Dict[str, Any], result_dir: Path) -> Dict[str, Any]:
    """nf-core/rnasplice parameter mapping."""
    return {
        "input":       arguments.get("bam_files") or arguments.get("input", ""),
        "contrasts":   arguments.get("contrasts", ""),
        "genome":      arguments.get("genome", "GRCh38"),
        "outdir":      str(result_dir),
    }


def _build_crispr_screen_params(arguments: Dict[str, Any], result_dir: Path) -> Dict[str, Any]:
    """nf-core/crisprseq parameter mapping."""
    return {
        "input":   arguments.get("count_table") or arguments.get("input", ""),
        "outdir":  str(result_dir),
    }


def _build_params(tool_name: str, arguments: Dict[str, Any], result_dir: Path) -> Dict[str, Any]:
    builders = {
        "chip_seq_analysis":      _build_chip_seq_params,
        "atac_seq_analysis":      _build_atac_seq_params,
        "genome_assembly":        _build_genome_assembly_params,
        "variant_calling":        _build_variant_calling_params,
        "metagenomics_16s":       _build_metagenomics_16s_params,
        "metagenomics_shotgun":   _build_metagenomics_shotgun_params,
        "rna_splicing_isoform":   _build_rna_splicing_params,
        "crispr_screen_analysis": _build_crispr_screen_params,
    }
    builder = builders.get(tool_name)
    if builder:
        return builder(arguments, result_dir)
    # Generic fallback: pass all non-internal arguments through
    _skip = {"session_id", "session_context", "command", "original_command", "_from_broker"}
    params: Dict[str, Any] = {k: v for k, v in arguments.items() if k not in _skip and v}
    params["outdir"] = str(result_dir)
    return params


# ---------------------------------------------------------------------------
# Core launcher
# ---------------------------------------------------------------------------

async def launch_pipeline(
    tool_name: str,
    arguments: Dict[str, Any],
    session_id: str,
) -> Dict[str, Any]:
    """Launch a Nextflow pipeline for *tool_name* and return immediately.

    Returns::

        {
            "job_id":   "<uuid>",
            "status":   "submitted",
            "tool":     "<tool_name>",
            "pipeline": "<pipeline>",
            "tier":     "local_docker",
            "stream_url": "/jobs/<job_id>/stream",
        }

    The pipeline runs asynchronously.  The caller should instruct the
    frontend to open ``EventSource("/jobs/<job_id>/stream")`` to receive
    progress events.
    """
    ok, err_msg = _prereqs_ok()
    if not ok:
        return {
            "status": "error",
            "tool": tool_name,
            "text": err_msg,
        }

    pipeline, profile = HELIX_TO_PIPELINE.get(tool_name, ("", "docker"))
    if not pipeline:
        return {
            "status": "error",
            "tool": tool_name,
            "text": f"No pipeline registered for tool '{tool_name}'.",
        }

    # Per-tool input validation before spinning up any resources
    _required: Dict[str, str] = {
        "chip_seq_analysis":      "treatment_bam",
        "atac_seq_analysis":      "bam_files",
        "genome_assembly":        "reads",
        "variant_calling":        "bam_files",
        "metagenomics_16s":       "fastq_dir",
        "metagenomics_shotgun":   "fastq_files",
        "rna_splicing_isoform":   "bam_files",
        "crispr_screen_analysis": "count_table",
    }
    _req_key = _required.get(tool_name)
    _alt_keys: Dict[str, list] = {
        "chip_seq_analysis": ["treatment_bam", "input_bam"],
    }
    _check_keys = _alt_keys.get(tool_name, [_req_key] if _req_key else [])
    if _req_key and not any(arguments.get(k) for k in _check_keys):
        return {
            "status": "needs_inputs",
            "tool": tool_name,
            "text": (
                f"Please provide the required input **{_req_key}** "
                f"(and optionally any other parameters) to launch the {tool_name} pipeline."
            ),
        }

    tier = choose_tier(arguments)
    if tier == "aws_batch":
        return {
            "status": "error",
            "tool": tool_name,
            "text": (
                "This dataset exceeds the local execution limits "
                f"({TIER3A_MAX_SAMPLES} samples / {TIER3A_MAX_GB} GB). "
                "AWS Batch execution (Tier 3b) is not yet enabled. "
                "Reduce the dataset size or contact your administrator to enable cloud execution."
            ),
        }

    job_id = str(uuid.uuid4())
    work_dir = WORK_BASE / job_id / "work"
    result_dir = RESULTS_BASE / job_id
    work_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)

    # Persist job record first so /jobs/{id} is queryable immediately
    from backend.job_manager import get_job_manager
    get_job_manager().create_nextflow_job(
        job_id=job_id,
        tool_name=tool_name,
        pipeline=pipeline,
        session_id=session_id,
        work_dir=str(work_dir),
        result_dir=str(result_dir),
    )

    # Launch subprocess in a background task — do not await
    asyncio.create_task(
        _run_nextflow(job_id, tool_name, pipeline, profile, arguments, work_dir, result_dir)
    )

    return {
        "job_id":     job_id,
        "status":     "submitted",
        "tool":       tool_name,
        "pipeline":   pipeline,
        "tier":       tier,
        "stream_url": f"/jobs/{job_id}/stream",
        "results_url": f"/jobs/{job_id}/results",
        "text": (
            f"Pipeline **{pipeline}** submitted (job `{job_id[:8]}…`). "
            f"Connect to `{f'/jobs/{job_id}/stream'}` for live progress."
        ),
    }


# ---------------------------------------------------------------------------
# Subprocess runner (runs in background task)
# ---------------------------------------------------------------------------

async def _run_nextflow(
    job_id: str,
    tool_name: str,
    pipeline: str,
    profile: str,
    arguments: Dict[str, Any],
    work_dir: Path,
    result_dir: Path,
) -> None:
    """Execute ``nextflow run`` and update job state via JobManager."""
    from backend.job_manager import get_job_manager
    from backend.nextflow_event_bus import get_event_bus

    jm = get_job_manager()
    bus = get_event_bus()
    log_path = work_dir.parent / "nextflow.log"

    try:
        params = _build_params(tool_name, arguments, result_dir)
        weblog_url = f"http://localhost:{HELIX_PORT}/internal/nextflow/events"

        # Resolve pipeline path: local .nf file or nf-core reference
        project_root = Path(__file__).resolve().parent.parent
        if pipeline.startswith("workflows/"):
            pipeline_arg = str(project_root / pipeline)
        else:
            pipeline_arg = pipeline

        # Build nextflow command
        # Nextflow run names must match ^[a-z][a-z\d_-]{0,79}$ — UUIDs can start
        # with a digit, so we prefix with "helix-" and truncate to 80 chars.
        nf_run_name = f"helix-{job_id}"[:80]

        cmd = [
            NEXTFLOW_BIN, "run", pipeline_arg,
            "-w", str(work_dir),
            "-name", nf_run_name,            # run name used for weblog correlation
            "-with-weblog", weblog_url,
            "--helix_job_id", job_id,        # param-level job_id for double-lookup
        ]
        # Use custom nextflow.config when present (enables Docker globally without a profile)
        # Only add -profile if the profile is not "docker" or "local" since those
        # are managed via nextflow.config rather than named profiles.
        if profile not in ("docker", "local", ""):
            cmd += ["-profile", profile]

        # Append pipeline parameters as --key value pairs
        for key, value in params.items():
            if value is not None and value != "":
                cmd += [f"--{key}", str(value)]

        # Use custom config if present
        if NF_CONFIG.exists():
            cmd += ["-c", str(NF_CONFIG)]

        logger.info("Launching Nextflow job %s: %s", job_id, " ".join(cmd))
        jm.set_nextflow_job_running(job_id, pid=0)  # pid updated below

        with open(log_path, "w") as log_fh:
            proc = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=log_fh,
                stderr=asyncio.subprocess.STDOUT,
                cwd=str(project_root),
            )
            jm.set_nextflow_job_running(job_id, pid=proc.pid)
            await proc.wait()

        if proc.returncode == 0:
            # Parse result files from result_dir
            result_files = _collect_results(result_dir)
            jm.set_nextflow_job_completed(job_id, result_dir=str(result_dir), result_files=result_files)
            await bus.publish(job_id, {
                "type": "completed",
                "job_id": job_id,
                "result_dir": str(result_dir),
                "files": result_files,
            })
            bus.clear_history(job_id)
        else:
            # Read last lines of log for error context
            try:
                with open(log_path) as f:
                    tail = "".join(f.readlines()[-20:])
            except Exception:
                tail = "(log unavailable)"
            error_msg = f"Nextflow exited with code {proc.returncode}.\n\nLog tail:\n{tail}"
            jm.set_nextflow_job_failed(job_id, error=error_msg)
            await bus.publish(job_id, {"type": "failed", "job_id": job_id, "error": error_msg})

    except Exception as exc:
        logger.exception("Unexpected error running Nextflow job %s", job_id)
        error_msg = f"Internal error: {exc}"
        try:
            jm.set_nextflow_job_failed(job_id, error=error_msg)
        except Exception:
            pass
        await bus.publish(job_id, {"type": "failed", "job_id": job_id, "error": error_msg})


def _collect_results(result_dir: Path) -> list:
    """Return a list of output file paths relative to *result_dir*."""
    if not result_dir.exists():
        return []
    files = []
    for f in result_dir.rglob("*"):
        if f.is_file() and not f.name.startswith("."):
            files.append(str(f.relative_to(result_dir)))
    return sorted(files)
