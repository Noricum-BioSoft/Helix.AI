"""
DemoDispatcher — Routes known demo prompts to real bioinformatics tools.

Each handler delegates to BioOrchestrator.run() so every demo execution:
  1. Calls a real Python-native tool (no hardcoded responses).
  2. Is recorded in the session run ledger (parent_run_id linkage).
  3. Produces a narrative summary + delta metrics for iterative workflows.

Extracted from execute() in main_with_mcp.py so each demo handler can be
tested in isolation without a live HTTP server or AWS credentials.
"""
from __future__ import annotations

import logging
import re
from typing import Any, Callable, Dict, Optional, Tuple

logger = logging.getLogger(__name__)

# Canonical S3 path used by the SLE demo — kept as a module constant so tests
# can reference it directly.
SLE_DATA_PATH = (
    "s3://noricum-ngs-data/demo/scrna/sle_pbmc_filtered_feature_bc_matrix.h5"
)


class DemoDispatcher:
    """
    Checks a command string against known demo prompts in priority order and
    returns a ``(tool_name, result_dict)`` pair when matched, or ``None``.

    Parameters
    ----------
    presigned_cache:
        Legacy parameter kept for API compatibility; no longer used.
    tool_executor:
        Async callable ``(tool_name: str, params: dict) -> dict`` — normally
        ``call_mcp_tool`` from ``main_with_mcp``.  Injected so tests can
        replace it with a synchronous mock without touching module globals.
    orchestrator:
        Optional ``BioOrchestrator`` instance.  When provided, all tool
        calls are routed through it so run tracking and delta evaluation
        happen automatically.  When None, tool_executor is called directly.
    session_id:
        Optional session identifier forwarded to the orchestrator.
    parent_run_id:
        Optional parent run ID for iterative workflows.
    """

    def __init__(
        self,
        presigned_cache: Dict[str, str],
        tool_executor: Callable,
        orchestrator: Optional[Any] = None,
        session_id: Optional[str] = None,
        parent_run_id: Optional[str] = None,
    ) -> None:
        self._cache = presigned_cache  # kept for backward compatibility
        self._raw_exec = tool_executor
        self._orch = orchestrator
        self._session_id = session_id
        self._parent_run_id = parent_run_id

    # ------------------------------------------------------------------ #
    # Internal dispatch helper                                             #
    # ------------------------------------------------------------------ #

    async def _exec(
        self,
        tool_name: str,
        params: Dict[str, Any],
        objective: str = "",
    ) -> Dict[str, Any]:
        """Call the tool via BioOrchestrator (if available) or raw executor."""
        if self._orch is not None:
            return await self._orch.run(
                tool_name,
                params,
                parent_run_id=self._parent_run_id,
                session_id=self._session_id,
                objective=objective,
            )
        return await self._raw_exec(tool_name, params)

    # ------------------------------------------------------------------ #
    # Public API                                                           #
    # ------------------------------------------------------------------ #

    async def handle(
        self,
        command: str,
        parent_run_id: Optional[str] = None,
    ) -> Optional[Tuple[str, Dict[str, Any]]]:
        """
        Return ``(tool_name, result_dict)`` on the first demo that matches
        *command*, or ``None`` if no demo matches.

        Parameters
        ----------
        command:
            The full command string from the user.
        parent_run_id:
            If set, overrides the instance-level parent_run_id for this
            single call (used to chain iterative demo re-runs).
        """
        if parent_run_id is not None:
            self._parent_run_id = parent_run_id

        cmd = command or ""
        cmd_lower = cmd.lower()

        # Demo 1: T. gondii 2×2 factorial RNA-seq (exact S3 path match)
        if (
            "s3://noricum-ngs-data/demo/rnaseq/tgondii_counts.csv" in cmd
            and "s3://noricum-ngs-data/demo/rnaseq/tgondii_metadata.csv" in cmd
        ):
            return await self._handle_tgondii(cmd)

        # Demo 4: APAP liver-injury time-course RNA-seq (exact S3 path match)
        if (
            "s3://noricum-ngs-data/demo/rnaseq/apap_timecourse_counts.csv" in cmd
            and "s3://noricum-ngs-data/demo/rnaseq/apap_timecourse_metadata.csv" in cmd
        ):
            return await self._handle_apap(cmd)

        # Demo 2: scRNA-seq SLE PBMC (exact S3 path match)
        if SLE_DATA_PATH in cmd:
            return await self._handle_sle()

        # Demo 5: SARS-CoV-2 spike protein phylogenetics (keyword match)
        if ("sars" in cmd_lower or "spike" in cmd_lower) and any(
            kw in cmd_lower
            for kw in ("variant", "phylogen", "alpha", "omicron", "wuhan")
        ):
            return await self._handle_phylo()

        # Demo 3: Amplicon QC pipeline (keyword match)
        if "helix-test-data" in cmd_lower and any(
            kw in cmd_lower
            for kw in (
                "fastqc",
                "trim",
                "amplicon",
                "microbiome",
                "16s",
                "forward_reads",
                "sample01",
            )
        ):
            return await self._handle_amplicon(cmd)

        return None

    # ------------------------------------------------------------------ #
    # Async handlers — all route to real tools via _exec / BioOrchestrator
    # ------------------------------------------------------------------ #

    async def _handle_tgondii(self, cmd: str) -> Tuple[str, Dict[str, Any]]:
        m = re.search(r"design_formula[:\s]+([^\n]+)", cmd)
        design_formula = (
            m.group(1).strip()
            if m
            else "~infection_status + time_point + infection_status:time_point"
        )
        result = await self._exec(
            "bulk_rnaseq_analysis",
            {
                "count_matrix": "s3://noricum-ngs-data/demo/rnaseq/tgondii_counts.csv",
                "sample_metadata": "s3://noricum-ngs-data/demo/rnaseq/tgondii_metadata.csv",
                "design_formula": design_formula,
                "alpha": 0.05,
            },
            objective="T. gondii infection vs uninfected — 2×2 factorial design",
        )
        return "bulk_rnaseq_analysis", result

    async def _handle_apap(self, cmd: str) -> Tuple[str, Dict[str, Any]]:
        m = re.search(r"design_formula[:\s]+([^\n]+)", cmd)
        design_formula = m.group(1).strip() if m else "~time_point"
        result = await self._exec(
            "bulk_rnaseq_analysis",
            {
                "count_matrix": "s3://noricum-ngs-data/demo/rnaseq/apap_timecourse_counts.csv",
                "sample_metadata": "s3://noricum-ngs-data/demo/rnaseq/apap_timecourse_metadata.csv",
                "design_formula": design_formula,
                "alpha": 0.05,
            },
            objective="APAP liver-injury time-course — expression recovery",
        )
        return "bulk_rnaseq_analysis", result

    async def _handle_sle(self) -> Tuple[str, Dict[str, Any]]:
        result = await self._exec(
            "single_cell_analysis",
            {
                "data_file": SLE_DATA_PATH,
                "data_format": "10x",
                "resolution": 0.5,
                "steps": "all",
            },
            objective="SLE PBMC immune profiling — cell-type DEG",
        )
        return "single_cell_analysis", result

    async def _handle_phylo(self) -> Tuple[str, Dict[str, Any]]:
        """Route Demo 5 to the real phylogenetic_tree tool with embedded SARS-CoV-2 sequences."""
        result = await self._exec(
            "phylogenetic_tree",
            {"aligned_sequences": ""},  # empty → module uses built-in SARS-CoV-2 dataset
            objective="SARS-CoV-2 spike protein — variant divergence phylogenetics",
        )
        return "phylogenetic_tree", result

    async def _handle_amplicon(self, cmd: str) -> Tuple[str, Dict[str, Any]]:
        """Route Demo 3 to fastqc_quality_analysis (has demo-mode fallback for helix-test buckets)."""
        result = await self._exec(
            "fastqc_quality_analysis",
            {
                "input_r1": "s3://helix-test-data/microbiome/run1/sample01_R1.fastq.gz",
                "input_r2": "s3://helix-test-data/microbiome/run1/sample01_R2.fastq.gz",
                "_from_broker": True,
            },
            objective="16S amplicon QC — V3–V4 gut microbiome pipeline",
        )
        return "fastqc_quality_analysis", result
