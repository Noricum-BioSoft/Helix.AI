"""
bio_pipeline — Bioinformatics iteration loop.

Mirrors the existing ds_pipeline architecture but for bioinformatics tools
(bulk RNA-seq, single-cell, phylogenetics, amplicon QC).

Usage:
    from bio_pipeline import BioOrchestrator
    orch = BioOrchestrator(dispatch_tool, history_manager, session_id)
    result = await orch.run(tool_name, params, parent_run_id=prior_run_id)
"""
from .bio_orchestrator import BioOrchestrator
from .bio_evaluator import BioEvaluator
from .bio_reviewer import BioReviewer
from .run_config import BioRunConfig

__all__ = ["BioOrchestrator", "BioEvaluator", "BioReviewer", "BioRunConfig"]
