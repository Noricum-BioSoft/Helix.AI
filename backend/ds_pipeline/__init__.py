"""
Data Science Pipeline — iterative plan→execute→evaluate→review loop.

Entry point:
    from backend.ds_pipeline.orchestrator import DataScienceOrchestrator, DSRunConfig
    orch = DataScienceOrchestrator(base_dir=".")
    run_data = orch.run(config)
"""
from backend.ds_pipeline.orchestrator import DataScienceOrchestrator, DSRunConfig

__all__ = ["DataScienceOrchestrator", "DSRunConfig"]
