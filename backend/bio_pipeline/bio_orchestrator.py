"""
BioOrchestrator — wraps any bioinformatics tool call with run tracking,
evaluation against a prior run, and narrative review.

Usage
-----
from bio_pipeline import BioOrchestrator

orchestrator = BioOrchestrator(
    tool_executor=dispatch_tool,          # async (tool_name, params) -> dict
    history_manager=history_mgr,          # HistoryManager instance (or None)
)

result = await orchestrator.run(
    tool_name="bulk_rnaseq_analysis",
    params={"count_matrix": "...", ...},
    parent_run_id="abc123",
    session_id="sess-xyz",
    objective="T. gondii infection vs control",
)
"""
from __future__ import annotations

import json
import logging
import os
import time
import uuid
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from .bio_evaluator import BioEvaluator
from .bio_reviewer import BioReviewer
from .run_config import BioRunConfig

logger = logging.getLogger(__name__)

# Artifacts are written to <repo-root>/artifacts/<run_id>/run.json
_ARTIFACTS_ROOT = Path(__file__).parent.parent.parent / "artifacts"


class BioOrchestrator:
    """
    Execute a bioinformatics tool, compare against a prior run,
    generate a narrative summary, and persist a run.json artifact.

    Parameters
    ----------
    tool_executor:
        Async callable ``(tool_name: str, params: dict) -> dict``.
        Typically ``dispatch_tool`` from main_with_mcp.
    history_manager:
        Optional HistoryManager instance.  When supplied, the run is
        registered via ``add_history_entry`` so it appears in the session ledger.
    """

    def __init__(
        self,
        tool_executor: Callable,
        history_manager: Optional[Any] = None,
    ) -> None:
        self._exec = tool_executor
        self._hm = history_manager
        self._evaluator = BioEvaluator()
        self._reviewer = BioReviewer()

    # ------------------------------------------------------------------ #
    # Public API                                                           #
    # ------------------------------------------------------------------ #

    async def run(
        self,
        tool_name: str,
        params: Dict[str, Any],
        *,
        parent_run_id: Optional[str] = None,
        session_id: Optional[str] = None,
        objective: str = "",
    ) -> Dict[str, Any]:
        """
        Execute *tool_name* with *params*, enrich with run metadata,
        compare to the parent run (if any), and return a unified result dict.

        Returns
        -------
        dict with keys:
          run_id, parent_run_id, tool_name, params, result,
          summary_text, delta, status, elapsed_s
        """
        run_id = str(uuid.uuid4())
        config = BioRunConfig(
            tool_name=tool_name,
            params=params,
            run_id=run_id,
            parent_run_id=parent_run_id,
            session_id=session_id,
            objective=objective,
        )

        t0 = time.time()
        try:
            result = await self._exec(tool_name, params)
        except Exception as exc:
            logger.exception("Tool %s failed", tool_name)
            result = {"status": "error", "message": str(exc)}

        elapsed = time.time() - t0

        # Load prior run for comparison
        prior_result: Optional[Dict] = None
        if parent_run_id:
            prior_result = self._load_run(parent_run_id)

        # Compute delta
        delta = self._evaluator.compare(tool_name, result, prior_result)

        # Generate narrative
        summary_text = self._reviewer.summarize(
            tool_name, result, delta=delta, objective=objective
        )

        # Persist artifact + analysis.py script
        self._persist(config, result, delta, summary_text, elapsed)

        # Determine script path so it can be registered in history
        _script_path = self._script_path_for(config)

        # Register with history manager
        if self._hm is not None and session_id:
            try:
                _artifacts_meta: List[Dict] = []
                if _script_path and _script_path.exists():
                    _artifacts_meta.append({
                        "type": "script",
                        "uri":  str(_script_path),
                        "title": "analysis.py",
                    })
                self._hm.add_history_entry(
                    session_id=session_id,
                    command=f"[{tool_name}] {objective or tool_name}",
                    tool=tool_name,
                    result=result,
                    metadata={
                        "run_id": run_id,
                        "parent_run_id": parent_run_id,
                        "summary": summary_text[:500],
                        "elapsed_s": round(elapsed, 2),
                        # Store tool params so iterative plot-update tools can
                        # re-run the same analysis with modified parameters.
                        "tool_args": params,
                        "produced_artifacts": _artifacts_meta,
                    },
                )
            except Exception as e:
                logger.warning("history_manager.add_history_entry failed: %s", e)

        return {
            "run_id": run_id,
            "parent_run_id": parent_run_id,
            "tool_name": tool_name,
            "params": params,
            "result": result,
            "summary_text": summary_text,
            "delta": delta,
            "status": result.get("status", "success"),
            "elapsed_s": round(elapsed, 2),
            # Flatten top-level convenience keys from the underlying result
            "text": result.get("text", summary_text),
            "visuals": result.get("visuals", []),
            "links": result.get("links", []),
            "visualization_type": result.get("visualization_type", "results_viewer"),
        }

    # ------------------------------------------------------------------ #
    # Artifact persistence                                                 #
    # ------------------------------------------------------------------ #

    def _persist(
        self,
        config: BioRunConfig,
        result: Dict,
        delta: Dict,
        summary_text: str,
        elapsed: float,
    ) -> None:
        try:
            run_dir = _ARTIFACTS_ROOT / config.run_id
            run_dir.mkdir(parents=True, exist_ok=True)
            run_json = {
                "config": config.to_dict(),
                "result_status": result.get("status", "unknown"),
                "summary_text": summary_text,
                "delta": delta,
                "elapsed_s": round(elapsed, 2),
                "n_visuals": len(result.get("visuals", [])),
            }
            (run_dir / "run.json").write_text(json.dumps(run_json, indent=2))
        except Exception as e:
            logger.debug("Failed to persist artifact for run %s: %s", config.run_id, e)

        # ── Save a re-runnable analysis.py so iterative tools can patch and re-execute ──
        self._save_script(config, result)

    def _script_path_for(self, config: BioRunConfig) -> Optional[Path]:
        """Return the expected analysis.py path for *config* (may not exist yet)."""
        from pathlib import Path as _Path
        sid = config.session_id or "default"
        sessions_root = _Path(__file__).parent.parent.parent / "sessions"
        return sessions_root / sid / "runs" / config.run_id / "analysis.py"

    def _save_script(self, config: BioRunConfig, result: Dict) -> None:
        """Generate and persist analysis.py under sessions/<sid>/runs/<run_id>/."""
        try:
            from backend import code_generator as _cg

            sid = config.session_id or "default"
            script_path = self._script_path_for(config)
            script_path.parent.mkdir(parents=True, exist_ok=True)

            script_text = _cg.generate(
                config.tool_name,
                config.params,
                run_id=config.run_id,
                session_id=sid,
                run_dir=str(script_path.parent),
            )
            if script_text is None:
                return

            script_path.write_text(script_text)
            logger.debug("Saved analysis script: %s", script_path)
        except Exception as exc:
            logger.debug("Script generation failed for run %s: %s", config.run_id, exc)

    def _load_run(self, run_id: str) -> Optional[Dict]:
        """Load the persisted result for *run_id*, or None if not found."""
        path = _ARTIFACTS_ROOT / run_id / "run.json"
        try:
            data = json.loads(path.read_text())
            return data
        except Exception:
            return None

    # ------------------------------------------------------------------ #
    # Convenience: re-run with parameter overrides                        #
    # ------------------------------------------------------------------ #

    async def rerun(
        self,
        prior_run_id: str,
        changes: Dict[str, Any],
        session_id: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Load a prior run's config, apply *changes* over its params, and re-run.

        Returns the same dict as :meth:`run`.
        """
        prior = self._load_run(prior_run_id)
        if prior is None:
            return {
                "status": "error",
                "message": f"Run {prior_run_id} not found in artifacts.",
            }
        config = BioRunConfig.from_dict(prior["config"])
        new_params = {**config.params, **changes}
        return await self.run(
            tool_name=config.tool_name,
            params=new_params,
            parent_run_id=prior_run_id,
            session_id=session_id or config.session_id,
            objective=config.objective,
        )

    async def diff_runs(self, run_id_a: str, run_id_b: str) -> Dict[str, Any]:
        """Return a structured diff of two runs."""
        run_a = self._load_run(run_id_a)
        run_b = self._load_run(run_id_b)
        if run_a is None or run_b is None:
            missing = run_id_a if run_a is None else run_id_b
            return {"status": "error", "message": f"Run {missing} not found."}

        config_a = run_a.get("config", {})
        config_b = run_b.get("config", {})

        param_diff: Dict[str, Any] = {}
        all_keys = set(config_a.get("params", {})) | set(config_b.get("params", {}))
        for k in all_keys:
            va = config_a.get("params", {}).get(k)
            vb = config_b.get("params", {}).get(k)
            if va != vb:
                param_diff[k] = {"run_a": va, "run_b": vb}

        narrative_a = run_a.get("summary_text", "")
        narrative_b = run_b.get("summary_text", "")

        return {
            "status": "success",
            "run_id_a": run_id_a,
            "run_id_b": run_id_b,
            "tool_name": config_a.get("tool_name"),
            "param_changes": param_diff,
            "delta_a_to_b": run_b.get("delta", {}),
            "narrative_a": narrative_a[:300],
            "narrative_b": narrative_b[:300],
        }
