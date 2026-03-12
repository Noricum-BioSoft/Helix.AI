"""
Unit tests for backend/bio_pipeline/.

Groups:
  A. BioRunConfig — serialisation round-trip
  B. BioEvaluator — per-tool delta metrics
  C. BioReviewer — narrative generation
  D. BioOrchestrator — run, rerun, diff_runs
  E. BioOrchestrator — parent_run_id linkage (iterative workflow)
"""
from __future__ import annotations

import asyncio
import json
import sys
import os
from pathlib import Path
from typing import Any, Dict

import pytest

# Ensure backend/ is on sys.path so `from bio_pipeline import ...` works.
_BACKEND = Path(__file__).parent.parent.parent.parent / "backend"
if str(_BACKEND) not in sys.path:
    sys.path.insert(0, str(_BACKEND))

from bio_pipeline import BioOrchestrator, BioEvaluator, BioReviewer, BioRunConfig


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

def _async(coro):
    return asyncio.get_event_loop().run_until_complete(coro)


def _make_rnaseq_result(n_sig: int = 10, top_gene: str = "Mx1") -> Dict[str, Any]:
    return {
        "status": "success",
        "mode": "real",
        "n_samples": 12,
        "n_genes_total": 1000,
        "summary": [
            {
                "contrast": "infection_status: infected vs uninfected",
                "total_genes": 1000,
                "significant": n_sig,
                "upregulated": n_sig // 2,
                "downregulated": n_sig - n_sig // 2,
            }
        ],
        "top_genes": f"\n**infection_status: infected vs uninfected** — top genes: {top_gene}, Stat1",
        "plots": {},
    }


def _make_scrna_result(n_clusters: int = 8) -> Dict[str, Any]:
    return {
        "status": "success",
        "mode": "synthetic",
        "n_cells": 600,
        "n_genes": 300,
        "n_clusters": n_clusters,
        "resolution": 0.5,
        "text": "scRNA analysis done",
        "plots": {},
        "markers": {"0": ["CD3D", "CD4"], "1": ["CD19", "MS4A1"]},
        "cell_type_composition": {"CD4 T cell": 150, "B cell": 108},
    }


def _make_phylo_result() -> Dict[str, Any]:
    return {
        "status": "success",
        "tree_newick": "(A:0.1,B:0.2);",
        "statistics": {
            "n_sequences": 5,
            "mean_pairwise_distance": 0.05,
            "mean_pairwise_identity_pct": 95.0,
        },
        "text": "Phylo analysis done",
        "plots": {},
    }


async def _fake_exec_success(tool_name: str, params: Dict) -> Dict:
    return {"status": "success", "tool_called": tool_name, "params": params}


async def _fake_exec_rnaseq(tool_name: str, params: Dict) -> Dict:
    return _make_rnaseq_result()


async def _fake_exec_error(tool_name: str, params: Dict) -> Dict:
    raise RuntimeError("Simulated tool failure")


# ---------------------------------------------------------------------------
# A. BioRunConfig
# ---------------------------------------------------------------------------

class TestBioRunConfig:
    def test_defaults(self):
        cfg = BioRunConfig(tool_name="bulk_rnaseq_analysis", params={"alpha": 0.05})
        assert cfg.tool_name == "bulk_rnaseq_analysis"
        assert cfg.run_id  # UUID generated
        assert cfg.parent_run_id is None

    def test_round_trip(self):
        cfg = BioRunConfig(
            tool_name="phylogenetic_tree",
            params={"aligned_sequences": ""},
            parent_run_id="parent-123",
            session_id="sess-xyz",
            objective="test run",
        )
        d = cfg.to_dict()
        restored = BioRunConfig.from_dict(d)
        assert restored.tool_name == cfg.tool_name
        assert restored.run_id == cfg.run_id
        assert restored.parent_run_id == cfg.parent_run_id
        assert restored.session_id == cfg.session_id
        assert restored.objective == cfg.objective


# ---------------------------------------------------------------------------
# B. BioEvaluator
# ---------------------------------------------------------------------------

class TestBioEvaluator:
    def setup_method(self):
        self.ev = BioEvaluator()

    def test_no_prior_returns_empty(self):
        delta = self.ev.compare("bulk_rnaseq_analysis", _make_rnaseq_result(), prior=None)
        assert delta == {}

    def test_rnaseq_delta_significant_genes(self):
        prior = _make_rnaseq_result(n_sig=10)
        current = _make_rnaseq_result(n_sig=15)
        delta = self.ev.compare("bulk_rnaseq_analysis", current, prior)
        assert delta["significant_genes_delta"] == 5
        assert delta["significant_genes_current"] == 15
        assert delta["significant_genes_prior"] == 10
        assert "narrative" in delta

    def test_rnaseq_top_gene_overlap(self):
        prior = {"result": {"top_genes": "\n**test** — top genes: Mx1, Stat1, Ifit1"}}
        current = {"result": {"top_genes": "\n**test** — top genes: Mx1, Stat1, NewGene"}}
        delta = self.ev.compare("bulk_rnaseq_analysis", current, prior)
        # mx1 and stat1 overlap
        assert delta["top_gene_overlap"] >= 2

    def test_scrna_delta_clusters(self):
        prior = _make_scrna_result(n_clusters=8)
        current = _make_scrna_result(n_clusters=12)
        delta = self.ev.compare("single_cell_analysis", current, prior)
        assert delta["cluster_count_delta"] == 4
        assert delta["cluster_count_current"] == 12

    def test_phylo_delta_distance(self):
        prior = _make_phylo_result()
        current = _make_phylo_result()
        current["statistics"]["mean_pairwise_distance"] = 0.08
        delta = self.ev.compare("phylogenetic_tree", current, prior)
        assert abs(delta["mean_distance_delta"] - 0.03) < 1e-9

    def test_unknown_tool_uses_generic(self):
        delta = self.ev.compare("unknown_tool", {"status": "success"}, {"status": "success"})
        assert "narrative" in delta


# ---------------------------------------------------------------------------
# C. BioReviewer
# ---------------------------------------------------------------------------

class TestBioReviewer:
    def setup_method(self):
        self.rv = BioReviewer()

    def test_rnaseq_review_contains_key_fields(self):
        text = self.rv.summarize("bulk_rnaseq_analysis", _make_rnaseq_result())
        assert "12" in text or "n_samples" in text.lower() or "bulk" in text.lower()
        assert "1000" in text or "1,000" in text

    def test_rnaseq_review_includes_delta_narrative(self):
        delta = {"narrative": "Genes went from 10 to 15."}
        text = self.rv.summarize("bulk_rnaseq_analysis", _make_rnaseq_result(), delta=delta)
        assert "Genes went from 10 to 15." in text
        assert "Prior Run" in text or "Comparison" in text

    def test_scrna_review(self):
        text = self.rv.summarize("single_cell_analysis", _make_scrna_result())
        assert "600" in text or "cells" in text.lower()

    def test_phylo_review(self):
        text = self.rv.summarize("phylogenetic_tree", _make_phylo_result())
        assert "5" in text or "sequence" in text.lower()

    def test_generic_fallback(self):
        text = self.rv.summarize("unknown_tool_xyz", {"status": "success"})
        assert "success" in text.lower()


# ---------------------------------------------------------------------------
# D. BioOrchestrator — basic run
# ---------------------------------------------------------------------------

class TestBioOrchestrator:
    def test_run_returns_run_id(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "bio_pipeline.bio_orchestrator._ARTIFACTS_ROOT", tmp_path
        )
        orch = BioOrchestrator(tool_executor=_fake_exec_rnaseq)
        result = _async(orch.run("bulk_rnaseq_analysis", {"alpha": 0.05}))
        assert result["run_id"]
        assert result["status"] == "success"

    def test_run_persists_artifact(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "bio_pipeline.bio_orchestrator._ARTIFACTS_ROOT", tmp_path
        )
        orch = BioOrchestrator(tool_executor=_fake_exec_rnaseq)
        result = _async(orch.run("bulk_rnaseq_analysis", {"alpha": 0.05}))
        run_dir = tmp_path / result["run_id"]
        assert (run_dir / "run.json").exists()
        data = json.loads((run_dir / "run.json").read_text())
        assert data["config"]["tool_name"] == "bulk_rnaseq_analysis"

    def test_run_includes_summary_text(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "bio_pipeline.bio_orchestrator._ARTIFACTS_ROOT", tmp_path
        )
        orch = BioOrchestrator(tool_executor=_fake_exec_rnaseq)
        result = _async(orch.run("bulk_rnaseq_analysis", {}))
        assert len(result["summary_text"]) > 20

    def test_run_tool_error_returns_error_status(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "bio_pipeline.bio_orchestrator._ARTIFACTS_ROOT", tmp_path
        )
        orch = BioOrchestrator(tool_executor=_fake_exec_error)
        result = _async(orch.run("bulk_rnaseq_analysis", {}))
        assert result["status"] == "error"
        assert result["result"]["status"] == "error"

    def test_run_propagates_parent_run_id(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "bio_pipeline.bio_orchestrator._ARTIFACTS_ROOT", tmp_path
        )
        orch = BioOrchestrator(tool_executor=_fake_exec_rnaseq)
        result = _async(orch.run("bulk_rnaseq_analysis", {}, parent_run_id="parent-abc"))
        assert result["parent_run_id"] == "parent-abc"
        data = json.loads((tmp_path / result["run_id"] / "run.json").read_text())
        assert data["config"]["parent_run_id"] == "parent-abc"


# ---------------------------------------------------------------------------
# E. BioOrchestrator — iterative workflows
# ---------------------------------------------------------------------------

class TestBioOrchestratorIterative:
    def test_rerun_with_changed_params(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "bio_pipeline.bio_orchestrator._ARTIFACTS_ROOT", tmp_path
        )
        orch = BioOrchestrator(tool_executor=_fake_exec_rnaseq)

        # First run
        run1 = _async(orch.run(
            "bulk_rnaseq_analysis",
            {"alpha": 0.05, "design_formula": "~condition"},
        ))
        run1_id = run1["run_id"]

        # Re-run with alpha=0.01
        run2 = _async(orch.rerun(run1_id, changes={"alpha": 0.01}))
        assert run2["parent_run_id"] == run1_id
        assert run2["params"]["alpha"] == 0.01
        assert run2["params"]["design_formula"] == "~condition"  # unchanged

    def test_rerun_missing_run_returns_error(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "bio_pipeline.bio_orchestrator._ARTIFACTS_ROOT", tmp_path
        )
        orch = BioOrchestrator(tool_executor=_fake_exec_rnaseq)
        result = _async(orch.rerun("non-existent-run-id", changes={}))
        assert result["status"] == "error"

    def test_diff_runs(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "bio_pipeline.bio_orchestrator._ARTIFACTS_ROOT", tmp_path
        )
        orch = BioOrchestrator(tool_executor=_fake_exec_rnaseq)
        run1 = _async(orch.run("bulk_rnaseq_analysis", {"alpha": 0.05}))
        run2 = _async(orch.run("bulk_rnaseq_analysis", {"alpha": 0.01}))
        diff = _async(orch.diff_runs(run1["run_id"], run2["run_id"]))
        assert diff["status"] == "success"
        assert diff["run_id_a"] == run1["run_id"]
        assert diff["run_id_b"] == run2["run_id"]
        assert "alpha" in diff["param_changes"]
        assert diff["param_changes"]["alpha"]["run_a"] == 0.05
        assert diff["param_changes"]["alpha"]["run_b"] == 0.01

    def test_diff_runs_missing_run_returns_error(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "bio_pipeline.bio_orchestrator._ARTIFACTS_ROOT", tmp_path
        )
        orch = BioOrchestrator(tool_executor=_fake_exec_rnaseq)
        run1 = _async(orch.run("bulk_rnaseq_analysis", {}))
        diff = _async(orch.diff_runs(run1["run_id"], "does-not-exist"))
        assert diff["status"] == "error"

    def test_three_generation_chain(self, tmp_path, monkeypatch):
        """Run → rerun → rerun builds a proper parent_run_id chain."""
        monkeypatch.setattr(
            "bio_pipeline.bio_orchestrator._ARTIFACTS_ROOT", tmp_path
        )
        orch = BioOrchestrator(tool_executor=_fake_exec_rnaseq)

        run1 = _async(orch.run("bulk_rnaseq_analysis", {"alpha": 0.05}))
        run2 = _async(orch.rerun(run1["run_id"], changes={"alpha": 0.01}))
        run3 = _async(orch.rerun(run2["run_id"], changes={"alpha": 0.001}))

        assert run3["parent_run_id"] == run2["run_id"]
        assert run2["parent_run_id"] == run1["run_id"]
        assert run1["parent_run_id"] is None

    def test_evaluator_invoked_on_rerun(self, tmp_path, monkeypatch):
        """When parent exists, result should have delta metrics."""
        monkeypatch.setattr(
            "bio_pipeline.bio_orchestrator._ARTIFACTS_ROOT", tmp_path
        )
        orch = BioOrchestrator(tool_executor=_fake_exec_rnaseq)
        run1 = _async(orch.run("bulk_rnaseq_analysis", {"alpha": 0.05}))
        run2 = _async(orch.rerun(run1["run_id"], changes={"alpha": 0.01}))
        assert run2["delta"] != {}
