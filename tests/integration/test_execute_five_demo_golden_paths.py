from __future__ import annotations

from typing import Any, Dict

import pytest
from fastapi.testclient import TestClient


@pytest.fixture(autouse=True)
def _isolate(tmp_path, monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")
    monkeypatch.setenv("HELIX_DEMO_MODE", "0")
    monkeypatch.setenv("HELIX_SANDBOX_HOST_FALLBACK", "1")

    from backend.history_manager import history_manager

    history_manager.storage_dir = tmp_path / "sessions"
    history_manager.storage_dir.mkdir(parents=True, exist_ok=True)
    history_manager.sessions = {}
    history_manager._sessions_loaded = True

    import backend.main as _mwm

    _mwm._daily_prompt_counters.clear()


@pytest.fixture()
def client() -> TestClient:
    from backend.main import app

    return TestClient(app)


@pytest.fixture()
def fake_execute_backend(monkeypatch):
    """
    Route-level integration harness:
    keep /execute orchestration real while making tool execution deterministic.
    """
    state: Dict[str, Dict[str, int]] = {}

    def _sess(params: Dict[str, Any]) -> str:
        return str(params.get("session_id") or "default")

    async def _fake_dispatch(tool: str, params: dict) -> dict:
        sid = _sess(params)
        s = state.setdefault(sid, {"analyses": 0, "reruns": 0, "diffs": 0, "patches": 0})

        if tool in {"bulk_rnaseq_analysis", "single_cell_analysis", "phylogenetic_tree", "fastqc_quality_analysis"}:
            if tool == "single_cell_analysis" and not params.get("data_file"):
                return {
                    "status": "needs_inputs",
                    "text": "single-cell analysis needs `data_file`.",
                    "needs_inputs": True,
                }
            s["analyses"] += 1
            run_id = f"{sid}-run-{s['analyses']}"
            return {
                "status": "success",
                "text": f"{tool} complete.",
                "run_id": run_id,
                "visualization_type": "results_viewer",
                "links": [{"label": "analysis.py", "url": f"/download/script?path=/tmp/{run_id}.py"}],
            }

        if tool == "bio_rerun":
            if s["analyses"] < 1:
                return {
                    "status": "needs_inputs",
                    "text": "No base run yet for rerun.",
                    "needs_inputs": True,
                }
            s["reruns"] += 1
            return {
                "status": "success",
                "text": "Re-run complete with requested changes.",
                "run_id": f"{sid}-rerun-{s['reruns']}",
                "visualization_type": "results_viewer",
            }

        if tool == "patch_and_rerun":
            if s["analyses"] < 1:
                return {
                    "status": "needs_inputs",
                    "text": "No patchable script in session history.",
                    "needs_inputs": True,
                }
            s["patches"] += 1
            return {
                "status": "success",
                "text": "## Analysis Updated\n\nPatched and re-executed successfully.",
                "run_id": f"{sid}-patch-{s['patches']}",
                "visualization_type": "results_viewer",
                "links": [{"label": "bundle.zip", "url": f"/download/bundle?session_id={sid}&run_id=patch"}],
            }

        if tool == "bio_diff_runs":
            if s["analyses"] + s["reruns"] + s["patches"] < 2:
                return {
                    "status": "needs_inputs",
                    "text": "Need at least two states to compare.",
                    "needs_inputs": True,
                }
            s["diffs"] += 1
            return {
                "status": "success",
                "text": "## Run Comparison\n\nComparison computed from historical states.",
                "visualization_type": "text",
                "result": {"source": "session_ledger_fallback", "diff_id": s["diffs"]},
            }

        if tool == "go_enrichment_analysis":
            genes = params.get("gene_list") if isinstance(params.get("gene_list"), list) else []
            if not genes:
                return {
                    "status": "needs_inputs",
                    "text": "GO enrichment recognized; provide a gene_list to execute.",
                    "needs_inputs": True,
                }
            return {
                "status": "success",
                "text": "GO enrichment complete.",
                "result": {"gene_count": len(genes)},
            }

        return {"status": "success", "text": f"mock {tool} complete"}

    monkeypatch.setattr("backend.main.dispatch_tool", _fake_dispatch)


def _post_execute(client: TestClient, command: str, session_id: str | None = None, execute_plan: bool = False) -> dict:
    payload: Dict[str, Any] = {"command": command, "execute_plan": execute_plan}
    if session_id:
        payload["session_id"] = session_id
    r = client.post("/execute", json=payload)
    assert r.status_code == 200, r.text
    out = r.json()
    assert out.get("session_id")
    assert isinstance(out.get("text", ""), str)
    return out


def test_execute_demo1_bulk_iterative_approval_and_rerun(client: TestClient, fake_execute_backend):
    first = _post_execute(client, "Analyze this RNA-seq dataset and tell me what is going on.")
    sid = first["session_id"]
    assert first["tool"] == "__plan__"
    assert first["status"] in {"workflow_planned", "needs_inputs"}
    if first["status"] == "needs_inputs":
        assert first.get("execute_ready") is False
        assert "all inputs are available" not in first["text"].lower()

    approved = _post_execute(client, "Approve.", sid)
    # `pipeline_executed` is the synchronous success status returned when
    # the planned pipeline has no remaining executable steps (see
    # CommandProcessor.execute_pipeline() in backend/agent.py).
    assert approved["status"] in {"success", "pipeline_submitted", "pipeline_executed"}
    assert "executed" in approved["text"].lower() or "submitted" in approved["text"].lower()

    rerun = _post_execute(client, "The PCA suggests 2 outlier samples. Exclude them and rerun.", sid)
    assert rerun["tool"] in {"bio_rerun", "__plan__"}
    assert rerun["text"].strip()


def test_execute_demo2_scrna_inputs_then_execution_and_update(client: TestClient, fake_execute_backend):
    first = _post_execute(client, "Run scRNA-seq analysis for PBMC disease/control.")
    sid = first["session_id"]
    assert first["tool"] == "single_cell_analysis"
    assert first["status"] == "needs_inputs"

    followup = _post_execute(
        client,
        "data_file: s3://demo/scrna/pbmc.h5 data_format: 10x resolution: 0.5 steps: all",
        sid,
    )
    assert followup["status"] in {"success", "workflow_planned"}

    update = _post_execute(client, "Now highlight IL6, TNF, CXCL8, and NFKB1 on the volcano plot.", sid)
    assert update["tool"] in {"patch_and_rerun", "bio_rerun", "__plan__"}
    assert update["text"].strip()


def test_execute_demo3_amplicon_fastqc_and_followup_outputs(client: TestClient, fake_execute_backend):
    first = _post_execute(
        client,
        "Run FastQC and amplicon QC on s3://bucket/test_R1.fq and s3://bucket/test_R2.fq",
    )
    sid = first["session_id"]
    assert first["tool"] == "fastqc_quality_analysis"
    assert first["status"] == "success"

    followup = _post_execute(client, "visualize the results of the latest job", sid)
    assert followup["text"].strip()


def test_execute_demo4_timecourse_bulk_and_historical_compare(client: TestClient, fake_execute_backend):
    first = _post_execute(client, "Run bulk RNA-seq time-course analysis for APAP recovery.")
    sid = first["session_id"]
    assert first["tool"] in {"bulk_rnaseq_analysis", "__plan__"}

    _post_execute(client, "Approve.", sid)
    _post_execute(client, "Run differential expression adjusting for sex and batch.", sid)
    _post_execute(client, "Approve.", sid)

    compare = _post_execute(
        client,
        "Show differences between my runs and summarize what changed.",
        sid,
    )
    assert compare["tool"] == "bio_diff_runs"
    assert compare["status"] in {"success", "needs_inputs"}
    assert compare["text"].strip()


def test_execute_demo4_workflow_design_intent_does_not_flatten_to_bulk_tool(client: TestClient, fake_execute_backend):
    first = _post_execute(
        client,
        (
            "Design the workflow before execution for this APAP bulk RNA-seq time-course study. "
            "Propose expected calculations, output artifacts, recommended plots, and QC checkpoints."
        ),
    )
    assert first["tool"] != "bulk_rnaseq_analysis"
    assert first["status"] in {"workflow_planned", "success", "needs_inputs"}
    assert first["text"].strip()
    assert "simple_operation" not in first["text"]
    assert "fastqc" not in first["text"].lower()


def test_execute_demo5_historical_recreation_executes_when_resolved(client: TestClient, fake_execute_backend):
    first = _post_execute(client, "Analyze this RNA-seq dataset and tell me what is going on.")
    sid = first["session_id"]
    _post_execute(client, "Approve.", sid)
    _post_execute(client, "Actually sample S08 is mislabeled. It should be control, not treated.", sid)
    _post_execute(client, "Approve the correction and rerun the comparison.", sid)

    recreate = _post_execute(
        client,
        "Recreate the figure set corresponding to the corrected metadata version before the fold-change bug fix.",
        sid,
    )
    assert recreate["tool"] == "bio_diff_runs"
    assert recreate["status"] in {"success", "needs_inputs"}
    assert recreate["status"] != "workflow_planned"
    assert recreate["text"].strip()
