"""
test_execute_dispatch.py — Unit tests for execute() dispatch logic.

Covers every dispatch branch in the execute() endpoint:

  Group A — Prologue (session auto-creation, ID reuse, rate limiting)
  Group B — S3 browse fast-path
  Group C — Demo fast-paths via DemoDispatcher (all 5 demos + history regression)
  Group D — CommandRouter fast-path (Phase 2c)
  Group E — BioAgent path (tool_mapped, full execution, timeout)
  Group F — Fallback path (Q&A intent, workflow plan IR, normal broker)
  Group G — Session-context side effects (mutation, alignment)
  Group H — Error handling and pure-unit helpers

Follows the pattern from test_iterative_workflows_local.py:
- ``monkeypatch.setenv("HELIX_MOCK_MODE", "1")`` — skip LLM entirely
- isolated ``history_manager.storage_dir`` + ``history_manager.sessions = {}``
- ``TestClient(app)`` for HTTP round-trips
- ``unittest.mock`` for heavy dependencies (agent, dispatch_tool, broker)
"""
from __future__ import annotations

import asyncio
import sys
from typing import Any, Dict, List
from unittest.mock import AsyncMock, MagicMock

import pytest
from fastapi.testclient import TestClient


# ── Shared fixtures ───────────────────────────────────────────────────────────


@pytest.fixture(autouse=True)
def _isolate(tmp_path, monkeypatch):
    """
    Baseline isolation for every test:
    - mock mode (no LLM)
    - isolated session storage
    - reset rate-limit counters
    """
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")
    monkeypatch.setenv("HELIX_DEMO_MODE", "0")
    monkeypatch.setenv("HELIX_SANDBOX_HOST_FALLBACK", "1")

    from backend.history_manager import history_manager

    history_manager.storage_dir = tmp_path / "sessions"
    history_manager.storage_dir.mkdir(parents=True, exist_ok=True)
    history_manager.sessions = {}
    history_manager._sessions_loaded = True

    import backend.main_with_mcp as _mwm

    _mwm._daily_prompt_counters.clear()


def _client() -> TestClient:
    from backend.main_with_mcp import app

    return TestClient(app)


def _get_runs(client: TestClient, session_id: str) -> List[Dict[str, Any]]:
    r = client.get(f"/session/{session_id}/runs")
    assert r.status_code == 200, r.text
    return r.json()["runs"]


# ── Group A: Prologue ─────────────────────────────────────────────────────────


class TestPrologue:
    def test_auto_creates_session_when_none_provided(self):
        """A request without session_id receives a fresh UUID in the response."""
        client = _client()
        r = client.post("/execute", json={"command": "Create demo plot"})
        assert r.status_code == 200
        session_id = r.json().get("session_id")
        assert session_id, "Expected a session_id in the response"
        assert len(session_id) >= 8

    def test_reuses_provided_session_id(self):
        """A request with a session_id always echoes the same ID back."""
        client = _client()
        r1 = client.post("/execute", json={"command": "Create demo plot"})
        session_id = r1.json()["session_id"]

        r2 = client.post(
            "/execute",
            json={"command": "Create demo plot", "session_id": session_id},
        )
        assert r2.json()["session_id"] == session_id

    def test_rate_limit_blocks_after_limit(self, monkeypatch):
        """
        Requests over the daily limit are rejected.

        The HTTPException(429) is caught by execute()'s outer try/except and
        converted to a structured error body (HTTP 200, success=False).
        """
        import backend.main_with_mcp as _mwm
        from backend.main_with_mcp import MAX_PROMPTS_PER_DAY, _today_iso

        session_id = "rate-limit-test-session"
        day = _today_iso()
        _mwm._daily_prompt_counters[f"session:{session_id}"] = {
            day: MAX_PROMPTS_PER_DAY
        }

        client = _client()
        r = client.post(
            "/execute",
            json={"command": "Create demo plot", "session_id": session_id},
        )
        payload = r.json()
        # The HTTPException(429) is caught by execute()'s broad outer except clause,
        # so it surfaces as HTTP 200 with success=False (the error field is empty
        # because Starlette's HTTPException does not populate args).
        assert r.status_code in (200, 429)
        if r.status_code == 200:
            assert payload.get("success") is False
            assert payload.get("session_id") == session_id
        # If FastAPI exception middleware intercepts first we get proper 429 — also fine.


# ── Group B: S3 browse fast-path ──────────────────────────────────────────────


class TestS3Browse:
    def test_s3_browse_routes_to_s3_browse_tool(self, monkeypatch):
        """
        'show s3://bucket/.../results.json' → s3_browse_results is invoked.
        History is NOT recorded (S3 browse skips the run ledger).
        """
        mock_result = {"status": "success", "text": "Files listed", "files": []}
        mock_tool = MagicMock()
        mock_tool.invoke = MagicMock(return_value=mock_result)
        monkeypatch.setattr("backend.agent_tools.s3_browse_results", mock_tool)

        client = _client()
        r = client.post(
            "/execute",
            json={
                "command": (
                    "show s3://bucket/prefix/ "
                    "s3://bucket/prefix/results.json"
                )
            },
        )
        assert r.status_code == 200
        mock_tool.invoke.assert_called_once()

        # S3 browse does NOT record a run ledger entry
        session_id = r.json().get("session_id")
        if session_id:
            runs = _get_runs(client, session_id)
            assert not any(run.get("tool") == "s3_browse_results" for run in runs)

    def test_s3_browse_guard_does_not_intercept_pipeline_commands(self, monkeypatch):
        """
        'run fastqc on s3://...' should NOT trigger the S3 browse fast-path
        because "run" is not a browse verb and the exec guard fires.
        """
        mock_tool = MagicMock()
        mock_tool.invoke = MagicMock(return_value={"status": "success"})
        monkeypatch.setattr("backend.agent_tools.s3_browse_results", mock_tool)

        # Phase 2c will handle the fastqc command; mock dispatch_tool for it
        async def _mock_call(tool, params):
            return {"status": "success", "text": f"mock {tool}"}

        monkeypatch.setattr("backend.main_with_mcp.dispatch_tool", _mock_call)

        client = _client()
        client.post(
            "/execute",
            json={"command": "run fastqc on s3://bucket/sample.fastq.gz"},
        )
        mock_tool.invoke.assert_not_called()


# ── Group C: Demo fast-paths ──────────────────────────────────────────────────

TGONDII_CMD = (
    "Analyze s3://noricum-ngs-data/demo/rnaseq/tgondii_counts.csv "
    "with sample metadata s3://noricum-ngs-data/demo/rnaseq/tgondii_metadata.csv "
    "design_formula: ~infection_status"
)
APAP_CMD = (
    "Analyze s3://noricum-ngs-data/demo/rnaseq/apap_timecourse_counts.csv "
    "with metadata s3://noricum-ngs-data/demo/rnaseq/apap_timecourse_metadata.csv "
    "design_formula: ~time_point"
)
SLE_CMD = (
    "Run scRNA-seq analysis on "
    "s3://noricum-ngs-data/demo/scrna/sle_pbmc_filtered_feature_bc_matrix.h5"
)
PHYLO_CMD = "Run SARS-CoV-2 spike protein variant phylogenetic analysis"
AMPLICON_CMD = (
    "Run FastQC and amplicon QC on "
    "s3://helix-test-data/sample01/forward_reads.fastq.gz"
)


@pytest.fixture()
def mock_tool_executor(monkeypatch):
    """Replace dispatch_tool so demo tests need no AWS credentials."""

    async def _mock(tool: str, params: dict) -> dict:
        if tool == "bulk_rnaseq_analysis":
            return {
                "status": "success",
                "text": "RNA-seq complete",
                "result": {"run_id": "rna-1"},
            }
        if tool == "single_cell_analysis":
            return {
                "status": "success",
                "text": "scRNA-seq complete",
                "result": {"run_id": "sc-1"},
            }
        return {"status": "success", "text": f"mock {tool}"}

    monkeypatch.setattr("backend.main_with_mcp.dispatch_tool", _mock)


class TestDemoFastPaths:
    def test_demo_tgondii_routes_correctly(self, mock_tool_executor):
        """T. gondii demo → tool == 'bulk_rnaseq_analysis'."""
        client = _client()
        r = client.post("/execute", json={"command": TGONDII_CMD})
        assert r.status_code == 200
        assert r.json().get("tool") == "bulk_rnaseq_analysis"

    def test_demo_apap_routes_correctly(self, mock_tool_executor):
        """APAP time-course demo → tool == 'bulk_rnaseq_analysis'."""
        client = _client()
        r = client.post("/execute", json={"command": APAP_CMD})
        assert r.status_code == 200
        assert r.json().get("tool") == "bulk_rnaseq_analysis"

    def test_demo_sle_routes_correctly(self, mock_tool_executor):
        """SLE PBMC scRNA-seq demo → tool == 'single_cell_analysis'."""
        client = _client()
        r = client.post("/execute", json={"command": SLE_CMD})
        assert r.status_code == 200
        assert r.json().get("tool") == "single_cell_analysis"

    def test_demo_phylo_routes_correctly(self, mock_tool_executor):
        """SARS-CoV-2 phylogenetics demo → tool == 'phylogenetic_tree'."""
        client = _client()
        r = client.post("/execute", json={"command": PHYLO_CMD})
        assert r.status_code == 200
        payload = r.json()
        assert payload.get("tool") == "phylogenetic_tree"

    def test_demo_amplicon_routes_correctly(self, mock_tool_executor):
        """Amplicon QC demo → tool == 'fastqc_quality_analysis'."""
        client = _client()
        r = client.post("/execute", json={"command": AMPLICON_CMD})
        assert r.status_code == 200
        assert r.json().get("tool") == "fastqc_quality_analysis"

    def test_all_demos_record_history_entry(self, mock_tool_executor):
        """
        Regression: every demo run must appear in the session run ledger.
        (Previously, demo fast-paths never called add_history_entry.)
        """
        client = _client()
        demo_cases = [
            (TGONDII_CMD,   "bulk_rnaseq_analysis"),
            (APAP_CMD,      "bulk_rnaseq_analysis"),
            (SLE_CMD,       "single_cell_analysis"),
            (PHYLO_CMD,     "phylogenetic_tree"),
            (AMPLICON_CMD,  "fastqc_quality_analysis"),
        ]
        for cmd, expected_tool in demo_cases:
            r = client.post("/execute", json={"command": cmd})
            assert r.status_code == 200, f"Demo '{cmd[:50]}' HTTP error: {r.text}"

            session_id = r.json().get("session_id")
            assert session_id, f"No session_id for '{cmd[:50]}'"

            runs = _get_runs(client, session_id)
            assert len(runs) == 1, (
                f"Expected exactly 1 run for demo '{cmd[:50]}', got {len(runs)}: {runs}"
            )
            assert runs[0]["tool"] == expected_tool, (
                f"Expected tool '{expected_tool}' but got '{runs[0]['tool']}' "
                f"for demo '{cmd[:50]}'"
            )


# ── Group D: CommandRouter fast-path (Phase 2c) ───────────────────────────────


class TestCommandRouterFastPath:
    def test_command_router_create_demo_plot(self):
        """'Create demo plot' → tool='local_demo_plot_script', history recorded."""
        client = _client()
        r = client.post("/execute", json={"command": "Create demo plot"})
        assert r.status_code == 200
        assert r.json().get("success") is True

        session_id = r.json()["session_id"]
        runs = _get_runs(client, session_id)
        assert len(runs) == 1
        assert runs[0]["tool"] == "local_demo_plot_script"

    def test_command_router_scale_change(self):
        """
        'change the plots from log to linear scale'
        → tool='patch_and_rerun' in run ledger.
        """
        client = _client()
        r1 = client.post("/execute", json={"command": "Create demo plot"})
        session_id = r1.json()["session_id"]

        r2 = client.post(
            "/execute",
            json={
                "command": "change the plots from log to linear scale",
                "session_id": session_id,
            },
        )
        assert r2.status_code == 200

        runs = _get_runs(client, session_id)
        edit_run = next(
            (run for run in runs if run.get("tool") == "patch_and_rerun"),
            None,
        )
        assert edit_run is not None, (
            f"Expected a patch_and_rerun run, got runs: {runs}"
        )

    def test_command_router_session_run_io_summary(self):
        """
        'What were the inputs/outputs of the first run?'
        → tool='session_run_io_summary'.
        """
        client = _client()
        r1 = client.post("/execute", json={"command": "Create demo plot"})
        session_id = r1.json()["session_id"]

        r2 = client.post(
            "/execute",
            json={
                "command": "What were the inputs/outputs of the first run?",
                "session_id": session_id,
            },
        )
        assert r2.status_code == 200
        payload = r2.json()
        assert payload.get("success") is True
        assert payload.get("tool") == "session_run_io_summary"

    def test_command_router_ds_run_analysis(self, tmp_path, monkeypatch):
        """'analyze my data file: <path>' → routes to ds_run_analysis."""
        csv_file = tmp_path / "sample.csv"
        csv_file.write_text("x,y\n1,2\n3,4\n5,6\n")

        async def _mock_call(tool: str, params: dict) -> dict:
            if tool == "ds_run_analysis":
                return {
                    "status": "success",
                    "text": "Analysis complete",
                    "result": {"run_id": "ds-1", "status": "success"},
                    "artifacts": [],
                }
            return {"status": "success", "text": f"mock {tool}"}

        monkeypatch.setattr("backend.main_with_mcp.dispatch_tool", _mock_call)

        client = _client()
        r = client.post(
            "/execute",
            json={"command": f"analyze my data file: {csv_file}"},
        )
        assert r.status_code == 200
        assert r.json().get("tool") == "ds_run_analysis"

    def test_command_router_bio_rerun_routes_correctly(self, monkeypatch):
        """
        'run again with alpha=0.01' → CommandRouter returns 'bio_rerun'
        with changes={'alpha': 0.01} and the /execute endpoint succeeds.
        """
        async def _mock_call(tool: str, params: dict) -> dict:
            if tool == "bio_rerun":
                return {
                    "status": "success",
                    "text": "Re-run complete with alpha=0.01",
                    "run_id": "new-run-id",
                    "parent_run_id": "prior-run-id",
                    "delta": {},
                }
            return {"status": "success", "text": f"mock {tool}"}

        monkeypatch.setattr("backend.main_with_mcp.dispatch_tool", _mock_call)

        client = _client()
        r = client.post(
            "/execute",
            json={"command": "run again with alpha=0.01"},
        )
        assert r.status_code == 200
        assert r.json().get("tool") == "bio_rerun"

    def test_command_router_bio_rerun_extracts_parameters(self):
        """CommandRouter.route_command parses alpha, resolution, design_formula."""
        from backend.command_router import CommandRouter

        router = CommandRouter()

        tool, params = router.route_command("re-run with alpha=0.001 and resolution=0.8", {})
        assert tool == "bio_rerun"
        assert params["changes"]["alpha"] == pytest.approx(0.001)
        assert params["changes"]["resolution"] == pytest.approx(0.8)

        tool2, params2 = router.route_command("change parameter design_formula=~batch+condition", {})
        assert tool2 == "bio_rerun"
        assert "design_formula" in params2["changes"]

    def test_command_router_bio_diff_runs_routes_correctly(self, monkeypatch):
        """
        'run compare results between runs' → CommandRouter returns 'bio_diff_runs'
        and the /execute endpoint succeeds.
        """
        async def _mock_call(tool: str, params: dict) -> dict:
            if tool == "bio_diff_runs":
                return {
                    "status": "success",
                    "text": "## Run Comparison\n\nNo changes detected.",
                    "result": {"param_changes": {}, "delta_a_to_b": {}},
                }
            return {"status": "success", "text": f"mock {tool}"}

        monkeypatch.setattr("backend.main_with_mcp.dispatch_tool", _mock_call)

        client = _client()
        r = client.post(
            "/execute",
            json={"command": "run compare results between runs"},
        )
        assert r.status_code == 200
        assert r.json().get("tool") == "bio_diff_runs"

    def test_command_router_bio_diff_runs_extracts_uuids(self):
        """CommandRouter.route_command extracts UUID run IDs when present."""
        from backend.command_router import CommandRouter

        router = CommandRouter()
        uuid_a = "a1b2c3d4-e5f6-7890-abcd-ef1234567890"
        uuid_b = "b2c3d4e5-f6a7-8901-bcde-f12345678901"

        tool, params = router.route_command(
            f"compare run {uuid_a} and {uuid_b}", {}
        )
        assert tool == "bio_diff_runs"
        assert params["run_id_a"] == uuid_a
        assert params["run_id_b"] == uuid_b

    def test_command_router_bio_diff_runs_defaults_to_latest_prior(self):
        """Without explicit UUIDs, defaults to latest/prior."""
        from backend.command_router import CommandRouter

        router = CommandRouter()
        tool, params = router.route_command("show differences between my runs", {})
        assert tool == "bio_diff_runs"
        assert params["run_id_a"] == "latest"
        assert params["run_id_b"] == "prior"


# ── Group E: BioAgent path ────────────────────────────────────────────────────


class TestAgentPath:
    """
    These tests enable the agent path (HELIX_MOCK_MODE=0) but inject a stub
    backend.agent module so no real LLM dependencies are imported.
    """

    @pytest.fixture(autouse=True)
    def _enable_agent(self, monkeypatch):
        monkeypatch.setenv("HELIX_MOCK_MODE", "0")

    @pytest.fixture()
    def stub_agent(self, monkeypatch):
        stub = MagicMock()
        monkeypatch.setitem(sys.modules, "backend.agent", stub)
        # CommandRouter returns "handle_natural_command" so Phase 2c falls through
        from backend.command_router import CommandRouter

        orig_route = CommandRouter.route_command
        monkeypatch.setattr(
            CommandRouter,
            "route_command",
            lambda self, cmd, ctx: ("handle_natural_command", {}),
        )
        return stub

    def test_agent_tool_mapped_executes_via_broker(self, monkeypatch, stub_agent):
        """
        Agent returns {status: tool_mapped, tool_name: ...}
        → broker executes that tool → response carries the tool name.
        """
        stub_agent.handle_command = AsyncMock(
            return_value={
                "status": "tool_mapped",
                "tool_name": "fastqc_quality_analysis",
                "parameters": {"input_r1": "s3://bucket/r1.fastq.gz"},
            }
        )

        mock_broker = MagicMock()
        mock_broker.execute_tool = AsyncMock(
            return_value={
                "status": "success",
                "text": "FastQC done",
                "result": {"run_id": "qc-1"},
            }
        )
        monkeypatch.setattr(
            "backend.main_with_mcp._get_execution_broker", lambda: mock_broker
        )

        client = _client()
        r = client.post("/execute", json={"command": "run custom fastqc analysis XYZ"})
        assert r.status_code == 200
        assert r.json().get("tool") == "fastqc_quality_analysis"

    def test_agent_full_execution_uses_agent_result(self, monkeypatch, stub_agent):
        """Agent returns a full result dict → response tool == 'agent'."""
        stub_agent.handle_command = AsyncMock(
            return_value={
                "status": "success",
                "text": "Pipeline complete: 3 steps executed.",
                "success": True,
            }
        )

        client = _client()
        r = client.post("/execute", json={"command": "run complex-pipeline-XYZ-99"})
        assert r.status_code == 200
        assert r.json().get("tool") == "agent"

    def test_agent_timeout_falls_through_to_fallback(self, monkeypatch, stub_agent):
        """
        Agent raises asyncio.TimeoutError → fallback path runs → HTTP 200,
        session_id present, no server error.
        """

        async def _slow(*args, **kwargs):
            raise asyncio.TimeoutError("deliberate timeout")

        stub_agent.handle_command = _slow

        mock_broker = MagicMock()
        mock_broker.execute_tool = AsyncMock(
            return_value={"status": "success", "text": "fallback done"}
        )
        monkeypatch.setattr(
            "backend.main_with_mcp._get_execution_broker", lambda: mock_broker
        )

        client = _client()
        r = client.post("/execute", json={"command": "run some-analysis-99"})
        assert r.status_code == 200
        assert r.json().get("session_id")


# ── Group F: Fallback path ────────────────────────────────────────────────────


def _disable_phase2c(monkeypatch):
    """
    Force Phase 2c to always fall through by making dispatch_tool raise.
    Phase 2c wraps the tool call in try/except, so this safely skips it.
    """

    async def _raise(tool: str, params: dict) -> dict:
        raise Exception("Phase2c intentionally disabled for fallback test")

    monkeypatch.setattr("backend.main_with_mcp.dispatch_tool", _raise)


class TestFallbackPath:
    """
    In HELIX_MOCK_MODE=1 the agent always raises, so after Phase 2c is
    disabled the fallback path runs.  Each test verifies a different branch.
    """

    @pytest.fixture(autouse=True)
    def _force_fallback(self, monkeypatch):
        _disable_phase2c(monkeypatch)

    def test_fallback_qa_intent_returns_safe_message(self, monkeypatch):
        """Q&A intent → fallback returns a human-readable 'Q&A intent' message."""
        mock_intent = MagicMock()
        mock_intent.intent = "question"
        mock_intent.reason = "user asked a question"
        monkeypatch.setattr(
            "backend.intent_classifier.classify_intent", lambda cmd: mock_intent
        )

        client = _client()
        r = client.post("/execute", json={"command": "What is RNA-seq analysis?"})
        assert r.status_code == 200
        payload = r.json()
        assert payload.get("success") is True
        text = payload.get("text") or ""
        assert any(kw in text.lower() for kw in ("q&a", "intent", "question", "agent"))

    def test_fallback_workflow_routes_to_plan_ir(self, monkeypatch):
        """
        Multi-step workflow command → broker is called with tool_name '__plan__'.
        """
        called_tools: list = []

        mock_broker = MagicMock()

        async def _capture(exec_req):
            called_tools.append(exec_req.tool_name)
            return {"status": "success", "text": "Plan executed", "result": {}}

        mock_broker.execute_tool = _capture
        monkeypatch.setattr(
            "backend.main_with_mcp._get_execution_broker", lambda: mock_broker
        )

        # Use a command with " then " that classify_intent won't block
        mock_intent = MagicMock()
        mock_intent.intent = "execute"
        monkeypatch.setattr(
            "backend.intent_classifier.classify_intent", lambda cmd: mock_intent
        )

        # Command has " then " → _looks_like_workflow True; no router keyword match
        cmd = "process-dataset-alpha-7 then validate-results-beta-7"
        client = _client()
        r = client.post("/execute", json={"command": cmd})
        assert r.status_code == 200
        assert "__plan__" in called_tools, (
            f"Expected __plan__ in called tools; got {called_tools}"
        )

    def test_fallback_normal_command_executes_via_broker(self, monkeypatch):
        """Fallback router executes a recognisable command via broker; history recorded."""
        mock_broker = MagicMock()
        mock_broker.execute_tool = AsyncMock(
            return_value={"status": "success", "text": "mutated", "result": {}}
        )
        monkeypatch.setattr(
            "backend.main_with_mcp._get_execution_broker", lambda: mock_broker
        )

        # Ensure fallback handles this as an "execute" intent
        mock_intent = MagicMock()
        mock_intent.intent = "execute"
        monkeypatch.setattr(
            "backend.intent_classifier.classify_intent", lambda cmd: mock_intent
        )

        client = _client()
        r = client.post(
            "/execute",
            json={"command": "mutate sequence ATCGATCG with substitution rate 0.01"},
        )
        assert r.status_code == 200
        payload = r.json()
        assert payload.get("success") is True
        session_id = payload.get("session_id")
        assert session_id

        runs = _get_runs(client, session_id)
        assert len(runs) >= 1


# ── Group G: Session-context side effects ─────────────────────────────────────


class TestSessionContextSideEffects:
    """
    _apply_session_context_side_effects is a pure function; test it directly
    without going through HTTP to keep these tests fast and deterministic.
    """

    def test_mutation_stores_variants_in_session_context(self):
        """mutate_sequence result with 'statistics.variants' → session populated."""
        from backend.history_manager import history_manager
        from backend.main_with_mcp import _apply_session_context_side_effects

        sid = "test-mutation-session"
        history_manager.sessions[sid] = {}

        result = {
            "status": "success",
            "statistics": {
                "variants": [
                    {"name": "mut1", "sequence": "ATCG"},
                    {"name": "mut2", "sequence": "TTCG"},
                ]
            },
        }
        _apply_session_context_side_effects(sid, "mutate_sequence", result)

        session = history_manager.sessions[sid]
        assert "mutated_sequences" in session
        assert len(session["mutated_sequences"]) == 2
        assert "mutation_results" in session
        assert session["mutation_results"] == session["mutated_sequences"]

    def test_alignment_stores_fasta_in_session_context(self):
        """sequence_alignment result → FASTA string stored in session."""
        from backend.history_manager import history_manager
        from backend.main_with_mcp import _apply_session_context_side_effects

        sid = "test-alignment-session"
        history_manager.sessions[sid] = {}

        result = {
            "status": "success",
            "alignment": [
                {"name": "seq1", "sequence": "ATCGATCG"},
                {"name": "seq2", "sequence": "TTCGATCG"},
            ],
        }
        _apply_session_context_side_effects(sid, "sequence_alignment", result)

        session = history_manager.sessions[sid]
        assert "aligned_sequences" in session
        fasta = session["aligned_sequences"]
        assert ">seq1" in fasta
        assert "ATCGATCG" in fasta
        assert ">seq2" in fasta

    def test_other_tools_do_not_pollute_session_context(self):
        """Unrelated tool results leave the session dict untouched."""
        from backend.history_manager import history_manager
        from backend.main_with_mcp import _apply_session_context_side_effects

        sid = "test-noop-session"
        history_manager.sessions[sid] = {}

        _apply_session_context_side_effects(
            sid, "fastqc_quality_analysis", {"status": "success"}
        )
        assert history_manager.sessions[sid] == {}

    def test_missing_session_id_does_not_crash(self):
        """If session_id is not in sessions, no exception is raised."""
        from backend.main_with_mcp import _apply_session_context_side_effects

        # Should silently do nothing — session doesn't exist
        _apply_session_context_side_effects(
            "no-such-session",
            "mutate_sequence",
            {"statistics": {"variants": [{"name": "x"}]}},
        )


# ── Group H: Error handling and pure-unit helpers ─────────────────────────────


class TestErrorHandling:
    def test_outer_exception_returns_structured_error(self, monkeypatch):
        """
        If session setup raises RuntimeError, execute() catches it and returns
        HTTP 200 with {success: False, error: '<message>', session_id: ...}.
        """

        def _boom(*args, **kwargs):
            raise RuntimeError("boom")

        monkeypatch.setattr(
            "backend.history_manager.history_manager.ensure_session_exists", _boom
        )

        client = _client()
        r = client.post(
            "/execute", json={"command": "anything", "session_id": "existing-sid"}
        )
        assert r.status_code == 200
        payload = r.json()
        assert payload.get("success") is False
        assert "boom" in payload.get("error", "")

    def test_is_success_classifies_all_variants(self):
        """_is_success covers all documented result shapes."""
        from backend.main_with_mcp import _is_success

        # Non-dict / unknown → success
        assert _is_success(None) is True
        assert _is_success("string result") is True
        assert _is_success(42) is True

        # Status-based
        assert _is_success({"status": "success"}) is True
        assert _is_success({"status": "completed"}) is True
        assert _is_success({"status": "error"}) is False
        assert _is_success({"status": "failed"}) is False
        assert _is_success({"status": "workflow_failed"}) is False

        # Explicit success field takes precedence
        assert _is_success({"success": True, "status": "error"}) is True
        assert _is_success({"success": False, "status": "success"}) is False

    def test_extract_metadata_top_level_fields(self):
        """_extract_metadata reads run_id/parent_run_id from the top-level dict."""
        from backend.main_with_mcp import _extract_metadata

        result = {
            "run_id": "r-top",
            "parent_run_id": "r-parent",
            "artifacts": [{"type": "csv", "uri": "s3://bucket/out.csv"}],
        }
        meta = _extract_metadata(result, tool_args={"alpha": 0.05})
        assert meta["run_id"] == "r-top"
        assert meta["parent_run_id"] == "r-parent"
        assert meta["tool_args"] == {"alpha": 0.05}
        assert meta["produced_artifacts"] == result["artifacts"]

    def test_extract_metadata_nested_result_subdict(self):
        """_extract_metadata falls back to nested 'result' sub-dict."""
        from backend.main_with_mcp import _extract_metadata

        result = {
            "status": "success",
            "result": {"run_id": "r-nested", "parent_run_id": "r-p-nested"},
        }
        meta = _extract_metadata(result)
        assert meta["run_id"] == "r-nested"
        assert meta["parent_run_id"] == "r-p-nested"

    def test_extract_metadata_non_dict_result(self):
        """_extract_metadata returns a valid dict even for non-dict results."""
        from backend.main_with_mcp import _extract_metadata

        meta = _extract_metadata("not a dict", tool_args={"x": 1})
        assert meta["inputs"] is None
        assert meta["produced_artifacts"] is None
        assert meta["tool_args"] == {"x": 1}
        assert meta["mcp_route"] == "/execute"
        # run_id key absent (non-dict path does not include it)
        assert "run_id" not in meta
