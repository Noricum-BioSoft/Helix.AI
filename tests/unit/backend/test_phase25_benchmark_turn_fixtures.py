from backend.artifact_resolver import resolve_semantic_reference
from backend.main import (
    _historical_recreation_ready_for_execution,
    _analyze_plan_execution,
    _is_approval_command,
    _preflight_tool_bindings,
    _requires_approval_semantics,
    _should_clear_pending_plan_after_execution,
    _should_stage_for_approval,
    build_standard_response,
)


def test_turn_03_approve_command_is_recognized():
    assert _is_approval_command("Approve.")


def test_turn_05_metadata_correction_stages_for_approval():
    prompt = "Actually sample S08 is mislabeled. It should be control, not treated."
    assert _requires_approval_semantics(prompt)
    assert _should_stage_for_approval("handle_natural_command", prompt, {}) is True


def test_turn_20_subset_request_stages_for_approval():
    prompt = "For the same dataset, now focus only on female samples and rerun the treated vs control analysis."
    assert _requires_approval_semantics(prompt)
    assert _should_stage_for_approval("handle_natural_command", prompt, {}) is True


def test_non_plan_binding_preflight_returns_needs_inputs():
    result = _preflight_tool_bindings("bulk_rnaseq_analysis", {"design_formula": "~condition"})
    assert isinstance(result, dict)
    assert result.get("status") == "needs_inputs"
    assert isinstance(result.get("binding_diagnostics", {}).get("issues"), list)


def test_patch_and_rerun_preflight_infers_change_request_from_command():
    result = _preflight_tool_bindings(
        "patch_and_rerun",
        {"session_id": "sid-1", "command": "Fix it and regenerate the table and plot."},
    )
    # Should no longer fail binding due to missing change_request
    assert result is None


def test_turn_18_reference_resolution_current_vs_first_deg():
    session = {
        "runs": [
            {"run_id": "run_v1", "tool": "bulk_rnaseq_analysis"},
            {"run_id": "run_v2", "tool": "bulk_rnaseq_analysis"},
        ],
        "artifacts": {},
    }
    first = resolve_semantic_reference(session, "first DEG results")
    current = resolve_semantic_reference(session, "current DEG results")
    assert first["status"] == "resolved" and first["target"]["run_id"] == "run_v1"
    assert current["status"] == "resolved" and current["target"]["run_id"] == "run_v2"


def test_turn_16_historical_reference_before_batch_exclusion():
    session = {
        "runs": [],
        "artifacts": [
            {
                "artifact_id": "a1",
                "title": "cleaned_dataset_v1",
                "state_tags": ["cleaned", "metadata_corrected"],
                "source_run_id": "run_clean",
            },
            {
                "artifact_id": "a2",
                "title": "cleaned_dataset_after_batch_exclusion_v2",
                "state_tags": ["cleaned", "batch_exclusion"],
                "source_run_id": "run_excluded",
            },
        ],
    }
    resolved = resolve_semantic_reference(session, "cleaned dataset from before batch exclusion")
    assert resolved["status"] == "resolved"
    assert resolved["target"]["artifact_id"] == "a1"


def test_turn_19_historical_reference_corrected_before_fold_change_fix():
    session = {
        "runs": [],
        "artifacts": [
            {
                "artifact_id": "fig_v2",
                "title": "figure_set_corrected_metadata_v2",
                "state_tags": ["corrected", "metadata"],
                "source_run_id": "run_v2",
            },
            {
                "artifact_id": "fig_v3",
                "title": "figure_set_corrected_metadata_after_fold_change_bug_fix_v3",
                "state_tags": ["corrected", "metadata", "fold_change_bug_fix"],
                "source_run_id": "run_v3",
            },
        ],
    }
    resolved = resolve_semantic_reference(
        session,
        "corrected metadata version before the fold-change bug fix",
    )
    assert resolved["status"] == "resolved"
    assert resolved["target"]["artifact_id"] == "fig_v2"


def test_approval_execution_marks_plan_as_executed_and_clearable():
    result = {
        "status": "success",
        "type": "execution_result",
        "result": {
            "status": "success",
            "type": "plan_result",
            "steps": [
                {
                    "id": "step1",
                    "tool_name": "bulk_rnaseq_analysis",
                    "result": {"status": "success", "text": "done"},
                }
            ],
        },
    }
    inner_steps = result["result"]["steps"]
    analyzed = _analyze_plan_execution(inner_steps)
    assert analyzed["executed"] is True
    assert analyzed["status"] == "success"
    assert _should_clear_pending_plan_after_execution(result) is True


def test_approval_needs_inputs_keeps_pending_plan_for_retry():
    result = {
        "status": "success",
        "type": "execution_result",
        "result": {
            "status": "success",
            "type": "plan_result",
            "steps": [
                {
                    "id": "step1",
                    "tool_name": "bulk_rnaseq_analysis",
                    "result": {"status": "needs_inputs", "needs_inputs": True},
                }
            ],
        },
    }
    assert _should_clear_pending_plan_after_execution(result) is False


def test_approval_async_submission_is_not_rendered_as_workflow_planned():
    rendered = build_standard_response(
        prompt="Approve.",
        tool="__plan__",
        result={"type": "job", "status": "submitted", "job_id": "job_123"},
        session_id="sid",
        mcp_route="/execute",
        success=True,
    )
    assert rendered["status"] == "pipeline_submitted"


def test_empty_success_output_is_upgraded_to_actionable_message():
    rendered = build_standard_response(
        prompt="Compare first and current results.",
        tool="bio_diff_runs",
        result={"status": "success", "text": ""},
        session_id="sid",
        mcp_route="/execute",
        success=True,
    )
    assert rendered["status"] == "success"
    assert isinstance(rendered.get("text"), str) and rendered["text"].strip()
    assert "historical state" in rendered["text"].lower() or "comparison" in rendered["text"].lower()


def test_empty_go_enrichment_output_is_actionable():
    rendered = build_standard_response(
        prompt="Run GO enrichment.",
        tool="go_enrichment_analysis",
        result={"status": "success", "text": ""},
        session_id="sid",
        mcp_route="/execute",
        success=True,
    )
    assert rendered["status"] == "success"
    assert "gene list" in rendered["text"].lower()


def test_empty_agent_output_is_actionable():
    rendered = build_standard_response(
        prompt="Now use the cleaned dataset from before batch exclusion.",
        tool="agent",
        result={"status": "success", "text": ""},
        session_id="sid",
        mcp_route="/execute",
        success=True,
    )
    assert rendered["status"] == "success"
    assert rendered["text"].strip()


def test_patch_and_rerun_runtime_import_failure_returns_structured_needs_inputs(monkeypatch, tmp_path):
    import importlib
    import sys
    import types
    from backend.history_manager import history_manager

    def _tool_decorator(fn):
        fn.func = fn
        return fn

    fake_langchain_core = types.ModuleType("langchain_core")
    fake_langchain_tools = types.ModuleType("langchain_core.tools")
    fake_langchain_tools.tool = _tool_decorator
    fake_langchain_core.tools = fake_langchain_tools
    monkeypatch.setitem(sys.modules, "langchain_core", fake_langchain_core)
    monkeypatch.setitem(sys.modules, "langchain_core.tools", fake_langchain_tools)
    import backend.agent_tools as _agent_tools
    importlib.reload(_agent_tools)

    session_id = "sid_patch_turn9"
    history_manager.ensure_session_exists(session_id)
    run_dir = tmp_path / "runs" / "run_base"
    run_dir.mkdir(parents=True, exist_ok=True)
    script_path = run_dir / "analysis.py"
    script_path.write_text("from backend import bulk_rnaseq as _mod\nprint('{}')\n")
    history_manager.add_run(
        session_id=session_id,
        command="base",
        tool="bulk_rnaseq_analysis",
        result={"status": "success", "text": "ok"},
        run_id="run_base",
        tool_args={
            "count_matrix": "counts.csv",
            "sample_metadata": "meta.csv",
            "design_formula": "~condition",
            "alpha": 0.05,
        },
        produced_artifacts=[{"type": "script", "uri": str(script_path), "title": "analysis.py"}],
    )

    monkeypatch.setattr("backend.agent_tools._llm_patch_script", lambda s, c: s)
    monkeypatch.setattr(
        "backend.script_executor.execute",
        lambda *args, **kwargs: {
            "status": "error",
            "error": "ImportError: cannot import name 'bulk_rnaseq' from 'backend'",
            "logs": "Traceback ... ImportError ...",
        },
    )

    result = _agent_tools.patch_and_rerun.func(
        session_id=session_id,
        change_request="Now highlight IL6, TNF, CXCL8, and NFKB1 on the volcano plot.",
        target_run="latest",
    )
    assert result["status"] == "needs_inputs"
    assert result.get("diagnostics", {}).get("issue") == "script_runtime_import_failure"
    assert isinstance(result.get("script_path"), str)


def test_patch_and_rerun_runtime_import_failure_heatmap_request_structured(monkeypatch, tmp_path):
    import importlib
    import sys
    import types
    from backend.history_manager import history_manager

    def _tool_decorator(fn):
        fn.func = fn
        return fn

    fake_langchain_core = types.ModuleType("langchain_core")
    fake_langchain_tools = types.ModuleType("langchain_core.tools")
    fake_langchain_tools.tool = _tool_decorator
    fake_langchain_core.tools = fake_langchain_tools
    monkeypatch.setitem(sys.modules, "langchain_core", fake_langchain_core)
    monkeypatch.setitem(sys.modules, "langchain_core.tools", fake_langchain_tools)
    import backend.agent_tools as _agent_tools
    importlib.reload(_agent_tools)

    session_id = "sid_patch_turn10"
    history_manager.ensure_session_exists(session_id)
    run_dir = tmp_path / "runs" / "run_base_2"
    run_dir.mkdir(parents=True, exist_ok=True)
    script_path = run_dir / "analysis.py"
    script_path.write_text("from backend import bulk_rnaseq as _mod\nprint('{}')\n")
    history_manager.add_run(
        session_id=session_id,
        command="base",
        tool="bulk_rnaseq_analysis",
        result={"status": "success", "text": "ok"},
        run_id="run_base_2",
        tool_args={
            "count_matrix": "counts.csv",
            "sample_metadata": "meta.csv",
            "design_formula": "~condition",
            "alpha": 0.05,
        },
        produced_artifacts=[{"type": "script", "uri": str(script_path), "title": "analysis.py"}],
    )

    monkeypatch.setattr("backend.agent_tools._llm_patch_script", lambda s, c: s)
    monkeypatch.setattr(
        "backend.script_executor.execute",
        lambda *args, **kwargs: {
            "status": "error",
            "error": "No module named backend.bulk_rnaseq",
            "logs": "Traceback ... No module named ...",
        },
    )

    result = _agent_tools.patch_and_rerun.func(
        session_id=session_id,
        change_request="Also give me a heatmap of the top 30 significant genes.",
        target_run="latest",
    )
    assert result["status"] == "needs_inputs"
    assert result.get("diagnostics", {}).get("issue") == "script_runtime_import_failure"
    assert "runtime import dependency issue" in result.get("text", "").lower()


def test_historical_recreation_with_resolved_state_skips_planning_gate():
    session = {
        "runs": [{"run_id": "run_v2", "tool": "bulk_rnaseq_analysis"}],
        "artifacts": [
            {
                "artifact_id": "fig_v2",
                "title": "figure_set_corrected_metadata_v2",
                "state_tags": ["corrected", "metadata", "pre_bugfix"],
                "source_run_id": "run_v2",
            }
        ],
    }
    cmd = "Recreate the figure set corresponding to the corrected metadata version before the fold-change bug fix."
    assert _historical_recreation_ready_for_execution(cmd, "bio_diff_runs", session) is True


def test_script_executor_retries_host_when_sandbox_import_fails(monkeypatch, tmp_path):
    from backend import script_executor as se

    script_path = tmp_path / "analysis.py"
    script_path.write_text("print('ok')\n")

    monkeypatch.setenv("HELIX_ANALYSIS_USE_SANDBOX", "true")
    monkeypatch.setattr(
        se,
        "_run_in_sandbox",
        lambda *args, **kwargs: {
            "status": "error",
            "error": "ImportError: cannot import name 'bulk_rnaseq' from backend",
            "logs": "Traceback ...",
            "execution_mode": "sandbox",
        },
    )
    monkeypatch.setattr(
        se,
        "_run_on_host",
        lambda *args, **kwargs: {
            "status": "success",
            "summary": {"status": "success"},
            "artifacts": [],
            "logs": "",
            "execution_mode": "host",
        },
    )

    result = se.execute(script_path)
    assert result["status"] == "success"
    assert result.get("execution_mode") == "host"


def test_patch_and_rerun_turn9_style_request_executes_successfully(monkeypatch, tmp_path):
    import importlib
    import sys
    import types
    from backend.history_manager import history_manager

    def _tool_decorator(fn):
        fn.func = fn
        return fn

    fake_langchain_core = types.ModuleType("langchain_core")
    fake_langchain_tools = types.ModuleType("langchain_core.tools")
    fake_langchain_tools.tool = _tool_decorator
    fake_langchain_core.tools = fake_langchain_tools
    monkeypatch.setitem(sys.modules, "langchain_core", fake_langchain_core)
    monkeypatch.setitem(sys.modules, "langchain_core.tools", fake_langchain_tools)
    import backend.agent_tools as _agent_tools
    importlib.reload(_agent_tools)

    session_id = "sid_patch_turn9_success"
    history_manager.ensure_session_exists(session_id)
    run_dir = tmp_path / "runs" / "run_base_success_9"
    run_dir.mkdir(parents=True, exist_ok=True)
    script_path = run_dir / "analysis.py"
    script_path.write_text("print('{}')\n")
    history_manager.add_run(
        session_id=session_id,
        command="base",
        tool="bulk_rnaseq_analysis",
        result={"status": "success", "text": "ok"},
        run_id="run_base_success_9",
        tool_args={"count_matrix": "counts.csv", "sample_metadata": "meta.csv"},
        produced_artifacts=[{"type": "script", "uri": str(script_path), "title": "analysis.py"}],
    )

    monkeypatch.setattr("backend.agent_tools._llm_patch_script", lambda s, c: s)
    monkeypatch.setattr(
        "backend.script_executor.execute",
        lambda *args, **kwargs: {
            "status": "success",
            "summary": {"plot_paths": {}, "status": "success"},
            "artifacts": [],
        },
    )

    result = _agent_tools.patch_and_rerun.func(
        session_id=session_id,
        change_request="Now highlight IL6, TNF, CXCL8, and NFKB1 on the volcano plot.",
        target_run="latest",
    )
    assert result["status"] == "success"
    assert "Analysis Updated" in result.get("text", "")


def test_patch_and_rerun_turn10_style_request_executes_successfully(monkeypatch, tmp_path):
    import importlib
    import sys
    import types
    from backend.history_manager import history_manager

    def _tool_decorator(fn):
        fn.func = fn
        return fn

    fake_langchain_core = types.ModuleType("langchain_core")
    fake_langchain_tools = types.ModuleType("langchain_core.tools")
    fake_langchain_tools.tool = _tool_decorator
    fake_langchain_core.tools = fake_langchain_tools
    monkeypatch.setitem(sys.modules, "langchain_core", fake_langchain_core)
    monkeypatch.setitem(sys.modules, "langchain_core.tools", fake_langchain_tools)
    import backend.agent_tools as _agent_tools
    importlib.reload(_agent_tools)

    session_id = "sid_patch_turn10_success"
    history_manager.ensure_session_exists(session_id)
    run_dir = tmp_path / "runs" / "run_base_success_10"
    run_dir.mkdir(parents=True, exist_ok=True)
    script_path = run_dir / "analysis.py"
    script_path.write_text("print('{}')\n")
    history_manager.add_run(
        session_id=session_id,
        command="base",
        tool="bulk_rnaseq_analysis",
        result={"status": "success", "text": "ok"},
        run_id="run_base_success_10",
        tool_args={"count_matrix": "counts.csv", "sample_metadata": "meta.csv"},
        produced_artifacts=[{"type": "script", "uri": str(script_path), "title": "analysis.py"}],
    )

    monkeypatch.setattr("backend.agent_tools._llm_patch_script", lambda s, c: s)
    monkeypatch.setattr(
        "backend.script_executor.execute",
        lambda *args, **kwargs: {
            "status": "success",
            "summary": {"plot_paths": {}, "status": "success"},
            "artifacts": [],
        },
    )

    result = _agent_tools.patch_and_rerun.func(
        session_id=session_id,
        change_request="Also give me a heatmap of the top 30 significant genes.",
        target_run="latest",
    )
    assert result["status"] == "success"
    assert "Analysis Updated" in result.get("text", "")
