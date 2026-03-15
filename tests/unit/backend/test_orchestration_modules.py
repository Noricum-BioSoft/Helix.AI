from backend.orchestration.approval_policy import (
    is_approval_command,
    requires_approval_semantics,
    should_stage_for_approval,
)
from backend.orchestration.execution_router import normalize_tool_selection
from backend.orchestration.artifact_resolver import resolve_semantic_reference
from backend.orchestration.state_selection import select_diff_anchors
from backend.orchestration.visualization_resolver import determine_visualization_type
from backend.command_router import CommandRouter
from backend.action_plan import infer_action_type, map_action_to_tool


def test_approval_policy_is_action_based():
    prompt = "Actually sample S08 is mislabeled. It should be control, not treated."
    assert is_approval_command("Approve.")
    assert requires_approval_semantics(prompt, action_type="correct_metadata")
    assert should_stage_for_approval("handle_natural_command", prompt, {}, action_type="correct_metadata")


def test_visualization_resolver_prefers_artifact_metadata():
    result = {"artifact_kind": "deg_table", "plot_family": "volcano"}
    vis = determine_visualization_type("handle_natural_command", result, "show plot")
    assert vis == "results_viewer"


def test_semantic_resolver_selectors():
    session = {
        "runs": [
            {"run_id": "run_v1", "tool": "bulk_rnaseq_analysis"},
            {"run_id": "run_v2", "tool": "bulk_rnaseq_analysis"},
        ],
        "artifacts": [
            {"artifact_id": "p1", "title": "pca_plot_v1"},
            {"artifact_id": "p2", "title": "pca_plot_v2"},
            {"artifact_id": "f1", "title": "filtered_dataset_before_batch_exclusion_v1"},
        ],
    }
    assert resolve_semantic_reference(session, "current DEG results")["target"]["run_id"] == "run_v2"
    assert resolve_semantic_reference(session, "earlier PCA")["target"]["artifact_id"] == "p1"
    assert resolve_semantic_reference(session, "prior filtered dataset")["target"]["artifact_id"] == "f1"


def test_approval_policy_does_not_stage_execute_viz_iterations():
    turn4 = "The PCA suggests 2 outlier samples. Exclude them and rerun."
    turn11 = "Actually make that top 50 genes, clustered by sample, and split by condition."
    assert should_stage_for_approval("bio_rerun", turn4, {}, action_type="subset_data") is False
    assert should_stage_for_approval("bulk_rnaseq_analysis", turn11, {}, action_type="generate_plot") is False


def test_execution_router_normalizes_wrong_single_cell_modality():
    session = {"runs": [{"run_id": "run_1", "tool": "bulk_rnaseq_analysis"}]}
    norm = normalize_tool_selection(
        command="Color the PCA by batch and sex and show both.",
        tool_name="single_cell_analysis",
        parameters={},
        session_context=session,
    )
    assert norm.tool_name in {"bio_rerun", "patch_and_rerun"}


def test_state_selection_diff_anchors_use_semantic_runs():
    session = {
        "runs": [
            {"run_id": "run_v1", "tool": "bulk_rnaseq_analysis"},
            {"run_id": "run_v2", "tool": "bulk_rnaseq_analysis"},
        ],
        "artifacts": {},
    }
    run_a, run_b, _diag = select_diff_anchors(
        session,
        "latest",
        "prior",
        "Compare current DEG results with first DEG results",
    )
    assert run_a == "run_v2"
    assert run_b == "run_v1"


def test_semantic_resolver_prefers_semantic_alias_and_structured_tags():
    session = {
        "runs": [],
        "artifacts": [
            {
                "artifact_id": "fig_old",
                "title": "figure_set_corrected_metadata_v2",
                "state_tags": ["corrected", "metadata"],
                "semantic_aliases": ["corrected_metadata_pre_bugfix"],
                "source_run_id": "run_v2",
            },
            {
                "artifact_id": "fig_new",
                "title": "figure_set_corrected_metadata_after_fold_change_bug_fix_v3",
                "state_tags": ["corrected", "metadata", "fold_change_bug_fix"],
                "source_run_id": "run_v3",
            },
        ],
    }
    via_alias = resolve_semantic_reference(session, "corrected_metadata_pre_bugfix")
    via_phrase = resolve_semantic_reference(
        session, "corrected metadata version before the fold-change bug fix"
    )
    assert via_alias["status"] == "resolved"
    assert via_alias["target"]["artifact_id"] == "fig_old"
    assert via_phrase["status"] == "resolved"
    assert via_phrase["target"]["artifact_id"] == "fig_old"


def test_execution_router_prefers_plot_update_tools_when_plot_artifacts_exist():
    session = {
        "runs": [{"run_id": "run_1", "tool": "bulk_rnaseq_analysis"}],
        "artifacts": {
            "a1": {"artifact_id": "a1", "artifact_kind": "pca_plot", "title": "PCA by condition"}
        },
    }
    norm = normalize_tool_selection(
        command="Update the volcano plot and highlight top 30 genes.",
        tool_name="bulk_rnaseq_analysis",
        parameters={},
        session_context=session,
    )
    assert norm.tool_name == "patch_and_rerun"


def test_state_selection_diff_anchors_resolve_before_state_from_metadata():
    session = {
        "runs": [
            {"run_id": "run_v2", "tool": "bulk_rnaseq_analysis"},
            {"run_id": "run_v3", "tool": "bulk_rnaseq_analysis"},
        ],
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
    run_a, run_b, diag = select_diff_anchors(
        session,
        "latest",
        "prior",
        "Recreate the figure set corresponding to the corrected metadata version before the fold-change bug fix.",
    )
    assert run_b == "run_v2"
    assert isinstance(diag.get("selectors"), list)


def test_state_selection_resolves_run_id_from_artifact_selector():
    session = {
        "runs": [{"run_id": "run_v1", "tool": "bulk_rnaseq_analysis"}],
        "artifacts": [
            {
                "artifact_id": "deg_v1",
                "title": "deg_results_v1",
                "artifact_kind": "deg_table",
                "source_run_id": "run_v1",
                "semantic_aliases": ["first_deg_results"],
            }
        ],
    }
    run_a, run_b, _diag = select_diff_anchors(
        session,
        "latest",
        "prior",
        "Compare current DEG results with first DEG results",
    )
    assert run_a == "run_v1"
    assert run_b == "run_v1"


def test_go_enrichment_routes_to_enrichment_tool_class():
    router = CommandRouter()
    tool, params = router.route_command(
        "Take the significantly upregulated genes and run GO enrichment.",
        {"session_id": "sid-1"},
    )
    assert tool == "go_enrichment_analysis"
    assert params.get("session_id") == "sid-1"
    action = infer_action_type("Take the significantly upregulated genes and run GO enrichment.", tool)
    assert action == "run_enrichment"
    assert map_action_to_tool(action, None, {"command": "run GO enrichment"}) == "go_enrichment_analysis"


def test_execution_router_normalizes_enrichment_intent_to_enrichment_tool():
    session = {"session_id": "sid-1", "runs": [{"run_id": "run_1", "tool": "bulk_rnaseq_analysis"}]}
    norm = normalize_tool_selection(
        command="Take the significantly upregulated genes and run GO enrichment.",
        tool_name="bulk_rnaseq_analysis",
        parameters={"session_id": "sid-1"},
        session_context=session,
    )
    assert norm.tool_name == "go_enrichment_analysis"
    assert norm.arguments.get("source_selector") == "current DEG results"


def test_command_router_routes_historical_recreation_to_diff_runs():
    router = CommandRouter()
    tool, params = router.route_command(
        "Recreate the figure set corresponding to the corrected metadata version before the fold-change bug fix.",
        {"session_id": "sid-1"},
    )
    assert tool == "bio_diff_runs"
    assert params.get("run_id_a") == "latest"
    assert params.get("run_id_b") == "prior"

