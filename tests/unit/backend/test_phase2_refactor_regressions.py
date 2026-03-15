from backend.action_plan import infer_action_type, map_action_to_tool
from backend.command_router import CommandRouter
from backend.plan_binding import validate_tool_bindings


def test_action_plan_infers_generic_metadata_correction_action():
    action = infer_action_type("Sample S08 is mislabeled and should be control.")
    assert action == "correct_metadata"


def test_action_plan_mapping_prefers_compare_versions():
    tool = map_action_to_tool("compare_versions", None, {})
    assert tool == "bio_diff_runs"


def test_router_routes_pca_color_update_to_rerun_path():
    router = CommandRouter()
    tool, params = router.route_command("Color the PCA by batch and sex and show both.", {})
    assert tool in {"bio_rerun", "patch_and_rerun", "bulk_rnaseq_analysis", "single_cell_analysis", "handle_natural_command"}
    # Regression expectation from Phase 2: avoid direct single-cell reroute for this prompt.
    assert tool != "single_cell_analysis"
    assert isinstance(params, dict)


def test_plan_binding_validator_flags_missing_required_bulk_inputs():
    issues = validate_tool_bindings("bulk_rnaseq_analysis", {"design_formula": "~condition"})
    assert issues
    keys = {i.get("input") for i in issues}
    assert "count_matrix" in keys and "sample_metadata" in keys

