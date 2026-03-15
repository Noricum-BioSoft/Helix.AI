from backend.artifact_resolver import resolve_semantic_reference
from backend.main_with_mcp import (
    _is_approval_command,
    _preflight_tool_bindings,
    _requires_approval_semantics,
    _should_stage_for_approval,
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
