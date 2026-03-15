from benchmarks.bio_benchmark_scorer import score_benchmark_run


def _spec_and_rows(turn_id: str, user: str, expected_behaviors, row):
    spec = {
        "evaluation_dimensions": [
            {"id": "artifact_persistence", "severity": "high"},
            {"id": "reproducibility", "severity": "high"},
            {"id": "workflow_selection", "severity": "high"},
        ],
        "scenarios": [
            {
                "id": "s1",
                "turns": [
                    {
                        "id": turn_id,
                        "user": user,
                        "expected_behaviors": expected_behaviors,
                        "focus_dimensions": ["artifact_persistence", "reproducibility"],
                    }
                ],
            }
        ],
    }
    run_data = {"session_id": "sid", "scenario_id": "s1", "rows": [row]}
    return spec, run_data


def test_unresolved_rerun_request_not_full_credit():
    spec, run_data = _spec_and_rows(
        "turn_04",
        "The PCA suggests 2 outlier samples. Exclude them and rerun.",
        ["rerun_affected_downstream_steps", "resolve_referenced_pca_artifact"],
        {
            "turn_id": "turn_04",
            "tool": "bio_rerun",
            "status": "success",
            "text": "I could not find a prior rerunnable analysis in this session yet. Please run the base analysis first.",
        },
    )
    scored = score_benchmark_run(run_data, spec)
    turn = scored["per_turn_scores"][0]
    assert turn["score"] == 0
    assert "fallback_instead_of_execution" in turn["failure_patterns"]
    assert "rerun_not_performed" in turn["failure_patterns"]


def test_enrichment_substituted_with_lookup_not_full_credit():
    spec, run_data = _spec_and_rows(
        "turn_12",
        "Take the significantly upregulated genes and run GO enrichment.",
        ["create_gene_list_and_go_enrichment_artifacts", "resolve_current_deg_table"],
        {
            "turn_id": "turn_12",
            "tool": "lookup_go_term",
            "status": "success",
            "text": "Lookup GO term GO:0008150",
        },
    )
    scored = score_benchmark_run(run_data, spec)
    turn = scored["per_turn_scores"][0]
    assert turn["score"] == 0
    assert "enrichment_substituted_with_lookup" in turn["failure_patterns"]


def test_unresolved_historical_comparison_not_full_credit():
    spec, run_data = _spec_and_rows(
        "turn_18",
        "Compare the current DEG results with the first DEG results.",
        ["resolve_current_and_first_deg_tables", "compare_significance_and_effect_direction"],
        {
            "turn_id": "turn_18",
            "tool": "bio_diff_runs",
            "status": "success",
            "text": "I could not resolve artifact-backed historical runs for comparison in this session yet.",
        },
    )
    scored = score_benchmark_run(run_data, spec)
    turn = scored["per_turn_scores"][0]
    assert turn["score"] == 0
    assert "unresolved_historical_reference" in turn["failure_patterns"]
    assert "comparison_not_performed" in turn["failure_patterns"]


def test_empty_success_output_not_full_credit():
    spec, run_data = _spec_and_rows(
        "turn_16",
        "Now use the cleaned dataset from before batch exclusion.",
        ["resolve_correct_historical_artifact_if_present"],
        {
            "turn_id": "turn_16",
            "tool": "handle_natural_command",
            "status": "success",
            "text": "",
        },
    )
    scored = score_benchmark_run(run_data, spec)
    turn = scored["per_turn_scores"][0]
    assert turn["score"] == 0
    assert "empty_output" in turn["failure_patterns"]


def test_plan_staged_when_execution_expected_not_full_credit():
    spec, run_data = _spec_and_rows(
        "turn_08",
        "Fix it and regenerate the table and plot.",
        ["rerun_only_affected_steps_if_possible", "clearly_state_what_changed"],
        {
            "turn_id": "turn_08",
            "tool": "__plan__",
            "status": "workflow_planned",
            "text": "## Pipeline Plan\n\n1. step1",
        },
    )
    scored = score_benchmark_run(run_data, spec)
    turn = scored["per_turn_scores"][0]
    assert turn["score"] == 0
    assert "plan_instead_of_execution" in turn["failure_patterns"]


def test_historical_figure_recreation_unresolved_not_full_credit():
    spec, run_data = _spec_and_rows(
        "turn_19",
        "Recreate the figure set corresponding to the corrected metadata version before the fold-change bug fix.",
        ["resolve_exact_historical_state", "retrieve_or_regenerate_correct_figure_set"],
        {
            "turn_id": "turn_19",
            "tool": "bio_diff_runs",
            "status": "success",
            "text": "Run at least two persisted analyses first, then retry comparison.",
        },
    )
    scored = score_benchmark_run(run_data, spec)
    turn = scored["per_turn_scores"][0]
    assert turn["score"] == 0
    assert "missing_artifact_resolution" in turn["failure_patterns"]
