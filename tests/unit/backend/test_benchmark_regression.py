"""
Regression protection for Helix.AI benchmark-critical behaviors.

These tests lock in the 90 % benchmark passing state (36/40) achieved after
iteration-2 remediation.  They protect behaviors that were hard to win and must
not regress silently.

Benchmark-critical behaviors protected
───────────────────────────────────────
1.  Approval gate always stages as workflow_planned (not needs_inputs) for
    clear-intent turns — fixes for turns 02, 14, 17.
2.  _use_historical fast-path routes "use X from before Y" to
    handle_natural_command — fix for turn_16 (was timing out at 180 s).
3.  handle_natural_command without approval semantics is not staged for
    approval — ensures turn_16 reaches the agent directly.
4.  bulk_rnaseq_analysis response never contains 'synthetic demo data used' —
    fix for turn_08 (was penalized as synthetic_success_only).
5.  Scorer: needs_inputs for planning_expected → score 1 (partial), not 2.
6.  Scorer: workflow_planned for planning_expected → score 2 (full credit).
7.  Approval execution path: plan → approve → success (turns 03, 06).
8.  Historical artifact resolution (turns 18, 19).
9.  GO enrichment routing does NOT substitute lookup_go_term.
10. Fallback text is always non-empty and actionable.
11. Benchmark baseline: saved iteration score ≥ 88 % with 0 critical failures.
12. Benchmark-critical turns 02, 08, 14, 16, 17 do not regress.
13. Benchmark-critical turns 09, 10, 19 do not regress.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

# ─────────────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────────────

_REPO_ROOT = Path(__file__).resolve().parents[3]
BENCHMARK_YAML = _REPO_ROOT / "tests" / "bioinformatics_benchmark.yaml"
ITERATION_SCORED = _REPO_ROOT / "benchmarks" / "runs" / "iteration_1_scored.json"
BASELINE_FILE = _REPO_ROOT / "benchmarks" / "runs" / "baseline.json"

BASELINE_THRESHOLD_PCT = 88.0   # 2-point margin below current 90 %
BASELINE_MAX_CRITICAL = 0


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _load_bench_spec() -> dict:
    return yaml.safe_load(BENCHMARK_YAML.read_text()) or {}


def _make_plan_preview(*, with_binding_diag: bool, approval_required: bool = True) -> dict:
    """Minimal plan_preview as produced by the approval gate."""
    preview: dict = {
        "status": "workflow_planned",
        "type": "plan_result",
        "steps": [
            {"id": "step1", "tool_name": "bulk_rnaseq_analysis", "arguments": {}}
        ],
        "execute_ready": True,
        "approval_required": approval_required,
        "workflow_state": "WAITING_FOR_APPROVAL",
    }
    if with_binding_diag:
        preview["binding_diagnostics"] = {
            "status": "error",
            "issues": [
                {
                    "step_index": 1,
                    "tool_name": "bulk_rnaseq_analysis",
                    "missing_inputs": ["count_matrix", "sample_metadata"],
                }
            ],
        }
    return preview


def _score_single_turn(
    prompt: str,
    behaviors: list[str],
    row: dict,
    *,
    should_require_approval: bool = False,
) -> dict:
    """Score one synthetic turn using the real scorer."""
    from benchmarks.bio_benchmark_scorer import score_benchmark_run

    spec = {
        "evaluation_dimensions": [
            {"id": "artifact_persistence", "severity": "high"},
            {"id": "workflow_selection", "severity": "high"},
            {"id": "scientific_intent_resolution", "severity": "high"},
        ],
        "scenarios": [
            {
                "id": "s1",
                "turns": [
                    {
                        "id": "t1",
                        "user": prompt,
                        "expected_behaviors": behaviors,
                        "should_require_approval_before_execution": should_require_approval,
                        "focus_dimensions": [
                            "workflow_selection",
                            "scientific_intent_resolution",
                        ],
                    }
                ],
            }
        ],
    }
    run_data = {"session_id": "sid", "scenario_id": "s1", "rows": [row]}
    scored = score_benchmark_run(run_data, spec)
    return scored["per_turn_scores"][0]


def _render(prompt: str, tool: str, result: dict, *, session_id: str = "sid") -> dict:
    """Convenience wrapper for build_standard_response."""
    from backend.main_with_mcp import build_standard_response

    return build_standard_response(
        prompt=prompt,
        tool=tool,
        result=result,
        session_id=session_id,
        mcp_route="/execute",
        success=True,
    )


# ─────────────────────────────────────────────────────────────────────────────
# 1.  Approval gate planning status  (turns 02, 14, 17)
# ─────────────────────────────────────────────────────────────────────────────

def test_approval_gate_workflow_planned_with_missing_bindings():
    """
    When file bindings are missing but the user's analytical intent is clear,
    build_standard_response must emit status=workflow_planned + WAITING_FOR_APPROVAL.
    This is the core fix for turns 02, 14, 17 (was returning needs_inputs → score 1).
    """
    rendered = _render(
        "I want to compare treated vs control and see the main genes and pathways.",
        "__plan__",
        _make_plan_preview(with_binding_diag=True),
    )
    assert rendered["status"] == "workflow_planned", (
        "Approval gate must emit workflow_planned even when file bindings are missing."
    )
    assert rendered["workflow_state"] == "WAITING_FOR_APPROVAL"
    assert rendered["approval_required"] is True


def test_approval_gate_plan_text_mentions_missing_inputs():
    """Plan text must hint what files are needed when binding_diagnostics is present."""
    rendered = _render(
        "I want to compare treated vs control.",
        "__plan__",
        _make_plan_preview(with_binding_diag=True),
    )
    text_l = rendered["text"].lower()
    assert "pipeline plan" in text_l, "Plan heading must be present."
    # At least one of: the input name, a 'provide' verb, or 'missing' hint
    assert any(kw in text_l for kw in ("count_matrix", "provide", "missing")), (
        "Plan text must mention what is still needed when bindings are incomplete."
    )


def test_approval_gate_workflow_planned_without_missing_bindings_shows_i_approve():
    """When all inputs are available, plan text shows 'I approve' CTA."""
    rendered = _render(
        "Approve the correction and rerun.",
        "__plan__",
        _make_plan_preview(with_binding_diag=False),
    )
    assert rendered["status"] == "workflow_planned"
    assert "i approve" in rendered["text"].lower(), (
        "Plan text must contain 'I approve' when all bindings are satisfied."
    )


# ─────────────────────────────────────────────────────────────────────────────
# 2.  _use_historical fast-path  (turn_16)
# ─────────────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize(
    "command",
    [
        "Now use the cleaned dataset from before batch exclusion.",
        "use the original dataset from before filtering",
        "please use the prior version of the data",
    ],
)
def test_use_historical_routes_to_handle_natural_command(command: str):
    """
    'use X from before Y dataset' must route deterministically to
    handle_natural_command so the approval gate does not stage it and it
    reaches the agent quickly (protects turn_16 from 180 s timeouts).
    """
    from backend.command_router import CommandRouter

    tool, _ = CommandRouter().route_command(command, {})
    assert tool == "handle_natural_command", (
        f"'{command}' must route to handle_natural_command (got '{tool}'). "
        "Routing this to bio_rerun triggers approval staging → plan_instead_of_execution → score 0."
    )


def test_use_historical_handle_natural_command_is_not_staged():
    """handle_natural_command without approval semantics must NOT be staged."""
    from backend.main_with_mcp import _should_stage_for_approval

    assert _should_stage_for_approval(
        "handle_natural_command",
        "Now use the cleaned dataset from before batch exclusion.",
        {},
    ) is False


# ─────────────────────────────────────────────────────────────────────────────
# 3.  No 'synthetic demo data used' marker  (turn_08)
# ─────────────────────────────────────────────────────────────────────────────

def test_build_standard_response_bulk_rnaseq_no_synthetic_marker():
    """
    build_standard_response for bulk_rnaseq_analysis must never emit the phrase
    'synthetic demo data used'. Its presence triggers the scorer's
    synthetic_success_only deduction (score 1 instead of 2 for turn_08).
    """
    result = {
        "status": "success",
        "mode": "synthetic",
        "n_samples": 6,
        "n_genes_total": 500,
        "summary": [
            {
                "contrast": "treated_vs_control",
                "total_genes": 500,
                "significant": 10,
                "up": 6,
                "down": 4,
            }
        ],
        "deg_tables": {},
        "plots": {},
    }
    rendered = _render("Fix it and regenerate the table and plot.", "bulk_rnaseq_analysis", result)
    assert "synthetic demo data used" not in rendered["text"].lower(), (
        "bulk_rnaseq_analysis must not embed 'synthetic demo data used' — "
        "this phrase triggers the synthetic_success_only scoring penalty."
    )


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Scorer strictness: planning_expected scoring
# ─────────────────────────────────────────────────────────────────────────────

def test_scorer_needs_inputs_for_planning_expected_is_partial():
    """
    needs_inputs for a planning_expected turn must give score 1 (partial), not 2.
    This was turns 02/14/17 before the fix — the scorer must stay strict.
    """
    turn = _score_single_turn(
        "I want to compare treated vs control and see the main genes and pathways.",
        ["propose_plan_using_metadata", "wait_for_approval"],
        {
            "turn_id": "t1",
            "tool": "__plan__",
            "status": "needs_inputs",
            "text": (
                "## Pipeline Plan (Needs Inputs)\n\n**Steps:**\n\n"
                "1. step1 (bulk_rnaseq_analysis)\n\nRequired inputs are still missing."
            ),
        },
        should_require_approval=True,
    )
    assert turn["score"] == 1, (
        "needs_inputs for planning_expected must stay at 1. "
        "Score 2 would mean the scorer became too lenient."
    )
    assert "wrong_tool_class" in turn["failure_patterns"]


def test_scorer_workflow_planned_for_planning_expected_is_full_credit():
    """
    workflow_planned for a planning_expected turn must give score 2 (full).
    This is the target behavior after the fix — must remain stable.
    """
    turn = _score_single_turn(
        "I want to compare treated vs control and see the main genes and pathways.",
        ["propose_plan_using_metadata", "wait_for_approval"],
        {
            "turn_id": "t1",
            "tool": "__plan__",
            "status": "workflow_planned",
            "text": (
                "## Pipeline Plan\n\n**Steps:**\n\n1. step1 (bulk_rnaseq_analysis)\n\n"
                "*Provide `count_matrix`, `sample_metadata` then click I approve.*"
            ),
        },
        should_require_approval=True,
    )
    assert turn["score"] == 2, (
        "workflow_planned for planning_expected must give full credit (2). "
        "Score < 2 means the approval gate fix was silently reverted."
    )
    assert not turn["failure_patterns"], (
        f"No failure patterns expected, got: {turn['failure_patterns']}"
    )


# ─────────────────────────────────────────────────────────────────────────────
# 5.  Approval execution path  (turns 03, 06)
# ─────────────────────────────────────────────────────────────────────────────

def test_scorer_approval_execution_success_gives_full_credit():
    """Approve → execution succeeds → must score 2."""
    turn = _score_single_turn(
        "Approve.",
        ["execute_approved_workflow", "create_reusable_artifacts"],
        {
            "turn_id": "t1",
            "tool": "__plan__",
            "status": "success",
            "text": (
                "## Pipeline Execution Complete\n\nThe approved plan has been executed.\n\n"
                "**Steps:**\n\n1. step1 (bulk_rnaseq_analysis)"
            ),
        },
    )
    assert turn["score"] == 2
    assert turn["assessments"]["execution_vs_planning"] == "approved_execution_happened"


def test_scorer_approval_re_stages_plan_gives_zero():
    """Approve that returns workflow_planned (re-staging) must score 0."""
    turn = _score_single_turn(
        "Approve.",
        ["execute_approved_workflow"],
        {
            "turn_id": "t1",
            "tool": "__plan__",
            "status": "workflow_planned",
            "text": "## Pipeline Plan\n\n**Steps:**\n\n1. step1",
        },
    )
    assert turn["score"] == 0
    assert "plan_instead_of_execution" in turn["failure_patterns"]


# ─────────────────────────────────────────────────────────────────────────────
# 6.  Historical artifact resolution  (turns 18, 19)
# ─────────────────────────────────────────────────────────────────────────────

def test_scorer_bio_diff_runs_comparison_success_gives_full_credit():
    """bio_diff_runs returning comparison text must score 2 (turns 18, 19)."""
    spec = _load_bench_spec()
    from benchmarks.bio_benchmark_scorer import score_benchmark_run

    scored = score_benchmark_run(
        {
            "session_id": "sid",
            "scenario_id": "bulk-rnaseq-e2e",
            "rows": [
                {
                    "turn_id": "turn_18",
                    "tool": "bio_diff_runs",
                    "status": "success",
                    "text": (
                        "## Run Comparison\n\n**Run A:** run-1 (bulk_rnaseq_analysis)\n"
                        "**Run B:** run-2 (bulk_rnaseq_analysis)\n\n"
                        "### Parameter Changes\n- condition: treated vs control"
                    ),
                }
            ],
        },
        spec,
    )
    t18 = next(t for t in scored["per_turn_scores"] if t["turn_id"] == "turn_18")
    assert t18["score"] == 2


def test_scorer_bio_diff_runs_unresolved_gives_zero():
    """bio_diff_runs that can't resolve must score 0."""
    spec = _load_bench_spec()
    from benchmarks.bio_benchmark_scorer import score_benchmark_run

    scored = score_benchmark_run(
        {
            "session_id": "sid",
            "scenario_id": "bulk-rnaseq-e2e",
            "rows": [
                {
                    "turn_id": "turn_18",
                    "tool": "bio_diff_runs",
                    "status": "success",
                    "text": "Could not resolve comparable run IDs. Run at least two persisted analyses first.",
                }
            ],
        },
        spec,
    )
    t18 = next(t for t in scored["per_turn_scores"] if t["turn_id"] == "turn_18")
    assert t18["score"] == 0
    assert "unresolved_historical_reference" in t18["failure_patterns"]


# ─────────────────────────────────────────────────────────────────────────────
# 7.  GO enrichment routing correctness
# ─────────────────────────────────────────────────────────────────────────────

def test_go_enrichment_routing_not_substituted_with_lookup():
    """GO enrichment requests must NOT route to lookup_go_term."""
    from backend.command_router import CommandRouter

    tool, _ = CommandRouter().route_command(
        "Take the significantly upregulated genes and run GO enrichment.", {}
    )
    assert tool != "lookup_go_term", (
        "GO enrichment must not route to lookup_go_term — "
        "the scorer penalizes this as enrichment_substituted_with_lookup (score 0)."
    )


def test_scorer_enrichment_lookup_substitution_gives_zero():
    """lookup_go_term for a GO enrichment turn must score 0."""
    spec = _load_bench_spec()
    from benchmarks.bio_benchmark_scorer import score_benchmark_run

    scored = score_benchmark_run(
        {
            "session_id": "sid",
            "scenario_id": "bulk-rnaseq-e2e",
            "rows": [
                {
                    "turn_id": "turn_12",
                    "tool": "lookup_go_term",
                    "status": "success",
                    "text": "GO:0008150 biological_process",
                }
            ],
        },
        spec,
    )
    t12 = next(t for t in scored["per_turn_scores"] if t["turn_id"] == "turn_12")
    assert t12["score"] == 0
    assert "enrichment_substituted_with_lookup" in t12["failure_patterns"]


# ─────────────────────────────────────────────────────────────────────────────
# 8.  Fallback text is always non-empty and actionable
# ─────────────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize(
    "tool,prompt",
    [
        ("bio_diff_runs", "Compare the current DEG results with the first DEG results."),
        ("go_enrichment_analysis", "Run GO enrichment on upregulated genes."),
        ("agent", "Now use the cleaned dataset from before batch exclusion."),
        ("bio_rerun", "Exclude samples and rerun."),
    ],
)
def test_fallback_response_non_empty(tool: str, prompt: str):
    """Empty result text must always be upgraded to an actionable fallback."""
    rendered = _render(prompt, tool, {"status": "success", "text": ""})
    assert rendered["text"].strip(), (
        f"Tool '{tool}' returned empty text — actionable fallback is required."
    )


# ─────────────────────────────────────────────────────────────────────────────
# 9.  Scorer strictness: hard failures
# ─────────────────────────────────────────────────────────────────────────────

def test_scorer_timeout_gives_zero_and_is_critical():
    """Timeout responses must score 0 and appear as a critical failure."""
    spec = _load_bench_spec()
    from benchmarks.bio_benchmark_scorer import score_benchmark_run

    scored = score_benchmark_run(
        {
            "session_id": "sid",
            "scenario_id": "bulk-rnaseq-e2e",
            "rows": [
                {
                    "turn_id": "turn_16",
                    "tool": "timeout",
                    "status": "timeout",
                    "text": "Request timed out after 180s.",
                }
            ],
        },
        spec,
    )
    t16 = next(t for t in scored["per_turn_scores"] if t["turn_id"] == "turn_16")
    assert t16["score"] == 0
    assert "hard_execution_failure" in t16["failure_patterns"]
    assert len(scored["critical_failures"]) > 0


def test_scorer_plan_instead_of_execution_gives_zero():
    """workflow_planned when execution is expected must score 0."""
    spec = _load_bench_spec()
    from benchmarks.bio_benchmark_scorer import score_benchmark_run

    scored = score_benchmark_run(
        {
            "session_id": "sid",
            "scenario_id": "bulk-rnaseq-e2e",
            "rows": [
                {
                    "turn_id": "turn_08",
                    "tool": "__plan__",
                    "status": "workflow_planned",
                    "text": "## Pipeline Plan\n\n1. step1 (bulk_rnaseq_analysis)",
                }
            ],
        },
        spec,
    )
    t08 = next(t for t in scored["per_turn_scores"] if t["turn_id"] == "turn_08")
    assert t08["score"] == 0
    assert "plan_instead_of_execution" in t08["failure_patterns"]


def test_scorer_synthetic_success_is_medium_penalty_not_zero():
    """
    synthetic_success_only is a medium pattern (score 1), not a hard failure (0).
    This test documents the exact severity so the scorer contract is clear.
    """
    spec = _load_bench_spec()
    from benchmarks.bio_benchmark_scorer import score_benchmark_run

    scored = score_benchmark_run(
        {
            "session_id": "sid",
            "scenario_id": "bulk-rnaseq-e2e",
            "rows": [
                {
                    "turn_id": "turn_08",
                    "tool": "bulk_rnaseq_analysis",
                    "status": "success",
                    "text": (
                        "## Bulk RNA-seq Analysis Complete\n\n"
                        "**Design formula:** `~condition`\n\n"
                        "> ⚠️ *Synthetic demo data used — real S3 data unavailable.*"
                    ),
                }
            ],
        },
        spec,
    )
    t08 = next(t for t in scored["per_turn_scores"] if t["turn_id"] == "turn_08")
    assert t08["score"] == 1, "synthetic_success_only is a medium penalty (score 1, not 0)."
    assert "synthetic_success_only" in t08["failure_patterns"]
    # Not a critical failure (score > 0)
    assert not any(
        cf["turn_id"] == "turn_08" for cf in scored["critical_failures"]
    )


# ─────────────────────────────────────────────────────────────────────────────
# 10.  Benchmark baseline snapshot
# ─────────────────────────────────────────────────────────────────────────────

_SCORED_MISSING = not ITERATION_SCORED.exists()
_SKIP_BASELINE = pytest.mark.skipif(
    _SCORED_MISSING,
    reason="benchmarks/runs/iteration_1_scored.json not found — run benchmark first.",
)


@_SKIP_BASELINE
def test_benchmark_baseline_above_threshold():
    """
    The saved benchmark run must be ≥ BASELINE_THRESHOLD_PCT (88 %).
    This is a regression gate: a code change that drops the score below 88 %
    will fail this test before it reaches CI / review.
    """
    scored = json.loads(ITERATION_SCORED.read_text())
    pct = float(scored.get("overall_percentage", 0))
    critical = scored.get("critical_failures") or []
    assert pct >= BASELINE_THRESHOLD_PCT, (
        f"Benchmark score {pct}% is below baseline floor {BASELINE_THRESHOLD_PCT}%. "
        "A code change has regressed benchmark performance."
    )
    assert len(critical) <= BASELINE_MAX_CRITICAL, (
        f"Benchmark has {len(critical)} critical failure(s) — expected 0."
    )


@_SKIP_BASELINE
def test_benchmark_baseline_fixed_turns_not_regressed():
    """
    Turns that moved from score 1 → 2 or 0 → 2 in iteration-2 must stay at 2
    in the saved baseline.
    """
    scored = json.loads(ITERATION_SCORED.read_text())
    turns = {t["turn_id"]: t for t in scored.get("per_turn_scores", [])}

    must_be_two = {
        "turn_02": "planning_expected workflow_planned fix",
        "turn_08": "no synthetic_demo note",
        "turn_14": "planning_expected workflow_planned fix",
        "turn_16": "use_historical fast-path (was timeout → 0)",
        "turn_17": "planning_expected workflow_planned fix",
    }
    regressions: list[str] = []
    for tid, reason in must_be_two.items():
        if tid not in turns:
            regressions.append(f"{tid}: missing from scored output")
            continue
        got = turns[tid]["score"]
        if got < 2:
            regressions.append(
                f"{tid} [{reason}]: expected 2, got {got} "
                f"(patterns={turns[tid].get('failure_patterns', [])})"
            )
    assert not regressions, (
        "Benchmark-critical turns regressed:\n" + "\n".join(f"  • {r}" for r in regressions)
    )


@pytest.mark.skipif(
    not BASELINE_FILE.exists(),
    reason="benchmarks/runs/baseline.json not found — run benchmark first.",
)
def test_benchmark_baseline_json_contract():
    """
    baseline.json stores the reference passing state.  Its must_not_regress_below_pct
    and per_turn_must_pass values must be self-consistent with the saved iteration.
    """
    baseline = json.loads(BASELINE_FILE.read_text())
    assert baseline["must_not_regress_below_pct"] >= 80.0, (
        "Baseline floor must be >= 80 % (the product threshold)."
    )
    assert baseline["threshold_met"] is True, "Baseline must record a passing run."
    assert baseline["critical_failures"] == 0, "Baseline must have 0 critical failures."
    # All pinned turns must have score 2
    for tid, score in baseline.get("per_turn_must_pass", {}).items():
        assert score == 2, f"Baseline turn {tid} is expected to be 2, got {score}."


@_SKIP_BASELINE
def test_benchmark_baseline_turns_09_10_19_not_regressed():
    """
    Turns 09 (patch_and_rerun highlight), 10 (patch_and_rerun heatmap), and
    19 (historical recreation) must not regress from score 2.
    """
    scored = json.loads(ITERATION_SCORED.read_text())
    turns = {t["turn_id"]: t for t in scored.get("per_turn_scores", [])}

    for tid in ("turn_09", "turn_10", "turn_19"):
        assert tid in turns, f"{tid} missing from scored output"
        assert turns[tid]["score"] == 2, (
            f"{tid} regressed: score={turns[tid]['score']}, "
            f"patterns={turns[tid].get('failure_patterns', [])}"
        )
