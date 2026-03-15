from pathlib import Path

from backend.orchestration.benchmark_controller import (
    BenchmarkLoopConfig,
    run_benchmark_iteration_loop,
)


def test_benchmark_controller_stops_on_threshold(tmp_path: Path):
    raws = [{"rows": []}, {"rows": []}]
    scores = [
        {"overall_percentage": 70.0, "critical_failures": [{"turn_id": "t1"}]},
        {"overall_percentage": 85.0, "critical_failures": []},
    ]
    state = {"i": 0}

    def _run_once():
        idx = state["i"]
        return raws[idx]

    def _score_once(_raw):
        idx = state["i"]
        state["i"] += 1
        return {
            "session_id": f"s{idx}",
            "scenario_id": "bench",
            "raw_score": 1 + idx,
            "max_score": 2,
            **scores[idx],
        }

    result = run_benchmark_iteration_loop(
        run_once=_run_once,
        score_once=_score_once,
        out_dir=tmp_path,
        config=BenchmarkLoopConfig(max_iters=5, threshold_pct=80.0, stall_iters=2),
    )
    assert result["stop_reason"] == "threshold_met_and_no_critical_failures"
    assert len(result["history"]) == 2
    assert (tmp_path / "iteration_2_report.md").exists()
    assert (tmp_path / "latest_iteration.json").exists()
    assert (tmp_path / "loop_summary.json").exists()


def test_benchmark_controller_stops_on_stall(tmp_path: Path):
    state = {"i": 0}

    def _run_once():
        return {"rows": []}

    def _score_once(_raw):
        state["i"] += 1
        return {
            "session_id": f"s{state['i']}",
            "scenario_id": "bench",
            "raw_score": 1,
            "max_score": 2,
            "overall_percentage": 60.0,
            "critical_failures": [{"turn_id": "t"}],
        }

    result = run_benchmark_iteration_loop(
        run_once=_run_once,
        score_once=_score_once,
        out_dir=tmp_path,
        config=BenchmarkLoopConfig(max_iters=6, threshold_pct=80.0, stall_delta_pct=0.1, stall_iters=2),
    )
    assert result["stop_reason"] == "progress_stalled"
    assert len(result["history"]) >= 3
    assert (tmp_path / "latest_iteration.json").exists()
    assert (tmp_path / "loop_summary.json").exists()

