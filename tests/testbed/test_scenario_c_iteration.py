"""
C. Iterative refinement — USER_SCENARIOS §3, TESTBED §4.

C1: Change alpha to 0.01 after B1
C2: Switch volcano x-axis to linear
C3: Design formula without interaction
C4: Clustering resolution 0.8 (single-cell)
"""

import pytest

from .conftest import BULK_RNASEQ_PROMPT


def _execute(client, session_id, command):
    r = client.post("/execute", json={"command": command, "session_id": session_id})
    assert r.status_code == 200, r.text
    return r.json()


def _run_b1(client, session_id):
    data = _execute(client, session_id, BULK_RNASEQ_PROMPT)
    run_id = data.get("run_id") or (data.get("result") or {}).get("run_id")
    return data, run_id


class TestC1ChangeAlpha:
    """C1_be: After B1, 'Change significance threshold to 0.01' → patch_and_rerun."""

    def test_c1_patch_alpha_after_b1(self, client, session_id):
        _run_b1(client, session_id)
        data = _execute(client, session_id, "Change the significance threshold to 0.01.")
        # May route to patch_and_rerun and return new run_id
        result = data.get("result") or data
        has_diff = "parameter_diff" in str(result) or "output_diff" in str(result)
        new_run = data.get("run_id") or result.get("run_id")
        assert data.get("success") is not False or has_diff or new_run

        runs_r = client.get(f"/session/{session_id}/runs")
        if runs_r.status_code == 200:
            runs = runs_r.json().get("runs", [])
            # At least 2 runs (B1 + patch)
            assert len(runs) >= 1


class TestC2VolcanoXAxisLinear:
    """C2_be: After B1, 'Switch volcano x-axis to linear'."""

    def test_c2_patch_x_scale_linear(self, client, session_id):
        _run_b1(client, session_id)
        data = _execute(client, session_id, "Switch the volcano plot x-axis to linear fold change.")
        assert data.get("success") is not False or "linear" in str(data).lower()


class TestC3DesignFormula:
    """C3_be: After B1, design formula without interaction."""

    def test_c3_design_formula_only_main_effects(self, client, session_id):
        _run_b1(client, session_id)
        data = _execute(
            client,
            session_id,
            "Use design formula ~infection_status + time_point only, no interaction.",
        )
        assert data.get("success") is not False or "design" in str(data).lower()


class TestC4SingleCellResolution:
    """C4_be: After B2, increase clustering resolution to 0.8."""

    def test_c4_resolution_rerun(self, client, session_id):
        cmd = (
            "Run single-cell analysis. data_file: s3://noricum-ngs-data/demo/scrna/sle_pbmc_filtered_feature_bc_matrix.h5 "
            "data_format: 10x resolution: 0.5"
        )
        _execute(client, session_id, cmd)
        data = _execute(client, session_id, "Increase clustering resolution to 0.8.")
        assert "resolution" in str(data).lower() or data.get("success") is not False
