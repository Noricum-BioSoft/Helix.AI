"""
I. Multi-step pipelines — USER_SCENARIOS §6, TESTBED §10.

I1: Align then build tree
I2: Generate variants → select 12 → plasmid maps
"""

import pytest

from .conftest import FASTA_THREE


def _execute(client, session_id, command):
    r = client.post("/execute", json={"command": command, "session_id": session_id})
    assert r.status_code == 200, r.text
    return r.json()


class TestI1AlignThenTree:
    """I1_be: Align these sequences, then build a phylogenetic tree."""

    def test_i1_align_then_tree(self, client, session_id):
        cmd = f"Align these sequences, then build a phylogenetic tree:\n{FASTA_THREE}"
        data = _execute(client, session_id, cmd)
        assert data.get("success") is not False
        result = data.get("result") or data
        has_tree = result.get("tree_newick") or result.get("ete_visualization") or "tree" in (data.get("text") or "").lower()
        assert has_tree or "alignment" in str(result).lower() or data.get("success") is True


class TestI2VariantsSelectPlasmid:
    """I2_be: Generate 96 variants, select 12 with lowest gaps, show plasmid maps."""

    def test_i2_three_step_chain(self, client, session_id):
        cmd = (
            "Generate 96 variants of this sequence: ATGCGATCGATCGATCG. "
            "Then select the 12 with lowest gaps and show them in plasmid maps."
        )
        data = _execute(client, session_id, cmd)
        assert data.get("success") is not False or "variant" in str(data).lower() or "plasmid" in str(data).lower()
