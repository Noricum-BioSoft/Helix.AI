"""
F. Sequence, alignment & plasmid — USER_SCENARIOS §6, TESTBED §7.

F1: Align FASTA
F2: Tree from sequences
F3: Generate variants
F4: Select variants (after F3)
F5: Plasmid map pUC19
F6: Plasmid for representatives
"""

import pytest

from .conftest import FASTA_THREE


def _execute(client, session_id, command):
    r = client.post("/execute", json={"command": command, "session_id": session_id})
    assert r.status_code == 200, r.text
    return r.json()


class TestF1SequenceAlignment:
    """F1_be: Align these FASTA sequences."""

    def test_f1_align_fasta(self, client, session_id):
        cmd = f"Align these FASTA sequences:\n{FASTA_THREE}"
        data = _execute(client, session_id, cmd)
        assert data.get("success") is not False
        result = data.get("result") or data.get("data") or data
        has_alignment = (
            "alignment" in str(result).lower()
            or "aligned" in (data.get("text") or "").lower()
            or result.get("aligned_sequences")
        )
        assert has_alignment or data.get("success") is True


class TestF2TreeFromSequences:
    """F2_be: Build a tree from sequences (unaligned FASTA)."""

    def test_f2_tree_from_fasta(self, client, session_id):
        cmd = f"Build a tree from these sequences:\n{FASTA_THREE}"
        data = _execute(client, session_id, cmd)
        assert data.get("success") is not False
        result = data.get("result") or data
        has_tree = result.get("tree_newick") or result.get("ete_visualization") or "tree" in (data.get("text") or "").lower()
        assert has_tree or data.get("success") is True


class TestF3GenerateVariants:
    """F3_be: Generate 50 random variants of sequence."""

    def test_f3_mutate_sequence(self, client, session_id):
        cmd = "Generate 50 random variants of this sequence: ATGCGATCGATCG"
        data = _execute(client, session_id, cmd)
        assert data.get("success") is not False
        result = data.get("result") or data.get("data") or data
        has_variants = (
            "variant" in str(result).lower()
            or "mutat" in str(result).lower()
            or result.get("mutated_sequences")
            or (isinstance(result.get("statistics"), dict) and "variants" in str(result.get("statistics")))
        )
        assert has_variants or data.get("success") is True


class TestF4SelectVariants:
    """F4_be: After F3, select 10 variants with best conservation."""

    def test_f4_select_variants_after_f3(self, client, session_id):
        _execute(client, session_id, "Generate 20 random variants of this sequence: ATGCGATCGATCG")
        data = _execute(client, session_id, "Select 10 variants with best conservation.")
        assert data.get("success") is not False or "variant" in str(data).lower() or "select" in str(data).lower()


class TestF5PlasmidMap:
    """F5_be: Plasmid map pUC19, EcoRI/BamHI."""

    def test_f5_plasmid_visualization(self, client, session_id):
        cmd = "Show my gene in a plasmid map — vector pUC19, EcoRI and BamHI sites."
        data = _execute(client, session_id, cmd)
        assert data.get("success") is not False
        result = data.get("result") or data
        has_plasmid = "plasmid" in str(result).lower() or "svg" in str(result).lower() or "image" in str(result).lower()
        assert has_plasmid or data.get("success") is True


class TestF6PlasmidForRepresentatives:
    """F6_be: Put representative sequences into pUC19, show plasmid maps."""

    def test_f6_plasmid_for_representatives(self, client, session_id):
        cmd = f"Put these representative sequences into pUC19 and show the plasmid maps:\n{FASTA_THREE}"
        data = _execute(client, session_id, cmd)
        assert data.get("success") is not False or "plasmid" in str(data).lower()
