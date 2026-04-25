from __future__ import annotations

from unittest.mock import patch

from backend.command_router import CommandRouter


def test_sequence_alignment_fallback_requires_sequence_cues() -> None:
    """Commands without actual nucleotide sequences should not route to sequence_alignment."""
    router = CommandRouter()
    tool, _params = router.route_command(
        "Use sequences from the previous analysis and summarize differences.",
        {},
    )
    assert tool != "sequence_alignment"


def test_sequence_alignment_fallback_allows_true_sequence_text() -> None:
    """Commands containing actual DNA/RNA sequences should route to sequence_alignment."""
    router = CommandRouter()
    tool, _params = router.route_command(
        "align ACGTACGT and ACGTTCGT to check mismatch positions",
        {},
    )
    assert tool == "sequence_alignment"


def test_route_command_with_shadow_attaches_metadata(monkeypatch) -> None:
    """route_command_with_shadow should always attach router_shadow metadata."""
    monkeypatch.setenv("HELIX_ROUTER_SHADOW_MODE", "1")

    router = CommandRouter()

    # Primary routing: _route_with_llm returns toolbox_inventory
    # Shadow call in route_command_with_shadow: also returns toolbox_inventory (same LLM)
    # → disagrees is False (both agree)
    with patch.object(
        CommandRouter,
        "_route_with_llm",
        return_value=("toolbox_inventory", {}),
    ):
        tool, params = router.route_command_with_shadow("list tools", {"session_id": "s1"})

    assert tool == "toolbox_inventory"
    assert "router_shadow" in params
    assert params["router_shadow"]["selected_tool"] == "toolbox_inventory"
    assert params["router_shadow"]["shadow_tool"] == "toolbox_inventory"
    assert params["router_shadow"]["disagrees"] is False


def test_route_command_with_shadow_detects_disagreement(monkeypatch) -> None:
    """shadow_tool differs from primary tool when LLM returns different results on repeated calls."""
    monkeypatch.setenv("HELIX_ROUTER_SHADOW_MODE", "1")

    router = CommandRouter()
    _calls = []

    def _alternating(self_or_cmd, cmd_or_ctx=None, ctx=None):
        # Handle both bound and unbound call signatures
        if cmd_or_ctx is None:
            return ("toolbox_inventory", {}) if not _calls else ("bulk_rnaseq_analysis", {})
        _calls.append(1)
        return ("bulk_rnaseq_analysis", {"session_id": "s1"})

    # Patch to return toolbox_inventory for primary, bulk_rnaseq_analysis for shadow
    with patch.object(
        CommandRouter,
        "_route_with_llm",
        side_effect=[
            ("toolbox_inventory", {}),      # primary call inside route_command
            ("bulk_rnaseq_analysis", {}),   # shadow call inside route_command_with_shadow
        ],
    ):
        tool, params = router.route_command_with_shadow("list tools", {"session_id": "s1"})

    assert tool == "toolbox_inventory"
    assert "router_shadow" in params
    assert params["router_shadow"]["shadow_tool"] == "bulk_rnaseq_analysis"
    assert params["router_shadow"]["disagrees"] is True
