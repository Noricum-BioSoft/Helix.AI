from __future__ import annotations

from backend.command_router import CommandRouter


def test_sequence_alignment_fallback_requires_sequence_cues() -> None:
    router = CommandRouter()
    tool, _params = router.route_command(
        "Use sequences from the previous analysis and summarize differences.",
        {},
    )
    assert tool != "sequence_alignment"


def test_sequence_alignment_fallback_allows_true_sequence_text() -> None:
    router = CommandRouter()
    tool, _params = router.route_command(
        "align ACGTACGT and ACGTTCGT to check mismatch positions",
        {},
    )
    assert tool == "sequence_alignment"


def test_route_command_with_shadow_attaches_disagreement(monkeypatch) -> None:
    monkeypatch.setenv("HELIX_ROUTER_SHADOW_MODE", "1")
    monkeypatch.setenv("HELIX_LLM_ROUTER_FIRST", "0")

    router = CommandRouter()
    monkeypatch.setattr(
        router,
        "_route_with_llm",
        lambda command, session_context: ("bulk_rnaseq_analysis", {"session_id": "s1"}),
    )

    tool, params = router.route_command_with_shadow("list tools", {"session_id": "s1"})
    assert tool == "toolbox_inventory"
    assert "router_shadow" in params
    assert params["router_shadow"]["selected_tool"] == "toolbox_inventory"
    assert params["router_shadow"]["shadow_tool"] == "bulk_rnaseq_analysis"
    assert params["router_shadow"]["disagrees"] is True
