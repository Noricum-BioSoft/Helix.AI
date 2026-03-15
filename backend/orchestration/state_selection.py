from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from backend.orchestration.artifact_resolver import resolve_semantic_reference


def _run_ids_from_session(session: Dict[str, Any]) -> List[str]:
    runs = session.get("runs", []) if isinstance(session, dict) else []
    out: List[str] = []
    for r in runs:
        if isinstance(r, dict) and isinstance(r.get("run_id"), str):
            out.append(r["run_id"])
    return out


def _resolve_run_id_from_selector(session: Dict[str, Any], selector: str) -> Optional[str]:
    resolved = resolve_semantic_reference(session, selector)
    if resolved.get("status") != "resolved":
        return None
    target = resolved.get("target") or {}
    run_id = target.get("run_id")
    if isinstance(run_id, str) and run_id:
        return run_id
    artifact_id = target.get("artifact_id")
    if isinstance(artifact_id, str) and artifact_id:
        artifacts = session.get("artifacts", {}) if isinstance(session, dict) else {}
        if isinstance(artifacts, dict):
            artifact = artifacts.get(artifact_id)
            if isinstance(artifact, dict):
                source_run_id = artifact.get("source_run_id") or artifact.get("run_id")
                if isinstance(source_run_id, str) and source_run_id:
                    return source_run_id
        elif isinstance(artifacts, list):
            for artifact in artifacts:
                if not isinstance(artifact, dict):
                    continue
                if artifact.get("artifact_id") != artifact_id:
                    continue
                source_run_id = artifact.get("source_run_id") or artifact.get("run_id")
                if isinstance(source_run_id, str) and source_run_id:
                    return source_run_id
    return None


def _infer_analysis_kind(command: str) -> str:
    cmd = (command or "").lower()
    if "deg" in cmd or "differential" in cmd or "significance" in cmd:
        return "deg results"
    if "pathway" in cmd or "enrichment" in cmd or "go term" in cmd:
        return "pathway analysis"
    if "pca" in cmd:
        return "pca"
    if any(k in cmd for k in ["figure", "heatmap", "volcano", "plot"]):
        return "figure set"
    return "results"


def _resolve_selector_with_diag(
    session: Dict[str, Any],
    selector: str,
    diag: Dict[str, Any],
    *,
    label: str,
) -> Optional[str]:
    resolved = resolve_semantic_reference(session, selector)
    diag.setdefault("selectors", []).append({"label": label, "selector": selector})
    if resolved.get("status") == "resolved":
        target = resolved.get("target") or {}
        run_id = target.get("run_id")
        if isinstance(run_id, str) and run_id:
            diag.setdefault("resolver_hits", []).append(
                {"label": label, "selector": selector, "run_id": run_id}
            )
            return run_id
    diagnostics = resolved.get("diagnostics")
    if diagnostics:
        diag.setdefault("resolver_misses", []).append(
            {"label": label, "selector": selector, "diagnostics": diagnostics}
        )
    return None


def select_rerun_anchor(session: Dict[str, Any], target_run: str, command: str) -> str:
    run_ids = _run_ids_from_session(session)
    if target_run not in {"latest", "prior"}:
        return target_run
    if not run_ids:
        return target_run

    cmd = (command or "").lower()
    kind = _infer_analysis_kind(command)
    if "current" in cmd or "latest" in cmd:
        rid = _resolve_run_id_from_selector(session, f"current {kind}")
        if rid:
            return rid
    if "prior" in cmd or "earlier" in cmd or "first" in cmd or "original" in cmd:
        rid = _resolve_run_id_from_selector(session, f"prior {kind}")
        if rid:
            return rid
    if target_run == "latest":
        return run_ids[-1]
    if len(run_ids) > 1:
        return run_ids[-2]
    return run_ids[-1]


def select_diff_anchors(
    session: Dict[str, Any],
    run_id_a: str,
    run_id_b: str,
    command: str,
) -> Tuple[str, str, Dict[str, Any]]:
    diag: Dict[str, Any] = {"selectors": [], "resolver_hits": [], "resolver_misses": []}
    a = run_id_a
    b = run_id_b
    cmd = (command or "").lower()
    kind = _infer_analysis_kind(command)

    if a in {"latest", ""}:
        for selector in [f"current {kind}", f"latest {kind}", command]:
            a_try = _resolve_selector_with_diag(session, selector, diag, label="run_a")
            if a_try:
                a = a_try
                break

    if b in {"prior", ""}:
        b_selectors = [f"prior {kind}", f"first {kind}", f"original {kind}"]
        if "before" in cmd:
            b_selectors.insert(0, command)
        for selector in b_selectors:
            b_try = _resolve_selector_with_diag(session, selector, diag, label="run_b")
            if b_try:
                b = b_try
                break

    run_ids = _run_ids_from_session(session)
    if a == "latest" and run_ids:
        a = run_ids[-1]
    if b == "prior" and len(run_ids) > 1:
        b = run_ids[-2]
    elif b == "prior" and run_ids:
        b = run_ids[-1]
    return a, b, diag

