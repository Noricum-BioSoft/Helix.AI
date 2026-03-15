from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple


def _version_from_title(title: str) -> Optional[int]:
    if not title:
        return None
    m = re.search(r"(?:^|[_\s-])v(\d+)(?:$|[_\s-])", title.lower())
    return int(m.group(1)) if m else None


def _normalize_key(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", (text or "").lower()).strip("_")


def _tokenize(text: str) -> List[str]:
    return [t for t in re.findall(r"[a-z0-9]+", (text or "").lower()) if len(t) > 2]


def _collect_artifacts(session: Dict[str, Any]) -> List[Dict[str, Any]]:
    artifacts_raw = session.get("artifacts", {}) if isinstance(session, dict) else {}
    if isinstance(artifacts_raw, dict):
        return [a for a in artifacts_raw.values() if isinstance(a, dict)]
    if isinstance(artifacts_raw, list):
        return [a for a in artifacts_raw if isinstance(a, dict)]
    return []


def _kind_tokens(kind: str) -> List[str]:
    return _tokenize(kind.replace("_", " "))


def build_alias_index(session: Dict[str, Any]) -> Dict[str, Any]:
    runs = session.get("runs", []) if isinstance(session, dict) else []
    artifacts = _collect_artifacts(session)

    aliases: Dict[str, Dict[str, Any]] = {}
    bulk_runs = [r for r in runs if isinstance(r, dict) and r.get("tool") == "bulk_rnaseq_analysis"]
    if bulk_runs:
        aliases["first_deg_results"] = {"run_id": bulk_runs[0].get("run_id"), "source": "runs"}
        aliases["latest_deg_results"] = {"run_id": bulk_runs[-1].get("run_id"), "source": "runs"}
        aliases["current_deg_results"] = {"run_id": bulk_runs[-1].get("run_id"), "source": "runs"}
        if len(bulk_runs) > 1:
            aliases["prior_deg_results"] = {"run_id": bulk_runs[-2].get("run_id"), "source": "runs"}
            aliases["original_deg_results"] = {"run_id": bulk_runs[0].get("run_id"), "source": "runs"}

    pathway_arts = [
        a
        for a in artifacts
        if "pathway" in str(a.get("title", "")).lower() or "enrichment" in str(a.get("title", "")).lower()
    ]
    if pathway_arts:
        aliases["first_pathway_analysis"] = {
            "artifact_id": pathway_arts[0].get("artifact_id"),
            "source": "artifacts",
        }
        aliases["latest_pathway_analysis"] = {
            "artifact_id": pathway_arts[-1].get("artifact_id"),
            "source": "artifacts",
        }

    pca_arts = [a for a in artifacts if "pca" in str(a.get("title", "")).lower()]
    if pca_arts:
        aliases["latest_pca"] = {"artifact_id": pca_arts[-1].get("artifact_id"), "source": "artifacts"}
        aliases["earlier_pca"] = {
            "artifact_id": (pca_arts[-2] if len(pca_arts) > 1 else pca_arts[0]).get("artifact_id"),
            "source": "artifacts",
        }

    filtered_arts = [
        a for a in artifacts if "filter" in str(a.get("title", "")).lower() or "clean" in str(a.get("title", "")).lower()
    ]
    if filtered_arts:
        aliases["prior_filtered_dataset"] = {
            "artifact_id": (filtered_arts[-2] if len(filtered_arts) > 1 else filtered_arts[0]).get("artifact_id"),
            "source": "artifacts",
        }

    for art in artifacts:
        title = str(art.get("title") or "").strip()
        artifact_id = art.get("artifact_id")
        if artifact_id and title:
            key = _normalize_key(title)
            if key and key not in aliases:
                aliases[key] = {"artifact_id": artifact_id, "source": "artifacts"}
        semantic_aliases = art.get("semantic_aliases") if isinstance(art.get("semantic_aliases"), list) else []
        for alias in semantic_aliases:
            alias_key = _normalize_key(str(alias))
            if alias_key and alias_key not in aliases and artifact_id:
                aliases[alias_key] = {"artifact_id": artifact_id, "source": "semantic_alias"}

    states: List[Dict[str, Any]] = []
    for idx, art in enumerate(artifacts):
        title = str(art.get("title") or "")
        artifact_kind = str(art.get("artifact_kind") or art.get("type") or "")
        alias_tokens = []
        if isinstance(art.get("semantic_aliases"), list):
            alias_tokens = [str(x).strip().lower() for x in art.get("semantic_aliases") if str(x).strip()]
        states.append(
            {
                "artifact_id": art.get("artifact_id"),
                "title": title,
                "artifact_kind": artifact_kind,
                "version": art.get("version") or _version_from_title(title),
                "analysis_branch": art.get("analysis_branch") or "main",
                "state_tags": art.get("state_tags") or [],
                "source_run_id": art.get("source_run_id") or art.get("run_id"),
                "parent_artifact_ids": art.get("parent_artifact_ids") or [],
                "semantic_aliases": alias_tokens,
                "created_at": art.get("created_at"),
                "ordinal": idx,
            }
        )
    return {"aliases": aliases, "historical_states": states}


def _temporal_rank(state: Dict[str, Any]) -> Tuple[int, int]:
    version = state.get("version")
    v = int(version) if isinstance(version, int) else -1
    return (v, int(state.get("ordinal") or 0))


def _resolve_historical_state(states: List[Dict[str, Any]], phrase: str) -> Dict[str, Any]:
    phrase_l = (phrase or "").strip().lower()
    if not phrase_l or not states:
        return {"status": "unresolved"}

    stop = {"the", "and", "from", "with", "that", "this", "version", "dataset", "results"}
    parts = phrase_l.split("before", 1)
    include_tokens = [t for t in _tokenize(parts[0]) if t not in stop]
    exclude_tokens = [t for t in _tokenize(parts[1]) if t not in stop] if len(parts) > 1 else []
    ranked: List[Dict[str, Any]] = []
    for state in states:
        title = str(state.get("title") or "")
        tags = " ".join(str(x) for x in (state.get("state_tags") or []))
        branch = str(state.get("analysis_branch") or "")
        kind = str(state.get("artifact_kind") or "")
        aliases = " ".join(str(x) for x in (state.get("semantic_aliases") or []))
        haystack = f"{title} {tags} {branch} {kind} {aliases}".lower()
        score = 0
        score += sum(1 for t in include_tokens if t in haystack)
        score -= sum(1 for t in exclude_tokens if t in haystack)
        score += sum(1 for t in _kind_tokens(kind) if t in include_tokens)
        if "corrected" in phrase_l and "correct" in haystack:
            score += 1
        if ("batch exclusion" in phrase_l or "exclude batch" in phrase_l) and "batch_exclusion" in haystack:
            score -= 2
        if "pre-bugfix" in phrase_l or "before bug fix" in phrase_l or "before the fold-change bug fix" in phrase_l:
            if "pre_bugfix" in haystack or "before_bug_fix" in haystack:
                score += 2
            elif "bug_fix" in haystack or "fold_change_bug_fix" in haystack:
                score -= 2
        if "latest" in phrase_l or "current" in phrase_l:
            v, ord_idx = _temporal_rank(state)
            score += max(v, 0) + ord_idx // 10
        if "first" in phrase_l or "original" in phrase_l or "initial" in phrase_l:
            v, ord_idx = _temporal_rank(state)
            if v >= 0:
                score -= v
            score -= ord_idx // 10
        if score > 0:
            ranked.append({"score": score, "state": state})

    if not ranked:
        return {"status": "unresolved"}
    ranked.sort(key=lambda x: (x["score"], _temporal_rank(x["state"])), reverse=True)
    top = ranked[0]
    top_candidates = [r for r in ranked if r["score"] == top["score"]]
    if len(top_candidates) > 1:
        return {
            "status": "ambiguous",
            "candidates": [
                {
                    "artifact_id": r["state"].get("artifact_id"),
                    "run_id": r["state"].get("source_run_id"),
                    "title": r["state"].get("title"),
                    "analysis_branch": r["state"].get("analysis_branch"),
                    "version": r["state"].get("version"),
                }
                for r in top_candidates[:3]
            ],
        }
    st = top["state"]
    return {
        "status": "resolved",
        "target": {
            "artifact_id": st.get("artifact_id"),
            "run_id": st.get("source_run_id"),
            "source": "historical_states",
            "analysis_branch": st.get("analysis_branch"),
            "version": st.get("version"),
        },
    }


def resolve_semantic_reference(session: Dict[str, Any], phrase: str) -> Dict[str, Any]:
    phrase_l = (phrase or "").strip().lower()
    idx = build_alias_index(session)
    aliases = idx.get("aliases", {})
    normalized_phrase = _normalize_key(phrase_l)
    if phrase_l in aliases:
        return {"status": "resolved", "selector": phrase, "target": aliases[phrase_l], "diagnostics": {"match": "exact_alias"}}
    if normalized_phrase in aliases:
        return {"status": "resolved", "selector": phrase, "target": aliases[normalized_phrase], "diagnostics": {"match": "normalized_alias"}}

    mapping = {
        "first deg results": "first_deg_results",
        "latest deg results": "latest_deg_results",
        "current deg results": "current_deg_results",
        "prior deg results": "prior_deg_results",
        "original deg results": "original_deg_results",
        "first pathway": "first_pathway_analysis",
        "first pathway analysis": "first_pathway_analysis",
        "latest pathway analysis": "latest_pathway_analysis",
        "earlier pca": "earlier_pca",
        "prior filtered dataset": "prior_filtered_dataset",
    }
    key = mapping.get(phrase_l)
    if key and key in aliases:
        return {"status": "resolved", "selector": phrase, "target": aliases[key], "diagnostics": {"match": "mapping"}}

    # Loose alias fallback: semantic alias/title keys may include additional tokens.
    phrase_tokens = set(_tokenize(phrase_l))
    if phrase_tokens:
        alias_hits: List[Tuple[str, Dict[str, Any]]] = []
        for alias_key, target in aliases.items():
            alias_tokens = set(_tokenize(alias_key))
            if phrase_tokens.issubset(alias_tokens) or alias_tokens.issubset(phrase_tokens):
                alias_hits.append((alias_key, target))
        if len(alias_hits) == 1:
            return {
                "status": "resolved",
                "selector": phrase,
                "target": alias_hits[0][1],
                "diagnostics": {"match": "loose_alias", "alias": alias_hits[0][0]},
            }
        if len(alias_hits) > 1:
            return {
                "status": "unresolved",
                "selector": phrase,
                "message": "Semantic reference is ambiguous. Please choose one candidate.",
                "diagnostics": {"type": "ambiguous_alias", "candidates": [a for a, _ in alias_hits[:5]]},
            }

    hist_target = _resolve_historical_state(idx.get("historical_states") or [], phrase_l)
    if hist_target.get("status") == "resolved":
        return {"status": "resolved", "selector": phrase, "target": hist_target.get("target"), "diagnostics": {"match": "historical_state"}}
    if hist_target.get("status") == "ambiguous":
        return {
            "status": "unresolved",
            "selector": phrase,
            "message": "Multiple historical states matched. Please clarify the branch/version.",
            "diagnostics": {"type": "ambiguous_historical_state", "candidates": hist_target.get("candidates")},
        }

    return {
        "status": "unresolved",
        "selector": phrase,
        "message": "No unique artifact/run matched this semantic reference in the current session.",
        "diagnostics": {"type": "no_match"},
    }

