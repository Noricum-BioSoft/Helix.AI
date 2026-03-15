from __future__ import annotations

import re
from typing import Any, Dict, List, Optional


def _version_from_title(title: str) -> Optional[int]:
    if not title:
        return None
    m = re.search(r"(?:^|[_\s-])v(\d+)(?:$|[_\s-])", title.lower())
    return int(m.group(1)) if m else None


def build_alias_index(session: Dict[str, Any]) -> Dict[str, Any]:
    runs = session.get("runs", []) if isinstance(session, dict) else []
    artifacts_raw = session.get("artifacts", {}) if isinstance(session, dict) else {}
    artifacts: List[Dict[str, Any]]
    if isinstance(artifacts_raw, dict):
        artifacts = [a for a in artifacts_raw.values() if isinstance(a, dict)]
    elif isinstance(artifacts_raw, list):
        artifacts = [a for a in artifacts_raw if isinstance(a, dict)]
    else:
        artifacts = []

    aliases: Dict[str, Dict[str, Any]] = {}
    deg_runs = [r for r in runs if isinstance(r, dict) and r.get("tool") == "bulk_rnaseq_analysis"]
    if deg_runs:
        aliases["first_deg_results"] = {"run_id": deg_runs[0].get("run_id"), "source": "runs"}
        aliases["latest_deg_results"] = {"run_id": deg_runs[-1].get("run_id"), "source": "runs"}

    # Pathway-ish alias fallback via GO lookup or enrichment-style titles.
    pathway_arts = [
        a for a in artifacts
        if "pathway" in str(a.get("title", "")).lower()
        or "enrichment" in str(a.get("title", "")).lower()
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

    # Create title-derived aliases (normalized snake keys).
    for art in artifacts:
        title = str(art.get("title") or "").strip()
        if not title:
            continue
        key = re.sub(r"[^a-z0-9]+", "_", title.lower()).strip("_")
        if key and key not in aliases:
            aliases[key] = {"artifact_id": art.get("artifact_id"), "source": "artifacts"}

    # Historical states by branch + version.
    states: List[Dict[str, Any]] = []
    for art in artifacts:
        title = str(art.get("title") or "")
        states.append(
            {
                "artifact_id": art.get("artifact_id"),
                "title": title,
                "artifact_kind": art.get("artifact_kind") or art.get("type"),
                "version": art.get("version") or _version_from_title(title),
                "analysis_branch": art.get("analysis_branch") or "main",
                "state_tags": art.get("state_tags") or [],
                "source_run_id": art.get("source_run_id") or art.get("run_id"),
            }
        )

    return {"aliases": aliases, "historical_states": states}


def _tokenize(text: str) -> List[str]:
    return [t for t in re.findall(r"[a-z0-9]+", (text or "").lower()) if len(t) > 2]


def _resolve_historical_state(states: List[Dict[str, Any]], phrase: str) -> Optional[Dict[str, Any]]:
    phrase_l = (phrase or "").strip().lower()
    if not phrase_l or not states:
        return None

    stop = {"the", "and", "from", "with", "that", "this", "version", "dataset", "results"}
    parts = phrase_l.split("before", 1)
    include_tokens = [t for t in _tokenize(parts[0]) if t not in stop]
    exclude_tokens = [t for t in _tokenize(parts[1]) if t not in stop] if len(parts) > 1 else []
    ranked: List[Dict[str, Any]] = []
    for state in states:
        title = str(state.get("title") or "")
        tags = " ".join(str(x) for x in (state.get("state_tags") or []))
        branch = str(state.get("analysis_branch") or "")
        haystack = f"{title} {tags} {branch}".lower()
        score = 0
        if include_tokens:
            score += sum(1 for t in include_tokens if t in haystack)
        if exclude_tokens:
            score -= 2 * sum(1 for t in exclude_tokens if t in haystack)
        if score > 0:
            ranked.append({"score": score, "state": state})

    if not ranked:
        return None
    ranked.sort(key=lambda x: x["score"], reverse=True)
    top = ranked[0]
    if len(ranked) > 1 and ranked[1]["score"] == top["score"]:
        return None
    st = top["state"]
    return {
        "artifact_id": st.get("artifact_id"),
        "run_id": st.get("source_run_id"),
        "source": "historical_states",
        "analysis_branch": st.get("analysis_branch"),
        "version": st.get("version"),
    }


def resolve_semantic_reference(session: Dict[str, Any], phrase: str) -> Dict[str, Any]:
    phrase_l = (phrase or "").strip().lower()
    idx = build_alias_index(session)
    aliases = idx.get("aliases", {})

    # Direct alias key hit.
    if phrase_l in aliases:
        return {"status": "resolved", "selector": phrase, "target": aliases[phrase_l]}

    # Common normalized references.
    mapping = {
        "first deg results": "first_deg_results",
        "current deg results": "latest_deg_results",
        "latest deg results": "latest_deg_results",
        "first pathway": "first_pathway_analysis",
        "first pathway analysis": "first_pathway_analysis",
        "latest pathway analysis": "latest_pathway_analysis",
    }
    key = mapping.get(phrase_l)
    if key and key in aliases:
        return {"status": "resolved", "selector": phrase, "target": aliases[key]}

    # Semantic historical-state lookup for references like:
    # "cleaned dataset from before batch exclusion" or
    # "corrected metadata version before the fold-change bug fix".
    hist_target = _resolve_historical_state(idx.get("historical_states") or [], phrase_l)
    if hist_target:
        return {"status": "resolved", "selector": phrase, "target": hist_target}

    return {
        "status": "unresolved",
        "selector": phrase,
        "message": "No unique artifact/run matched this semantic reference in the current session.",
    }

