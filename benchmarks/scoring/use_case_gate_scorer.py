from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Tuple


@dataclass
class GateThresholds:
    min_total_score: float = 0.80
    require_no_forbidden_route: bool = True
    require_no_forbidden_fallback: bool = True
    require_provenance_complete: bool = True
    require_replay_pass: bool = True


def _to_set(values: Any) -> set[str]:
    if not values:
        return set()
    if isinstance(values, list):
        return {str(v) for v in values}
    return {str(values)}


def _check_routing(case_def: Dict[str, Any], run_obs: Dict[str, Any]) -> Tuple[bool, List[str]]:
    reasons: List[str] = []
    expected_route = str(case_def.get("expected_route") or "").strip()
    observed_route = str(run_obs.get("observed_route") or "").strip()
    forbidden_routes = _to_set(case_def.get("forbidden_routes"))

    if expected_route and observed_route != expected_route:
        reasons.append(f"observed_route '{observed_route}' != expected_route '{expected_route}'")

    if observed_route in forbidden_routes:
        reasons.append(f"observed_route '{observed_route}' is forbidden")

    ok = len(reasons) == 0
    return ok, reasons


def _check_fallbacks(case_def: Dict[str, Any], run_obs: Dict[str, Any]) -> Tuple[bool, List[str]]:
    reasons: List[str] = []
    forbidden = _to_set(case_def.get("forbidden_fallbacks"))
    observed = _to_set(run_obs.get("fallback_targets"))

    bad = sorted(observed.intersection(forbidden))
    if bad:
        reasons.append("forbidden_fallbacks_invoked=" + ",".join(bad))

    return len(reasons) == 0, reasons


def _check_required_steps(case_def: Dict[str, Any], run_obs: Dict[str, Any]) -> Tuple[bool, List[str]]:
    reasons: List[str] = []
    required = _to_set(case_def.get("required_steps"))
    observed = _to_set(run_obs.get("observed_steps"))

    missing = sorted(required.difference(observed))
    if missing:
        reasons.append("missing_required_steps=" + ",".join(missing))

    return len(reasons) == 0, reasons


def _check_required_outputs(case_def: Dict[str, Any], run_obs: Dict[str, Any]) -> Tuple[bool, List[str]]:
    reasons: List[str] = []
    required = _to_set(case_def.get("required_outputs"))
    observed = _to_set(run_obs.get("observed_output_types"))

    missing = sorted(required.difference(observed))
    if missing:
        reasons.append("missing_required_outputs=" + ",".join(missing))

    return len(reasons) == 0, reasons


def _check_safety(case_def: Dict[str, Any], run_obs: Dict[str, Any]) -> Tuple[bool, List[str]]:
    reasons: List[str] = []
    required = _to_set(case_def.get("safety_requirements"))
    observed = _to_set(run_obs.get("observed_safety_flags"))

    missing = sorted(required.difference(observed))
    if missing:
        reasons.append("missing_safety_requirements=" + ",".join(missing))

    return len(reasons) == 0, reasons


def _check_provenance(case_def: Dict[str, Any], run_obs: Dict[str, Any]) -> Tuple[bool, List[str]]:
    reasons: List[str] = []
    required = _to_set(case_def.get("provenance_requirements"))
    observed = _to_set(run_obs.get("observed_provenance_signals"))

    missing = sorted(required.difference(observed))
    if missing:
        reasons.append("missing_provenance_requirements=" + ",".join(missing))

    return len(reasons) == 0, reasons


def _check_replay(run_obs: Dict[str, Any]) -> Tuple[bool, List[str]]:
    reasons: List[str] = []
    replay_status = str(run_obs.get("replay_status") or "")
    if replay_status not in {"replay_match", "replay_within_tolerance"}:
        reasons.append(f"replay_status='{replay_status}'")
    return len(reasons) == 0, reasons


def score_case_gates(
    case_def: Dict[str, Any],
    run_obs: Dict[str, Any],
    *,
    thresholds: GateThresholds | None = None,
) -> Dict[str, Any]:
    cfg = thresholds or GateThresholds()

    checks = {
        "routing": _check_routing(case_def, run_obs),
        "fallbacks": _check_fallbacks(case_def, run_obs),
        "required_steps": _check_required_steps(case_def, run_obs),
        "required_outputs": _check_required_outputs(case_def, run_obs),
        "safety": _check_safety(case_def, run_obs),
        "provenance": _check_provenance(case_def, run_obs),
        "replay": _check_replay(run_obs),
    }

    failed_checks: Dict[str, List[str]] = {
        name: reasons for name, (ok, reasons) in checks.items() if not ok
    }
    passed = len(failed_checks) == 0

    # Weighted score for non-binary quality trend tracking.
    # Hard-gate conditions are still enforced independently via thresholds.
    weight = {
        "routing": 0.25,
        "fallbacks": 0.20,
        "required_steps": 0.15,
        "required_outputs": 0.10,
        "safety": 0.10,
        "provenance": 0.10,
        "replay": 0.10,
    }
    total_score = 0.0
    for name, (ok, _) in checks.items():
        total_score += weight[name] if ok else 0.0

    gate_reasons: List[str] = []
    if cfg.require_no_forbidden_route and "routing" in failed_checks:
        gate_reasons.append("forbidden_or_wrong_route")
    if cfg.require_no_forbidden_fallback and "fallbacks" in failed_checks:
        gate_reasons.append("forbidden_fallback")
    if cfg.require_provenance_complete and "provenance" in failed_checks:
        gate_reasons.append("provenance_incomplete")
    if cfg.require_replay_pass and "replay" in failed_checks:
        gate_reasons.append("replay_failed")
    if total_score < cfg.min_total_score:
        gate_reasons.append("below_min_total_score")

    gate_pass = len(gate_reasons) == 0

    return {
        "case_id": case_def.get("id"),
        "passed_all_checks": passed,
        "total_score": round(total_score, 3),
        "checks": {k: {"ok": v[0], "reasons": v[1]} for k, v in checks.items()},
        "failed_checks": failed_checks,
        "gate_pass": gate_pass,
        "gate_fail_reasons": gate_reasons,
        "thresholds": {
            "min_total_score": cfg.min_total_score,
            "require_no_forbidden_route": cfg.require_no_forbidden_route,
            "require_no_forbidden_fallback": cfg.require_no_forbidden_fallback,
            "require_provenance_complete": cfg.require_provenance_complete,
            "require_replay_pass": cfg.require_replay_pass,
        },
    }

