from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict

import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
POLICY_PATH = REPO_ROOT / "benchmarks" / "performance_cost_policy.yaml"
RELEASE_THRESHOLDS_PATH = REPO_ROOT / "benchmarks" / "release_thresholds.yaml"
ENV_CAPABILITIES_PATH = REPO_ROOT / "backend" / "config" / "environment_capabilities.yaml"
HEURISTICS_PATH = REPO_ROOT / "backend" / "config" / "cost_heuristics.yaml"


def _load_yaml(path: Path) -> Dict[str, Any]:
    return yaml.safe_load(path.read_text()) or {}


def _require_keys(obj: Dict[str, Any], keys: set[str], context: str) -> None:
    missing = sorted(keys.difference(obj.keys()))
    if missing:
        raise ValueError(f"{context} missing keys: {missing}")


def _validate_latency_targets(policy: Dict[str, Any]) -> None:
    latency = policy.get("latency_targets_seconds") or {}
    by_family = latency.get("by_capability_family") or {}
    if not by_family:
        raise ValueError("latency_targets_seconds.by_capability_family is empty")
    for family, tiers in by_family.items():
        if not isinstance(tiers, dict):
            raise ValueError(f"{family} latency tiers must be an object")
        for tier_name, limits in tiers.items():
            _require_keys(limits, {"p50_max", "p95_max"}, f"{family}.{tier_name}")
            p50 = float(limits["p50_max"])
            p95 = float(limits["p95_max"])
            if p50 <= 0 or p95 <= 0:
                raise ValueError(f"{family}.{tier_name} latency thresholds must be positive")
            if p95 < p50:
                raise ValueError(f"{family}.{tier_name} p95_max must be >= p50_max")


def _validate_budget_guardrails(policy: Dict[str, Any], heuristics: Dict[str, Any]) -> None:
    guardrails = policy.get("budget_guardrails_usd") or {}
    _require_keys(guardrails, {"per_run_max", "per_session_max"}, "budget_guardrails_usd")
    per_run = guardrails.get("per_run_max") or {}
    per_session = guardrails.get("per_session_max") or {}
    _require_keys(per_run, {"Local", "EC2", "EMR", "Batch", "Lambda"}, "per_run_max")
    _require_keys(per_session, {"Local", "EC2", "EMR", "Batch", "Lambda"}, "per_session_max")

    for env in per_run:
        if float(per_run[env]) < 0 or float(per_session[env]) < 0:
            raise ValueError(f"Budget guardrails must be non-negative for {env}")
        if float(per_session[env]) < float(per_run[env]):
            raise ValueError(f"per_session_max must be >= per_run_max for {env}")

    # Sanity check against cost heuristic upper bounds.
    heuristics_rows = list((heuristics.get("heuristics") or []))
    for row in heuristics_rows:
        env = str(row.get("environment") or "")
        rng = row.get("cost_range_usd") or [0.0, 0.0]
        if env in per_run and isinstance(rng, list) and len(rng) == 2:
            observed_max = float(rng[1])
            if observed_max > float(per_session[env]):
                raise ValueError(
                    f"Cost heuristic max {observed_max} for {env} exceeds per_session_max {per_session[env]}"
                )


def _validate_capacity_assumptions(policy: Dict[str, Any], env_capabilities: Dict[str, Any]) -> None:
    capacity = policy.get("capacity_assumptions") or {}
    _require_keys(capacity, {"sandbox_tier", "async_tier"}, "capacity_assumptions")
    sandbox = capacity.get("sandbox_tier") or {}
    async_tier = capacity.get("async_tier") or {}
    _require_keys(
        sandbox,
        {"max_concurrent_jobs", "queue_backlog_warning_threshold", "queue_backlog_block_threshold"},
        "capacity_assumptions.sandbox_tier",
    )
    _require_keys(
        async_tier,
        {"max_concurrent_jobs", "queue_backlog_warning_threshold", "queue_backlog_block_threshold"},
        "capacity_assumptions.async_tier",
    )
    if int(sandbox["queue_backlog_block_threshold"]) <= int(sandbox["queue_backlog_warning_threshold"]):
        raise ValueError("sandbox_tier block threshold must exceed warning threshold")
    if int(async_tier["queue_backlog_block_threshold"]) <= int(async_tier["queue_backlog_warning_threshold"]):
        raise ValueError("async_tier block threshold must exceed warning threshold")

    # Capacity sanity: at least one cloud environment marked available for async tier.
    envs = list((env_capabilities.get("environments") or []))
    available_cloud = [
        e for e in envs if str(e.get("name")) in {"EC2", "EMR", "Batch", "Lambda"} and bool(e.get("available"))
    ]
    if not available_cloud:
        raise ValueError("No available cloud execution environments in environment_capabilities.yaml")


def _validate_release_threshold_alignment(policy: Dict[str, Any], release_thresholds: Dict[str, Any]) -> None:
    regression_policy = release_thresholds.get("performance_regression_policy") or {}
    required_source = regression_policy.get("policy_source")
    if required_source != "benchmarks/performance_cost_policy.yaml":
        raise ValueError("release_thresholds.performance_regression_policy.policy_source mismatch")
    if not bool(regression_policy.get("require_capability_family_latency_targets")):
        raise ValueError("release_thresholds must require capability family latency targets")
    if not bool(regression_policy.get("require_budget_guardrails")):
        raise ValueError("release_thresholds must require budget guardrails")
    if not bool(regression_policy.get("require_capacity_assumptions")):
        raise ValueError("release_thresholds must require capacity assumptions")

    latency = release_thresholds.get("latency") or {}
    _require_keys(latency, {"workflow_median_seconds_max", "workflow_p95_seconds_max"}, "release_thresholds.latency")
    cost = release_thresholds.get("cost") or {}
    _require_keys(cost, {"per_run_usd_max", "per_session_usd_max"}, "release_thresholds.cost")

    if float(cost["per_session_usd_max"]) < float(cost["per_run_usd_max"]):
        raise ValueError("release_thresholds.cost.per_session_usd_max must be >= per_run_usd_max")


def validate() -> int:
    policy = _load_yaml(POLICY_PATH)
    release_thresholds = _load_yaml(RELEASE_THRESHOLDS_PATH)
    env_capabilities = _load_yaml(ENV_CAPABILITIES_PATH)
    heuristics = _load_yaml(HEURISTICS_PATH)

    _validate_latency_targets(policy)
    _validate_budget_guardrails(policy, heuristics)
    _validate_capacity_assumptions(policy, env_capabilities)
    _validate_release_threshold_alignment(policy, release_thresholds)

    print("Performance/cost policy validation passed.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(validate())
    except Exception as exc:  # noqa: BLE001
        print(f"Performance/cost policy validation failed: {exc}")
        sys.exit(1)
