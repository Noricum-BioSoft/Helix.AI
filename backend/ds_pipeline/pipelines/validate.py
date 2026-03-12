"""
Validation pipeline step: runs all policy checks against the loaded DataFrame.

Returns a structured validation_results dict and raises on hard stops.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from backend.ds_pipeline.policies import PolicyResult, format_policy_summary, has_hard_stop, run_all_checks


class ValidationError(Exception):
    """Raised when a hard-stop policy violation is found."""


def validate(
    df: Any,
    *,
    target_col: Optional[str] = None,
    time_col: Optional[str] = None,
    entity_col: Optional[str] = None,
    split_strategy: str = "random",
    missingness_threshold: float = 0.5,
    error_threshold: float = 0.9,
    raise_on_error: bool = True,
) -> Dict[str, Any]:
    """
    Run all policy checks.

    Returns a dict with 'results' (list of PolicyResult dicts), 'passed' bool,
    and 'summary' string.  If raise_on_error=True and a hard stop is found,
    raises ValidationError.
    """
    policy_results: List[PolicyResult] = run_all_checks(
        df,
        target_col=target_col,
        time_col=time_col,
        entity_col=entity_col,
        split_strategy=split_strategy,
        missingness_threshold=missingness_threshold,
        error_threshold=error_threshold,
    )

    errors = [r for r in policy_results if r.is_error]
    warnings = [r for r in policy_results if not r.is_error]
    passed = len(errors) == 0

    result: Dict[str, Any] = {
        "passed": passed,
        "n_errors": len(errors),
        "n_warnings": len(warnings),
        "summary": format_policy_summary(policy_results),
        "results": [r.to_dict() for r in policy_results],
    }

    if not passed and raise_on_error:
        raise ValidationError(f"Validation failed: {result['summary']}")

    return result
