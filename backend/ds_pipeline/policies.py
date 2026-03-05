"""
Data quality and safety policy checks.

Each check returns a list of PolicyResult objects, which are either warnings
(logged and reported) or errors (hard stops that prevent the pipeline from continuing).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional


class PolicyLevel(str, Enum):
    WARNING = "warning"
    ERROR = "error"  # hard stop — pipeline will not proceed


@dataclass
class PolicyResult:
    level: PolicyLevel
    check: str
    message: str
    details: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_error(self) -> bool:
        return self.level == PolicyLevel.ERROR

    def to_dict(self) -> Dict[str, Any]:
        return {
            "level": self.level.value,
            "check": self.check,
            "message": self.message,
            "details": self.details,
        }


# ── Individual checks ─────────────────────────────────────────────────────────

def check_missingness(
    df: Any,
    threshold: float = 0.5,
    error_threshold: float = 0.9,
) -> List[PolicyResult]:
    """Warn (≥threshold) or hard-stop (≥error_threshold) on high missingness."""
    results: List[PolicyResult] = []
    missing_rates = df.isnull().mean()
    for col, rate in missing_rates.items():
        rate = float(rate)
        if rate >= error_threshold:
            results.append(PolicyResult(
                level=PolicyLevel.ERROR,
                check="missingness",
                message=(
                    f"Column '{col}' is {rate:.0%} missing — exceeds hard-stop "
                    f"threshold ({error_threshold:.0%})."
                ),
                details={"column": col, "missing_rate": rate},
            ))
        elif rate >= threshold:
            results.append(PolicyResult(
                level=PolicyLevel.WARNING,
                check="missingness",
                message=f"Column '{col}' is {rate:.0%} missing (threshold: {threshold:.0%}).",
                details={"column": col, "missing_rate": rate},
            ))
    return results


def check_duplicates(df: Any) -> List[PolicyResult]:
    """Warn when duplicate rows are detected."""
    n_dups = int(df.duplicated().sum())
    if n_dups > 0:
        pct = n_dups / max(len(df), 1)
        return [PolicyResult(
            level=PolicyLevel.WARNING,
            check="duplicates",
            message=f"Found {n_dups} duplicate rows ({pct:.1%} of dataset).",
            details={"n_duplicates": n_dups, "n_rows": len(df)},
        )]
    return []


_PII_PATTERNS: List[tuple[str, str]] = [
    (r"email", "email-like column name"),
    (r"phone|mobile|tel$", "phone-like column name"),
    (r"^(first_?name|last_?name|full_?name|given_?name|surname)$", "name-like column name"),
    (r"ssn|social_security", "SSN-like column name"),
    (r"passport|license_?number", "ID document-like column name"),
    (r"ip_?addr(ess)?$", "IP address-like column name"),
    (r"(birth|dob|date_of_birth)", "date-of-birth-like column name"),
]


def check_pii(df: Any) -> List[PolicyResult]:
    """Warn when column names suggest the presence of PII."""
    results: List[PolicyResult] = []
    for col in df.columns:
        col_lower = col.lower()
        for pattern, reason in _PII_PATTERNS:
            if re.search(pattern, col_lower):
                results.append(PolicyResult(
                    level=PolicyLevel.WARNING,
                    check="pii",
                    message=(
                        f"Column '{col}' may contain PII ({reason}). "
                        "Review before sharing artifacts."
                    ),
                    details={"column": col, "reason": reason},
                ))
                break
    return results


def check_leakage(
    df: Any,
    target_col: Optional[str],
    time_col: Optional[str] = None,
    entity_col: Optional[str] = None,
    split_strategy: str = "random",
) -> List[PolicyResult]:
    """Heuristic data leakage checks: target name in features, temporal split, entity leakage."""
    results: List[PolicyResult] = []
    if not target_col:
        return results

    target_lower = target_col.lower()

    # Target name appears as substring of another column name
    for col in df.columns:
        if col == target_col:
            continue
        if target_lower in col.lower():
            results.append(PolicyResult(
                level=PolicyLevel.WARNING,
                check="leakage",
                message=(
                    f"Feature column '{col}' contains the target name '{target_col}' — "
                    "possible data leakage."
                ),
                details={"column": col, "target": target_col},
            ))

    # Time column present but random split used
    if time_col and split_strategy == "random":
        results.append(PolicyResult(
            level=PolicyLevel.WARNING,
            check="leakage",
            message=(
                f"A time column '{time_col}' exists but a random split is configured. "
                "Consider a temporal split to avoid future data leakage."
            ),
            details={"time_col": time_col, "split_strategy": split_strategy},
        ))

    # Entity column present — warn about entity leakage
    if entity_col:
        results.append(PolicyResult(
            level=PolicyLevel.WARNING,
            check="leakage",
            message=(
                f"An entity column '{entity_col}' is present. Ensure the same entity "
                "does not appear in both train and test sets (entity leakage)."
            ),
            details={"entity_col": entity_col},
        ))

    return results


def run_all_checks(
    df: Any,
    target_col: Optional[str] = None,
    time_col: Optional[str] = None,
    entity_col: Optional[str] = None,
    split_strategy: str = "random",
    missingness_threshold: float = 0.5,
    error_threshold: float = 0.9,
) -> List[PolicyResult]:
    """Run all policy checks and return the combined list of results."""
    results: List[PolicyResult] = []
    results.extend(check_missingness(df, missingness_threshold, error_threshold))
    results.extend(check_duplicates(df))
    results.extend(check_pii(df))
    results.extend(check_leakage(df, target_col, time_col, entity_col, split_strategy))
    return results


def has_hard_stop(results: List[PolicyResult]) -> bool:
    return any(r.is_error for r in results)


def format_policy_summary(results: List[PolicyResult]) -> str:
    if not results:
        return "All policy checks passed."
    lines = [f"Policy check results ({len(results)} issue(s)):"]
    for r in results:
        icon = "🚫" if r.is_error else "⚠️"
        lines.append(f"  {icon} [{r.check}] {r.message}")
    return "\n".join(lines)
