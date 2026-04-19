"""Generate executable benchmark scenario files from the use-case catalog.

Reads `benchmarks/use_case_catalog.yaml`, maps each use case to a normalized
capability family, and emits a scenario YAML file under
`benchmarks/cases/<group>/<case-id>.yaml` suitable for the gate scorer at
`benchmarks/scoring/use_case_gate_scorer.py`.

Existing scenario files are preserved unless `--overwrite` is passed.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Dict, List

import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
CATALOG_PATH = REPO_ROOT / "benchmarks" / "use_case_catalog.yaml"
CASES_DIR = REPO_ROOT / "benchmarks" / "cases"


FAMILY_TO_GROUP: Dict[str, str] = {
    "tabular_target_ranking": "tabular",
    "survival_biomarker_modeling": "tabular",
    "bulk_rnaseq_de": "transcriptomics",
    "scrna_annotation": "transcriptomics",
    "cell_cell_communication": "transcriptomics",
    "variant_annotation": "genomics",
    "epigenomics_peak_analysis": "epigenomics",
    "network_regulatory_analysis": "transcriptomics",
}

FAMILY_TO_ROUTE: Dict[str, str] = {
    "tabular_target_ranking": "tabular_analysis",
    "survival_biomarker_modeling": "tabular_analysis",
    "bulk_rnaseq_de": "bulk_rnaseq_analysis",
    "scrna_annotation": "single_cell_analysis",
    "cell_cell_communication": "single_cell_analysis",
    "variant_annotation": "variant_annotation",
    "epigenomics_peak_analysis": "epigenomics_analysis",
    "network_regulatory_analysis": "network_analysis",
}

FORBIDDEN_BY_FAMILY: Dict[str, Dict[str, List[str]]] = {
    "tabular_target_ranking": {
        "routes": ["sequence_alignment", "single_cell_analysis", "bulk_rnaseq_analysis"],
        "fallbacks": ["sequence_alignment", "generic_rephrase_and_retry"],
    },
    "survival_biomarker_modeling": {
        "routes": ["sequence_alignment", "single_cell_analysis"],
        "fallbacks": ["sequence_alignment", "generic_rephrase_and_retry"],
    },
    "bulk_rnaseq_de": {
        "routes": ["sequence_alignment", "single_cell_analysis"],
        "fallbacks": ["sequence_alignment"],
    },
    "scrna_annotation": {
        "routes": ["sequence_alignment", "bulk_rnaseq_analysis"],
        "fallbacks": ["sequence_alignment"],
    },
    "cell_cell_communication": {
        "routes": ["sequence_alignment", "bulk_rnaseq_analysis"],
        "fallbacks": ["sequence_alignment"],
    },
    "variant_annotation": {
        "routes": ["sequence_alignment", "single_cell_analysis", "bulk_rnaseq_analysis"],
        "fallbacks": ["sequence_alignment", "generic_rephrase_and_retry"],
    },
    "epigenomics_peak_analysis": {
        "routes": ["single_cell_analysis", "bulk_rnaseq_analysis"],
        "fallbacks": ["sequence_alignment"],
    },
    "network_regulatory_analysis": {
        "routes": ["sequence_alignment", "single_cell_analysis"],
        "fallbacks": ["sequence_alignment"],
    },
}

STEPS_BY_FAMILY: Dict[str, List[str]] = {
    "tabular_target_ranking": [
        "ingest_tabular",
        "detect_numeric_columns",
        "derive_column_ratio",
        "rank_and_select",
        "explain_selection",
    ],
    "survival_biomarker_modeling": [
        "ingest_tabular",
        "validate_schema",
        "feature_selection",
        "model_fit_or_score",
        "summarize_and_explain",
    ],
    "bulk_rnaseq_de": [
        "ingest_tabular",
        "validate_metadata",
        "run_qc",
        "run_differential_expression",
        "run_pathway_enrichment",
    ],
    "scrna_annotation": [
        "ingest_scrna",
        "run_qc",
        "cluster",
        "annotate_cell_types",
        "summarize_and_explain",
    ],
    "cell_cell_communication": [
        "ingest_scrna",
        "run_qc",
        "infer_interactions",
        "summarize_and_explain",
    ],
    "variant_annotation": [
        "ingest_variants",
        "annotate",
        "prioritize",
        "summarize_and_explain",
    ],
    "epigenomics_peak_analysis": [
        "ingest_reads_or_peaks",
        "run_qc",
        "call_or_compare_peaks",
        "summarize_and_explain",
    ],
    "network_regulatory_analysis": [
        "ingest_tabular",
        "build_network",
        "analyze_modules",
        "summarize_and_explain",
    ],
}

SAFETY_REQS_DEFAULT = ["explicit_assumptions", "uncertainty_note"]
SAFETY_REQS_BY_RISK: Dict[str, List[str]] = {
    "low": ["explicit_assumptions"],
    "medium": ["explicit_assumptions", "uncertainty_note"],
    "high": [
        "explicit_assumptions",
        "uncertainty_note",
        "clinical_interpretation_disclaimer",
    ],
}

PROVENANCE_DEFAULT = [
    "manifest_v1_present",
    "lineage_edges_present",
    "fallback_events_present",
    "route_decision_trace_present",
]


def _prompt_for_case(case: Dict[str, Any]) -> str:
    name = str(case.get("source_use_case_name") or case.get("id"))
    family = str(case.get("normalized_primary_family") or "unknown")
    return (
        f"Perform the '{name}' workflow for the provided inputs. "
        f"Use the appropriate {family.replace('_', ' ')} capability and return "
        "a scientifically defensible result with explicit assumptions."
    )


def _scenario_from_case(case: Dict[str, Any]) -> Dict[str, Any]:
    family = str(case.get("normalized_primary_family") or "")
    route = FAMILY_TO_ROUTE.get(family, "generic_analysis")
    forbidden = FORBIDDEN_BY_FAMILY.get(family, {"routes": [], "fallbacks": []})
    steps = STEPS_BY_FAMILY.get(family, ["ingest", "analyze", "summarize_and_explain"])
    risk = str(case.get("risk_class") or "medium")
    outputs = list(case.get("expected_output_types") or ["report", "explanation"])

    scenario: Dict[str, Any] = {
        "id": f"case-{case.get('id')}",
        "title": str(case.get("source_use_case_name") or case.get("id")),
        "source_refs": [
            {
                "source": case.get("source"),
                "source_use_case_id": case.get("id"),
                "source_url": case.get("source_url"),
            }
        ],
        "family": family,
        "difficulty": case.get("difficulty", "intermediate"),
        "risk_class": risk,
        "input_fixtures": [
            {
                "fixture_id": f"{case.get('id')}-input",
                "path": f"fixtures/{FAMILY_TO_GROUP.get(family, 'misc')}/{case.get('id')}.placeholder",
                "expected_format": "TBD",
                "checksum_sha256": "TBD",
            }
        ],
        "prompt": _prompt_for_case(case),
        "expected_route": route,
        "forbidden_routes": forbidden["routes"],
        "forbidden_fallbacks": forbidden["fallbacks"],
        "required_steps": steps,
        "required_outputs": outputs,
        "safety_requirements": SAFETY_REQS_BY_RISK.get(risk, SAFETY_REQS_DEFAULT),
        "provenance_requirements": PROVENANCE_DEFAULT,
        "scoring": {
            "weights": {
                "routing_correctness": 0.25,
                "workflow_correctness": 0.20,
                "output_quality": 0.20,
                "safety_caveats": 0.10,
                "provenance_completeness": 0.15,
                "replay_reproducibility": 0.10,
            }
        },
    }
    return scenario


def _target_path_for_case(case: Dict[str, Any]) -> Path:
    family = str(case.get("normalized_primary_family") or "")
    group = FAMILY_TO_GROUP.get(family, "misc")
    return CASES_DIR / group / f"{case.get('id')}.yaml"


def generate(overwrite: bool = False) -> Dict[str, Any]:
    catalog = yaml.safe_load(CATALOG_PATH.read_text()) or {}
    cases: List[Dict[str, Any]] = list(catalog.get("use_cases") or [])

    summary = {"written": [], "skipped_existing": [], "errors": []}
    for case in cases:
        try:
            target = _target_path_for_case(case)
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists() and not overwrite:
                summary["skipped_existing"].append(str(target.relative_to(REPO_ROOT)))
                continue
            scenario = _scenario_from_case(case)
            target.write_text(yaml.safe_dump(scenario, sort_keys=False))
            summary["written"].append(str(target.relative_to(REPO_ROOT)))
        except Exception as exc:  # noqa: BLE001 - surface generation errors
            summary["errors"].append({"case_id": case.get("id"), "error": str(exc)})
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing scenario files")
    args = parser.parse_args()

    summary = generate(overwrite=args.overwrite)

    print(f"written: {len(summary['written'])}")
    for p in summary["written"]:
        print(f"  + {p}")
    print(f"skipped_existing: {len(summary['skipped_existing'])}")
    for p in summary["skipped_existing"]:
        print(f"  = {p}")
    if summary["errors"]:
        print(f"errors: {len(summary['errors'])}")
        for e in summary["errors"]:
            print(f"  ! {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
