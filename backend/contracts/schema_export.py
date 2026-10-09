"""
Export platform contracts as JSON Schema.

DataWeaver validates Context API request bodies against these snapshots
(copied to ``DataWeaver.AI/backend/app/schemas/contracts/``) until a shared
contracts package exists. ``python -m backend.contracts.schema_export`` writes
them to ``shared/schemas/contracts/``; a unit test asserts the committed
snapshot is current.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, Type

from pydantic import BaseModel

from backend.contracts.evidence_assessment import EvidenceAssessment, NextDecision, Observation
from backend.contracts.execution_intent import ExecutionIntent
from backend.contracts.execution_recommendation import ExecutionRecommendation
from backend.contracts.execution_request import ExecutionRequest, ExecutionRun
from backend.contracts.human_approval import HumanApproval
from backend.contracts.provenance import ProvenanceRecord
from backend.contracts.provider_authorization import ProviderAuthorization
from backend.contracts.scientific_objective import ScientificObjective
from backend.contracts.scientific_plan import ScientificPlan
from backend.contracts.security_assessment import SecurityAssessment
from backend.contracts.security_review import SecurityReview
from shared.capability_registry import CapabilityDescriptor

SCHEMA_DIR = Path(__file__).resolve().parents[2] / "shared" / "schemas" / "contracts"

CONTRACTS: Dict[str, Type[BaseModel]] = {
    "scientific_objective": ScientificObjective,
    "scientific_plan": ScientificPlan,
    "security_assessment": SecurityAssessment,
    "security_review": SecurityReview,
    "execution_recommendation": ExecutionRecommendation,
    "provider_authorization": ProviderAuthorization,
    "execution_intent": ExecutionIntent,
    "human_approval": HumanApproval,
    "execution_request": ExecutionRequest,
    "execution_run": ExecutionRun,
    "observation": Observation,
    "next_decision": NextDecision,
    "evidence_assessment": EvidenceAssessment,
    "provenance_record": ProvenanceRecord,
    "capability_descriptor": CapabilityDescriptor,
}


def render(name: str) -> str:
    schema = CONTRACTS[name].model_json_schema()
    schema["$id"] = f"https://noricum.bio/schemas/contracts/{name}.json"
    return json.dumps(schema, indent=2, sort_keys=True) + "\n"


def export_all(target: Path = SCHEMA_DIR) -> Dict[str, Path]:
    target.mkdir(parents=True, exist_ok=True)
    written = {}
    for name in CONTRACTS:
        path = target / f"{name}.json"
        path.write_text(render(name), encoding="utf-8")
        written[name] = path
    return written


def stale(target: Path = SCHEMA_DIR) -> Dict[str, str]:
    """Return {name: reason} for every schema whose committed snapshot differs from the model."""
    problems = {}
    for name in CONTRACTS:
        path = target / f"{name}.json"
        if not path.exists():
            problems[name] = "missing"
        elif path.read_text(encoding="utf-8") != render(name):
            problems[name] = "outdated"
    return problems


if __name__ == "__main__":
    if "--check" in sys.argv:
        bad = stale()
        if bad:
            print("stale contract schemas:", bad)
            sys.exit(1)
        print("contract schemas up to date")
    else:
        for name, path in export_all().items():
            print(f"wrote {path}")
