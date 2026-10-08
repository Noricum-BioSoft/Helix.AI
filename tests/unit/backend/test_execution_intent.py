"""Phase 1.3 — intent builder: hash moves with provider/inputs/params; frozen; manifest covers local files."""

from __future__ import annotations

import copy
from types import SimpleNamespace

import pytest
from pydantic import ValidationError

from backend.orchestration.intent_builder import (
    build_intent,
    collect_input_manifest,
    compute_input_manifest_hash,
    legacy_provider_id,
)
from backend.orchestration.objective_builder import build_objective
from backend.orchestration.plan_staging import build_scientific_plan


@pytest.fixture(autouse=True)
def _mock_mode(monkeypatch):
    monkeypatch.setenv("HELIX_MOCK_MODE", "1")


def _plan(fastq="s3://bucket/reads.fastq.gz", threads=4):
    plan_dict = {
        "version": "v1",
        "steps": [
            {
                "id": "step1",
                "tool_name": "fastqc_quality_analysis",
                "arguments": {"fastq_file": fastq, "threads": threads},
            }
        ],
    }
    return build_scientific_plan(plan_dict, build_objective("x"), "x")


def test_intent_hash_changes_with_provider():
    plan = _plan()
    local = build_intent(plan, None)
    emr = build_intent(plan, SimpleNamespace(infrastructure="EMR"))
    assert local.provider_id == "legacy:Local" and emr.provider_id == "legacy:EMR"
    assert local.execution_intent_hash != emr.execution_intent_hash
    assert local.plan_hash == emr.plan_hash  # same plan, different intent


def test_intent_hash_changes_with_inputs_and_params():
    base = build_intent(_plan(), None)
    other_input = build_intent(_plan(fastq="s3://bucket/other.fastq.gz"), None)
    other_param = build_intent(_plan(threads=8), None)
    assert base.input_manifest_hash != other_input.input_manifest_hash
    assert base.execution_parameters_hash != other_param.execution_parameters_hash
    assert len({base.execution_intent_hash, other_input.execution_intent_hash, other_param.execution_intent_hash}) == 3


def test_intent_is_frozen_and_hash_is_validated():
    intent = build_intent(_plan(), None)
    with pytest.raises(ValidationError):
        intent.provider_id = "legacy:EC2"  # type: ignore[misc]
    tampered = {**intent.model_dump(), "provider_id": "legacy:EC2"}  # keeps old hash
    with pytest.raises(ValidationError):
        type(intent)(**tampered)


def test_input_manifest_hashes_local_files(tmp_path):
    f = tmp_path / "reads.fastq"
    f.write_bytes(b"@r1\nACGT\n+\nIIII\n")
    plan = _plan(fastq=str(f))
    manifest = collect_input_manifest(plan)
    assert len(manifest) == 1
    assert manifest[0]["size"] == f.stat().st_size
    assert manifest[0]["content_hash"] and len(manifest[0]["content_hash"]) == 64
    h1 = compute_input_manifest_hash(manifest)
    f.write_bytes(b"@r1\nACGA\n+\nIIII\n")
    h2 = compute_input_manifest_hash(collect_input_manifest(plan))
    assert h1 != h2  # same path, different content → different intent


def test_input_manifest_is_order_independent():
    plan_dict = {
        "steps": [
            {"id": "s1", "tool_name": "t", "arguments": {"input_a": "s3://b/a", "input_b": "s3://b/b"}},
        ]
    }
    swapped = copy.deepcopy(plan_dict)
    swapped["steps"][0]["arguments"] = {"input_b": "s3://b/b", "input_a": "s3://b/a"}
    obj = build_objective("x")
    a = build_scientific_plan(plan_dict, obj, "x")
    b = build_scientific_plan(swapped, obj, "x")
    assert compute_input_manifest_hash(collect_input_manifest(a)) == compute_input_manifest_hash(collect_input_manifest(b))


def test_legacy_provider_id_defaults_to_local():
    assert legacy_provider_id(None) == "legacy:Local"
    assert legacy_provider_id("Batch") == "legacy:Batch"
