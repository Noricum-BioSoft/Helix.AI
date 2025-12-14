import os
from unittest.mock import patch

import pytest


@pytest.mark.asyncio
async def test_policy_threshold_env_var_promotes_to_async(monkeypatch):
    from backend.execution_broker import ExecutionBroker, ExecutionRequest

    monkeypatch.setenv("HELIX_ASYNC_BYTES_THRESHOLD", str(10))

    async def dummy_executor(tool_name: str, arguments: dict):
        return {"status": "success", "text": "ok", "echo": arguments}

    broker = ExecutionBroker(tool_executor=dummy_executor)

    # fake input asset discovery via local path size mock
    with patch.object(broker, "_try_get_local_size", return_value=11):
        out = await broker.execute_tool(
            ExecutionRequest(
                tool_name="sequence_alignment",
                arguments={"path": "/tmp/bigfile.bin"},
                session_id="s",
                original_command="x",
                session_context={},
            )
        )

    assert out["type"] == "execution_result"
    assert out["routing"]["threshold_bytes"] == 10
    assert out["routing"]["estimated_bytes"] == 11
    assert out["routing"]["mode"] == "async"


@pytest.mark.asyncio
async def test_tool_override_force_async(monkeypatch):
    from backend.execution_broker import ExecutionBroker, ExecutionRequest

    # Large threshold; override should still force async for fastqc
    monkeypatch.setenv("HELIX_ASYNC_BYTES_THRESHOLD", str(10_000_000_000))

    async def dummy_executor(tool_name: str, arguments: dict):
        return {"status": "success", "text": "ok"}

    broker = ExecutionBroker(tool_executor=dummy_executor)

    async def fake_submit(req):
        return {"type": "job", "status": "submitted", "job_id": "job-123"}

    monkeypatch.setattr(broker, "_submit_fastqc_job", fake_submit)

    out = await broker.execute_tool(
        ExecutionRequest(
            tool_name="fastqc_quality_analysis",
            arguments={"input_r1": "s3://b/k1", "input_r2": "s3://b/k2", "session_id": "s"},
            session_id="s",
            original_command="fastqc",
            session_context={},
        )
    )

    assert out["mode"] == "async"
    assert out["routing"]["override"] == "force_async"
    assert out["result"]["job_id"] == "job-123"


def test_input_discovery_includes_session_dataset_references():
    from backend.execution_broker import ExecutionBroker

    async def dummy_executor(tool_name: str, arguments: dict):
        return {"status": "success"}

    broker = ExecutionBroker(tool_executor=dummy_executor)

    session_context = {
        "metadata": {
            "dataset_references": [
                {"s3_bucket": "bucket", "s3_key": "datasets/x/file.fq", "size": 123, "dataset_id": "x"},
            ]
        }
    }
    inputs = broker._discover_inputs(arguments={}, session_context=session_context)
    assert any(i.uri == "s3://bucket/datasets/x/file.fq" and i.size_bytes == 123 for i in inputs)


@pytest.mark.asyncio
async def test_sync_execution_still_calls_tool_executor():
    from backend.execution_broker import ExecutionBroker, ExecutionRequest

    called = {}

    async def dummy_executor(tool_name: str, arguments: dict):
        called["tool_name"] = tool_name
        return {"status": "success", "text": "done", "value": 1}

    broker = ExecutionBroker(tool_executor=dummy_executor)
    out = await broker.execute_tool(
        ExecutionRequest(
            tool_name="fetch_ncbi_sequence",
            arguments={"accession": "NC_000001.11"},
            session_id="s",
            original_command="fetch",
            session_context={},
        )
    )
    assert called["tool_name"] == "fetch_ncbi_sequence"
    assert out["mode"] == "sync"
    assert out["status"] == "success"
    assert out["result"]["value"] == 1


