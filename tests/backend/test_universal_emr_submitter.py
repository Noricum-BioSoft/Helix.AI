import os
from unittest.mock import Mock, patch

import pytest


def test_submit_universal_emr_job_creates_job_record(tmp_path, monkeypatch):
    from backend.job_manager import JobManager, STATUS_SUBMITTED

    monkeypatch.setenv("EMR_CLUSTER_ID", "j-TEST123456789")
    monkeypatch.setenv("AWS_REGION", "us-east-1")
    monkeypatch.setenv("S3_DATASET_BUCKET", "noricum-ngs-data")
    monkeypatch.setenv("S3_SCRIPT_BUCKET", "noricum-ngs-data")

    jm = JobManager()
    jm.project_root = tmp_path  # type: ignore[attr-defined]

    # Provide a minimal tools directory for bundling
    tools_dir = tmp_path / "tools"
    tools_dir.mkdir(parents=True)
    (tools_dir / "mutations.py").write_text("def run_mutation_raw(sequence: str, num_variants: int=96):\n    return {'status':'success'}\n")

    # Provide submit script
    script_dir = tmp_path / "scripts" / "emr"
    script_dir.mkdir(parents=True)
    submit_script = script_dir / "submit-universal-job.sh"
    submit_script.write_text("#!/bin/bash\necho 'Step ID: s-ABC123456'\n")
    submit_script.chmod(0o755)

    with patch.object(jm, "_check_cluster_state", return_value="WAITING"), patch(
        "subprocess.run",
        return_value=Mock(returncode=0, stdout="Step ID: s-ABC123456\n", stderr=""),
    ) as sp:
        job_id = jm.submit_universal_emr_job(
            tool_name="mutate_sequence",
            tool_args={"sequence": "ATGC", "num_variants": 2},
            session_id="sess-1",
        )

    assert job_id in jm.jobs
    job = jm.jobs[job_id]
    assert job["status"] == STATUS_SUBMITTED
    assert job["tool_name"] == "mutate_sequence"
    assert job["infra"]["service"] == "emr"
    assert job["step_id"] == "s-ABC123456"


def test_execution_broker_async_submits_universal_job(monkeypatch):
    import asyncio

    from backend.execution_broker import ExecutionBroker, ExecutionRequest

    monkeypatch.setenv("HELIX_ASYNC_BYTES_THRESHOLD", "1")  # force async based on size

    async def dummy_executor(tool_name: str, arguments: dict):
        return {"status": "success", "text": "ok"}

    broker = ExecutionBroker(tool_executor=dummy_executor)

    async def fake_submit(req):
        return {"type": "job", "status": "submitted", "job_id": "job-xyz"}

    monkeypatch.setattr(broker, "_submit_universal_emr_job", fake_submit)

    # Force size estimation above threshold by providing a local path and mocking its size.
    with patch.object(broker, "_try_get_local_size", return_value=2):
        args = {"sequence": "ATGC", "num_variants": 2, "session_id": "s", "path": "/tmp/x.bin"}

        out = asyncio.get_event_loop().run_until_complete(
            broker.execute_tool(
                ExecutionRequest(
                    tool_name="mutate_sequence",
                    arguments=args,
                    session_id="s",
                    original_command="mutate",
                    session_context={},
                )
            )
        )
    assert out["mode"] == "async"
    assert out["result"]["job_id"] == "job-xyz"


