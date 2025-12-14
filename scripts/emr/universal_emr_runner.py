#!/usr/bin/env python3
"""
Universal EMR runner (Phase 2).

Responsibilities:
- Read a payload either from S3 (s3://bucket/key) or local file
- Fetch referenced S3 inputs to local disk (best-effort)
- Execute a single step (v1): run a known Helix tool (from tools bundle) OR run arbitrary python code (from payload)
- Write results.json + logs to S3 output prefix

This runner is intentionally self-contained and uses the AWS CLI for S3 IO on EMR.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple


S3_URI_RE = re.compile(r"^s3://([^/]+)/(.+)$")


def _now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _run(cmd: list[str], *, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, check=check, capture_output=True, text=True)


def _aws_s3_cp(src: str, dst: str) -> None:
    _run(["aws", "s3", "cp", src, dst], check=True)


def _aws_s3_sync(src: str, dst: str) -> None:
    _run(["aws", "s3", "sync", src, dst], check=True)


def _download_s3_to_file(uri: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    _aws_s3_cp(uri, str(path))


def _upload_file_to_s3(path: Path, uri: str) -> None:
    _aws_s3_cp(str(path), uri)


def _looks_like_s3_uri(value: Any) -> bool:
    return isinstance(value, str) and value.startswith("s3://")


def _iter_strings(obj: Any) -> Iterable[str]:
    if isinstance(obj, str):
        yield obj
    elif isinstance(obj, dict):
        for v in obj.values():
            yield from _iter_strings(v)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            yield from _iter_strings(v)


def _materialize_s3_inputs(arguments: Dict[str, Any], workdir: Path) -> Dict[str, Any]:
    """
    Download any S3 URIs found in arguments to local files and rewrite the argument values
    to local paths. This is intentionally simple: it preserves only basenames.
    """
    rewritten = json.loads(json.dumps(arguments))  # deep copy (json-safe)

    def rewrite(obj: Any) -> Any:
        if _looks_like_s3_uri(obj):
            basename = Path(obj.split("/", 3)[-1]).name
            local_path = workdir / "inputs" / basename
            try:
                _download_s3_to_file(obj, local_path)
                return str(local_path)
            except Exception:
                # Best-effort: keep original URI if download fails.
                return obj
        if isinstance(obj, dict):
            return {k: rewrite(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [rewrite(v) for v in obj]
        return obj

    return rewrite(rewritten)


def _ensure_tools_bundle(tools_bundle_s3: Optional[str], workdir: Path) -> Optional[Path]:
    """
    Download and unzip a tools bundle (zip) to workdir/tools_bundle and add it to sys.path.
    Returns the extracted path or None.
    """
    if not tools_bundle_s3:
        return None
    extract_dir = workdir / "tools_bundle"
    extract_dir.mkdir(parents=True, exist_ok=True)
    zip_path = workdir / "tools_bundle.zip"
    _download_s3_to_file(tools_bundle_s3, zip_path)
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(extract_dir)
    # Add extracted root to path so `import mutations`, etc. works
    sys.path.insert(0, str(extract_dir))
    return extract_dir


def _dispatch_tool(tool_name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
    """
    Minimal tool dispatcher for EMR v1.
    We intentionally support a small set of pure-Python tools that are safe on EMR.
    """
    if tool_name == "mutate_sequence":
        import mutations

        return mutations.run_mutation_raw(
            sequence=arguments.get("sequence", ""),
            num_variants=int(arguments.get("num_variants", 96)),
        )

    if tool_name == "sequence_alignment":
        import alignment

        return alignment.run_alignment(arguments.get("sequences", ""))

    if tool_name == "read_trimming":
        import read_trimming

        return read_trimming.run_read_trimming_raw(
            reads=arguments.get("reads", ""),
            adapter=arguments.get("adapter"),
            quality_threshold=int(arguments.get("quality_threshold", 20)),
        )

    if tool_name == "read_merging":
        import read_merging

        return read_merging.run_read_merging_raw(
            forward_reads=arguments.get("forward_reads", ""),
            reverse_reads=arguments.get("reverse_reads", ""),
            min_overlap=int(arguments.get("min_overlap", 12)),
        )

    if tool_name == "quality_assessment":
        import quality_assessment

        return quality_assessment.run_quality_assessment_raw(arguments.get("sequences", ""))

    raise ValueError(f"Unsupported tool for EMR universal runner: {tool_name}")


def _execute_python_code(code: str, workdir: Path) -> Dict[str, Any]:
    """
    Execute arbitrary python code in a subprocess.
    Convention: if the process prints JSON to stdout, we parse it as the result.
    Otherwise we return stdout/stderr/returncode.
    """
    script_path = workdir / "payload_code.py"
    script_path.write_text(code)
    proc = subprocess.run(
        [sys.executable, str(script_path)],
        cwd=str(workdir),
        capture_output=True,
        text=True,
    )

    stdout = (proc.stdout or "").strip()
    stderr = (proc.stderr or "").strip()
    result: Dict[str, Any] = {
        "returncode": proc.returncode,
        "stdout": stdout,
        "stderr": stderr,
    }
    # Best-effort parse JSON stdout
    if stdout:
        try:
            parsed = json.loads(stdout)
            if isinstance(parsed, dict):
                result["parsed"] = parsed
        except Exception:
            pass
    return result


def _load_payload(args: argparse.Namespace, workdir: Path) -> Dict[str, Any]:
    if args.payload_s3:
        payload_path = workdir / "payload.json"
        _download_s3_to_file(args.payload_s3, payload_path)
        return json.loads(payload_path.read_text())
    if args.payload_file:
        return json.loads(Path(args.payload_file).read_text())
    raise ValueError("Must provide --payload-s3 or --payload-file")


def main() -> int:
    parser = argparse.ArgumentParser(description="Universal EMR runner (Helix.AI Phase 2)")
    parser.add_argument("--payload-s3", help="S3 URI to payload JSON")
    parser.add_argument("--payload-file", help="Local path to payload JSON")
    parser.add_argument("--output-s3", required=True, help="S3 prefix to write outputs (ends with /)")
    parser.add_argument("--job-id", default=os.getenv("EMR_STEP_ID") or "", help="Job/step id for logging")
    args = parser.parse_args()

    started = time.time()
    job_id = args.job_id or f"job-{int(started)}"

    workdir = Path(os.getenv("HELIX_EMR_WORKDIR") or f"/tmp/helix-emr-{job_id}")
    workdir.mkdir(parents=True, exist_ok=True)
    log_path = workdir / "runner.log"

    # Tee stdout/stderr into a log file
    class Tee:
        def __init__(self, stream, file):
            self.stream = stream
            self.file = file

        def write(self, data):
            self.stream.write(data)
            self.file.write(data)

        def flush(self):
            self.stream.flush()
            self.file.flush()

    with log_path.open("a") as lf:
        sys.stdout = Tee(sys.stdout, lf)  # type: ignore[assignment]
        sys.stderr = Tee(sys.stderr, lf)  # type: ignore[assignment]

        print("==========================================")
        print("Helix.AI Universal EMR Runner")
        print("==========================================")
        print(f"job_id: {job_id}")
        print(f"started_at: {_now_iso()}")
        print(f"workdir: {workdir}")
        print(f"output_s3: {args.output_s3}")

        status = "success"
        error: Optional[str] = None
        tool_output: Any = None
        payload: Dict[str, Any] = {}

        try:
            payload = _load_payload(args, workdir)

            tool_name = payload.get("tool_name") or payload.get("tool") or ""
            tool_args = payload.get("arguments") or payload.get("args") or {}
            tools_bundle_s3 = payload.get("tools_bundle_s3")
            python_code = payload.get("python_code")
            python_code_s3 = payload.get("python_code_s3")

            # Ensure tools bundle if provided
            _ensure_tools_bundle(tools_bundle_s3, workdir)

            # Download S3 inputs and rewrite args (best-effort)
            if isinstance(tool_args, dict):
                tool_args = _materialize_s3_inputs(tool_args, workdir)

            if python_code_s3 and isinstance(python_code_s3, str):
                code_path = workdir / "python_code.py"
                _download_s3_to_file(python_code_s3, code_path)
                python_code = code_path.read_text()

            if python_code and isinstance(python_code, str):
                tool_output = _execute_python_code(python_code, workdir)
            else:
                if not tool_name:
                    raise ValueError("payload missing tool_name (or python_code/python_code_s3)")
                tool_output = _dispatch_tool(tool_name, tool_args if isinstance(tool_args, dict) else {})

        except Exception as e:
            status = "error"
            error = str(e)
            tool_output = {"status": "error", "error": error}

        finished = time.time()
        results = {
            "job_id": job_id,
            "status": status,
            "error": error,
            "started_at": started,
            "finished_at": finished,
            "duration_s": round(finished - started, 3),
            "payload": payload,
            "result": tool_output,
        }

        results_path = workdir / "results.json"
        results_path.write_text(json.dumps(results, indent=2, default=str))

        # Upload results + logs
        output_prefix = args.output_s3 if args.output_s3.endswith("/") else args.output_s3 + "/"
        try:
            _upload_file_to_s3(results_path, f"{output_prefix}results.json")
            _upload_file_to_s3(log_path, f"{output_prefix}runner.log")
        except Exception as e:
            print(f"WARNING: failed to upload outputs to S3: {e}")
            # still exit non-zero if job failed

        return 0 if status == "success" else 2


if __name__ == "__main__":
    raise SystemExit(main())


