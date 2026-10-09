#!/usr/bin/env python3
"""
LabOS — a thin stand-in for “another system that calls Helix.”

This is a minimal caller that uses Helix as an API:
  create session → submit intent → (optional) approve plan → poll runs → download bundle

It is intentionally small: Helix is an execution plane other systems call — not
a chat UI.

Usage
-----
  export HELIX_BASE_URL=http://localhost:8001   # or your beta API URL

  # Dry-run the amplicon multi-step plan (stops at the plan approval gate)
  python examples/labos_client.py plan

  # Full path: plan → approve → wait for runs → download bundle
  python examples/labos_client.py run

  # Resume an existing session (second-caller beat)
  python examples/labos_client.py iterate --session-id <id> \\
      --command "Rerun with a stricter quality threshold"

  # Health check
  python examples/labos_client.py health

Environment
-----------
  HELIX_BASE_URL   Base URL of the Helix API (default: http://localhost:8001)
  HELIX_TIMEOUT_S  HTTP timeout seconds (default: 300)
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional

# ── Defaults ──────────────────────────────────────────────────────────────────

DEFAULT_BASE = os.environ.get("HELIX_BASE_URL", "http://localhost:8001").rstrip("/")
DEFAULT_TIMEOUT = float(os.environ.get("HELIX_TIMEOUT_S", "300"))

# Multi-step amplicon prompt — matches frontend demoScenarios amplicon-qc-pipeline.
# Reliable for Plan → Approve staging in beta.
AMPLICON_INTENT = """You are processing a 16S rRNA amplicon sequencing dataset from a gut microbiome study. Raw paired-end FASTQ files are on S3 and need a full preprocessing pipeline before downstream diversity analysis.

Dataset
  Illumina MiSeq 2×250 bp paired-end reads; V3–V4 hypervariable region.
  Forward reads: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq
  Reverse reads: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq
  Output prefix:  s3://noricum-ngs-data/test-output/amplicon-demo/

Pipeline Steps
  1. Run FastQC quality assessment on both raw R1 and R2 files.
  2. Trim adapter sequences (CTGTCTCTTATACACATCT) and low-quality bases
     (Phred < 20) from both ends; minimum read length 150 bp.
  3. Merge overlapping paired-end reads with minimum overlap of 20 bp.
  4. Generate a quality report summarizing read counts before and after each step.

Desired Outputs
  - FastQC HTML reports for raw R1 and R2.
  - Trimmed FASTQ files saved to the output S3 prefix.
  - Merged FASTA file of consensus amplicon sequences.
  - Quality summary CSV (sample, raw_reads, post_trim_reads, merged_reads, merge_rate).
"""

APPROVE_COMMAND = "approve"


# ── HTTP helpers (stdlib only — no extra deps for the demo prop) ───────────────

class HelixClient:
    def __init__(self, base_url: str = DEFAULT_BASE, timeout: float = DEFAULT_TIMEOUT):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout

    def _request(
        self,
        method: str,
        path: str,
        body: Optional[Dict[str, Any]] = None,
        raw: bool = False,
    ) -> Any:
        url = f"{self.base_url}{path}"
        data = None
        headers = {"Accept": "application/json"}
        if body is not None:
            data = json.dumps(body).encode("utf-8")
            headers["Content-Type"] = "application/json"
        req = urllib.request.Request(url, data=data, headers=headers, method=method)
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                payload = resp.read()
                if raw:
                    return payload, dict(resp.headers)
                if not payload:
                    return {}
                return json.loads(payload.decode("utf-8"))
        except urllib.error.HTTPError as e:
            detail = e.read().decode("utf-8", errors="replace")
            raise SystemExit(f"HTTP {e.code} {method} {path}: {detail}") from e
        except urllib.error.URLError as e:
            raise SystemExit(
                f"Cannot reach Helix at {self.base_url} ({e.reason}). "
                f"Set HELIX_BASE_URL or start the backend."
            ) from e

    def health(self) -> Dict[str, Any]:
        return self._request("GET", "/health")

    def create_session(self, user_id: str = "labos-demo") -> str:
        out = self._request("POST", "/session/create", {"user_id": user_id})
        sid = out.get("session_id")
        if not sid:
            raise SystemExit(f"No session_id in response: {out}")
        return sid

    def execute(self, session_id: str, command: str, execute_plan: bool = False) -> Dict[str, Any]:
        return self._request(
            "POST",
            "/execute",
            {
                "session_id": session_id,
                "command": command,
                "execute_plan": execute_plan,
            },
        )

    def list_runs(self, session_id: str) -> Dict[str, Any]:
        return self._request("GET", f"/session/{session_id}/runs")

    def lineage(self, session_id: str) -> Dict[str, Any]:
        return self._request("GET", f"/session/{session_id}/lineage")

    def download_bundle(
        self,
        session_id: str,
        run_id: Optional[str] = None,
        dest: Optional[Path] = None,
    ) -> Path:
        qs = {"session_id": session_id}
        if run_id:
            qs["run_id"] = run_id
        path = "/download/bundle?" + urllib.parse.urlencode(qs)
        payload, headers = self._request("GET", path, raw=True)
        filename = "helix-bundle.zip"
        cd = headers.get("Content-Disposition") or headers.get("content-disposition") or ""
        if "filename=" in cd:
            filename = cd.split("filename=")[-1].strip().strip('"')
        out = dest or Path.cwd() / filename
        out.write_bytes(payload)
        return out


# ── Pretty printing ────────────────────────────────────────────────────────────

def _banner(title: str) -> None:
    line = "─" * 60
    print(f"\n{line}")
    print(f"  LabOS → Helix  ·  {title}")
    print(line)


def _pp(label: str, obj: Any, max_chars: int = 2400) -> None:
    print(f"\n› {label}")
    text = json.dumps(obj, indent=2, default=str)
    if len(text) > max_chars:
        text = text[:max_chars] + "\n  … (truncated for demo)"
    print(text)


def _plan_summary(response: Dict[str, Any]) -> None:
    """Extract a short human-readable plan outline."""
    tool = response.get("tool") or response.get("tool_name")
    approval = response.get("approval_required")
    execute_ready = response.get("execute_ready")
    status = response.get("status")
    print(f"\n› gate  tool={tool!r}  status={status!r}  "
          f"approval_required={approval}  execute_ready={execute_ready}")

    result = response.get("result") or response.get("raw_result") or {}
    if isinstance(result, dict) and result.get("type") == "execution_result":
        result = result.get("result") or {}
    steps = None
    if isinstance(result, dict):
        if isinstance(result.get("steps"), list):
            steps = result["steps"]
        elif isinstance(result.get("result"), dict) and isinstance(result["result"].get("steps"), list):
            steps = result["result"]["steps"]
        elif isinstance(result.get("plan"), dict) and isinstance(result["plan"].get("steps"), list):
            steps = result["plan"]["steps"]

    if steps:
        print("› plan steps")
        for i, step in enumerate(steps, 1):
            if isinstance(step, dict):
                name = step.get("name") or step.get("tool") or step.get("id") or step.get("action") or "step"
                print(f"    {i}. {name}")
            else:
                print(f"    {i}. {step}")


def _needs_approval(response: Dict[str, Any]) -> bool:
    if response.get("approval_required"):
        return True
    tool = response.get("tool") or response.get("tool_name")
    if tool == "__plan__":
        # Plan responses are staged even when the flag is nested oddly
        if response.get("execute_ready") or response.get("status") in {
            "workflow_planned",
            "waiting_for_approval",
        }:
            return True
    return False


# ── Commands ──────────────────────────────────────────────────────────────────

def cmd_health(client: HelixClient) -> None:
    _banner("health")
    _pp("GET /health", client.health())


def cmd_plan(client: HelixClient, command: str) -> str:
    """Create a session, submit the intent, and show the plan approval gate."""
    _banner("1 · create session")
    print(f"› POST {client.base_url}/session/create")
    sid = client.create_session()
    print(f"› session_id = {sid}")

    _banner("2 · submit intent (caller → Helix)")
    print(f"› POST {client.base_url}/execute")
    print(f"› command preview: {command[:120].replace(chr(10), ' ')}…")
    resp = client.execute(sid, command)
    _plan_summary(resp)
    _pp("response", resp)

    if _needs_approval(resp):
        print("\n✓ Plan gate engaged — Helix is waiting for approval.")
        print(f"  Next:  python examples/labos_client.py approve --session-id {sid}")
        print(f"  Or:    python examples/labos_client.py run   # does plan+approve+bundle")
    else:
        print("\n⚠ No approval gate on this response (tool may have executed directly).")
        print("  Prefer the default amplicon multi-step intent.")

    print(f"\nSESSION_ID={sid}")
    return sid


def cmd_approve(client: HelixClient, session_id: str) -> Dict[str, Any]:
    _banner("3 · approve plan")
    print(f"› POST /execute  command={APPROVE_COMMAND!r}  session_id={session_id}")
    resp = client.execute(session_id, APPROVE_COMMAND)
    _plan_summary(resp)
    _pp("response", resp)
    return resp


def cmd_runs(client: HelixClient, session_id: str) -> List[Dict[str, Any]]:
    _banner("4 · runs + lineage")
    runs = client.list_runs(session_id)
    _pp("GET /session/{id}/runs", runs)
    try:
        lin = client.lineage(session_id)
        _pp("GET /session/{id}/lineage", lin, max_chars=1200)
    except SystemExit as e:
        print(f"› lineage skipped: {e}")
    return runs.get("runs") or []


def cmd_bundle(client: HelixClient, session_id: str, run_id: Optional[str], dest: Path) -> Path:
    _banner("5 · download reproducibility bundle")
    print(f"› GET /download/bundle?session_id={session_id}"
          + (f"&run_id={run_id}" if run_id else ""))
    path = client.download_bundle(session_id, run_id=run_id, dest=dest)
    print(f"✓ wrote {path}  ({path.stat().st_size:,} bytes)")
    print("  Bundle contents: README.md · run_manifest.json · analysis.py · plots/")
    return path


def cmd_run(
    client: HelixClient,
    command: str,
    out_dir: Path,
    poll_s: float,
    max_wait_s: float,
) -> None:
    """Full path for a single take: plan → approve → poll → bundle."""
    sid = cmd_plan(client, command)

    if True:  # always attempt approve; harmless if nothing pending
        time.sleep(0.5)
        cmd_approve(client, sid)

    _banner("4 · poll session runs")
    deadline = time.time() + max_wait_s
    runs: List[Dict[str, Any]] = []
    while time.time() < deadline:
        runs = client.list_runs(sid).get("runs") or []
        print(f"› runs so far: {len(runs)}")
        for r in runs[-5:]:
            print(f"    - {r.get('run_id')}  tool={r.get('tool')}  cmd={(r.get('command') or '')[:60]}")
        # Bundle needs at least one scriptable run; don't wait forever on long jobs
        if runs:
            break
        time.sleep(poll_s)
    else:
        print("⚠ Timed out waiting for runs — downloading bundle anyway (may 404).")

    run_id = runs[-1]["run_id"] if runs else None
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / f"helix-bundle-{sid[:8]}.zip"
    try:
        cmd_bundle(client, sid, run_id, dest)
    except SystemExit as e:
        print(f"⚠ Bundle not ready yet: {e}")
        print("  Re-try later:")
        print(f"    python examples/labos_client.py bundle --session-id {sid}")

    print(f"\nSESSION_ID={sid}")
    print("Done. Use this session_id for the ‘second caller’ beat:")
    print(f"  python examples/labos_client.py iterate --session-id {sid}")


def cmd_iterate(client: HelixClient, session_id: str, command: str) -> None:
    _banner("second caller · same session")
    print(f"› session_id = {session_id}")
    print(f"› command = {command}")
    resp = client.execute(session_id, command)
    _plan_summary(resp)
    _pp("response", resp)
    cmd_runs(client, session_id)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        description="LabOS — thin Helix API client",
    )
    parser.add_argument(
        "--base-url",
        default=DEFAULT_BASE,
        help=f"Helix API base URL (default: {DEFAULT_BASE})",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    sub.add_parser("health", help="GET /health")

    p_plan = sub.add_parser("plan", help="Create session + submit intent; stop at plan gate")
    p_plan.add_argument("--command", default=AMPLICON_INTENT, help="Intent text (default: amplicon demo)")

    p_approve = sub.add_parser("approve", help="Approve a pending plan")
    p_approve.add_argument("--session-id", required=True)

    p_runs = sub.add_parser("runs", help="List runs + lineage")
    p_runs.add_argument("--session-id", required=True)

    p_bundle = sub.add_parser("bundle", help="Download reproducibility ZIP")
    p_bundle.add_argument("--session-id", required=True)
    p_bundle.add_argument("--run-id", default=None)
    p_bundle.add_argument("--out", type=Path, default=Path("helix-bundle.zip"))

    p_run = sub.add_parser("run", help="Full demo path: plan → approve → poll → bundle")
    p_run.add_argument("--command", default=AMPLICON_INTENT)
    p_run.add_argument("--out-dir", type=Path, default=Path("tmp/labos-demo"))
    p_run.add_argument("--poll", type=float, default=3.0, help="Seconds between run polls")
    p_run.add_argument("--max-wait", type=float, default=120.0, help="Max seconds waiting for runs")

    p_iter = sub.add_parser("iterate", help="Second-caller beat against an existing session")
    p_iter.add_argument("--session-id", required=True)
    p_iter.add_argument(
        "--command",
        default="Rerun the quality report with a stricter Phred threshold of 25",
    )

    args = parser.parse_args(argv)
    client = HelixClient(base_url=args.base_url)

    print(f"LabOS client  →  {client.base_url}")

    if args.cmd == "health":
        cmd_health(client)
    elif args.cmd == "plan":
        cmd_plan(client, args.command)
    elif args.cmd == "approve":
        cmd_approve(client, args.session_id)
    elif args.cmd == "runs":
        cmd_runs(client, args.session_id)
    elif args.cmd == "bundle":
        cmd_bundle(client, args.session_id, args.run_id, args.out)
    elif args.cmd == "run":
        cmd_run(client, args.command, args.out_dir, args.poll, args.max_wait)
    elif args.cmd == "iterate":
        cmd_iterate(client, args.session_id, args.command)
    else:
        parser.error(f"unknown command {args.cmd}")


if __name__ == "__main__":
    main()
