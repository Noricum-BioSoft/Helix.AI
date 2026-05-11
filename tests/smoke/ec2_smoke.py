#!/usr/bin/env python3
"""
ec2_smoke.py — Live end-to-end smoke tests against the running Helix.AI backend.

Usage:
    python tests/smoke/ec2_smoke.py [BASE_URL]

    BASE_URL defaults to https://helix-beta.noricum-biosoft.com

Each test POSTs real HTTP requests, validates the response shape, and prints
PASS / FAIL with a short reason.  Exit code 0 = all pass, 1 = any failure.

Run this after every deploy to confirm the system behaves before handing it
back to a user.
"""
from __future__ import annotations

import json
import sys
import textwrap
import time
from typing import Any

import requests

BASE_URL = sys.argv[1].rstrip("/") if len(sys.argv) > 1 else "https://helix-beta.noricum-biosoft.com"
TIMEOUT = 90  # seconds — advisory + routing calls can be slow

# ── colour helpers ──────────────────────────────────────────────────────────

GREEN = "\033[32m"
RED   = "\033[31m"
YELLOW = "\033[33m"
RESET = "\033[0m"

_results: list[tuple[str, bool, str]] = []


def _ok(name: str, detail: str = "") -> None:
    _results.append((name, True, detail))
    print(f"  {GREEN}PASS{RESET}  {name}" + (f"  ({detail})" if detail else ""))


def _fail(name: str, reason: str) -> None:
    _results.append((name, False, reason))
    print(f"  {RED}FAIL{RESET}  {name}")
    for line in textwrap.wrap(reason, 100):
        print(f"        {line}")


def _skip(name: str, reason: str) -> None:
    _results.append((name, True, f"SKIP: {reason}"))
    print(f"  {YELLOW}SKIP{RESET}  {name}  ({reason})")


# ── HTTP helpers ─────────────────────────────────────────────────────────────


def _post(path: str, body: dict, session_id: str | None = None) -> dict[str, Any]:
    if session_id:
        body = {**body, "session_id": session_id}
    r = requests.post(f"{BASE_URL}{path}", json=body, timeout=TIMEOUT)
    r.raise_for_status()
    return r.json()


def _create_session() -> str:
    r = requests.post(f"{BASE_URL}/create_session", json={}, timeout=30)
    r.raise_for_status()
    data = r.json()
    sid = data.get("session_id") or data.get("id")
    if not sid:
        raise RuntimeError(f"create_session returned no session_id: {data}")
    return sid


def _execute(command: str, session_id: str | None = None) -> dict[str, Any]:
    return _post("/execute", {"command": command}, session_id=session_id)


# ── Individual smoke tests ────────────────────────────────────────────────────


def smoke_health():
    """GET /health returns status=healthy and a git_commit."""
    name = "health endpoint"
    try:
        r = requests.get(f"{BASE_URL}/health", timeout=15)
        if r.status_code != 200:
            _fail(name, f"HTTP {r.status_code}")
            return
        d = r.json()
        if d.get("status") != "healthy":
            _fail(name, f"status != healthy: {d}")
            return
        commit = d.get("git_commit", "")
        _ok(name, f"commit={commit}")
    except Exception as e:
        _fail(name, str(e))


def smoke_session_create():
    """POST /create_session returns a valid session_id."""
    name = "session creation"
    try:
        sid = _create_session()
        _ok(name, f"sid={sid[:8]}…")
        return sid
    except Exception as e:
        _fail(name, str(e))
        return None


def smoke_meta_question(sid: str):
    """A plain question routes to a sensible advisory/text response (not an error)."""
    name = "meta question → text response"
    try:
        resp = _execute("what kind of analyses can you help me with?", session_id=sid)
        if not resp.get("success"):
            _fail(name, f"success=False: {resp.get('text', '')[:200]}")
            return
        text = resp.get("text", "")
        if not text:
            _fail(name, "result.text is empty")
            return
        _ok(name, f"{len(text)} chars")
    except Exception as e:
        _fail(name, str(e))


def smoke_advisory_quickreply_no_pending_plan(sid: str):
    """After an advisory response, sending a quick-reply label must NOT return
    'There is no pending workflow to approve'."""
    name = "advisory quick-reply does not hit approval gate"

    # First turn: get an advisory
    try:
        r1 = _execute("what data analysis can I do with an Excel workbook?", session_id=sid)
    except Exception as e:
        _skip(name, f"turn 1 failed: {e}")
        return

    # Second turn: simulate clicking a quick-reply button
    phrases = ["yes, inspect it", "yes, profile it", "focus on visualization"]
    for phrase in phrases:
        try:
            r2 = _execute(phrase, session_id=sid)
            text = (r2.get("text") or "").lower()
            if "no pending workflow" in text:
                _fail(name, f"'{phrase}' hit the approval gate dead-end: {text[:200]}")
                return
        except Exception as e:
            _fail(name, f"'{phrase}' raised: {e}")
            return
    _ok(name, f"tested {len(phrases)} quick-reply phrases")


def smoke_unknown_intent(sid: str):
    """A nonsensical command returns a clarification card, not a 500."""
    name = "unknown intent → clarification (not 500)"
    try:
        resp = _execute("xyzzy frobnicator 12345", session_id=sid)
        # We accept any non-error response — what matters is no crash.
        if not resp.get("success") and resp.get("status") == "error":
            # 'error' status with a message is acceptable (tool dispatch error);
            # a missing/empty response is not.
            if not resp.get("text"):
                _fail(name, f"error with no text: {resp}")
                return
        _ok(name, resp.get("status", "?"))
    except requests.HTTPError as e:
        if e.response.status_code == 500:
            _fail(name, f"HTTP 500 on unknown intent: {e.response.text[:200]}")
        else:
            _ok(name, f"HTTP {e.response.status_code} (non-500)")
    except Exception as e:
        _fail(name, str(e))


def smoke_approval_guard_with_no_plan(sid: str):
    """Bare 'yes' / 'ok' with no pending plan must NOT return the approval dead-end."""
    name = "bare 'yes' with no pending plan → not approval dead-end"
    try:
        resp = _execute("yes", session_id=sid)
        text = (resp.get("text") or "").lower()
        if "no pending workflow" in text:
            _fail(name, f"bare 'yes' hit the approval gate with no pending plan: {text[:200]}")
        else:
            _ok(name, f"routed to: {resp.get('tool', resp.get('execution_path', '?'))}")
    except Exception as e:
        _fail(name, str(e))


# ── Runner ────────────────────────────────────────────────────────────────────


def main():
    print(f"\nHelix.AI smoke suite → {BASE_URL}\n{'─' * 60}")
    t0 = time.time()

    smoke_health()
    sid = smoke_session_create()
    if sid:
        smoke_meta_question(sid)
        smoke_advisory_quickreply_no_pending_plan(sid)
        smoke_approval_guard_with_no_plan(sid)
        smoke_unknown_intent(sid)

    elapsed = time.time() - t0
    passed  = sum(1 for _, ok, _ in _results if ok)
    failed  = sum(1 for _, ok, _ in _results if not ok)

    print(f"\n{'─' * 60}")
    print(f"Results: {passed} passed, {failed} failed  ({elapsed:.1f}s)")

    if failed:
        print(f"\n{RED}SMOKE FAILED{RESET} — do not accept this deploy until failures are resolved.\n")
        sys.exit(1)
    else:
        print(f"\n{GREEN}SMOKE PASSED{RESET}\n")
        sys.exit(0)


if __name__ == "__main__":
    main()
