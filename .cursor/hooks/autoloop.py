import json
import os
from pathlib import Path

MAX_ITERATIONS = int(os.getenv("CURSOR_AUTOLOOP_MAX", "8"))

ARTIFACT_DIR = Path("artifacts")
STATE_FILE = ARTIFACT_DIR / "autoloop_state.json"
FAILURE_SUMMARY = ARTIFACT_DIR / "failure_summaries" / "latest.md"
RELEASE_FILE = ARTIFACT_DIR / "release_readiness.json"
WORKLOG = ARTIFACT_DIR / "worklog.md"


def load_state():
    if STATE_FILE.exists():
        try:
            return json.loads(STATE_FILE.read_text())
        except Exception:
            pass
    return {"iterations": 0}


def save_state(state):
    STATE_FILE.parent.mkdir(parents=True, exist_ok=True)
    STATE_FILE.write_text(json.dumps(state, indent=2))


def release_ready():
    if not RELEASE_FILE.exists():
        return False, "release_readiness.json not found"
    try:
        data = json.loads(RELEASE_FILE.read_text())
    except Exception as exc:
        return False, f"invalid release_readiness.json: {exc}"

    ready = bool(data.get("release_ready", False))
    blockers = data.get("blockers", [])
    if ready:
        return True, "release thresholds satisfied"
    if blockers:
        return False, f"blockers: {blockers[:3]}"
    return False, "release_ready is false"


def latest_failure_hint():
    if FAILURE_SUMMARY.exists():
        try:
            text = FAILURE_SUMMARY.read_text().strip()
            if text:
                return text[:1200]
        except Exception:
            pass
    return "Read the latest failing test output and fix the highest-impact root cause."


def append_worklog(msg):
    WORKLOG.parent.mkdir(parents=True, exist_ok=True)
    with WORKLOG.open("a", encoding="utf-8") as f:
        f.write(msg.rstrip() + "\n")


def main():
    _ = input() if not os.isatty(0) else ""  # hook payload, unused but consumed

    state = load_state()
    state["iterations"] += 1
    save_state(state)

    ready, reason = release_ready()
    if ready:
        append_worklog(f"- Autoloop stopped: {reason}")
        print(json.dumps({}))
        return

    if state["iterations"] >= MAX_ITERATIONS:
        append_worklog(f"- Autoloop stopped after max iterations: {reason}")
        print(json.dumps({}))
        return

    hint = latest_failure_hint()
    append_worklog(
        f"- Autoloop iteration {state['iterations']}/{MAX_ITERATIONS}: continuing because {reason}"
    )

    message = (
        f"[Autoloop iteration {state['iterations']}/{MAX_ITERATIONS}] "
        f"Tests or release gates are not yet satisfied. "
        f"Fix the highest-impact blocker by root cause, update tests, rerun the smallest "
        f"relevant verification suite, and refresh artifacts. Context: {hint}"
    )

    print(json.dumps({"followup_message": message}))


if __name__ == "__main__":
    main()