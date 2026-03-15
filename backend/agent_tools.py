# /backend/agent_tools.py
# All @tool decorated functions for the Helix.AI agent

import logging
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from langchain_core.tools import tool

logger = logging.getLogger(__name__)

# Ensure tools directory is in Python path for imports
# This allows importing modules from the tools/ directory
_TOOLS_DIR = Path(__file__).resolve().parent.parent / "tools"
if str(_TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(_TOOLS_DIR))


@tool
def toolbox_inventory() -> Dict:
    """
    List all tools Helix.AI has access to (registered tools, discovered @tool functions, and local/EC2 CLI tools).
    
    **IMPORTANT: Only use this tool when the user explicitly asks about available tools, capabilities, or what you can do.**
    Examples: "What tools do you have?", "Show me your capabilities", "List available tools"
    
    **DO NOT use this tool as a fallback when you don't recognize a command or operation.**
    If the user requests a bioinformatics operation (e.g., "merge reads", "align sequences"), and no matching tool exists,
    do NOT use this tool - instead, let the system fall through to the tool-generator-agent which will dynamically generate the appropriate tool.
    """
    from backend.tool_inventory import build_toolbox_inventory, format_toolbox_inventory_markdown
    inv = build_toolbox_inventory()
    return {
        "text": format_toolbox_inventory_markdown(inv),
        "input": {},
        "output": inv,
        "plot": {},
    }


def _get_session_local_dir(session_id: str) -> Path:
    """
    Best-effort local session directory resolver.

    Local-only feature (used for iterative workflows and artifacts). In cloud mode,
    artifacts may be stored in S3 instead.
    """
    try:
        from backend.history_manager import history_manager

        session = history_manager.get_session(session_id) if session_id else None
        local_path = None
        if isinstance(session, dict):
            local_path = (session.get("metadata") or {}).get("local_path")
        if local_path:
            p = Path(local_path)
        else:
            p = Path("sessions") / session_id
        p.mkdir(parents=True, exist_ok=True)
        return p
    except Exception:
        p = Path("sessions") / (session_id or "unknown")
        p.mkdir(parents=True, exist_ok=True)
        return p


def _new_local_run_id() -> str:
    import uuid

    return f"run_{uuid.uuid4().hex[:12]}"


@tool
def local_demo_scatter_plot(
    session_id: str,
    x_scale: str = "log",
    n_points: int = 200,
    seed: int = 7,
    title: str = "Demo scatter",
) -> Dict:
    """
    Local-only demo tool that produces:
    - a CSV data artifact
    - a parameterized visualization spec
    - (if matplotlib available) a rendered PNG

    Intended to prove iterative workflows:
    1) produce an output
    2) modify the output (e.g., log->linear) without re-running upstream analysis
    """
    import json
    import math
    import random

    run_id = _new_local_run_id()
    session_dir = _get_session_local_dir(session_id)
    run_dir = session_dir / "runs" / run_id
    artifacts_dir = run_dir / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    # Generate synthetic data with a heavy-tailed x distribution
    random.seed(int(seed))
    xs = [10 ** random.uniform(-2, 3) for _ in range(max(5, int(n_points)))]
    ys = [math.log10(x) + random.gauss(0, 0.25) for x in xs]

    data_csv_path = artifacts_dir / "demo_data.csv"
    with open(data_csv_path, "w") as f:
        f.write("x,y\n")
        for x, y in zip(xs, ys):
            f.write(f"{x},{y}\n")

    viz_spec = {
        "viz_type": "scatter",
        "title": title,
        "x_field": "x",
        "y_field": "y",
        "x_scale": (x_scale or "log").lower(),
        "data_uri": str(data_csv_path),
    }

    viz_spec_path = artifacts_dir / "viz_spec.json"
    with open(viz_spec_path, "w") as f:
        json.dump(viz_spec, f, indent=2)

    plot_png_path = artifacts_dir / "plot.png"
    rendered = False
    render_error = None
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(7, 4))
        ax = fig.add_subplot(1, 1, 1)
        ax.scatter(xs, ys, s=14, alpha=0.8)
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        if viz_spec["x_scale"] == "log":
            ax.set_xscale("log")
        fig.tight_layout()
        fig.savefig(plot_png_path, dpi=150)
        plt.close(fig)
        rendered = True
    except Exception as e:
        render_error = str(e)

    artifacts: List[Dict[str, Any]] = [
        {
            "type": "csv",
            "title": "demo_data",
            "uri": str(data_csv_path),
            "format": "csv",
            "params": {"n_points": int(n_points), "seed": int(seed)},
        },
        {
            "type": "viz_spec",
            "title": "demo_scatter_spec",
            "uri": str(viz_spec_path),
            "format": "json",
            "params": {"x_scale": viz_spec["x_scale"]},
        },
    ]
    if rendered:
        artifacts.append(
            {
                "type": "plot",
                "title": "demo_scatter",
                "uri": str(plot_png_path),
                "format": "png",
                "params": {"x_scale": viz_spec["x_scale"]},
                "extra": {"viz_spec_uri": str(viz_spec_path), "data_uri": str(data_csv_path)},
            }
        )

    text = (
        f"Created demo scatter plot (local). x_scale={viz_spec['x_scale']}. "
        f"Artifacts saved under {artifacts_dir}."
    )
    if not rendered:
        text += f" (PNG render skipped: {render_error})"

    return {
        "status": "success",
        "text": text,
        "run_id": run_id,
        "viz_spec": viz_spec,
        "artifacts": artifacts,
    }


def _apply_unified_diff(original_text: str, diff_text: str) -> str:
    """
    Apply a unified diff (single-file) to a text blob.
    Supports hunks with context/add/remove lines.
    """
    import re

    base_lines = original_text.splitlines(keepends=True)
    diff_lines = diff_text.splitlines(keepends=True)

    # Skip leading headers until first hunk
    i = 0
    while i < len(diff_lines) and not diff_lines[i].startswith("@@"):
        i += 1

    out: List[str] = []
    base_idx = 0

    hunk_re = re.compile(r"^@@\s+-(\d+)(?:,(\d+))?\s+\+(\d+)(?:,(\d+))?\s+@@")

    while i < len(diff_lines):
        line = diff_lines[i]
        if not line.startswith("@@"):
            i += 1
            continue

        m = hunk_re.match(line.rstrip("\n"))
        if not m:
            raise ValueError("Invalid unified diff hunk header")
        old_start = int(m.group(1))

        # Copy unchanged prefix
        target_base_idx = max(0, old_start - 1)
        if target_base_idx < base_idx:
            raise ValueError("Overlapping or out-of-order hunks")
        out.extend(base_lines[base_idx:target_base_idx])
        base_idx = target_base_idx

        i += 1
        # Process hunk body
        while i < len(diff_lines) and not diff_lines[i].startswith("@@"):
            h = diff_lines[i]
            if h.startswith("\\"):
                # "\ No newline at end of file" – ignore
                i += 1
                continue
            if not h:
                i += 1
                continue
            tag = h[0]
            content = h[1:]
            if tag == " ":
                if base_idx >= len(base_lines) or base_lines[base_idx] != content:
                    raise ValueError("Context mismatch while applying diff")
                out.append(content)
                base_idx += 1
            elif tag == "-":
                if base_idx >= len(base_lines) or base_lines[base_idx] != content:
                    raise ValueError("Delete mismatch while applying diff")
                base_idx += 1
            elif tag == "+":
                out.append(content)
            else:
                # Unexpected line; ignore defensively
                pass
            i += 1

    out.extend(base_lines[base_idx:])
    return "".join(out)


def _read_code_fence(text: str) -> Optional[str]:
    """
    Extract the first triple-backtick fenced block content, if present.
    """
    import re

    m = re.search(r"```(?:python|diff)?\s*([\s\S]*?)\s*```", text or "", flags=re.IGNORECASE)
    if not m:
        return None
    return m.group(1)


@tool
def local_demo_plot_script(
    session_id: str,
    x_scale: str = "log",
    n_points: int = 200,
    seed: int = 7,
    title: str = "Demo scatter (script)",
) -> Dict:
    """
    Local-only demo that models a data scientist workflow:
    - write a runnable `workspace/script.py`
    - write a patchable `workspace/viz_spec.json`
    - execute the script (sandboxed) to produce artifacts
    """
    import json
    from backend.sandbox_executor import get_sandbox_executor

    run_id = _new_local_run_id()
    session_dir = _get_session_local_dir(session_id)
    run_dir = session_dir / "runs" / run_id
    workspace_dir = run_dir / "workspace"
    artifacts_dir = run_dir / "artifacts"
    workspace_dir.mkdir(parents=True, exist_ok=True)
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    data_csv_path = artifacts_dir / "demo_data.csv"
    plot_png_path = artifacts_dir / "plot.png"
    logs_path = artifacts_dir / "execution.log"

    viz_spec = {
        "viz_type": "scatter",
        "title": title,
        "x_field": "x",
        "y_field": "y",
        "x_scale": (x_scale or "log").lower(),
        "x_label": "x",
        "y_label": "y",
        "data_uri": str(data_csv_path),
        "plot_uri": str(plot_png_path),
    }
    viz_spec_path = workspace_dir / "viz_spec.json"
    with open(viz_spec_path, "w") as f:
        json.dump(viz_spec, f, indent=2)

    script_path = workspace_dir / "script.py"
    script = f"""\
import json
import math
import os
import random

def main():
    spec_path = os.environ.get("HELIX_VIZ_SPEC_PATH") or "workspace/viz_spec.json"
    if not os.path.exists(spec_path):
        spec_path = "workspace/viz_spec.json"
    out_dir = os.environ.get("HELIX_OUTPUT_DIR") or ""
    if not out_dir or not os.path.isdir(out_dir):
        out_dir = os.environ.get("HELIX_OUTPUT_DIR_HOST") or "."

    with open(spec_path, "r") as f:
        spec = json.load(f)

    n_points = int(os.environ.get("HELIX_N_POINTS", "{int(n_points)}"))
    seed = int(os.environ.get("HELIX_SEED", "{int(seed)}"))
    random.seed(seed)
    xs = [10 ** random.uniform(-2, 3) for _ in range(max(5, n_points))]
    ys = [math.log10(x) + random.gauss(0, 0.25) for x in xs]

    data_path = os.path.join(out_dir, "demo_data.csv")
    with open(data_path, "w") as f:
        f.write("x,y\\n")
        for x, y in zip(xs, ys):
            f.write(f"{{x}},{{y}}\\n")

    # Render plot if matplotlib is available
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(7, 4))
        ax = fig.add_subplot(1, 1, 1)
        ax.scatter(xs, ys, s=14, alpha=0.8)
        ax.set_title(spec.get("title") or "Demo plot")
        ax.set_xlabel(spec.get("x_label") or spec.get("x_field") or "x")
        ax.set_ylabel(spec.get("y_label") or spec.get("y_field") or "y")
        if (spec.get("x_scale") or "").lower() == "log":
            ax.set_xscale("log")
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "plot.png"), dpi=150)
        plt.close(fig)
        print("PLOT_OK")
    except Exception as e:
        print("PLOT_SKIPPED:", str(e))

if __name__ == "__main__":
    main()
"""
    with open(script_path, "w") as f:
        f.write(script)

    exec_env = {
        "HELIX_VIZ_SPEC_PATH": "/sandbox/work/workspace/viz_spec.json",
        "HELIX_OUTPUT_DIR": "/sandbox/output",
        "HELIX_OUTPUT_DIR_HOST": str(artifacts_dir),
        "HELIX_N_POINTS": str(int(n_points)),
        "HELIX_SEED": str(int(seed)),
    }

    executor = get_sandbox_executor()
    res = executor.execute_command(
        ["python3", "workspace/script.py"],
        working_dir=str(run_dir),
        output_dir=str(artifacts_dir),
        timeout=120,
        env_vars=exec_env,
    )

    with open(logs_path, "w") as f:
        f.write("STDOUT:\n")
        f.write(res.stdout or "")
        f.write("\n\nSTDERR:\n")
        f.write(res.stderr or "")

    artifacts: List[Dict[str, Any]] = [
        {"type": "code", "title": "script", "uri": str(script_path), "format": "py"},
        {"type": "viz_spec", "title": "viz_spec", "uri": str(viz_spec_path), "format": "json", "params": {"x_scale": viz_spec["x_scale"]}},
        {"type": "log", "title": "execution_log", "uri": str(logs_path), "format": "text"},
        {"type": "csv", "title": "demo_data", "uri": str(data_csv_path), "format": "csv"},
    ]
    if plot_png_path.exists():
        artifacts.append({"type": "plot", "title": "demo_scatter", "uri": str(plot_png_path), "format": "png"})

    status = "success" if res.success else "error"
    text = f"Created demo plot script and executed it. exit_code={res.exit_code}."
    if not res.success:
        text += " (Script execution failed; see execution_log.)"

    return {
        "status": status,
        "text": text,
        "run_id": run_id,
        "artifacts": artifacts,
        "viz_spec": viz_spec,
        "execution": {
            "success": res.success,
            "exit_code": res.exit_code,
            "execution_time": res.execution_time,
        },
    }


@tool
def local_edit_and_rerun_script(
    session_id: str,
    target_run: str = "latest",
    code_patch: Optional[str] = None,
    new_code: Optional[str] = None,
) -> Dict:
    """
    Local-only: edit a prior run's `workspace/script.py` (diff or replacement) and re-run it.
    """
    import json
    from backend.history_manager import history_manager
    from backend.sandbox_executor import get_sandbox_executor

    session = history_manager.get_session(session_id)
    if not session:
        return {"status": "error", "text": f"Session not found: {session_id}", "error": "SESSION_NOT_FOUND"}

    # Resolve a script-producing run (prefer those with code/viz_spec artifacts)
    target = None
    if (target_run or "").strip().lower() in {"latest", "last", "current"}:
        runs = history_manager.list_runs(session_id)
        for r in reversed(runs):
            if not isinstance(r, dict):
                continue
            arts = r.get("produced_artifacts") or []
            # For script edits, require a prior run that actually produced code.
            if any(isinstance(a, dict) and a.get("type") == "code" for a in arts):
                target = r
                break
        if target is None:
            target = history_manager.resolve_run_reference(session_id, target_run)
    else:
        target = history_manager.resolve_run_reference(session_id, target_run)
    if not target:
        return {"status": "error", "text": "No prior runs in session.", "error": "NO_RUNS"}

    parent_run_id = target.get("run_id")
    if not parent_run_id:
        return {"status": "error", "text": "Target run missing run_id.", "error": "RUN_ID_MISSING"}

    script_uri = None
    spec_uri = None
    for art in (target.get("produced_artifacts") or []):
        if not isinstance(art, dict):
            continue
        if art.get("type") == "code" and isinstance(art.get("uri"), str) and art["uri"].endswith(".py"):
            script_uri = art["uri"]
        if art.get("type") == "viz_spec" and isinstance(art.get("uri"), str) and art["uri"].endswith(".json"):
            spec_uri = art["uri"]

    if not script_uri:
        session_dir = _get_session_local_dir(session_id)
        cand = session_dir / "runs" / parent_run_id / "workspace" / "script.py"
        if cand.exists():
            script_uri = str(cand)
    if not spec_uri:
        session_dir = _get_session_local_dir(session_id)
        cand = session_dir / "runs" / parent_run_id / "workspace" / "viz_spec.json"
        if cand.exists():
            spec_uri = str(cand)

    if not script_uri:
        return {"status": "error", "text": "Could not locate a script to edit.", "error": "SCRIPT_NOT_FOUND"}

    # Allow code input as fenced block in the string
    if isinstance(new_code, str):
        fenced = _read_code_fence(new_code)
        if fenced:
            new_code = fenced
    if isinstance(code_patch, str):
        fenced = _read_code_fence(code_patch)
        if fenced:
            code_patch = fenced

    if not (new_code or code_patch):
        return {"status": "error", "text": "Provide either new_code or code_patch.", "error": "NO_EDIT_PROVIDED"}

    with open(script_uri, "r") as f:
        base_code = f.read()

    updated_code = base_code
    if new_code:
        updated_code = new_code
    elif code_patch:
        try:
            updated_code = _apply_unified_diff(base_code, code_patch)
        except Exception as e:
            return {"status": "error", "text": f"Failed to apply patch: {e}", "error": "PATCH_APPLY_FAILED"}

    new_run_id = _new_local_run_id()
    session_dir = _get_session_local_dir(session_id)
    run_dir = session_dir / "runs" / new_run_id
    workspace_dir = run_dir / "workspace"
    artifacts_dir = run_dir / "artifacts"
    workspace_dir.mkdir(parents=True, exist_ok=True)
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    # Copy viz_spec forward (best-effort)
    viz_spec = None
    if spec_uri:
        try:
            with open(spec_uri, "r") as f:
                viz_spec = json.load(f)
            with open(workspace_dir / "viz_spec.json", "w") as f:
                json.dump(viz_spec, f, indent=2)
        except Exception:
            viz_spec = None

    new_script_path = workspace_dir / "script.py"
    with open(new_script_path, "w") as f:
        f.write(updated_code)

    logs_path = artifacts_dir / "execution.log"
    data_csv_path = artifacts_dir / "demo_data.csv"
    plot_png_path = artifacts_dir / "plot.png"

    exec_env = {
        "HELIX_VIZ_SPEC_PATH": "/sandbox/work/workspace/viz_spec.json",
        "HELIX_OUTPUT_DIR": "/sandbox/output",
        "HELIX_OUTPUT_DIR_HOST": str(artifacts_dir),
    }

    executor = get_sandbox_executor()
    res = executor.execute_command(
        ["python3", "workspace/script.py"],
        working_dir=str(run_dir),
        output_dir=str(artifacts_dir),
        timeout=120,
        env_vars=exec_env,
    )

    with open(logs_path, "w") as f:
        f.write("STDOUT:\n")
        f.write(res.stdout or "")
        f.write("\n\nSTDERR:\n")
        f.write(res.stderr or "")

    artifacts: List[Dict[str, Any]] = [
        {"type": "code", "title": "script", "uri": str(new_script_path), "format": "py"},
        {"type": "log", "title": "execution_log", "uri": str(logs_path), "format": "text"},
    ]
    if (workspace_dir / "viz_spec.json").exists():
        artifacts.append({"type": "viz_spec", "title": "viz_spec", "uri": str(workspace_dir / "viz_spec.json"), "format": "json"})
    if data_csv_path.exists():
        artifacts.append({"type": "csv", "title": "demo_data", "uri": str(data_csv_path), "format": "csv"})
    if plot_png_path.exists():
        artifacts.append({"type": "plot", "title": "demo_scatter", "uri": str(plot_png_path), "format": "png"})

    status = "success" if res.success else "error"
    text = f"Edited script and re-ran it. exit_code={res.exit_code}."
    if not res.success:
        text += " (Execution failed; see execution_log.)"

    return {
        "status": status,
        "text": text,
        "run_id": new_run_id,
        "parent_run_id": parent_run_id,
        "artifacts": artifacts,
        "viz_spec": viz_spec,
        "execution": {
            "success": res.success,
            "exit_code": res.exit_code,
            "execution_time": res.execution_time,
        },
    }


@tool
def local_update_scatter_x_scale(
    session_id: str,
    x_scale: str,
    target_run: str = "latest",
) -> Dict:
    """
    Local-only: update the latest demo scatter plot's x-axis scale (log <-> linear)
    and re-render a new artifact without regenerating upstream data.
    """
    import json

    from backend.history_manager import history_manager

    session = history_manager.get_session(session_id)
    if not session:
        return {"status": "error", "text": f"Session not found: {session_id}", "error": "SESSION_NOT_FOUND"}

    # Find target run (by run_id, "first", "latest", or iteration index).
    # For "latest", prefer the most recent run that actually produced a viz spec/plot,
    # not necessarily the most recent run overall (which might be a Q&A query).
    target = None
    if (target_run or "").strip().lower() in {"latest", "last", "current"}:
        runs = history_manager.list_runs(session_id)
        for r in reversed(runs):
            if not isinstance(r, dict):
                continue
            arts = r.get("produced_artifacts") or []
            if any(isinstance(a, dict) and a.get("type") in {"viz_spec", "plot"} for a in arts):
                target = r
                break
        if target is None:
            target = history_manager.resolve_run_reference(session_id, target_run)
    else:
        target = history_manager.resolve_run_reference(session_id, target_run)
    if not target:
        return {"status": "error", "text": "No prior runs in session.", "error": "NO_RUNS"}

    parent_run_id = target.get("run_id")
    # Find a viz spec URI (prefer artifact registry from produced_artifacts)
    spec_uri = None
    data_uri = None
    for art in (target.get("produced_artifacts") or []):
        if not isinstance(art, dict):
            continue
        if art.get("type") == "viz_spec" and isinstance(art.get("uri"), str):
            spec_uri = art["uri"]
        if art.get("type") == "csv" and isinstance(art.get("uri"), str) and "demo_data" in art["uri"]:
            data_uri = art["uri"]
    # Fallback: try to locate spec file under the run directory
    if not spec_uri and parent_run_id:
        session_dir = _get_session_local_dir(session_id)
        cand = session_dir / "runs" / parent_run_id / "artifacts" / "viz_spec.json"
        if cand.exists():
            spec_uri = str(cand)

    if not spec_uri:
        # Fallback: re-run the last analysis with the updated x_scale.
        return _rerun_with_plot_patch(
            session_id, parent_run_id, {"x_scale": (x_scale or "log2")}
        )

    # Load spec + data
    with open(spec_uri, "r") as f:
        viz_spec = json.load(f)

    if not data_uri:
        data_uri = viz_spec.get("data_uri")
    if not data_uri:
        return {"status": "error", "text": "Could not locate data for visualization.", "error": "DATA_NOT_FOUND"}

    # Read CSV (simple parser; avoid pandas dependency)
    xs: List[float] = []
    ys: List[float] = []
    with open(str(data_uri), "r") as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split(",")
            if len(parts) != 2:
                continue
            try:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
            except Exception:
                continue

    new_run_id = _new_local_run_id()
    session_dir = _get_session_local_dir(session_id)
    run_dir = session_dir / "runs" / new_run_id
    artifacts_dir = run_dir / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    new_spec = dict(viz_spec)
    new_spec["x_scale"] = (x_scale or "linear").lower()
    new_spec["derived_from_run_id"] = parent_run_id
    new_spec["derived_from_spec_uri"] = spec_uri

    new_spec_path = artifacts_dir / "viz_spec.json"
    with open(new_spec_path, "w") as f:
        json.dump(new_spec, f, indent=2)

    plot_png_path = artifacts_dir / "plot.png"
    rendered = False
    render_error = None
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(7, 4))
        ax = fig.add_subplot(1, 1, 1)
        ax.scatter(xs, ys, s=14, alpha=0.8)
        ax.set_title(new_spec.get("title") or "Demo scatter")
        ax.set_xlabel(new_spec.get("x_field") or "x")
        ax.set_ylabel(new_spec.get("y_field") or "y")
        if new_spec["x_scale"] == "log":
            ax.set_xscale("log")
        fig.tight_layout()
        fig.savefig(plot_png_path, dpi=150)
        plt.close(fig)
        rendered = True
    except Exception as e:
        render_error = str(e)

    artifacts: List[Dict[str, Any]] = [
        {
            "type": "viz_spec",
            "title": "demo_scatter_spec",
            "uri": str(new_spec_path),
            "format": "json",
            "params": {"x_scale": new_spec["x_scale"]},
            "derived_from": None,
        }
    ]
    if rendered:
        artifacts.append(
            {
                "type": "plot",
                "title": "demo_scatter",
                "uri": str(plot_png_path),
                "format": "png",
                "params": {"x_scale": new_spec["x_scale"]},
                "extra": {"viz_spec_uri": str(new_spec_path), "data_uri": str(data_uri)},
            }
        )

    text = f"Updated demo scatter x_scale to {new_spec['x_scale']} (local)."
    if not rendered:
        text += f" (PNG render skipped: {render_error})"

    return {
        "status": "success",
        "text": text,
        "run_id": new_run_id,
        "parent_run_id": parent_run_id,
        "viz_spec": new_spec,
        "artifacts": artifacts,
    }


def _resolve_latest_viz_run(session_id: str, target_run: str) -> Optional[Dict[str, Any]]:
    from backend.history_manager import history_manager

    if (target_run or "").strip().lower() in {"latest", "last", "current"}:
        runs = history_manager.list_runs(session_id)
        for r in reversed(runs):
            if not isinstance(r, dict):
                continue
            arts = r.get("produced_artifacts") or []
            if any(isinstance(a, dict) and a.get("type") in {"viz_spec", "plot", "image"} for a in arts):
                return r
        return history_manager.resolve_run_reference(session_id, target_run)
    return history_manager.resolve_run_reference(session_id, target_run)


def _deep_merge(base: Dict[str, Any], patch: Dict[str, Any]) -> Dict[str, Any]:
    out = dict(base)
    for k, v in patch.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)  # type: ignore[arg-type]
        else:
            out[k] = v
    return out


def _render_viz_if_possible(viz_spec: Dict[str, Any], xs: Optional[List[float]], ys: Optional[List[float]], out_png_path) -> Tuple[bool, Optional[str]]:
    """
    Best-effort local renderer. Supports the demo scatter spec today.
    Returns (rendered, error).
    """
    try:
        viz_type = (viz_spec.get("viz_type") or "scatter").lower()
        if viz_type not in {"scatter"}:
            return False, f"Unsupported viz_type for local render: {viz_type}"

        if xs is None or ys is None:
            return False, "Missing data for rendering"

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(7, 4))
        ax = fig.add_subplot(1, 1, 1)
        ax.scatter(xs, ys, s=14, alpha=0.8)
        ax.set_title(viz_spec.get("title") or "Visualization")
        ax.set_xlabel(viz_spec.get("x_label") or viz_spec.get("x_field") or "x")
        ax.set_ylabel(viz_spec.get("y_label") or viz_spec.get("y_field") or "y")
        if (viz_spec.get("x_scale") or "").lower() == "log":
            ax.set_xscale("log")
        fig.tight_layout()
        fig.savefig(out_png_path, dpi=150)
        plt.close(fig)
        return True, None
    except Exception as e:
        return False, str(e)


# ---------------------------------------------------------------------------
# Helper: re-execute the last analytical run with an updated plot parameter.
# Used as a fallback when local_edit_visualization cannot locate a viz_spec
# (e.g. the last run used inline base64 plots via BioOrchestrator).
# ---------------------------------------------------------------------------
_PLOT_PARAM_TOOLS = {"bulk_rnaseq_analysis", "single_cell_analysis"}
_SCALE_KEYWORDS = {"x_scale", "y_scale", "scale", "fold_change_scale", "axis_scale"}


def _rerun_with_plot_patch(
    session_id: str,
    parent_run_id: Optional[str],
    patch_obj: Dict[str, Any],
) -> Dict[str, Any]:
    """Re-run the last relevant analytical tool with updated plot parameters.

    Returns a result dict compatible with the standard response builder.
    """
    from backend.history_manager import history_manager

    runs = history_manager.list_runs(session_id)
    target_run: Optional[Dict[str, Any]] = None
    for r in reversed(runs):
        if not isinstance(r, dict):
            continue
        if r.get("tool") in _PLOT_PARAM_TOOLS:
            target_run = r
            break

    if not target_run:
        return {
            "status": "error",
            "text": (
                "No prior analysis run found in this session. "
                "Please run an analysis first, then request plot changes."
            ),
            "error": "NO_ANALYSIS_RUN",
            "parent_run_id": parent_run_id,
        }

    tool_name = target_run["tool"]
    tool_args: Dict[str, Any] = dict(target_run.get("tool_args") or {})
    prev_run_id = target_run.get("run_id")

    # Normalise the scale value from the patch to the canonical key expected
    # by the tool (e.g. {"x_scale": "linear"} or {"scale": "log"}).
    x_scale = "log2"
    for key in ("x_scale", "scale", "fold_change_scale", "axis_scale"):
        if key in patch_obj:
            raw = str(patch_obj[key]).lower()
            x_scale = "linear" if raw.startswith("lin") else "log2"
            break

    tool_args["x_scale"] = x_scale

    # Re-execute the tool (only bulk RNA-seq supported for now).
    try:
        if tool_name == "bulk_rnaseq_analysis":
            from backend import bulk_rnaseq as _brna
            result = _brna.run_deseq2_analysis(
                count_matrix=tool_args.get("count_matrix", ""),
                sample_metadata=tool_args.get("sample_metadata", ""),
                design_formula=tool_args.get("design_formula", "~condition"),
                alpha=float(tool_args.get("alpha", 0.05)),
                x_scale=x_scale,
            )
        else:
            return {
                "status": "error",
                "text": f"On-the-fly re-render is not yet supported for {tool_name}.",
                "error": "UNSUPPORTED_TOOL_FOR_RERENDER",
                "parent_run_id": parent_run_id,
            }
    except Exception as exc:
        return {
            "status": "error",
            "text": f"Re-run failed: {exc}",
            "error": "RERUN_FAILED",
            "parent_run_id": parent_run_id,
        }

    # Build visuals list from the new plots
    plots = result.get("plots", {})
    visuals = []
    for plot_key, b64 in plots.items():
        if not b64:
            continue
        label = plot_key.replace("volcano_", "Volcano: ").replace("_", " ").replace("__", " — ").title()
        if "pca" in plot_key.lower():
            label = "PCA: Sample Overview"
        visuals.append({"type": "image_b64", "data": b64, "title": label})

    scale_label = "linear fold change" if x_scale == "linear" else "log₂ fold change"
    return {
        "status": "success",
        "text": (
            f"## Volcano Plots Updated\n\n"
            f"Re-rendered **{len(visuals)} plots** with x-axis set to **{scale_label}**.\n\n"
            f"All other analysis parameters (samples, contrasts, significance thresholds) are unchanged."
        ),
        "visuals": visuals,
        "visualization_type": "results_viewer",
        "parent_run_id": prev_run_id,
        "tool_name": tool_name,
        "x_scale": x_scale,
    }


@tool
def local_edit_visualization(
    session_id: str,
    patch: Any,
    target_run: str = "latest",
    allow_data_uri_change: bool = False,
) -> Dict:
    """
    Local-only generalized visualization edit tool.

    Applies a patch (dict or JSON string) to a prior run's viz_spec.json, writes a new viz_spec,
    optionally re-renders a PNG (when supported), and creates a new run linked via parent_run_id.
    """
    import json
    from pathlib import Path

    from backend.history_manager import history_manager

    session = history_manager.get_session(session_id)
    if not session:
        return {"status": "error", "text": f"Session not found: {session_id}", "error": "SESSION_NOT_FOUND"}

    target = _resolve_latest_viz_run(session_id, target_run)
    if not target:
        return {"status": "error", "text": "No prior visualization runs in session.", "error": "NO_RUNS"}

    parent_run_id = target.get("run_id")
    if not parent_run_id:
        return {"status": "error", "text": "Target run missing run_id.", "error": "RUN_ID_MISSING"}

    # Parse patch
    patch_obj: Dict[str, Any]
    if isinstance(patch, str):
        try:
            patch_obj = json.loads(patch)
        except Exception:
            return {"status": "error", "text": "Patch must be valid JSON when provided as a string.", "error": "INVALID_PATCH_JSON"}
    elif isinstance(patch, dict):
        patch_obj = patch
    else:
        return {"status": "error", "text": "Patch must be an object (dict) or JSON string.", "error": "INVALID_PATCH_TYPE"}

    # Locate viz spec + (optional) csv data
    spec_uri = None
    data_uri = None
    spec_artifact_id = None
    for art in (target.get("produced_artifacts") or []):
        if not isinstance(art, dict):
            continue
        if art.get("type") == "viz_spec" and isinstance(art.get("uri"), str):
            spec_uri = art["uri"]
            spec_artifact_id = art.get("artifact_id")
        if art.get("type") == "csv" and isinstance(art.get("uri"), str) and "demo_data" in art["uri"]:
            data_uri = art["uri"]

    if not spec_uri:
        session_dir = _get_session_local_dir(session_id)
        cand = session_dir / "runs" / parent_run_id / "artifacts" / "viz_spec.json"
        if cand.exists():
            spec_uri = str(cand)

    if not spec_uri:
        # ── Fallback: no viz_spec on disk (e.g. inline base64 runs from BioOrchestrator).
        # If the patch is purely a plot-parameter change (x_scale, y_scale, color_scheme…)
        # we can re-execute the last analytical run with the updated parameter.
        return _rerun_with_plot_patch(session_id, parent_run_id, patch_obj)

    with open(spec_uri, "r") as f:
        viz_spec = json.load(f)

    # Guard against arbitrary file reads by patching data_uri unless explicitly allowed.
    if not allow_data_uri_change and "data_uri" in patch_obj:
        patch_obj = dict(patch_obj)
        patch_obj.pop("data_uri", None)

    # Determine data uri (allow spec default)
    if not data_uri:
        data_uri = viz_spec.get("data_uri")

    # If data_uri exists, ensure it's within the session directory when present.
    xs: Optional[List[float]] = None
    ys: Optional[List[float]] = None
    if isinstance(data_uri, str) and data_uri:
        try:
            session_root = _get_session_local_dir(session_id).resolve()
            data_path = Path(data_uri).resolve()
            if session_root not in data_path.parents and data_path != session_root:
                return {"status": "error", "text": "Data URI is outside the session directory.", "error": "DATA_URI_OUT_OF_BOUNDS"}

            # Read CSV (simple parser; avoid pandas dependency)
            xs = []
            ys = []
            with open(str(data_path), "r") as f:
                _ = f.readline()
                for line in f:
                    parts = line.strip().split(",")
                    if len(parts) != 2:
                        continue
                    try:
                        xs.append(float(parts[0]))
                        ys.append(float(parts[1]))
                    except Exception:
                        continue
        except Exception:
            xs, ys = None, None

    new_run_id = _new_local_run_id()
    session_dir = _get_session_local_dir(session_id)
    run_dir = session_dir / "runs" / new_run_id
    artifacts_dir = run_dir / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    new_spec = _deep_merge(viz_spec if isinstance(viz_spec, dict) else {}, patch_obj)
    new_spec["derived_from_run_id"] = parent_run_id
    new_spec["derived_from_spec_uri"] = spec_uri
    if spec_artifact_id:
        new_spec["derived_from_spec_artifact_id"] = spec_artifact_id

    new_spec_path = artifacts_dir / "viz_spec.json"
    with open(new_spec_path, "w") as f:
        json.dump(new_spec, f, indent=2)

    plot_png_path = artifacts_dir / "plot.png"
    rendered, render_error = _render_viz_if_possible(new_spec, xs, ys, plot_png_path)

    artifacts: List[Dict[str, Any]] = [
        {
            "type": "viz_spec",
            "title": (new_spec.get("title") or "viz_spec"),
            "uri": str(new_spec_path),
            "format": "json",
            "params": {"patch_keys": sorted(list(patch_obj.keys()))},
            "derived_from": spec_artifact_id,
        }
    ]
    if rendered:
        artifacts.append(
            {
                "type": "plot",
                "title": (new_spec.get("title") or "visualization"),
                "uri": str(plot_png_path),
                "format": "png",
                "params": {"patch_keys": sorted(list(patch_obj.keys()))},
                "extra": {"viz_spec_uri": str(new_spec_path), "data_uri": str(data_uri) if data_uri else None},
            }
        )

    text = f"Updated visualization (local). Applied patch keys: {', '.join(sorted(patch_obj.keys())) or '(none)'}."
    if not rendered and render_error:
        text += f" (PNG render skipped: {render_error})"

    return {
        "status": "success",
        "text": text,
        "run_id": new_run_id,
        "parent_run_id": parent_run_id,
        "viz_spec": new_spec,
        "artifacts": artifacts,
    }


@tool
def session_run_io_summary(session_id: str, run_ref: str = "first") -> Dict:
    """
    Deterministic, read-only helper for iterative Q&A:
    "What were the inputs/outputs of the first run?"
    """
    from backend.history_manager import history_manager

    run = history_manager.resolve_run_reference(session_id, run_ref)
    if not run:
        return {"status": "error", "text": "No runs found.", "error": "NO_RUNS"}

    artifacts = run.get("produced_artifacts") or []
    lines = []
    lines.append("### Run summary")
    lines.append(f"- **iteration_index**: {run.get('iteration_index')}")
    lines.append(f"- **run_id**: `{run.get('run_id')}`")
    lines.append(f"- **tool**: `{run.get('tool')}`")
    lines.append("")
    lines.append("### Inputs")
    inps = run.get("inputs") or []
    if inps:
        for i in inps:
            if isinstance(i, dict):
                lines.append(f"- `{i.get('uri')}`")
    else:
        lines.append("- (no structured inputs recorded)")
    lines.append("")
    lines.append("### Outputs / Artifacts")
    if artifacts:
        for a in artifacts:
            if isinstance(a, dict):
                lines.append(f"- **{a.get('type')}** `{a.get('uri')}`")
    else:
        lines.append("- (no artifacts recorded)")

    return {
        "status": "success",
        "text": "\n".join(lines),
        "run": run,
    }


@tool
def s3_browse_results(
    prefix: str,
    show: Optional[str] = None,
    recursive: bool = True,
    max_keys: int = 200,
    mode: str = "display",
) -> Dict:
    """
    List objects under an S3 prefix and (optionally) fetch and display a JSON file (e.g. results.json).

    Use this when the user asks to display/show/list/view results stored at an S3 location.
    """
    import json
    import os
    from datetime import datetime, timezone

    try:
        import boto3
    except Exception as e:
        return {
            "status": "error",
            "success": False,
            "message": "S3 browsing unavailable (boto3 not installed).",
            "text": f"S3 browsing unavailable: {e}",
            "result": {},
        }

    def parse_s3_uri(uri: str) -> tuple[str, str]:
        if not isinstance(uri, str) or not uri.startswith("s3://"):
            raise ValueError(f"Not an S3 URI: {uri}")
        path = uri[5:]
        parts = path.split("/", 1)
        bucket = parts[0]
        key = parts[1] if len(parts) > 1 else ""
        return bucket, key

    if not prefix or not isinstance(prefix, str):
        return {
            "status": "error",
            "success": False,
            "message": "prefix is required",
            "text": "prefix is required",
            "result": {},
        }

    # Normalize prefix to s3://bucket/prefix/
    if not prefix.startswith("s3://"):
        return {
            "status": "error",
            "success": False,
            "message": "prefix must be an s3:// URI",
            "text": "prefix must be an s3:// URI",
            "result": {},
        }
    if not prefix.endswith("/"):
        prefix = prefix + "/"

    bucket, key_prefix = parse_s3_uri(prefix)
    try:
        from botocore.config import Config
        s3 = boto3.client(
            "s3",
            config=Config(
                connect_timeout=5,
                read_timeout=20,
                retries={"max_attempts": 3, "mode": "standard"},
            ),
        )
    except Exception:
        s3 = boto3.client("s3")

    # List objects (best-effort; bounded by max_keys)
    objects: list[dict] = []
    continuation: Optional[str] = None
    remaining = max(1, int(max_keys))
    while remaining > 0:
        kwargs = {"Bucket": bucket, "Prefix": key_prefix, "MaxKeys": min(1000, remaining)}
        if continuation:
            kwargs["ContinuationToken"] = continuation
        resp = s3.list_objects_v2(**kwargs)
        for obj in resp.get("Contents", []) or []:
            objects.append(
                {
                    "key": obj.get("Key"),
                    "size": obj.get("Size"),
                    "last_modified": obj.get("LastModified").isoformat() if obj.get("LastModified") else None,
                }
            )
        remaining = max_keys - len(objects)
        if not resp.get("IsTruncated"):
            break
        continuation = resp.get("NextContinuationToken")
        if not continuation:
            break

    # Prefer showing results.json if the caller didn't specify a file.
    show_uri = show
    if not show_uri:
        candidate = prefix + "results.json"
        show_uri = candidate

    results_json = None
    results_json_error = None
    if show_uri:
        try:
            show_bucket, show_key = parse_s3_uri(show_uri)
            body = s3.get_object(Bucket=show_bucket, Key=show_key)["Body"].read()
            # Guardrail: avoid extremely large responses in text
            if len(body) > 2_000_000:
                results_json_error = f"{show_uri} is too large to display inline ({len(body)} bytes)"
            else:
                text = body.decode("utf-8", errors="replace")
                try:
                    results_json = json.loads(text)
                except Exception:
                    results_json = {"_raw_text": text}
        except Exception as e:
            results_json_error = str(e)

    # Build presigned URLs for common artifacts (HTML reports + JSON).
    # These are used by the frontend to render/embed instead of merely listing.
    def presign_get(bucket_name: str, key_name: str, expires_seconds: int = 3600) -> Optional[str]:
        try:
            return s3.generate_presigned_url(
                "get_object",
                Params={"Bucket": bucket_name, "Key": key_name},
                ExpiresIn=expires_seconds,
            )
        except Exception:
            return None

    def s3_uri_for_key(bucket_name: str, key_name: str) -> str:
        return f"s3://{bucket_name}/{key_name}"

    artifacts = []
    for o in objects:
        k = o.get("key")
        if not isinstance(k, str):
            continue
        base = k.rsplit("/", 1)[-1].lower()
        if base.endswith(".html") or base.endswith(".json"):
            artifacts.append(o)

    def pick_artifact(name: str) -> Optional[dict]:
        for o in artifacts:
            k = o.get("key") or ""
            if k.rsplit("/", 1)[-1] == name:
                return o
        return None

    main_html = pick_artifact("fastqc_results.html")
    if not main_html:
        # fallback: first HTML file if any
        main_html = next((o for o in artifacts if isinstance(o.get("key"), str) and o["key"].lower().endswith(".html")), None)

    links: list[dict] = []
    for o in artifacts:
        k = o.get("key")
        if not isinstance(k, str):
            continue
        url = presign_get(bucket, k)
        if not url:
            continue
        label = k.rsplit("/", 1)[-1]
        kind = "html" if label.lower().endswith(".html") else "json"
        links.append(
            {
                "label": label,
                "kind": kind,
                "s3_uri": s3_uri_for_key(bucket, k),
                "url": url,
                "size": o.get("size"),
                "last_modified": o.get("last_modified"),
            }
        )

    visuals: list[dict] = []
    if isinstance(main_html, dict) and isinstance(main_html.get("key"), str):
        main_url = presign_get(bucket, main_html["key"])
        if main_url:
            visuals.append(
                {
                    "type": "iframe",
                    "title": "FastQC report",
                    "url": main_url,
                    "s3_uri": s3_uri_for_key(bucket, main_html["key"]),
                }
            )

    # Build a human-readable text response (markdown)
    now = datetime.now(timezone.utc).isoformat()
    mode_norm = (mode or "display").strip().lower()
    want_listing = mode_norm in {"list", "listing"}

    lines: list[str] = []
    lines.append("### S3 results")
    lines.append(f"- **prefix**: `{prefix}`")
    lines.append(f"- **artifacts_found**: {len(artifacts)} (from {len(objects)} objects, max {max_keys})")

    # Add a short interpreted summary for known result shapes
    if isinstance(results_json, dict) and "basic_statistics" in results_json:
        bs = results_json.get("basic_statistics") or {}
        r1 = bs.get("r1") or {}
        r2 = bs.get("r2") or {}
        lines.append("")
        lines.append("### Summary")
        lines.append(f"- **R1 reads**: {r1.get('total_sequences')}")
        lines.append(f"- **R2 reads**: {r2.get('total_sequences')}")
        lines.append(f"- **R1 avg length**: {r1.get('average_length')}")
        lines.append(f"- **R2 avg length**: {r2.get('average_length')}")

    # Display intent: highlight the rendered report
    if visuals:
        lines.append("")
        lines.append("### Report")
        lines.append("- Rendering the main report in the UI below (FastQC HTML).")

    # Provide links for the frontend (and markdown fallback)
    if links:
        lines.append("")
        lines.append("### Artifacts")
        # Prefer main report first, then charts, then JSON
        def sort_key(l: dict) -> tuple[int, str]:
            label = (l.get("label") or "").lower()
            if label == "fastqc_results.html":
                return (0, label)
            if label.endswith(".html"):
                return (1, label)
            if label.endswith("results.json"):
                return (2, label)
            return (3, label)
        for l in sorted(links, key=sort_key):
            lines.append(f"- **{l.get('label')}**: `{l.get('s3_uri')}`")

    if want_listing and objects:
        lines.append("")
        lines.append("### Objects (listing)")
        for o in objects[: min(len(objects), 80)]:
            k = o.get("key")
            sz = o.get("size")
            lines.append(f"- `{k}` ({sz} bytes)")
        if len(objects) > 80:
            lines.append(f"- ... ({len(objects) - 80} more)")

    lines.append("")
    lines.append("### results.json")
    lines.append(f"- **requested**: `{show_uri}`")
    if results_json_error:
        lines.append(f"- **error**: {results_json_error}")
    elif results_json is None:
        lines.append("- **status**: not found")
    else:
        # In display mode, avoid spamming huge JSON in markdown; keep a short excerpt.
        limit = 5000 if not want_listing else 20000
        rendered = json.dumps(results_json, indent=2)[:limit]
        lines.append("")
        lines.append("```json")
        lines.append(rendered)
        lines.append("```")

    return {
        "status": "success",
        "success": True,
        "message": "S3 results fetched",
        "text": "\n".join(lines),
        "links": links,
        "visuals": visuals,
        "result": {
            "prefix": prefix,
            "objects": objects,
            "results_json_uri": show_uri,
            "results_json": results_json,
            "results_json_error": results_json_error,
            "links": links,
            "visuals": visuals,
            "mode": mode_norm,
            "timestamp": now,
        },
    }


@tool
def sequence_alignment(sequences: str) -> Dict:
    """Performs a sequence alignment on a given set of sequences.
    
    Use this tool ONLY when the user explicitly asks to align sequences without
    building a phylogenetic tree. If the user asks to visualize or create a
    phylogenetic tree, use phylogenetic_tree instead (which handles alignment internally).
    """
    
    # Import the alignment function
    from alignment import run_alignment_tool
    
    result = run_alignment_tool(sequences)
    
    return {
        "text": result.get("text", "Sequences aligned successfully."),
        "input": sequences,
        "output": result.get("alignment", []),
        "statistics": result.get("statistics", {}),
        "plot": {
            "data": [{"x": [1, 2, 3], "y": [3, 3, 3], "type": "bar"}],
            "layout": {"title": "Alignment Visualization"},
        },
    }


@tool
def mutate_sequence(sequence: str, num_variants: int = 96) -> Dict:
    """Mutates a given sequence and returns the specified number of variants. Use this when you need to CREATE new sequence variants from an input sequence."""
    
    # Import the mutation function
    from mutations import run_mutation_raw
    
    result = run_mutation_raw(sequence, num_variants)
    
    return {
        "text": result.get("text", "Sequence mutated successfully."),
        "input": {
            "sequence": sequence,
            "variants": num_variants,
        },
        "output": result.get("statistics", {}),
        "plot": result.get("plot", {
            "data": [{"x": [1, 2, 3], "y": [1, 1, 1], "type": "bar"}],
            "layout": {"title": "Mutation Visualization"},
        }),
    }


@tool
def dna_vendor_research(command: str, sequence_length: int = None, quantity: str = "standard") -> Dict:
    """Research DNA synthesis vendors and testing options for experimental validation. Use this when the user wants to ORDER sequences, find VENDORS, or research TESTING options. Keywords: order, vendor, synthesis, test, assay, expression, function, binding."""
    
    # Import the research function
    from dna_vendor_research import run_dna_vendor_research_raw
    
    result = run_dna_vendor_research_raw(command, sequence_length, quantity)
    
    return {
        "text": f"DNA vendor research completed: {result.get('message', 'Research done')}",
        "input": {
            "command": command,
            "sequence_length": sequence_length,
            "quantity": quantity
        },
        "output": result,
        "plot": {
            "data": [{"x": ["Vendors", "Testing Options"], "y": [result.get('total_vendors', 0), result.get('total_testing_options', 0)], "type": "bar"}],
            "layout": {"title": "DNA Vendor Research Results"},
        },
    }


@tool
def phylogenetic_tree(aligned_sequences: str) -> Dict:
    """Create phylogenetic tree visualization from sequences.
    
    This tool accepts either aligned or unaligned sequences in FASTA format.
    If sequences are unaligned, it will automatically align them first, then
    build the phylogenetic tree and create visualizations.
    
    Use this tool when the user asks to:
    - Visualize a phylogenetic tree
    - Create a phylogenetic tree
    - Build an evolutionary tree
    - Show phylogenetic relationships
    
    Do NOT use sequence_alignment for this - phylogenetic_tree handles alignment internally.
    """
    
    logger.debug("phylogenetic_tree called, sequences length=%d", len(aligned_sequences or ""))

    # Import the phylogenetic tree function
    from phylogenetic_tree import run_phylogenetic_tree

    result = run_phylogenetic_tree(aligned_sequences)
    logger.debug("phylogenetic_tree: has_newick=%s", bool(isinstance(result, dict) and result.get("tree_newick")))
    
    # Return all fields from result, especially tree_newick and ete_visualization for frontend
    if isinstance(result, dict):
        return {
            "text": result.get("text", "Phylogenetic tree created successfully."),
            "input": aligned_sequences,
            "output": result.get("statistics", {}),
            "plot": result.get("plot", {
                "data": [{"x": [1, 2, 3], "y": [1, 1, 1], "type": "bar"}],
                "layout": {"title": "Phylogenetic Tree"},
            }),
            # Pass through tree data for frontend visualization
            "tree_newick": result.get("tree_newick"),
            "ete_visualization": result.get("ete_visualization"),
            "aligned_sequences": result.get("aligned_sequences"),
            "clustering_result": result.get("clustering_result"),
            "clustered_visualization": result.get("clustered_visualization"),
            "statistics": result.get("statistics", {}),
        }
    else:
        return {
            "text": "Phylogenetic tree created successfully.",
            "input": aligned_sequences,
            "output": {},
            "plot": {
                "data": [{"x": [1, 2, 3], "y": [1, 1, 1], "type": "bar"}],
                "layout": {"title": "Phylogenetic Tree"},
            },
        }


@tool
def sequence_selection(aligned_sequences: str, selection_type: str = "random", num_sequences: int = 1) -> Dict:
    """Select sequences from aligned sequences based on various criteria (random, best_conservation, lowest_gaps, highest_gc, longest, shortest)."""
    
    logger.debug("sequence_selection: type=%s n=%s", selection_type, num_sequences)

    # Import the sequence selection function
    from sequence_selection import run_sequence_selection_raw

    result = run_sequence_selection_raw(aligned_sequences, selection_type, num_sequences)
    
    return {
        "text": result.get("text", "Sequence selection completed successfully."),
        "input": {
            "aligned_sequences": aligned_sequences,
            "selection_type": selection_type,
            "num_sequences": num_sequences
        },
        "output": result.get("selected_sequences", []),
        "selection_criteria": result.get("selection_criteria", {}),
    }


@tool
def synthesis_submission(sequences: str, vendor_preference: Optional[str] = None, quantity: str = "standard", delivery_time: str = "standard") -> Dict:
    """Submit sequences for DNA synthesis and get pricing quote. Quantity options: standard, large, custom. Delivery options: rush, standard, economy."""
    
    # Import the synthesis submission function
    from synthesis_submission import run_synthesis_submission_raw
    
    result = run_synthesis_submission_raw(sequences, vendor_preference, quantity, delivery_time)
    
    return {
        "text": result.get("text", "Synthesis submission completed successfully."),
        "input": {
            "sequences": sequences,
            "vendor_preference": vendor_preference,
            "quantity": quantity,
            "delivery_time": delivery_time
        },
        "output": result.get("quote", {}),
        "validation_results": result.get("validation_results", {}),
    }


@tool
def create_session() -> Dict:
    """Create a new session for tracking user interactions. Use this when the user wants to start a new session or create a new experiment."""
    
    import uuid
    
    session_id = str(uuid.uuid4())
    
    return {
        "text": f"Session created successfully with ID: {session_id}",
        "input": {},
        "output": {"session_id": session_id},
        "plot": {}
    }


@tool
def plasmid_visualization(vector_name: str = None, cloning_sites: str = "", insert_sequence: str = "", 
                          full_plasmid_sequence: str = None, insert_position: int = None) -> Dict:
    """Generate plasmid visualization data.
    
    Supports two modes:
    1. Full plasmid: When full_plasmid_sequence is provided, visualize the complete plasmid
    2. Insert mode: When vector_name and insert_sequence are provided, create plasmid with insert
    """
    
    # Import the plasmid visualization function
    from plasmid_visualizer import run_plasmid_visualization_raw
    
    result = run_plasmid_visualization_raw(
        vector_name=vector_name,
        cloning_sites=cloning_sites,
        insert_sequence=insert_sequence,
        full_plasmid_sequence=full_plasmid_sequence,
        insert_position=insert_position
    )
    
    return {
        "text": result.get("text", "Plasmid visualization created successfully."),
        "input": {
            "vector_name": vector_name,
            "cloning_sites": cloning_sites,
            "insert_sequence": insert_sequence,
            "full_plasmid_sequence": full_plasmid_sequence,
            "insert_position": insert_position
        },
        "output": result.get("plasmid_data", {}),
        "visualization_type": result.get("visualization_type", "circular_plasmid"),
    }

@tool
def plasmid_for_representatives(representatives: List[str], aligned_sequences: str, vector_name: str = "pUC19", cloning_sites: str = "EcoRI, BamHI, HindIII") -> Dict:
    """Create plasmid visualizations for representative sequences from clustering analysis."""
    
    # Import the plasmid visualization function
    from plasmid_visualizer import create_plasmid_for_representatives
    
    result = create_plasmid_for_representatives(representatives, aligned_sequences, vector_name, cloning_sites)
    
    return {
        "text": result.get("text", "Plasmid visualizations for representatives created successfully."),
        "plasmid_results": result.get("plasmid_results", []),
        "total_representatives": result.get("total_representatives", 0),
        "vector_name": result.get("vector_name", vector_name),
        "cloning_sites": result.get("cloning_sites", cloning_sites),
    }


@tool
def single_cell_analysis(
    data_file: Optional[str] = None,
    data_format: str = "10x",
    steps: str = "all",
    question: Optional[str] = None
) -> Dict:
    """Answer questions about single-cell RNA-seq (scRNA-seq) sequencing and data analysis, OR perform actual analysis on data.
    
    USE THIS TOOL for ANY question about:
    - Single-cell sequencing methods and workflows
    - scRNA-seq data analysis
    - Seurat analysis pipelines
    - Cell type identification and annotation
    - Marker gene discovery
    - Differential expression in single-cell data
    - Pathway analysis for single-cell data
    - Batch correction methods
    - Single-cell data formats (10x, H5, Seurat objects)
    
    This tool provides comprehensive information about single-cell analysis capabilities and can also perform actual analysis when data is provided.
    
    Args:
        data_file: Path to input data file (10x format directory, H5, CSV, or Seurat RDS). Leave empty for informational questions.
        data_format: Format of input data ('10x', 'h5', 'csv', 'seurat'). Default: '10x'. Only needed if data_file is provided.
        steps: Comma-separated list of analysis steps: 'preprocessing', 'markers', 'differential', 'pathways', 'annotation', 'batch_correction', or 'all'. Default: 'all'. Only needed if data_file is provided.
        question: The user's question about single-cell analysis. Use this parameter to capture what the user is asking about.
    
    Returns:
        Dictionary containing detailed information about single-cell analysis capabilities, or analysis results if data_file is provided.
    """
    
    # Import the single-cell analysis function
    from single_cell_analysis import analyze_single_cell_data
    
    # If no data file is provided, return informational response about single-cell analysis
    if not data_file:
        # Build a comprehensive response about single-cell analysis
        intro = "I can help you with single-cell RNA-seq (scRNA-seq) analysis! "
        if question:
            intro = f"Regarding your question about single-cell sequencing and data analysis: "
        
        info_text = f"""{intro}Here's what I can help you with:

**Single-Cell Analysis Capabilities:**
- **Preprocessing**: Quality control, normalization, and filtering of single-cell data
- **Marker Gene Identification**: Find genes that distinguish cell clusters using tools like Seurat
- **Differential Expression Analysis**: Identify genes differentially expressed between cell types or conditions
- **Pathway Enrichment Analysis**: Discover biological pathways enriched in specific cell types
- **Cell-Type Annotation**: Automatically annotate cell types using reference datasets
- **Batch Correction**: Correct for batch effects in multi-sample datasets

**Supported Data Formats:**
- 10x Genomics format (standard output from Cell Ranger)
- H5/H5AD files (AnnData format)
- CSV files (count matrices)
- Seurat RDS objects

**Analysis Workflows:**
I use the scPipeline R package (built on Seurat) to perform comprehensive single-cell analysis. You can run individual steps or a complete end-to-end analysis.

**What I can do:**
- Answer questions about single-cell sequencing methods, workflows, and interpretation
- Explain Seurat analysis pipelines and best practices
- Help with cell type identification and marker gene discovery
- Perform actual analysis when you provide single-cell data
- Guide you through differential expression and pathway analysis

**To get started with analysis:**
- Upload your single-cell data (10x format directory, H5 file, or CSV)
- Ask me to analyze it, or specify which steps you want (e.g., "find marker genes", "perform differential expression analysis")

Feel free to ask me any specific questions about single-cell sequencing and data analysis!"""
        
        return {
            "text": info_text,
            "input": {"question": question or "general inquiry about single-cell analysis"},
            "output": {
                "capabilities": [
                    "Preprocessing and quality control",
                    "Marker gene identification",
                    "Differential expression analysis",
                    "Pathway enrichment analysis",
                    "Cell-type annotation",
                    "Batch correction"
                ],
                "supported_formats": ["10x", "h5", "h5ad", "csv", "seurat"],
                "analysis_package": "scPipeline (R) / Seurat",
                "can_perform_analysis": True,
                "can_answer_questions": True
            },
            "plot": {}
        }
    
    # If data file is provided, perform actual analysis
    steps_list = [s.strip() for s in steps.split(",")] if isinstance(steps, str) else steps
    
    result = analyze_single_cell_data(
        data_file=data_file,
        data_format=data_format,
        steps=steps_list
    )
    
    return {
        "text": result.get("text", f"Single-cell analysis completed. Steps: {', '.join(steps_list)}"),
        "input": {
            "data_file": data_file,
            "data_format": data_format,
            "steps": steps_list
        },
        "output": result,
        "plot": result.get("plot", {})
    }


@tool
def fetch_ncbi_sequence(
    accession: str,
    database: str = "nucleotide"
) -> Dict:
    """Fetch DNA/RNA/protein sequence from NCBI databases by accession number.
    
    USE THIS TOOL when users ask to:
    - Fetch sequences from NCBI
    - Get sequences by accession number
    - Retrieve reference sequences
    - Download GenBank/RefSeq sequences
    
    Examples:
    - "fetch sequence NC_000001.11"
    - "get protein sequence NP_000483.1"
    - "download BRCA1 gene sequence"
    
    Args:
        accession: NCBI accession number (e.g., "NC_000001.11", "NM_000492.3", "NP_000483.1")
        database: Database to search ('nucleotide' for DNA/RNA, 'protein' for proteins). Default: 'nucleotide'
    
    Returns:
        Dictionary containing sequence data and metadata.
    """

    # Deterministic offline behavior for CI/unit tests.
    if os.getenv("HELIX_MOCK_MODE") == "1":
        mock_sequence = "ATGCGTACGTAGCTAGCTAG"
        return {
            "status": "success",
            "accession": accession,
            "text": f"Fetched sequence {accession} (mock, {len(mock_sequence)} bp)",
            "input": {"accession": accession, "database": database},
            "output": {
                "status": "success",
                "accession": accession,
                "sequence": mock_sequence,
                "description": "Mock NCBI sequence (HELIX_MOCK_MODE=1)",
                "length": len(mock_sequence),
                "database": database,
            },
            "plot": {},
        }
    
    # Import the NCBI tool function
    from ncbi_tools import fetch_sequence_from_ncbi
    
    result = fetch_sequence_from_ncbi(accession, database)
    
    if result.get("status") == "success":
        # For very large sequences (like full chromosomes), we don't need the full sequence
        # in the LLM call - just metadata and a small preview
        full_sequence = result.get("sequence", "")
        sequence_length = len(full_sequence)
        
        # For LLM: only include metadata and a tiny preview (20 bases) for validation
        # The full sequence is available in the response for frontend display
        if sequence_length > 20:
            # For large sequences, just show first 20 bases + metadata
            sequence_preview = full_sequence[:20] + f"... (full length: {sequence_length:,} bp)"
        else:
            sequence_preview = full_sequence
        
        return {
            "status": "success",
            "accession": accession,
            "text": f"Fetched sequence {accession} ({result.get('length', 0):,} bp): {result.get('description', '')}",
            "input": {"accession": accession, "database": database},
            "output": {
                "status": "success",
                "accession": accession,
                "sequence": sequence_preview,  # Minimal preview for LLM (20 bases + metadata)
                # Don't include full_sequence here - it causes LLM timeouts and huge responses
                # Full sequence is stored in tool execution result, not in LLM context
                "description": result.get("description"),
                "length": result.get("length"),
                "database": result.get("database"),
                "note": f"Full sequence ({sequence_length:,} bp) fetched and stored, preview shown above" if sequence_length > 20 else None
            },
            "plot": {}
        }
    else:
        return {
            "status": "error",
            "accession": accession,
            "text": f"Error fetching sequence {accession}: {result.get('error', 'Unknown error')}",
            "input": {"accession": accession, "database": database},
            "output": result,
            "plot": {}
        }


@tool
def query_uniprot(
    query: str,
    format: str = "fasta",
    limit: int = 10
) -> Dict:
    """Query UniProt protein database for sequences and metadata."""
    
    from uniprot_tools import query_uniprot as query_uniprot_api
    
    result = query_uniprot_api(query, format=format, limit=limit)
    
    if result.get("status") == "success":
        summary = f"{result.get('count', 0)} result(s)" if "count" in result else "Query executed"
        return {
            "text": f"UniProt query for '{query}' completed: {summary}",
            "input": {"query": query, "format": format, "limit": limit},
            "output": result,
            "plot": {}
        }
    else:
        return {
            "text": f"Error querying UniProt for '{query}': {result.get('error', 'Unknown error')}",
            "input": {"query": query, "format": format, "limit": limit},
            "output": result,
            "plot": {}
        }


@tool
def lookup_go_term(go_id: str) -> Dict:
    """Lookup Gene Ontology (GO) term details by ID."""
    
    from go_tools import lookup_go_term as lookup_go
    
    result = lookup_go(go_id)
    
    if result.get("status") == "success":
        return {
            "text": f"GO term {go_id}: {result.get('name', '')}",
            "input": {"go_id": go_id},
            "output": result,
            "plot": {}
        }
    else:
        return {
            "text": f"Error looking up GO term {go_id}: {result.get('error', 'Unknown error')}",
            "input": {"go_id": go_id},
            "output": result,
            "plot": {}
        }


@tool
def bulk_rnaseq_analysis(
    count_matrix: str,
    sample_metadata: str,
    design_formula: str = "~condition",
    alpha: float = 0.05
) -> Dict:
    """Run bulk RNA-seq differential expression analysis using DESeq2."""
    
    from backend.bulk_rnaseq import run_deseq2_analysis
    
    if not count_matrix or not sample_metadata:
        return {
            "text": "Bulk RNA-seq requires count matrix and sample metadata CSV paths.",
            "input": {
                "count_matrix": count_matrix,
                "sample_metadata": sample_metadata,
                "design_formula": design_formula,
                "alpha": alpha
            },
            "output": {"status": "error", "message": "Missing input files"},
            "plot": {}
        }
    
    result = run_deseq2_analysis(
        count_matrix=count_matrix,
        sample_metadata=sample_metadata,
        design_formula=design_formula,
        alpha=alpha
    )
    
    status = result.get("status", "success")
    summary = result.get("summary", [])
    # summary is a list of per-contrast dicts; sum up significant genes across all contrasts
    if isinstance(summary, list):
        total_sig = sum(c.get("significant", 0) for c in summary if isinstance(c, dict))
        text_summary = f"DESeq2 complete: {total_sig} significant genes across {len(summary)} contrast(s)"
    elif isinstance(summary, dict):
        text_summary = f"DESeq2 complete: {summary.get('significant_genes', 'n/a')} significant genes"
    else:
        text_summary = result.get("message", "DESeq2 analysis completed")
    
    return {
        "text": text_summary if status == "success" else f"DESeq2 error: {result.get('message', 'Unknown error')}",
        "input": {
            "count_matrix": count_matrix,
            "sample_metadata": sample_metadata,
            "design_formula": design_formula,
            "alpha": alpha
        },
        "output": result,
        "plot": {}
    }


@tool
def read_merging(
    forward_reads: str,
    reverse_reads: str,
    min_overlap: int = 12,
    output: Optional[str] = None
) -> Dict:
    """
    Merge paired-end reads (R1 and R2) into consensus sequences.
    
    Use this tool when the user wants to merge forward and reverse paired-end sequencing reads.
    Supports both S3 file paths and FASTQ content strings.
    
    Args:
        forward_reads: S3 path to forward/read 1 FASTQ file (e.g., "s3://bucket/path/R1.fq") or FASTQ content string
        reverse_reads: S3 path to reverse/read 2 FASTQ file (e.g., "s3://bucket/path/R2.fq") or FASTQ content string
        min_overlap: Minimum overlap length required to merge reads (default: 12)
        output: Optional S3 path for output merged FASTQ file (e.g., "s3://bucket/path/merged.fq")
    
    Returns:
        Dictionary containing merge status and summary metrics.
    """
    # NOTE: This function is only used by the agent for tool mapping/recognition.
    # Actual execution happens in main_with_mcp.py dispatch_tool() which directly
    # calls the read_merging module functions (merge_reads_from_s3 for S3, run_read_merging_raw for content).
    # The code below is never executed but serves as documentation of the tool's behavior.
    import read_merging
    
    # Check if inputs are S3 paths
    is_s3_path = forward_reads.startswith("s3://") or reverse_reads.startswith("s3://")
    
    if is_s3_path:
        # For S3 paths, use merge_reads_from_s3 (requires output path)
        if not output:
            # Try to infer output path from input paths
            if forward_reads.startswith("s3://"):
                # Use same directory as R1, with _merged suffix
                if forward_reads.endswith(".fq") or forward_reads.endswith(".fastq"):
                    output = forward_reads.rsplit(".", 1)[0] + "_merged.fq"
                else:
                    output = forward_reads + "_merged.fq"
            else:
                output = reverse_reads.rsplit(".", 1)[0] + "_merged.fq" if reverse_reads.endswith(".fq") or reverse_reads.endswith(".fastq") else reverse_reads + "_merged.fq"
        
        result = read_merging.merge_reads_from_s3(
            r1_path=forward_reads,
            r2_path=reverse_reads,
            output_path=output,
            min_overlap=min_overlap
        )
    else:
        # For FASTQ content strings, use run_read_merging_raw
        result = read_merging.run_read_merging_raw(
            forward_reads=forward_reads,
            reverse_reads=reverse_reads,
            min_overlap=min_overlap
        )
    
    return {
        "text": result.get("text", "Read merging completed successfully."),
        "input": {
            "forward_reads": forward_reads,
            "reverse_reads": reverse_reads,
            "min_overlap": min_overlap,
            "output": output
        },
        "output": result,
        "summary": result.get("summary", {}),
        "plot": {}
    }


@tool
def fastqc_quality_analysis(
    input_r1: str,
    input_r2: str,
    output: Optional[str] = None,
    _from_broker: bool = False,
    **kwargs
) -> Dict:
    """
    Run FastQC quality control analysis on paired-end FASTQ files.
    
    Automatically routes based on file size:
    - Small files (<100MB): Local execution (fast, 1-2 minutes)
    - Large files (>100MB): EMR execution (10-30 minutes)
    
    IMPORTANT: For paired-end reads, you must provide TWO DIFFERENT files:
    - input_r1: The FORWARD/R1/READ1 file (often contains _R1, _1, or mate_R1 in filename)
    - input_r2: The REVERSE/R2/READ2 file (often contains _R2, _2, or mate_R2 in filename)
    
    Make sure input_r1 and input_r2 are DIFFERENT files - do not use the same file for both parameters!
    
    Args:
        input_r1: S3 URI for the forward/read 1 FASTQ file. Must be different from input_r2.
                  Example: "s3://my-bucket/data/sample_R1.fastq" or "s3://my-bucket/data/mate_R1.fq"
        input_r2: S3 URI for the reverse/read 2 FASTQ file. Must be different from input_r1.
                  Example: "s3://my-bucket/data/sample_R2.fastq" or "s3://my-bucket/data/mate_R2.fq"
        output: Optional S3 URI for output directory (defaults to same directory as input files with /fastqc-results suffix)
        _from_broker: Internal flag indicating call is from ExecutionBroker (handles routing)
        **kwargs: Additional context (session_id, original_command, etc.)
    
    Returns:
        For async (EMR): Dict with job_id
        For sync (local): Dict with immediate results
    """
    # When called from ExecutionBroker with sync mode, execute locally
    if _from_broker:
        logger.info(f"✅ FastQC LOCAL execution mode - processing small files locally")

        # Local execution for small files (sandboxed or direct)
        import os
        import tempfile
        import boto3
        from pathlib import Path
        
        try:
            # Try sandbox execution first (Docker container with all tools)
            use_sandbox = os.getenv("HELIX_USE_SANDBOX", "true").lower() == "true"
            
            if use_sandbox:
                try:
                    from backend.sandbox_executor import get_sandbox_executor
                    logger.info(f"🐳 Using sandboxed execution (Docker)")
                    
                    s3_client = boto3.client('s3')
                    
                    # Parse S3 URIs
                    def parse_s3_uri(uri):
                        parts = uri.replace("s3://", "").split("/", 1)
                        return parts[0], parts[1] if len(parts) > 1 else ""
                    
                    r1_bucket, r1_key = parse_s3_uri(input_r1)
                    r2_bucket, r2_key = parse_s3_uri(input_r2)
                    
                    # Create temp directory for downloads
                    with tempfile.TemporaryDirectory() as tmpdir:
                        tmpdir_path = Path(tmpdir)
                        
                        # Download files
                        logger.info(f"📥 Downloading {input_r1}...")
                        r1_local = tmpdir_path / Path(r1_key).name
                        s3_client.download_file(r1_bucket, r1_key, str(r1_local))
                        
                        logger.info(f"📥 Downloading {input_r2}...")
                        r2_local = tmpdir_path / Path(r2_key).name
                        s3_client.download_file(r2_bucket, r2_key, str(r2_local))
                        
                        # Create output directory
                        output_dir = tmpdir_path / "fastqc_results"
                        output_dir.mkdir()
                        
                        # Execute FastQC in sandbox
                        executor = get_sandbox_executor()
                        result = executor.execute_tool(
                            tool="fastqc",
                            args=[Path(r1_key).name, Path(r2_key).name, "-o", "/sandbox/output", "-t", "2"],
                            input_files={
                                Path(r1_key).name: str(r1_local),
                                Path(r2_key).name: str(r2_local)
                            },
                            working_dir=str(tmpdir_path),
                            output_dir=str(output_dir),
                            max_memory="2g",
                            max_cpus=2,
                            timeout=300
                        )
                        
                        if not result.success:
                            raise RuntimeError(f"FastQC sandbox execution failed: {result.error_message}")
                        
                        # Upload results to S3
                        if not output:
                            # Default: same directory as input with /fastqc-results suffix
                            output = f"s3://{r1_bucket}/{Path(r1_key).parent}/fastqc-results/"
                        
                        out_bucket, out_prefix = parse_s3_uri(output)
                        
                        logger.info(f"📤 Uploading results to {output}...")
                        for result_file in output_dir.glob("*"):
                            s3_key = f"{out_prefix.rstrip('/')}/{result_file.name}"
                            s3_client.upload_file(str(result_file), out_bucket, s3_key)
                        
                        logger.info(f"✅ FastQC sandbox execution completed in {result.execution_time:.2f}s")
                        
                        return {
                            "type": "local_execution",
                            "status": "completed",
                            "message": f"FastQC completed successfully in sandboxed environment ({result.execution_time:.2f}s)",
                            "input_r1": input_r1,
                            "input_r2": input_r2,
                            "output": output,
                            "execution_mode": "sandbox",
                            "execution_time": result.execution_time,
                            "results_available": True
                        }
                
                except (ImportError, RuntimeError) as e:
                    logger.warning(f"⚠️  Sandbox execution not available: {e}. Falling back to direct execution.")
                    use_sandbox = False
            
            # Fallback: Direct host execution (if sandbox disabled or unavailable)
            if not use_sandbox:
                import subprocess
                import time
                logger.info(f"🖥️  Using direct host execution (no sandbox)")
                
                s3_client = boto3.client('s3')
                
                # Parse S3 URIs
                def parse_s3_uri(uri):
                    parts = uri.replace("s3://", "").split("/", 1)
                    return parts[0], parts[1] if len(parts) > 1 else ""
                
                r1_bucket, r1_key = parse_s3_uri(input_r1)
                r2_bucket, r2_key = parse_s3_uri(input_r2)
                
                # Create temp directory
                with tempfile.TemporaryDirectory() as tmpdir:
                    tmpdir_path = Path(tmpdir)
                    
                    # Download files
                    logger.info(f"📥 Downloading {input_r1}...")
                    r1_local = tmpdir_path / Path(r1_key).name
                    s3_client.download_file(r1_bucket, r1_key, str(r1_local))
                    
                    logger.info(f"📥 Downloading {input_r2}...")
                    r2_local = tmpdir_path / Path(r2_key).name
                    s3_client.download_file(r2_bucket, r2_key, str(r2_local))
                    
                    # Create output directory
                    output_dir = tmpdir_path / "fastqc_results"
                    output_dir.mkdir()
                    
                    # Run FastQC directly on host
                    logger.info(f"🔬 Running FastQC on host...")
                    cmd = ["fastqc", str(r1_local), str(r2_local), "-o", str(output_dir), "-t", "2"]

                    # FastQC uses Java. On minimal/cloud hosts, Java can block on low entropy.
                    # Ensure fast startup and headless mode.
                    env = os.environ.copy()
                    java_opts = env.get("JAVA_TOOL_OPTIONS", "")
                    safe_opts = "-Djava.awt.headless=true -Djava.security.egd=file:/dev/./urandom"
                    env["JAVA_TOOL_OPTIONS"] = (java_opts + " " + safe_opts).strip() if java_opts else safe_opts

                    # FastQC occasionally fails to exit cleanly on headless hosts even after writing outputs.
                    # We run it as a subprocess, then:
                    # - Prefer a normal clean exit.
                    # - If it runs "too long" but expected output files are present and stable, we terminate
                    #   the process and proceed with uploading the outputs.
                    soft_exit_seconds = int(os.getenv("HELIX_FASTQC_SOFT_EXIT_SECONDS", "90"))
                    hard_timeout_seconds = int(os.getenv("HELIX_FASTQC_HARD_TIMEOUT_SECONDS", "900"))
                    stable_seconds = int(os.getenv("HELIX_FASTQC_STABLE_OUTPUT_SECONDS", "5"))

                    expected = [
                        output_dir / f"{Path(r1_key).stem}_fastqc.zip",
                        output_dir / f"{Path(r2_key).stem}_fastqc.zip",
                    ]

                    def _stable_outputs_present() -> bool:
                        try:
                            for p in expected:
                                if not p.exists():
                                    return False
                            sizes1 = [p.stat().st_size for p in expected]
                            time.sleep(stable_seconds)
                            sizes2 = [p.stat().st_size for p in expected]
                            return sizes1 == sizes2 and all(s > 0 for s in sizes2)
                        except Exception:
                            return False

                    proc = subprocess.Popen(
                        cmd,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        env=env,
                    )
                    start = time.time()
                    terminated_for_outputs = False
                    stdout_tail = ""
                    stderr_tail = ""

                    while True:
                        rc = proc.poll()
                        elapsed = time.time() - start

                        if rc is not None:
                            # Drain pipes (avoid huge memory; keep tail only)
                            try:
                                out, err = proc.communicate(timeout=1)
                            except Exception:
                                out, err = ("", "")
                            stdout_tail = (out or "")[-2000:]
                            stderr_tail = (err or "")[-2000:]
                            if rc != 0:
                                raise RuntimeError(f"FastQC failed (rc={rc}). stderr_tail={stderr_tail!r}")
                            break

                        if elapsed >= hard_timeout_seconds:
                            proc.kill()
                            try:
                                out, err = proc.communicate(timeout=2)
                            except Exception:
                                out, err = ("", "")
                            stdout_tail = (out or "")[-2000:]
                            stderr_tail = (err or "")[-2000:]
                            raise RuntimeError(
                                f"FastQC hard-timeout after {hard_timeout_seconds}s. "
                                f"stdout_tail={stdout_tail!r} stderr_tail={stderr_tail!r}"
                            )

                        # If outputs are already present and stable, don't block user on a stuck JVM.
                        if elapsed >= soft_exit_seconds and _stable_outputs_present():
                            terminated_for_outputs = True
                            proc.terminate()
                            try:
                                proc.wait(timeout=10)
                            except Exception:
                                proc.kill()
                            try:
                                out, err = proc.communicate(timeout=2)
                            except Exception:
                                out, err = ("", "")
                            stdout_tail = (out or "")[-2000:]
                            stderr_tail = (err or "")[-2000:]
                            logger.warning(
                                "⚠️  FastQC produced expected outputs but did not exit in time; "
                                "terminated process to continue. "
                                f"stdout_tail={stdout_tail[:300]!r} stderr_tail={stderr_tail[:300]!r}"
                            )
                            break

                        time.sleep(1.0)
                    
                    # Upload results to S3
                    if not output:
                        # Default: same directory as input with /fastqc-results suffix
                        output = f"s3://{r1_bucket}/{Path(r1_key).parent}/fastqc-results/"
                    
                    out_bucket, out_prefix = parse_s3_uri(output)
                    
                    logger.info(f"📤 Uploading results to {output}...")
                    for result_file in output_dir.glob("*"):
                        s3_key = f"{out_prefix.rstrip('/')}/{result_file.name}"
                        s3_client.upload_file(str(result_file), out_bucket, s3_key)
                    
                    logger.info(
                        "✅ FastQC host execution completed successfully"
                        + (" (process terminated after outputs stabilized)" if terminated_for_outputs else "")
                    )
                    
                    return {
                        "type": "local_execution",
                        "status": "completed",
                        "message": "FastQC completed successfully (direct host execution for small files)"
                        + ("; process was terminated after outputs stabilized" if terminated_for_outputs else ""),
                        "input_r1": input_r1,
                        "input_r2": input_r2,
                        "output": output,
                        "execution_mode": "host",
                        "results_available": True
                    }
                
        except FileNotFoundError as e:
            # FastQC not installed - fall back to EMR
            logger.warning(f"⚠️  FastQC not found on host - falling back to EMR execution: {e}")
            # Fall through to EMR path below
        except Exception as e:
            logger.error(f"❌ Local FastQC execution failed: {e}")
            allow_emr_fallback = os.getenv("HELIX_FASTQC_ALLOW_EMR_FALLBACK_ON_LOCAL_FAILURE", "false").lower() == "true"
            if allow_emr_fallback:
                logger.warning("⚠️  Falling back to EMR execution (HELIX_FASTQC_ALLOW_EMR_FALLBACK_ON_LOCAL_FAILURE=true)")
                # Fall through to EMR path below
            else:
                return {
                    "type": "error",
                    "status": "error",
                    "message": f"Local FastQC execution failed (EMR fallback disabled): {str(e)}",
                    "error": str(e),
                }
    
    # EMR path: Large files or fallback from local execution
    from backend.job_manager import get_job_manager
    
    job_manager = get_job_manager()
    
    try:
        job_id = job_manager.submit_fastqc_job(
            r1_path=input_r1,
            r2_path=input_r2,
            output_path=output
        )
        
        return {
            "type": "job",
            "job_id": job_id,
            "status": "submitted",
            "message": "FastQC job submitted to EMR. Processing will take 10-30 minutes.",
            "input_r1": input_r1,
            "input_r2": input_r2,
            "output": output,
            "execution_mode": "emr"
        }
    except Exception as e:
        logger.error(f"FastQC job submission failed: {e}")
        return {
            "type": "error",
            "status": "error",
            "message": f"Failed to submit FastQC job: {str(e)}",
            "error": str(e)
        }


# ── Data Science Pipeline tools ───────────────────────────────────────────────

@tool
def ds_run_analysis(
    session_id: str,
    data_path: str,
    target_col: str = "",
    task_type: str = "auto",
    objective: str = "Analyze dataset",
    hypothesis: str = "",
    changes: str = "Initial run",
    random_seed: int = 42,
    time_col: str = "",
    entity_col: str = "",
) -> Dict:
    """
    Run the full iterative data science pipeline on a CSV dataset.

    Executes the plan→execute→evaluate→review loop:
      data audit → cleaning → EDA → baseline model → evaluation → report → next-step planning

    Parameters
    ----------
    session_id  : Helix.AI session for artifact tracking
    data_path   : path to the CSV file
    target_col  : target/label column name (empty = EDA-only, no model)
    task_type   : "auto" | "classification" | "regression"
    objective   : what you are trying to achieve this iteration
    hypothesis  : what you expect to find / change
    changes     : what is different from the previous run
    random_seed : for reproducibility
    time_col    : optional time column name (triggers temporal-split warnings)
    entity_col  : optional entity-ID column (triggers entity-leakage warnings)
    """
    from backend.ds_pipeline.orchestrator import DataScienceOrchestrator, DSRunConfig
    from backend.ds_pipeline.reviewer import review as make_review

    session_dir = _get_session_local_dir(session_id) if session_id else Path(".")
    config = DSRunConfig(
        data_path=data_path,
        target_col=target_col or None,
        task_type=task_type,
        objective=objective,
        hypothesis=hypothesis,
        changes=changes,
        random_seed=random_seed,
        time_col=time_col or None,
        entity_col=entity_col or None,
    )

    # Resolve the previous DS run in this session so iterations are linked.
    # This must happen before the orchestrator writes the new run.
    prev_run_id: Optional[str] = None
    if session_id:
        try:
            from backend.history_manager import history_manager as _hm
            _hm.ensure_session_exists(session_id)
            _runs = _hm.sessions.get(session_id, {}).get("runs", [])
            for _r in reversed(_runs):
                if _r.get("tool") == "ds_run_analysis" and _r.get("run_id"):
                    prev_run_id = _r["run_id"]
                    break
        except Exception:
            pass

    # Pass session_id=None to suppress _register_with_helix: the Helix API
    # (Phase2c / dispatch_tool) will call add_history_entry on our behalf, so
    # letting the orchestrator also call it would create duplicate run ledger
    # entries.  The orchestrator still writes artifacts/{run_id}/ and
    # experiments/experiment_log.csv via its RunStore regardless of session_id.
    orch = DataScienceOrchestrator(base_dir=session_dir, session_id=None)
    run_data = orch.run(config, parent_run_id=prev_run_id)

    # Build a structured artifacts list that add_history_entry expects:
    # [{type, title, uri, format}, ...]
    produced_artifacts: list = []
    raw_artifacts: dict = run_data.get("artifacts") or {}
    report_path = raw_artifacts.get("report_md")
    if report_path and Path(report_path).exists():
        produced_artifacts.append({
            "type": "report",
            "title": "ds_report",
            "uri": report_path,
            "format": "md",
        })
    for fig_path in run_data.get("figure_paths") or []:
        if Path(fig_path).exists():
            produced_artifacts.append({
                "type": "plot",
                "title": Path(fig_path).stem,
                "uri": fig_path,
                "format": "png",
            })
    run_json_path = str(session_dir / "artifacts" / run_data["run_id"] / "run.json")
    if Path(run_json_path).exists():
        produced_artifacts.append({
            "type": "run_json",
            "title": "run_metadata",
            "uri": run_json_path,
            "format": "json",
        })

    summary = make_review(run_data)
    return {
        "status": "success",
        "text": summary,
        "run_id": run_data["run_id"],
        "parent_run_id": prev_run_id,
        "metrics": run_data.get("metrics"),
        "decision": run_data.get("decision"),
        "next_steps": run_data.get("next_steps"),
        "artifacts": raw_artifacts,
        "produced_artifacts": produced_artifacts,
        "visualization_type": "ds_report",
        "report_md": raw_artifacts.get("report_md"),
    }


@tool
def ds_reproduce_run(
    session_id: str,
    run_id: str,
) -> Dict:
    """
    Re-run a prior data science iteration using its exact recorded configuration.

    Verifies that the data hash matches and produces a new run_id for traceability.
    """
    from backend.ds_pipeline.orchestrator import DataScienceOrchestrator
    from backend.ds_pipeline.reviewer import review as make_review

    session_dir = _get_session_local_dir(session_id) if session_id else Path(".")
    orch = DataScienceOrchestrator(base_dir=session_dir, session_id=session_id)
    run_data = orch.reproduce(run_id)

    summary = make_review(run_data)
    return {
        "status": "success",
        "text": summary,
        "run_id": run_data["run_id"],
        "parent_run_id": run_id,
        "metrics": run_data.get("metrics"),
        "decision": run_data.get("decision"),
        "visualization_type": "ds_report",
        "report_md": run_data.get("artifacts", {}).get("report_md"),
    }


@tool
def ds_diff_runs(
    session_id: str,
    run_id_a: str,
    run_id_b: str,
) -> Dict:
    """
    Compare two data science runs: configs, metrics, and next-step suggestions.
    """
    from backend.ds_pipeline.orchestrator import DataScienceOrchestrator

    session_dir = _get_session_local_dir(session_id) if session_id else Path(".")
    orch = DataScienceOrchestrator(base_dir=session_dir, session_id=session_id)
    diff = orch.diff(run_id_a, run_id_b)

    lines = [f"## Run Diff: {run_id_a} vs {run_id_b}\n"]
    if diff["config_diff"]:
        lines.append("### Config Changes")
        for k, v in diff["config_diff"].items():
            lines.append(f"- **{k}**: `{v['a']}` → `{v['b']}`")
    else:
        lines.append("Config is identical.")
    if diff["metric_diff"]:
        lines.append("\n### Metric Changes")
        for k, v in diff["metric_diff"].items():
            try:
                delta = float(v["b"]) - float(v["a"])
                lines.append(f"- **{k}**: {v['a']} → {v['b']} (Δ {delta:+.4f})")
            except (TypeError, ValueError):
                lines.append(f"- **{k}**: {v['a']} → {v['b']}")
    else:
        lines.append("\nMetrics are identical.")

    return {
        "status": "success",
        "text": "\n".join(lines),
        "diff": diff,
        "visualization_type": "markdown",
    }


@tool
def ds_list_runs(session_id: str) -> Dict:
    """
    List all data science runs in the current session with their key metrics and decisions.
    """
    from backend.ds_pipeline.run_store import RunStore

    session_dir = _get_session_local_dir(session_id) if session_id else Path(".")
    store = RunStore(base_dir=session_dir)
    log = store.read_experiment_log()

    if not log:
        return {
            "status": "success",
            "text": "No data science runs found in this session.",
            "runs": [],
        }

    lines = ["## Data Science Run History\n", "| Run ID | Timestamp | Objective | Decision | Key Metric |",
             "|---|---|---|---|---|"]
    for row in log:
        metric_str = ""
        for k in ["metric_roc_auc", "metric_f1", "metric_accuracy", "metric_r2", "metric_rmse"]:
            if row.get(k):
                name = k.replace("metric_", "")
                metric_str = f"{name}={float(row[k]):.4f}"
                break
        lines.append(
            f"| {row.get('run_id','')} | {row.get('timestamp','')[:16]} "
            f"| {row.get('objective','')[:40]} | {row.get('decision','')} | {metric_str} |"
        )

    return {
        "status": "success",
        "text": "\n".join(lines),
        "runs": log,
        "visualization_type": "markdown",
    }


@tool
def patch_and_rerun(
    session_id: str,
    change_request: str,
    target_run: str = "latest",
) -> Dict:
    """Apply a natural-language change to the most recent analysis script and re-execute it.

    This tool:
      1. Finds the ``analysis.py`` saved for the last scriptable tool run.
      2. Applies *change_request* as a parameter patch (for simple changes) or
         delegates to the LLM agent for complex code edits.
      3. Executes the new script in a subprocess.
      4. Returns the new results plus a human-readable diff of what changed.

    Parameters
    ----------
    session_id      : Active session ID.
    change_request  : Plain-English description of the change,
                      e.g. ``"change alpha to 0.01"`` or
                      ``"use linear fold change on the x-axis"``.
    target_run      : ``"latest"`` (default) or a specific run_id.
    """
    from backend.history_manager import history_manager
    from backend import script_executor as _se

    def _normalize_backend_imports(script_text: str) -> str:
        """
        Make generated scripts resilient to package-export differences.
        """
        import re as _re
        return _re.sub(
            r"^\s*from\s+backend\s+import\s+([a-zA-Z_][a-zA-Z0-9_]*)\s+as\s+([a-zA-Z_][a-zA-Z0-9_]*)\s*$",
            r"import backend.\1 as \2",
            script_text,
            flags=_re.MULTILINE,
        )

    # ── 1. Find the script to patch ──────────────────────────────────────────
    runs = history_manager.list_runs(session_id)
    target: Optional[Dict] = None
    for r in reversed(runs):
        if not isinstance(r, dict):
            continue
        arts = r.get("produced_artifacts") or []
        if any(
            isinstance(a, dict) and a.get("type") == "script" and a.get("uri", "").endswith("analysis.py")
            for a in arts
        ):
            if target_run in {"latest", "last", "current"} or r.get("run_id") == target_run:
                target = r
                break

    if not target:
        return {
            "status": "needs_inputs",
            "text": (
                "No patchable analysis script is currently linked to this session history. "
                "Provide a specific `target_run` that has a saved script artifact, or request "
                "a deterministic rerun action (for example, update PCA coloring via `bio_rerun`)."
            ),
            "error": "NO_SCRIPT_FOUND",
            "diagnostics": {
                "issue": "missing_script_artifact",
                "target_run": target_run,
                "available_run_ids": [r.get("run_id") for r in runs if isinstance(r, dict)][-5:],
            },
        }

    script_art = next(
        a for a in (target.get("produced_artifacts") or [])
        if isinstance(a, dict) and a.get("type") == "script"
    )
    script_path = Path(script_art["uri"])
    if not script_path.exists():
        return {
            "status": "error",
            "text": f"Script file not found on disk: {script_path}",
            "error": "SCRIPT_FILE_MISSING",
        }

    old_script = script_path.read_text()
    old_run_id = target.get("run_id", "")

    # ── 2. Parse change_request into a parameter patch ────────────────────────
    patch = _parse_change_request(change_request)

    if patch:
        # Simple deterministic patch
        new_script = _se.patch_parameters(old_script, patch)
        patch_description = ", ".join(f"{k}={v!r}" for k, v in patch.items())
    else:
        # Fallback: use the LLM to edit the Parameters block
        new_script = _llm_patch_script(old_script, change_request)
        patch_description = f"(LLM-applied) {change_request}"

    new_script = _normalize_backend_imports(new_script)

    # ── 3. Write new script in a new run directory ────────────────────────────
    import time, re as _re
    new_run_id = f"run-{int(time.time()*1000)}-patch"
    old_run_dir = script_path.parent
    # Keep new run under sessions/<sid>/runs/ alongside the original
    new_run_dir = old_run_dir.parent / new_run_id
    new_run_dir.mkdir(parents=True, exist_ok=True)

    # Update RUN_DIR inside the script — preserve Path() wrapper
    new_script = _re.sub(
        r"^RUN_DIR\s*=.*$",
        f"RUN_DIR          = Path({repr(str(new_run_dir))})",
        new_script,
        flags=_re.MULTILINE,
    )

    new_script_path = new_run_dir / "analysis.py"
    new_script_path.write_text(new_script)

    # ── 4. Execute ────────────────────────────────────────────────────────────
    repo_root = str(Path(__file__).resolve().parent.parent)
    py_path = os.pathsep.join([repo_root, os.getenv("PYTHONPATH", "")]).strip(os.pathsep)
    exec_result = _se.execute(new_script_path, env={"PYTHONPATH": py_path})

    if exec_result.get("status") == "error":
        _err = str(exec_result.get("error") or "")
        _logs = str(exec_result.get("logs") or "")
        _err_l = f"{_err}\n{_logs}".lower()
        if "importerror" in _err_l or "cannot import name" in _err_l or "no module named" in _err_l:
            return {
                "status": "needs_inputs",
                "text": (
                    "Patch execution reached a runtime import dependency issue in the isolated script runner. "
                    "The change request was captured, but this run could not be executed in script mode. "
                    "Retry with an explicit rerun anchor or concrete analysis inputs."
                ),
                "diagnostics": {
                    "issue": "script_runtime_import_failure",
                    "change_request": change_request,
                    "target_run": target.get("run_id"),
                    "script_path": str(new_script_path),
                    "execution_error": _err,
                    "execution_logs_tail": _logs[-1200:],
                },
                "script_path": str(new_script_path),
                "run_id": new_run_id,
                "parent_run_id": old_run_id,
                "links": [
                    {"label": "analysis.py", "url": f"/download/script?path={str(new_script_path)}"},
                ],
            }
        return {
            "status": "error",
            "text": (
                f"Script execution failed after applying: {patch_description}\n\n"
                f"**Error:** {exec_result.get('error', 'unknown')}\n\n"
                f"```\n{exec_result.get('logs', '')[-800:]}\n```"
            ),
            "error": exec_result.get("error"),
        }

    # ── 5. Compute diffs ──────────────────────────────────────────────────────
    param_changes = _se.parameter_diff(old_script, new_script)
    file_diff     = _se.output_diff(old_run_dir, new_run_dir)

    # ── 6. Build visuals from new plot files ──────────────────────────────────
    visuals = []
    summary = exec_result.get("summary", {})
    for plot_key, path_str in (summary.get("plot_paths") or {}).items():
        try:
            import base64
            b64 = base64.b64encode(Path(path_str).read_bytes()).decode()
            label = plot_key.replace("volcano_", "Volcano: ").replace("_", " ").replace("__", " — ").title()
            if "pca" in plot_key.lower():
                label = "PCA — Sample Overview"
            visuals.append({"type": "image_b64", "data": b64, "title": label})
        except Exception:
            pass

    # ── 7. Persist new run in history ─────────────────────────────────────────
    try:
        history_manager.add_run(
            session_id,
            command=f"[patch_and_rerun] {change_request}",
            tool=target.get("tool", "unknown"),
            result=exec_result.get("summary", {}),
            run_id=new_run_id,
            parent_run_id=old_run_id,
            produced_artifacts=[
                {"type": "script", "uri": str(new_script_path), "title": "analysis.py"},
                *[
                    {"type": a["type"], "uri": a["uri"], "title": a["name"]}
                    for a in exec_result.get("artifacts", [])
                    if a.get("type") in ("plot", "table")
                ],
            ],
        )
    except Exception:
        pass

    # ── 8. Assemble response ──────────────────────────────────────────────────
    changes_md = "\n".join(f"- `{c}`" for c in param_changes) or "- (no parameter changes detected)"
    n_new = len(file_diff.get("new_files", []))
    n_changed = len(file_diff.get("changed_files", []))

    text = (
        f"## Analysis Updated\n\n"
        f"Applied: **{change_request}**\n\n"
        f"### What Changed\n{changes_md}\n\n"
        f"### Outputs\n"
        f"- {len(visuals)} plot(s) re-rendered\n"
        f"- {n_changed} file(s) updated, {n_new} new file(s)\n"
        f"- New run: `{new_run_id}` (parent: `{old_run_id[:8] if old_run_id else '?'}…`)\n"
    )

    script_download_url = f"/download/script?path={str(new_script_path)}"
    bundle_url = f"/download/bundle?session_id={session_id}&run_id={new_run_id}"
    return {
        "status":               "success",
        "visualization_type":   "results_viewer",
        "text":                 text,
        "visuals":              visuals,
        "links": [
            {"label": "analysis.py", "url": script_download_url},
            {"label": "bundle.zip",  "url": bundle_url},
        ],
        "parent_run_id":        old_run_id,
        "run_id":               new_run_id,
        "script_path":          str(new_script_path),
        "script_download_url":  script_download_url,
        "bundle_url":           bundle_url,
        "param_changes":        param_changes,
        "file_diff":            file_diff,
    }


# ── helpers for patch_and_rerun ───────────────────────────────────────────────

# Map of keyword patterns → parameter patches. Order matters (first match wins).
_CHANGE_PATTERNS: List[Tuple[List[str], Dict[str, Any]]] = [
    # x-axis / fold change scale
    (["linear fold", "linear scale", "linear x", "x.*linear", "fold change.*linear",
      "x-axis.*linear", "x axis.*linear"],
     {"X_SCALE": "linear"}),
    (["log.*fold", "log2.*fold", "log scale", "log x", "x.*log", "logarithmic"],
     {"X_SCALE": "log2"}),
    # alpha / significance threshold
    (["alpha.*0\\.001", "p.*<.*0\\.001", "significance.*0\\.001"],
     {"ALPHA": 0.001}),
    (["alpha.*0\\.01", "p.*<.*0\\.01", "significance.*0\\.01"],
     {"ALPHA": 0.01}),
    (["alpha.*0\\.1", "p.*<.*0\\.1"],
     {"ALPHA": 0.1}),
    # clustering resolution
    (["resolution.*0\\.3"],  {"RESOLUTION": 0.3}),
    (["resolution.*0\\.5"],  {"RESOLUTION": 0.5}),
    (["resolution.*0\\.8"],  {"RESOLUTION": 0.8}),
    (["resolution.*1\\.0"],  {"RESOLUTION": 1.0}),
    (["more clusters", "higher resolution"], {"RESOLUTION": 0.8}),
    (["fewer clusters", "lower resolution"],  {"RESOLUTION": 0.3}),
]


def _parse_change_request(request: str) -> Dict[str, Any]:
    """Return a simple parameter patch dict, or {} if no pattern matched."""
    import re
    req_lower = request.lower()
    for patterns, patch in _CHANGE_PATTERNS:
        for pat in patterns:
            if re.search(pat, req_lower):
                return dict(patch)

    # Try to extract alpha from free text: "change alpha to 0.005"
    m = re.search(r"\balpha\s*(?:to|=|:)?\s*(0\.\d+)", req_lower)
    if m:
        try:
            return {"ALPHA": float(m.group(1))}
        except ValueError:
            pass

    return {}


def _llm_patch_script(script_text: str, change_request: str) -> str:
    """Ask the LLM to modify the Parameters block.  Falls back to original if unavailable."""
    try:
        from backend.agent import get_llm
        llm = get_llm()
        prompt = (
            "You are editing a Python analysis script. "
            "Modify ONLY the '# ── Parameters ──' block to implement the following change:\n"
            f"{change_request}\n\n"
            "Return the COMPLETE modified script. Do not change anything outside the Parameters block.\n\n"
            f"```python\n{script_text}\n```"
        )
        response = llm.invoke(prompt)
        text = response.content if hasattr(response, "content") else str(response)
        # Extract code block if present
        import re
        m = re.search(r"```python\n(.*?)```", text, re.DOTALL)
        if m:
            return m.group(1)
        return text
    except Exception:
        return script_text  # unchanged; caller will still save + execute


@tool
def bio_rerun(
    session_id: str = "",
    changes: Optional[Dict] = None,
    target_run: str = "latest",
) -> Dict:
    """
    Re-run the most recent bioinformatics analysis with optional parameter overrides.

    Use this tool when the user asks to:
    - Re-run an analysis with different parameters (e.g., "run again with alpha=0.01")
    - Change a parameter and re-run (e.g., "change resolution to 0.8 and re-run")
    - Repeat the last analysis with modifications

    Args:
        session_id: Session identifier.
        changes: Dict of parameter overrides, e.g. {"alpha": 0.01} or {"resolution": 0.8}.
        target_run: Run ID to re-run, or "latest" to use the most recent run.

    Returns:
        Dict with status, text summary, run_id, parent_run_id, delta metrics.
    """
    # NOTE: Actual execution is routed through dispatch_tool("bio_rerun", ...)
    # in main_with_mcp.py which calls BioOrchestrator.rerun().
    # This @tool definition exists for agent routing and documentation.
    return {
        "status": "routed",
        "text": f"Re-running analysis with changes: {changes or {}}",
        "session_id": session_id,
        "changes": changes or {},
        "target_run": target_run,
    }


@tool
def bio_diff_runs(
    session_id: str = "",
    run_id_a: str = "latest",
    run_id_b: str = "prior",
) -> Dict:
    """
    Compare two bioinformatics runs and return a structured diff.

    Use this tool when the user asks to:
    - Compare two runs (e.g., "compare run X and run Y")
    - Show what changed between runs (e.g., "what changed?")
    - Diff two analyses (e.g., "diff runs abc and def")

    Args:
        session_id: Session identifier.
        run_id_a: First run ID (or "latest" for most recent).
        run_id_b: Second run ID (or "prior" for second most recent).

    Returns:
        Dict with status, text summary of changes, parameter diffs, metric deltas.
    """
    # NOTE: Actual execution is routed through dispatch_tool("bio_diff_runs", ...)
    # in main_with_mcp.py which calls BioOrchestrator.diff_runs().
    # This @tool definition exists for agent routing and documentation.
    return {
        "status": "routed",
        "text": f"Comparing runs {run_id_a} and {run_id_b}",
        "session_id": session_id,
        "run_id_a": run_id_a,
        "run_id_b": run_id_b,
    }
