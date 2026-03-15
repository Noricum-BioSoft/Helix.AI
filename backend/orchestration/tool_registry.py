from __future__ import annotations

from pathlib import Path
from typing import Any, Awaitable, Callable, Dict, Tuple


async def _call_tool_obj(tool_obj: Any, arguments: Dict[str, Any]) -> Dict[str, Any]:
    if tool_obj is None:
        raise ValueError("Tool object is missing")
    if hasattr(tool_obj, "invoke"):
        return tool_obj.invoke(arguments)
    if hasattr(tool_obj, "func"):
        return tool_obj.func(**arguments)
    if callable(tool_obj):
        return tool_obj(**arguments)
    raise ValueError("Tool object is not callable")


async def _handle_toolbox_inventory(arguments: Dict[str, Any]) -> Dict[str, Any]:
    from tool_inventory import build_toolbox_inventory, format_toolbox_inventory_markdown

    inv = build_toolbox_inventory()
    return {"status": "success", "text": format_toolbox_inventory_markdown(inv), "result": inv}


async def _handle_visualize_job_results(arguments: Dict[str, Any]) -> Dict[str, Any]:
    from backend.job_manager import get_job_manager
    import re as _re

    jm = get_job_manager()
    requested_job_id = arguments.get("job_id")
    original_command = str(arguments.get("original_command") or "")
    session_id = arguments.get("session_id")

    if not requested_job_id and original_command:
        m = _re.search(
            r"\b[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}\b",
            original_command.lower(),
        )
        if m:
            requested_job_id = m.group(0)

    if not requested_job_id:
        candidates = [
            job for job in jm.jobs.values() if isinstance(job, dict) and job.get("status") == "completed"
        ]
        if session_id:
            scoped = [j for j in candidates if j.get("session_id") == session_id]
            if scoped:
                candidates = scoped
        candidates.sort(key=lambda j: str(j.get("updated_at") or j.get("completed_at") or ""), reverse=True)
        if candidates:
            requested_job_id = candidates[0].get("job_id")
    if not requested_job_id:
        return {
            "status": "error",
            "visualization_type": "text",
            "text": "No job ID was provided and no completed jobs were found to visualize.",
            "result": {},
        }

    try:
        job = jm.get_job_status(requested_job_id)
    except Exception as e:
        return {
            "status": "error",
            "visualization_type": "text",
            "text": f"Could not find job `{requested_job_id}`: {e}",
            "result": {"job_id": requested_job_id},
        }
    if job.get("status") != "completed":
        return {
            "status": "success",
            "visualization_type": "text",
            "text": f"Job `{requested_job_id}` is currently `{job.get('status')}`. Please wait for completion.",
            "result": {"job_id": requested_job_id, "job_status": job.get("status"), "job": job},
        }

    results_info = jm.get_job_results(requested_job_id)
    tool_result = results_info.get("result") if isinstance(results_info, dict) else {}
    if not isinstance(tool_result, dict):
        tool_result = {"value": tool_result}
    visuals = tool_result.get("visuals", []) if isinstance(tool_result.get("visuals"), list) else []
    links = []
    if results_info.get("session_html_path"):
        links.append({"title": "Session HTML Results", "url": str(results_info.get("session_html_path"))})
    if results_info.get("results_path"):
        links.append({"title": "Results JSON", "url": str(results_info.get("results_path"))})
    return {
        "status": "success",
        "visualization_type": "results_viewer" if (visuals or tool_result) else "text",
        "text": tool_result.get("text") or f"Loaded visualization payload for job `{requested_job_id}`.",
        "visuals": visuals,
        "links": links,
        "result": {"job_id": requested_job_id, "job": job, "results": results_info, "tool_result": tool_result},
    }


async def dispatch_via_registry(
    tool_name: str,
    arguments: Dict[str, Any],
    *,
    needs_inputs_builder: Callable[[str, Dict[str, Any]], Dict[str, Any]],
) -> Tuple[bool, Dict[str, Any]]:
    """
    Registry-first dispatch for common tools. Returns (handled, result).
    """
    if arguments.get("needs_inputs"):
        return True, needs_inputs_builder(tool_name, arguments)

    if tool_name == "toolbox_inventory":
        return True, await _handle_toolbox_inventory(arguments)
    if tool_name == "visualize_job_results":
        return True, await _handle_visualize_job_results(arguments)

    if tool_name in {"ds_run_analysis", "ds_reproduce_run", "ds_diff_runs", "ds_list_runs"}:
        from backend import agent_tools as _agent_tools

        return True, await _call_tool_obj(getattr(_agent_tools, tool_name, None), arguments)

    if tool_name in {
        "local_demo_scatter_plot",
        "local_demo_plot_script",
        "local_update_scatter_x_scale",
        "local_edit_visualization",
        "local_edit_and_rerun_script",
        "session_run_io_summary",
        "patch_and_rerun",
        "s3_browse_results",
    }:
        from backend import agent_tools as _agent_tools

        return True, await _call_tool_obj(getattr(_agent_tools, tool_name, None), arguments)

    return False, {}

