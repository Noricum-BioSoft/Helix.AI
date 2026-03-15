from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

from backend.action_plan import infer_action_type, map_action_to_tool
from backend.orchestration.approval_policy import READ_ONLY_ROUTER_TOOLS
from backend.plan_binding import validate_tool_bindings


@dataclass
class ValidationResult:
    ok: bool
    tool_name: str
    arguments: Dict[str, Any]
    result: Optional[Dict[str, Any]] = None


def normalize_tool_selection(
    *,
    command: str,
    tool_name: str,
    parameters: Dict[str, Any],
    session_context: Dict[str, Any],
) -> ValidationResult:
    """
    Normalize tool choice using action semantics and available artifact context.
    Prevents modality drift (e.g. PCA recolor routed to single-cell executor).
    """
    tool = tool_name
    params = dict(parameters or {})
    action = infer_action_type(command, tool_name)
    runs = (session_context or {}).get("runs") or []
    artifacts_raw = (session_context or {}).get("artifacts") or {}
    artifacts = list(artifacts_raw.values()) if isinstance(artifacts_raw, dict) else (
        artifacts_raw if isinstance(artifacts_raw, list) else []
    )
    has_bulk_context = any(isinstance(r, dict) and r.get("tool") == "bulk_rnaseq_analysis" for r in runs)
    has_plot_artifacts = any(
        isinstance(a, dict)
        and any(
            token in str(a.get("title", "")).lower() or token in str(a.get("artifact_kind", "")).lower()
            for token in ("pca", "volcano", "heatmap", "plot", "figure")
        )
        for a in artifacts
    )
    cmd = (command or "").lower()

    if tool_name in {"single_cell_analysis", "handle_natural_command"} and has_bulk_context:
        if action in {"generate_plot", "compare_versions", "rerun_downstream_steps"}:
            mapped = map_action_to_tool(
                action,
                None,
                {**params, "command": command, "original_command": command},
            )
            if mapped:
                tool = mapped
    if action == "run_enrichment" and tool_name in {
        "lookup_go_term",
        "bulk_rnaseq_analysis",
        "single_cell_analysis",
        "handle_natural_command",
    }:
        tool = "go_enrichment_analysis"
        params.setdefault("source_selector", "current DEG results")
        params.setdefault("session_id", (session_context or {}).get("session_id", ""))
        params.setdefault("command", command)
    if action == "generate_plot" and has_plot_artifacts:
        # Prefer artifact-aware update tools for plot iterations, independent of assay tool family.
        if "pca" in cmd and ("batch" in cmd or "sex" in cmd):
            tool = "bio_rerun"
        elif tool_name in {"single_cell_analysis", "bulk_rnaseq_analysis", "handle_natural_command"}:
            tool = "patch_and_rerun"
    return ValidationResult(ok=True, tool_name=tool, arguments=params)


def preflight_tool_bindings(
    tool_name: str,
    parameters: Optional[Dict[str, Any]],
    needs_inputs_builder,
) -> Optional[Dict[str, Any]]:
    """
    Validate non-plan tool arguments before execution.
    Returns a needs-inputs response payload when required bindings are missing.
    """
    if not tool_name or tool_name in READ_ONLY_ROUTER_TOOLS:
        return None
    params = parameters if isinstance(parameters, dict) else {}
    if tool_name == "patch_and_rerun":
        # Agent tool-mapped flows may omit `change_request`; infer it from command text.
        inferred_change = str(
            params.get("change_request")
            or params.get("original_command")
            or params.get("command")
            or ""
        ).strip()
        if inferred_change:
            params["change_request"] = inferred_change
    issues = validate_tool_bindings(tool_name, params)
    if not issues:
        return None
    payload = needs_inputs_builder(tool_name, {**params, "needs_inputs": True})
    payload["binding_diagnostics"] = {
        "status": "error",
        "error_type": "artifact_binding_failed",
        "issues": issues,
    }
    return payload


def normalize_fastqc_parameters(
    *,
    tool_name: str,
    arguments: Dict[str, Any],
    command: str,
    session_context: Dict[str, Any],
) -> ValidationResult:
    """
    Apply session-aware FastQC input extraction and input-shape validation.
    """
    if tool_name != "fastqc_quality_analysis":
        return ValidationResult(ok=True, tool_name=tool_name, arguments=arguments)

    params = dict(arguments or {})
    if not params.get("input_r1") or not params.get("input_r2"):
        from backend.session_param_extractor import get_fastqc_inputs_from_session

        sess_r1, sess_r2, sess_out, sess_err = get_fastqc_inputs_from_session(session_context or {}, command)
        if sess_err:
            return ValidationResult(
                ok=False,
                tool_name=tool_name,
                arguments=params,
                result={"status": "success", "text": sess_err, "validation_message": True},
            )
        if sess_r1 and sess_r2:
            params["input_r1"] = params.get("input_r1") or sess_r1
            params["input_r2"] = params.get("input_r2") or sess_r2
            if sess_out and not params.get("output"):
                params["output"] = sess_out

    from backend.input_validation import validate_fastqc_inputs

    val = validate_fastqc_inputs(params.get("input_r1", ""), params.get("input_r2", ""), command)
    if not val.valid:
        msg = val.message or "Invalid FastQC inputs."
        if val.suggestion:
            msg += "\n\n**What to do:**\n" + val.suggestion
        return ValidationResult(
            ok=False,
            tool_name=tool_name,
            arguments=params,
            result={"status": "success", "text": msg, "validation_message": True},
        )

    return ValidationResult(ok=True, tool_name=tool_name, arguments=params)

