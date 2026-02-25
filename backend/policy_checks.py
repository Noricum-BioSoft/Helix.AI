"""
Policy Checks Module

Enforces "must not" rules from agents/agent-responsibilities.md at the contract level.
Each agent has clear boundaries on what they may NOT do. This module provides
validation functions that raise PolicyViolationError when those boundaries are crossed.

Usage:
    from backend.policy_checks import (
        check_planner_output,
        check_plan_not_mutated,
        check_codegen_output,
        check_visualizer_output,
    )
    
    # After planner produces a plan
    check_planner_output(plan)
    
    # Before/after infrastructure agent runs
    before_steps = serialize_plan_steps(plan.steps)
    # ... infra agent runs ...
    check_plan_not_mutated(before_steps, plan.steps, "InfrastructureExpert")
    
    # After code generator runs
    check_codegen_output(original_plan, execution_spec)
    
    # After visualizer runs
    check_visualizer_output(visualization_artifacts)
"""

import json
import logging
from typing import Any, Dict, List, Optional, Union
from pydantic import BaseModel

logger = logging.getLogger(__name__)


class PolicyViolationError(Exception):
    """
    Raised when an agent violates its responsibilities contract.
    
    Per agent-responsibilities.md, each agent has clear "must not" rules.
    This exception is raised when those rules are violated, with a message
    explaining what was found and how to fix it.
    """
    pass


# =============================================================================
# Infrastructure Keywords Detection
# =============================================================================

# Infrastructure-related keywords that should NOT appear in planner output
INFRA_KEYWORDS = {
    # Environment types
    "ec2", "emr", "batch", "lambda", "local",
    "aws", "cloud", "cluster",
    
    # Instance/resource types
    "instance", "instance_type", "instance-type",
    "m5.xlarge", "m5.2xlarge", "m5.4xlarge", "m5.8xlarge",
    "r5.xlarge", "r5.2xlarge", "r5.4xlarge", "r5.8xlarge",
    "c5.xlarge", "c5.2xlarge", "c5.4xlarge", "c5.9xlarge",
    
    # EMR-specific
    "emr_steps", "emr-steps", "spark", "hadoop", "yarn",
    "master_node", "core_node", "task_node",
    
    # Container/deployment
    "container", "docker", "image", "ecs", "fargate",
    
    # Infrastructure config
    "vpc", "subnet", "security_group", "iam_role",
    "ebs", "volume", "storage_class",
    
    # Execution environment hints
    "executor", "driver", "worker", "node_count",
    "cpu", "memory", "vcpu", "ram", "disk",
}


def _scan_for_infra_keywords(data: Any, path: str = "") -> List[str]:
    """
    Recursively scan a data structure for infrastructure keywords.
    
    Args:
        data: Data to scan (dict, list, str, or primitive)
        path: Current path in the data structure (for error reporting)
        
    Returns:
        List of violation messages describing what was found and where
    """
    violations = []
    
    if isinstance(data, dict):
        for key, value in data.items():
            current_path = f"{path}.{key}" if path else key
            
            # Check key name
            key_lower = key.lower()
            for keyword in INFRA_KEYWORDS:
                if keyword in key_lower:
                    violations.append(
                        f"Infrastructure keyword '{keyword}' found in field name: {current_path}"
                    )
            
            # Recurse into value
            violations.extend(_scan_for_infra_keywords(value, current_path))
    
    elif isinstance(data, list):
        for i, item in enumerate(data):
            current_path = f"{path}[{i}]"
            violations.extend(_scan_for_infra_keywords(item, current_path))
    
    elif isinstance(data, str):
        # Check string values (case-insensitive)
        data_lower = data.lower()
        for keyword in INFRA_KEYWORDS:
            if keyword in data_lower and len(keyword) > 3:  # Skip very short keywords to reduce false positives
                violations.append(
                    f"Infrastructure keyword '{keyword}' found in string at {path}: '{data[:100]}...'"
                )
    
    return violations


# =============================================================================
# 1. Planner Output Validation
# =============================================================================

def check_planner_output(plan: Union[Dict, BaseModel]) -> None:
    """
    Validate that Bioinformatics Executor (Planner) output does NOT include
    infrastructure decisions.
    
    Per agent-responsibilities.md:
        "Bioinformatics Executor (Planner) - Owner of: creating a workflow plan
        
        Must not:
        - choose infrastructure (no EC2/EMR/Batch/Lambda decisions)
        - generate new tools/code (unless explicitly delegated to Code Generator)
        - execute jobs, upload data, or perform side effects"
    
    Checks:
    - Plan must not include environment selection (EC2/EMR/Batch/Lambda/Local)
    - Plan must not include instance types, cluster configs, EMR steps
    - Plan must not include container/deployment details
    
    Args:
        plan: WorkflowPlan or dict representation of a plan
        
    Raises:
        PolicyViolationError: If infrastructure-related fields/keywords are found
    """
    # Convert Pydantic model to dict if needed
    if isinstance(plan, BaseModel):
        plan_dict = plan.model_dump() if hasattr(plan, 'model_dump') else plan.dict()
    else:
        plan_dict = plan
    
    # Scan for infrastructure keywords
    violations = _scan_for_infra_keywords(plan_dict, path="plan")
    
    # Check for explicit environment field
    if "environment" in plan_dict:
        violations.append(
            "Field 'environment' found in plan. "
            "Planner must not choose infrastructure; that's the Infrastructure Expert's role."
        )
    
    if "recommended_environment" in plan_dict:
        violations.append(
            "Field 'recommended_environment' found in plan. "
            "Planner must not choose infrastructure; that's the Infrastructure Expert's role."
        )
    
    if "target" in plan_dict and isinstance(plan_dict["target"], str):
        target_lower = plan_dict["target"].lower()
        if target_lower in {"ec2", "emr", "batch", "lambda", "local"}:
            violations.append(
                f"Execution target '{plan_dict['target']}' found in plan. "
                "Planner must not choose infrastructure; that's the Infrastructure Expert's role."
            )
    
    # If violations found, raise error with actionable message
    if violations:
        message = (
            "PolicyViolation: Bioinformatics Executor (Planner) must not choose infrastructure.\n\n"
            "The Planner's role is to create a workflow plan (steps, tools, inputs/outputs) "
            "without making infrastructure decisions. Infrastructure selection is owned by the "
            "Infrastructure Expert agent.\n\n"
            "Violations found:\n" + "\n".join(f"  - {v}" for v in violations) +
            "\n\nHow to fix:\n"
            "  - Remove infrastructure-related fields from the plan\n"
            "  - Keep the plan focused on WHAT to do (steps/tools), not WHERE to do it (infra)\n"
            "  - Let the Infrastructure Expert decide environment, instance types, etc."
        )
        raise PolicyViolationError(message)
    
    logger.debug("[PolicyCheck] ✓ Planner output validated: no infrastructure decisions found")


# =============================================================================
# 2. Plan Mutation Detection
# =============================================================================

def serialize_plan_steps(steps: Union[List[Dict], List[BaseModel]]) -> str:
    """
    Create a stable canonical serialization of plan steps for comparison.
    
    Args:
        steps: List of plan steps (dicts or Pydantic models)
        
    Returns:
        JSON string with deterministic ordering
    """
    # Convert Pydantic models to dicts
    if steps and isinstance(steps[0], BaseModel):
        steps_dicts = [step.model_dump() if hasattr(step, 'model_dump') else step.dict() for step in steps]
    else:
        steps_dicts = steps
    
    # Serialize with sorted keys for stability
    return json.dumps(steps_dicts, sort_keys=True, indent=2)


def check_plan_not_mutated(
    before_steps: Union[str, List],
    after_steps: Union[str, List],
    actor: str
) -> None:
    """
    Validate that an agent did NOT mutate the workflow plan steps.
    
    Per agent-responsibilities.md:
        "3) No-Mutation Rule
        Agents do not mutate upstream contracts. They may:
        - add annotations (warnings, assumptions, constraints)
        - produce a new downstream contract"
    
    Specifically for Infrastructure Expert:
        "Must not:
        - rewrite workflow steps"
    
    Args:
        before_steps: Plan steps before agent ran (serialized string or list)
        after_steps: Plan steps after agent ran (serialized string or list)
        actor: Name of agent that ran (for error message)
        
    Raises:
        PolicyViolationError: If steps were modified
    """
    # Normalize inputs to strings
    if not isinstance(before_steps, str):
        before_str = serialize_plan_steps(before_steps)
    else:
        before_str = before_steps
    
    if not isinstance(after_steps, str):
        after_str = serialize_plan_steps(after_steps)
    else:
        after_str = after_steps
    
    # Compare
    if before_str != after_str:
        # Compute a simple diff for error message
        before_lines = before_str.split('\n')
        after_lines = after_str.split('\n')
        
        diff_lines = []
        for i, (b, a) in enumerate(zip(before_lines, after_lines)):
            if b != a:
                diff_lines.append(f"  Line {i}: '{b}' → '{a}'")
        
        # Include extra lines if lengths differ
        if len(before_lines) != len(after_lines):
            diff_lines.append(f"  Length changed: {len(before_lines)} → {len(after_lines)} lines")
        
        message = (
            f"PolicyViolation: {actor} must not mutate workflow plan steps.\n\n"
            f"Per the No-Mutation Rule in agent-responsibilities.md, agents may only:\n"
            f"  - Add annotations (warnings, assumptions, constraints)\n"
            f"  - Produce a new downstream contract\n\n"
            f"But {actor} modified the plan steps:\n" +
            "\n".join(diff_lines[:10]) +  # Limit diff to first 10 lines
            ("\n  ... (more differences)" if len(diff_lines) > 10 else "") +
            "\n\nHow to fix:\n"
            f"  - {actor} should return its output separately (e.g., InfraDecision)\n"
            f"  - {actor} should not modify the input WorkflowPlan.steps\n"
            f"  - Use constraints/advice fields to suggest non-invasive changes"
        )
        raise PolicyViolationError(message)
    
    logger.debug(f"[PolicyCheck] ✓ Plan steps unchanged after {actor}")


# =============================================================================
# 3. Code Generator Output Validation
# =============================================================================

def check_codegen_output(
    original_plan: Union[Dict, BaseModel],
    execution_spec: Union[Dict, BaseModel, None]
) -> None:
    """
    Validate that Code Generator output does NOT change scientific intent.
    
    Per agent-responsibilities.md:
        "Code Generator (Tooling/Packaging Agent) - Owner of: filling tool gaps
        and producing a runnable execution specification
        
        Must not:
        - change scientific intent of the workflow
        - choose infrastructure
        - submit jobs or execute anything"
    
    Checks:
    - If codegen returns a plan, verify steps match original (allow wrapper steps)
    - Execution spec should not include infrastructure decisions (checked separately)
    
    Args:
        original_plan: Original WorkflowPlan from Planner
        execution_spec: ExecutionSpec or dict returned by CodeGen (may be None)
        
    Raises:
        PolicyViolationError: If scientific intent was changed
    """
    if execution_spec is None:
        logger.debug("[PolicyCheck] ✓ CodeGen output is None, nothing to validate")
        return
    
    # Convert to dicts
    if isinstance(original_plan, BaseModel):
        original_dict = original_plan.model_dump() if hasattr(original_plan, 'model_dump') else original_plan.dict()
    else:
        original_dict = original_plan
    
    if isinstance(execution_spec, BaseModel):
        spec_dict = execution_spec.model_dump() if hasattr(execution_spec, 'model_dump') else execution_spec.dict()
    else:
        spec_dict = execution_spec
    
    # Check if execution_spec has 'steps' field (should not rewrite plan)
    if "steps" in spec_dict:
        # This is suspicious - ExecutionSpec should not include steps
        # Compare to original plan steps
        original_steps = original_dict.get("steps", [])
        spec_steps = spec_dict.get("steps", [])
        
        # Serialize for comparison
        original_str = serialize_plan_steps(original_steps)
        spec_str = serialize_plan_steps(spec_steps)
        
        if original_str != spec_str:
            message = (
                "PolicyViolation: Code Generator must not change scientific intent.\n\n"
                "The Code Generator's role is to produce a runnable ExecutionSpec "
                "(container, commands, entrypoints) without rewriting the workflow steps. "
                "Scientific intent is owned by the Planner.\n\n"
                "Detected step rewriting in ExecutionSpec.\n\n"
                "How to fix:\n"
                "  - CodeGen should output ExecutionSpec (container/commands)\n"
                "  - CodeGen should not rewrite WorkflowPlan.steps\n"
                "  - If wrapper steps are needed, add them as staging/setup steps "
                "    with clear provenance, not by replacing scientific steps"
            )
            raise PolicyViolationError(message)
    
    # Check for infrastructure decisions in execution_spec (should be pre-decided by InfraExpert)
    # Allow 'target' field since that's expected in ExecutionSpec (from InfraDecision)
    # But disallow instance types, cluster configs, etc.
    violations = []
    for field in ["instance_type", "cluster_config", "emr_steps", "master_node", "core_node"]:
        if field in spec_dict:
            violations.append(
                f"Field '{field}' found in ExecutionSpec. "
                "CodeGen must not choose infrastructure; that's the Infrastructure Expert's role."
            )
    
    if violations:
        message = (
            "PolicyViolation: Code Generator must not choose infrastructure.\n\n"
            "Violations found:\n" + "\n".join(f"  - {v}" for v in violations) +
            "\n\nHow to fix:\n"
            "  - Remove infrastructure decision fields from ExecutionSpec\n"
            "  - Use the target from InfraDecision (pre-decided by Infrastructure Expert)"
        )
        raise PolicyViolationError(message)
    
    logger.debug("[PolicyCheck] ✓ CodeGen output validated: no scientific intent changes")


# =============================================================================
# 4. Visualizer Output Validation
# =============================================================================

def check_visualizer_output(output: Union[Dict, BaseModel, None]) -> None:
    """
    Validate that Data Visualizer output does NOT introduce workflow steps or
    infrastructure changes.
    
    Per agent-responsibilities.md:
        "Data Visualizer - Owner of: visualization/reporting of results
        
        Must not:
        - redesign workflow
        - choose infrastructure
        - generate execution code
        - execute jobs"
    
    Checks:
    - Output must not include 'steps' or 'workflow' fields
    - Output must not include infrastructure fields
    - Output should be VisualizationArtifacts (plots/tables/reports)
    
    Args:
        output: VisualizationArtifacts or dict returned by Visualizer
        
    Raises:
        PolicyViolationError: If workflow or infra changes are found
    """
    if output is None:
        logger.debug("[PolicyCheck] ✓ Visualizer output is None, nothing to validate")
        return
    
    # Convert to dict
    if isinstance(output, BaseModel):
        output_dict = output.model_dump() if hasattr(output, 'model_dump') else output.dict()
    else:
        output_dict = output
    
    violations = []
    
    # Check for workflow-related fields
    if "steps" in output_dict:
        violations.append(
            "Field 'steps' found in Visualizer output. "
            "Visualizer must not redesign workflow; it should only consume ExecutionResult."
        )
    
    if "workflow" in output_dict and isinstance(output_dict["workflow"], dict):
        if "steps" in output_dict["workflow"]:
            violations.append(
                "Field 'workflow.steps' found in Visualizer output. "
                "Visualizer must not redesign workflow; it should only consume ExecutionResult."
            )
    
    # Check for infrastructure fields
    for field in ["environment", "recommended_environment", "target", "infrastructure"]:
        if field in output_dict:
            violations.append(
                f"Field '{field}' found in Visualizer output. "
                "Visualizer must not choose infrastructure; it should only visualize results."
            )
    
    # Check for execution code
    for field in ["command", "commands", "script", "entrypoint"]:
        if field in output_dict and output_dict[field]:
            violations.append(
                f"Field '{field}' found in Visualizer output. "
                "Visualizer must not generate execution code; that's the Code Generator's role."
            )
    
    if violations:
        message = (
            "PolicyViolation: Data Visualizer must not introduce workflow steps or infrastructure changes.\n\n"
            "Violations found:\n" + "\n".join(f"  - {v}" for v in violations) +
            "\n\nHow to fix:\n"
            "  - Visualizer should only consume ExecutionResult (artifacts, logs, metrics)\n"
            "  - Visualizer should produce VisualizationArtifacts (plots/tables/reports)\n"
            "  - Remove workflow and infrastructure fields from output"
        )
        raise PolicyViolationError(message)
    
    logger.debug("[PolicyCheck] ✓ Visualizer output validated: no workflow/infra changes")


# =============================================================================
# Convenience function for full plan validation
# =============================================================================

def validate_workflow_contracts(
    planner_output: Union[Dict, BaseModel],
    infra_input_steps: Union[List, str],
    infra_output_steps: Union[List, str],
    codegen_execution_spec: Optional[Union[Dict, BaseModel]] = None,
    visualizer_output: Optional[Union[Dict, BaseModel]] = None,
) -> None:
    """
    Convenience function to validate all workflow contracts in one call.
    
    This is useful for end-to-end validation in orchestrator code.
    
    Args:
        planner_output: WorkflowPlan from Planner
        infra_input_steps: Plan steps before Infrastructure Expert ran
        infra_output_steps: Plan steps after Infrastructure Expert ran
        codegen_execution_spec: Optional ExecutionSpec from CodeGen
        visualizer_output: Optional VisualizationArtifacts from Visualizer
        
    Raises:
        PolicyViolationError: If any validation fails
    """
    check_planner_output(planner_output)
    check_plan_not_mutated(infra_input_steps, infra_output_steps, "InfrastructureExpert")
    
    if codegen_execution_spec is not None:
        check_codegen_output(planner_output, codegen_execution_spec)
    
    if visualizer_output is not None:
        check_visualizer_output(visualizer_output)
    
    logger.info("[PolicyCheck] ✓ All workflow contracts validated")
