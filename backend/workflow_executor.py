"""
Workflow Executor - Executes multi-step workflows with tool chaining.

This module coordinates the execution of multi-step bioinformatics workflows,
handling tool availability checking, dynamic code generation, and step-to-step
data passing.

Architecture:
    User Command
        ↓
    workflow_planner.plan_workflow() → WorkflowPlan
        ↓
    WorkflowExecutor.execute_workflow()
        ↓
    For each step:
        - Check if tool exists
        - Generate tool if missing (via Code Generator Agent)
        - Execute tool
        - Pass outputs to next step
        ↓
    Return complete workflow results
"""

import logging
import asyncio
import time
from typing import Dict, Any, List, Optional
from pathlib import Path

from backend.contracts.workflow_plan import WorkflowPlan, OperationSpec
from backend.workflow_planner_agent import plan_workflow

logger = logging.getLogger(__name__)


class WorkflowExecutionResult:
    """Result of a workflow execution."""
    
    def __init__(self):
        self.workflow_id: Optional[str] = None
        self.status: str = "pending"  # pending, running, completed, failed
        self.start_time: float = time.time()
        self.end_time: Optional[float] = None
        self.steps_completed: int = 0
        self.total_steps: int = 0
        self.step_results: List[Dict[str, Any]] = []
        self.errors: List[str] = []
        self.outputs: Dict[str, Any] = {}
        
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "workflow_id": self.workflow_id,
            "status": self.status,
            "start_time": self.start_time,
            "end_time": self.end_time,
            "duration_seconds": (self.end_time - self.start_time) if self.end_time else None,
            "steps_completed": self.steps_completed,
            "total_steps": self.total_steps,
            "step_results": self.step_results,
            "errors": self.errors,
            "outputs": self.outputs
        }


class WorkflowExecutor:
    """
    Executes multi-step bioinformatics workflows.
    
    Handles:
    - Tool availability checking
    - Dynamic code generation for missing tools
    - Sequential step execution
    - Output-to-input data passing
    - Error recovery and partial results
    
    Usage:
        executor = WorkflowExecutor()
        result = await executor.execute_workflow(workflow_plan, session_id)
    """
    
    def __init__(self):
        """Initialize workflow executor."""
        self.tool_registry: Dict[str, bool] = {}  # Cache of tool availability
        
    def _is_multi_step_command(self, command: str) -> bool:
        """
        Detect if command requires multi-step EXECUTION workflow.

        Only triggers on explicit pipeline/execution indicators, NOT on
        analytical objectives that happen to be written as numbered lists
        (e.g. "1. Perform QC  2. Normalize  3. Cluster" in a scRNA-seq
        study design — those are objectives, not execution steps).
        """
        command_lower = command.lower()

        # Strong, unambiguous sequential-execution indicators
        strong_indicators = [
            '→', '->',
            'pipeline steps',           # Explicit "Pipeline Steps" section header
            'preprocessing pipeline',
            'full preprocessing',
            'then trim', 'then merge', 'then run', 'then align',
            'then perform fastqc', 'then perform qc',
            'followed by',
            'run fastqc', 'run trimming', 'run merging',
        ]

        if any(kw in command_lower for kw in strong_indicators):
            return True

        # Explicit numbered *execution* steps — only when they appear in a
        # "Step N:" / "step N." pattern AND are paired with pipeline verbs.
        import re as _re
        has_numbered_steps = bool(_re.search(r'\bstep\s+\d+\s*[:.)]', command_lower))
        pipeline_verbs = ['run ', 'trim', 'merg', 'align', 'assemble', 'filter', 'qc ', 'map ']
        if has_numbered_steps and any(v in command_lower for v in pipeline_verbs):
            return True

        # Multiple specific bioinformatics *execution* tools mentioned together
        tool_keywords = [
            'fastqc', 'star ', 'hisat', ' bwa ', 'featurecounts',
            'trimmomatic', 'trim galore', 'cutadapt', 'flash ', 'panda',
        ]
        tools_mentioned = sum(1 for tool in tool_keywords if tool in command_lower)
        if tools_mentioned >= 2:
            return True

        return False
    
    async def execute_workflow(
        self,
        workflow_plan: WorkflowPlan,
        session_id: str,
        session_context: Optional[Dict] = None
    ) -> WorkflowExecutionResult:
        """
        Execute a multi-step workflow.
        
        Args:
            workflow_plan: Workflow plan with operations to execute
            session_id: Session ID for tracking
            session_context: Optional session context
        
        Returns:
            WorkflowExecutionResult with complete execution trace
        """
        result = WorkflowExecutionResult()
        result.workflow_id = f"wf_{session_id}_{int(time.time())}"
        result.total_steps = len(workflow_plan.operations)
        result.status = "running"
        
        logger.info(
            f"[WorkflowExecutor] Starting workflow {result.workflow_id}: "
            f"{len(workflow_plan.operations)} steps"
        )
        
        previous_outputs: Dict[str, Any] = {}
        
        # Store workflow plan in context for first step parameter extraction
        execution_context = dict(session_context or {})
        execution_context["workflow_plan"] = workflow_plan
        
        try:
            for step_idx, operation in enumerate(workflow_plan.operations, 1):
                logger.info(
                    f"[WorkflowExecutor] Step {step_idx}/{result.total_steps}: "
                    f"{operation.operation_name} ({operation.tool_name})"
                )
                
                # Execute step
                step_result = await self._execute_step(
                    operation=operation,
                    step_index=step_idx,
                    previous_outputs=previous_outputs,
                    session_id=session_id,
                    session_context=execution_context
                )
                
                # Record step result
                result.step_results.append(step_result)
                result.steps_completed = step_idx
                
                # Check if step failed
                if step_result.get("status") != "completed":
                    error_msg = step_result.get("error", "Unknown error")
                    result.errors.append(
                        f"Step {step_idx} ({operation.operation_name}) failed: {error_msg}"
                    )
                    result.status = "failed"
                    logger.error(f"[WorkflowExecutor] Step {step_idx} failed: {error_msg}")
                    break
                
                # Extract outputs for next step
                if "outputs" in step_result:
                    previous_outputs.update(step_result["outputs"])
                
                logger.info(f"[WorkflowExecutor] Step {step_idx} completed successfully")
            
            # Mark as completed if all steps succeeded
            if result.steps_completed == result.total_steps and result.status == "running":
                result.status = "completed"
                result.outputs = previous_outputs
                logger.info(f"[WorkflowExecutor] Workflow {result.workflow_id} completed successfully")
            
        except Exception as e:
            logger.error(f"[WorkflowExecutor] Workflow {result.workflow_id} failed: {e}", exc_info=True)
            result.status = "failed"
            result.errors.append(f"Workflow execution error: {str(e)}")
        
        finally:
            result.end_time = time.time()
        
        return result
    
    async def _execute_step(
        self,
        operation: OperationSpec,
        step_index: int,
        previous_outputs: Dict[str, Any],
        session_id: str,
        session_context: Dict
    ) -> Dict[str, Any]:
        """
        Execute a single workflow step.
        
        Returns:
            Dict with:
                - status: "completed" or "failed"
                - tool_name: Name of tool executed
                - execution_mode: "local" or "emr"
                - outputs: Dict of output files/data
                - error: Error message if failed
        """
        step_result = {
            "step_index": step_index,
            "operation_name": operation.operation_name,
            "tool_name": operation.tool_name,
            "status": "pending",
            "start_time": time.time()
        }
        
        try:
            # Step 1: Check if tool exists
            tool_available = await self._check_tool_availability(operation.tool_name)
            
            # Step 2: Prepare execution parameters
            execution_params = self._prepare_execution_params(
                operation,
                previous_outputs,
                session_context
            )
            
            # Step 3: Execute tool (generate if missing)
            if not tool_available:
                logger.info(f"[WorkflowExecutor] Tool '{operation.tool_name}' not found, generating and executing...")
                step_result["tool_generated"] = True
                
                # Generate and execute in one step
                execution_result = await self._generate_and_execute_tool(
                    operation=operation,
                    parameters=execution_params,
                    session_id=session_id,
                    session_context=session_context
                )
            else:
                logger.info(f"[WorkflowExecutor] Tool '{operation.tool_name}' found, executing...")
                step_result["tool_generated"] = False
                
                # Execute existing tool
                execution_result = await self._execute_tool(
                    tool_name=operation.tool_name,
                    parameters=execution_params,
                    session_id=session_id
                )
            
            # Step 4: Extract outputs
            step_result["status"] = execution_result.get("status", "failed")
            step_result["execution_mode"] = execution_result.get("execution_mode", "unknown")
            step_result["message"] = execution_result.get("message")
            step_result["outputs"] = execution_result.get("outputs", {})
            if "raw_result" in execution_result:
                step_result["raw_result"] = execution_result.get("raw_result")
            step_result["end_time"] = time.time()
            step_result["duration_seconds"] = step_result["end_time"] - step_result["start_time"]
            
            if step_result["status"] != "completed":
                step_result["error"] = execution_result.get("message", "Execution failed")
            
            return step_result
            
        except Exception as e:
            logger.error(f"[WorkflowExecutor] Step {step_index} error: {e}", exc_info=True)
            step_result["status"] = "failed"
            step_result["error"] = str(e)
            step_result["end_time"] = time.time()
            return step_result
    
    async def _check_tool_availability(self, tool_name: str) -> bool:
        """
        Check if a tool is available in the system.
        
        Returns:
            True if tool exists, False otherwise
        """
        # Check cache first
        if tool_name in self.tool_registry:
            return self.tool_registry[tool_name]
        
        # Check if tool module exists in agent_tools
        try:
            from backend import agent_tools
            
            # Common tool name patterns
            possible_names = [
                tool_name.lower(),
                tool_name.lower().replace('-', '_'),
                tool_name.lower().replace(' ', '_'),
                f"{tool_name.lower()}_analysis",
                f"{tool_name.lower()}_tool"
            ]
            
            available = any(hasattr(agent_tools, name) for name in possible_names)
            
            # Also check for FastQC as a special case
            if tool_name.lower() == "fastqc":
                available = hasattr(agent_tools, "fastqc_quality_analysis")
            
            self.tool_registry[tool_name] = available
            return available
            
        except Exception as e:
            logger.warning(f"[WorkflowExecutor] Error checking tool '{tool_name}': {e}")
            self.tool_registry[tool_name] = False
            return False
    
    async def _generate_and_execute_tool(
        self,
        operation: OperationSpec,
        parameters: Dict[str, Any],
        session_id: str,
        session_context: Dict
    ) -> Dict[str, Any]:
        """
        Generate and execute a missing tool using the Code Generator Agent.
        
        This method combines generation and execution in one step, which is more
        efficient than generating separately and then executing.
        
        Args:
            operation: Operation specification
            parameters: Execution parameters (including inputs from previous steps)
            session_id: Session ID
            session_context: Session context with file info
        
        Returns:
            Dict with execution result (compatible with _execute_tool format)
        """
        try:
            from backend.tool_generator_agent import generate_and_execute_tool
            
            logger.info(
                f"[WorkflowExecutor] 🔧 Generating and executing tool '{operation.tool_name}' "
                f"for operation '{operation.operation_name}'"
            )
            
            # Build command for tool generator
            tool_command = f"Run {operation.tool_name} for {operation.operation_name}"
            
            # Add input files if available
            input_files = []
            if "input_r1" in parameters:
                input_files.append(parameters["input_r1"])
            if "input_r2" in parameters:
                input_files.append(parameters["input_r2"])
            if "input_files" in parameters:
                if isinstance(parameters["input_files"], list):
                    input_files.extend(parameters["input_files"])
                else:
                    input_files.append(parameters["input_files"])
            if "input_bam" in parameters:
                input_files.append(parameters["input_bam"])
            
            if input_files:
                tool_command += f". Input files: {', '.join(input_files)}"
            
            # Add parameter details
            if operation.parameters:
                param_desc = ", ".join(f"{k}={v}" for k, v in operation.parameters.items() if v is not None)
                if param_desc:
                    tool_command += f". Parameters: {param_desc}"
            
            # Add output location if specified
            if "output" in parameters:
                tool_command += f". Output to: {parameters['output']}"
            
            # Create user request
            user_request = f"Generate and run {operation.tool_name} for {operation.operation_name}"
            
            # Call tool generator with actual inputs/outputs
            logger.info(f"[WorkflowExecutor] Calling tool generator: {tool_command}")
            output_uris = []
            if "output" in parameters and isinstance(parameters["output"], str):
                output_uris.append(parameters["output"])
            if "output_path" in parameters and isinstance(parameters["output_path"], str):
                output_uris.append(parameters["output_path"])

            generation_result = await generate_and_execute_tool(
                command=tool_command,
                user_request=user_request,
                session_id=session_id,
                inputs=[{"uri": f, "size_bytes": None} for f in input_files],
                outputs=[{"uri": o} for o in output_uris],
                session_context=session_context
            )
            
            # Convert generation result to execution result format
            if generation_result.get("status") == "success":
                logger.info(f"[WorkflowExecutor] ✅ Tool '{operation.tool_name}' generated and executed successfully")
                
                return {
                    "status": "completed",
                    "execution_mode": "generated_tool",
                    "message": "Tool generated and executed successfully",
                    "outputs": self._extract_outputs(generation_result.get("result", {})),
                    "raw_result": generation_result
                }
            else:
                error_msg = generation_result.get("error") or generation_result.get("explanation", "Unknown error")
                logger.warning(f"[WorkflowExecutor] ❌ Tool generation/execution failed: {error_msg}")
                
                # Log full execution result for debugging
                exec_result = generation_result.get("execution_result", {})
                if exec_result.get("stderr_full"):
                    logger.error(f"[WorkflowExecutor] Full stderr from generated tool:\n{exec_result['stderr_full']}")
                
                return {
                    "status": "failed",
                    "execution_mode": "generated_tool",
                    "message": f"Tool generation/execution failed: {error_msg}",
                    "outputs": {},
                    "raw_result": generation_result
                }
            
        except Exception as e:
            logger.error(f"[WorkflowExecutor] Tool generation/execution exception: {e}", exc_info=True)
            return {
                "status": "failed",
                "execution_mode": "generated_tool",
                "message": f"Tool generation/execution exception: {str(e)}",
                "outputs": {}
            }
    
    async def _generate_tool(
        self,
        operation: OperationSpec,
        session_context: Dict
    ) -> Dict[str, Any]:
        """
        DEPRECATED: Use _generate_and_execute_tool instead.
        
        This method is kept for backward compatibility but is no longer used.
        """
        return {
            "success": False,
            "error": "This method is deprecated, use _generate_and_execute_tool instead"
        }
    
    def _prepare_execution_params(
        self,
        operation: OperationSpec,
        previous_outputs: Dict[str, Any],
        session_context: Dict
    ) -> Dict[str, Any]:
        """
        Prepare execution parameters for a step.
        
        Combines:
        - Operation parameters
        - Outputs from previous steps
        - Input files from workflow plan (for first step)
        - Session context
        
        Returns:
            Dict of execution parameters
        """
        params = dict(operation.parameters)  # Copy operation parameters
        
        def _filter_fastq_inputs(files: List[str]) -> List[str]:
            fastq_exts = (".fq", ".fastq", ".fq.gz", ".fastq.gz")
            fastqs = [f for f in files if isinstance(f, str) and f.lower().endswith(fastq_exts)]
            candidates = fastqs if fastqs else files
            if len(candidates) >= 2:
                r1 = next((f for f in candidates if "r1" in f.lower()), None)
                r2 = next((f for f in candidates if "r2" in f.lower()), None)
                if r1 and r2:
                    ordered = [r1, r2] + [f for f in candidates if f not in {r1, r2}]
                    return ordered
            return candidates

        # For first step, get inputs from workflow plan
        if not previous_outputs and "workflow_plan" in session_context:
            workflow_plan = session_context["workflow_plan"]
            if workflow_plan.data_inputs:
                input_files = _filter_fastq_inputs([inp.uri for inp in workflow_plan.data_inputs])
                
                # Map to tool-specific parameter names
                if operation.tool_name.lower() in ["fastqc", "trimmomatic", "quality_control"]:
                    # FastQC/QC tools expect input_r1 and input_r2
                    if len(input_files) >= 1:
                        params["input_r1"] = input_files[0]
                    if len(input_files) >= 2:
                        params["input_r2"] = input_files[1]
                else:
                    # Generic mapping
                    params["input_files"] = input_files
        
        # Add previous outputs as inputs
        if previous_outputs:
            # Map common output types to input parameters
            if "bam_file" in previous_outputs:
                params["input_bam"] = previous_outputs["bam_file"]
            
            if "fastq_files" in previous_outputs:
                params["input_files"] = previous_outputs["fastq_files"]
            
            if "alignment_file" in previous_outputs:
                params["alignment"] = previous_outputs["alignment_file"]

        # If we still don't have inputs for tools that need raw reads, fall back to workflow inputs
        if "workflow_plan" in session_context:
            workflow_plan = session_context["workflow_plan"]
            if workflow_plan.data_inputs:
                input_files = _filter_fastq_inputs([inp.uri for inp in workflow_plan.data_inputs])
                if operation.tool_name.lower() in ["trimmomatic", "star", "featurecounts"] and not params.get("input_r1"):
                    if len(input_files) >= 1:
                        params["input_r1"] = input_files[0]
                    if len(input_files) >= 2:
                        params["input_r2"] = input_files[1]

        # Add output targets from command, if available
        if "workflow_command" in session_context and "output" not in params:
            try:
                from backend.tool_generator_agent import _discover_outputs_from_args

                command_text = session_context.get("workflow_command") or ""
                outputs = _discover_outputs_from_args({}, command_text)
                output_uris = [o.get("uri") for o in outputs if isinstance(o, dict) and o.get("uri")]

                if output_uris:
                    tool = (operation.tool_name or "").lower()

                    def pick_output(match_terms: List[str]) -> Optional[str]:
                        for uri in output_uris:
                            if any(term in uri.lower() for term in match_terms):
                                return uri
                        return None

                    def derive_prefix(suffix: str) -> Optional[str]:
                        base_uri = output_uris[0]
                        if base_uri.startswith("s3://"):
                            parts = base_uri.rstrip("/").rsplit("/", 1)
                            if len(parts) == 2:
                                return f"{parts[0]}/{suffix}/"
                        return None

                    selected = None
                    if tool in ["fastqc"]:
                        selected = pick_output(["fastqc", "qc"])
                    elif tool in ["trimmomatic", "trimming"]:
                        selected = pick_output(["trim", "trimmomatic"]) or derive_prefix("trimmed")
                    elif tool in ["star", "alignment"]:
                        selected = pick_output(["aligned", "alignment", "alignments", "bam"])
                    elif tool in ["featurecounts"]:
                        selected = pick_output(["count"])

                    params["output"] = selected or output_uris[0]
            except Exception as e:
                logger.debug(f"[WorkflowExecutor] Output discovery skipped: {e}")
        
        # Add session context info
        if "uploaded_files" in session_context:
            params["session_files"] = session_context["uploaded_files"]
        
        return params
    
    async def _execute_tool(
        self,
        tool_name: str,
        parameters: Dict[str, Any],
        session_id: str
    ) -> Dict[str, Any]:
        """
        Execute a tool using the execution broker.
        
        Args:
            tool_name: Name of tool to execute
            parameters: Tool parameters
            session_id: Session ID
        
        Returns:
            Dict with execution result
        """
        try:
            # Map tool names to actual function names in agent_tools
            tool_mapping = {
                "fastqc": "fastqc_quality_analysis",
                "STAR": "star_alignment",  # Will need generation
                "featureCounts": "feature_counts",  # Will need generation
                "trimmomatic": "trimmomatic_trimming"
            }
            
            actual_tool_name = tool_mapping.get(tool_name, tool_name)
            
            # Use execution broker
            from backend import agent_tools
            
            if hasattr(agent_tools, actual_tool_name):
                tool_func = getattr(agent_tools, actual_tool_name)
                
                # Check if it's a LangChain Tool (has invoke/run methods)
                if hasattr(tool_func, 'invoke'):
                    # LangChain Tool - use invoke() method
                    logger.info(f"[WorkflowExecutor] Executing LangChain tool: {actual_tool_name}")
                    result = tool_func.invoke(parameters)
                elif hasattr(tool_func, 'run'):
                    # LangChain Tool - use run() method (older API)
                    logger.info(f"[WorkflowExecutor] Executing LangChain tool (run): {actual_tool_name}")
                    result = tool_func.run(parameters)
                elif callable(tool_func):
                    # Regular Python function
                    logger.info(f"[WorkflowExecutor] Executing regular function: {actual_tool_name}")
                    result = tool_func(**parameters)
                else:
                    return {
                        "status": "failed",
                        "message": f"Tool '{actual_tool_name}' is not callable",
                        "outputs": {}
                    }
                
                # Handle result
                if isinstance(result, dict):
                    # Respect tool-level errors. Many tools return {"status": "error", ...}
                    tool_status = str(result.get("status", "")).lower()
                    tool_success = result.get("success", True)
                    is_error = tool_status in {"error", "failed", "failure"} or tool_success is False

                    return {
                        "status": "failed" if is_error else "completed",
                        "execution_mode": result.get("execution_mode", "local"),
                        "message": result.get("text") or result.get("message") or ("Tool failed" if is_error else "Tool executed successfully"),
                        "outputs": self._extract_outputs(result),
                        "raw_result": result
                    }
                else:
                    # Result is not a dict (might be string or other type)
                    return {
                        "status": "completed",
                        "execution_mode": "local",
                        "message": str(result),
                        "outputs": {"result": result},
                        "raw_result": {"result": result}
                    }
            else:
                return {
                    "status": "failed",
                    "message": f"Tool '{actual_tool_name}' not found in agent_tools",
                    "outputs": {}
                }
                
        except Exception as e:
            logger.error(f"[WorkflowExecutor] Tool execution error: {e}", exc_info=True)
            return {
                "status": "failed",
                "message": str(e),
                "outputs": {}
            }
    
    def _extract_outputs(self, tool_result: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract outputs from tool result for passing to next step.
        
        Args:
            tool_result: Raw tool execution result
        
        Returns:
            Dict of output files/data
        """
        outputs = {}
        
        # Extract common output types
        if "output" in tool_result:
            outputs["output_path"] = tool_result["output"]
        
        if "output_bam" in tool_result:
            outputs["bam_file"] = tool_result["output_bam"]
        
        if "counts_file" in tool_result:
            outputs["counts_file"] = tool_result["counts_file"]
        
        if "results_path" in tool_result:
            outputs["results_path"] = tool_result["results_path"]
        
        return outputs


# Global workflow executor instance
_workflow_executor = None

def get_workflow_executor() -> WorkflowExecutor:
    """Get global workflow executor instance (singleton)."""
    global _workflow_executor
    if _workflow_executor is None:
        _workflow_executor = WorkflowExecutor()
    return _workflow_executor
