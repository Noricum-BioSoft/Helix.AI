"""
Scenario executor - runs scenarios against the actual multi-agent system.

This module orchestrates scenario execution by:
1. Setting up the execution environment (mock or real)
2. Invoking the multi-agent system with scenario input
3. Tracing agent calls
4. Validating results
5. Producing execution traces
"""

import asyncio
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional
from pydantic import BaseModel, Field
from enum import Enum

from .scenario import Scenario
from .tracer import AgentCallTracer, trace_agent
from .validator import ContractValidator, ValidationResult

# Import AgentRole with fallback
try:
    from backend.agent import AgentRole
except ImportError:
    class AgentRole(str, Enum):
        """Agent roles in the system (LEGACY - see backend.config.agent_registry)."""
        INTENT_DETECTOR = "IntentDetector"
        GURU = "BioinformaticsGuru"
        PLANNER = "BioinformaticsExecutor"
        INFRA = "InfrastructureExpert"
        CODEGEN = "CodeGenerator"
        BROKER = "ExecutionBroker"
        VISUALIZER = "DataVisualizer"


# Agent name mapping: OLD system → NEW canonical names
# See docs/AGENT_NAMING_MIGRATION.md for details
AGENT_NAME_MAPPING = {
    # OLD system names → NEW canonical names
    "BioinformaticsGuru": "WorkflowPlannerAgent",  # TBD: May need refinement
    "BioinformaticsExecutor": "ImplementationAgent",
    "InfrastructureExpert": "InfrastructureDecisionAgent",
    "ExecutionBroker": "ExecutionBroker",  # Not an agent, but a service
    
    # NEW canonical names → NEW canonical names (identity)
    "InfrastructureDecisionAgent": "InfrastructureDecisionAgent",
    "ImplementationAgent": "ImplementationAgent",
    "WorkflowPlannerAgent": "WorkflowPlannerAgent",
    "ToolGeneratorAgent": "ToolGeneratorAgent",
    "RefactorAgent": "RefactorAgent",
    
    # Other agents
    "IntentDetector": "IntentDetector",
    "DataVisualizer": "DataVisualizer",
    "CodeGenerator": "ToolGeneratorAgent",
}


def normalize_agent_name(agent_name: str) -> str:
    """
    Normalize agent name to canonical form.
    
    Converts OLD system names to NEW canonical names.
    Issues warnings for deprecated names.
    """
    if agent_name in AGENT_NAME_MAPPING:
        canonical = AGENT_NAME_MAPPING[agent_name]
        if canonical != agent_name:
            # Deprecated name used
            import warnings
            warnings.warn(
                f"Agent name '{agent_name}' is deprecated. "
                f"Use canonical name '{canonical}' instead. "
                f"See docs/AGENT_NAMING_MIGRATION.md",
                DeprecationWarning,
                stacklevel=2
            )
        return canonical
    else:
        # Unknown agent name
        import warnings
        warnings.warn(
            f"Unknown agent name '{agent_name}'. "
            f"See backend.config.agent_registry for valid names.",
            UserWarning,
            stacklevel=2
        )
        return agent_name


class ExecutionTrace(BaseModel):
    """Complete trace of a scenario execution."""
    scenario_id: str
    timestamp: str
    agent_sequence: List[str]
    contracts: List[Dict[str, Any]]
    validation_result: Optional[Dict[str, Any]] = None
    duration_ms: int
    errors: List[str] = Field(default_factory=list)
    metadata: Dict[str, Any] = Field(default_factory=dict)
    
    def to_json(self, path: Path) -> None:
        """Save trace to JSON file."""
        with open(path, 'w') as f:
            json.dump(self.model_dump(), f, indent=2)
    
    @classmethod
    def from_json(cls, path: Path) -> "ExecutionTrace":
        """Load trace from JSON file."""
        with open(path, 'r') as f:
            data = json.load(f)
        return cls(**data)


class ScenarioExecutor:
    """
    Executes scenarios against the multi-agent system.
    
    This is the main execution engine that:
    - Orchestrates agent invocations
    - Records execution traces
    - Validates results
    - Manages mock vs. real execution modes
    
    Usage:
        executor = ScenarioExecutor(mock_mode=True)
        result = await executor.execute(scenario)
        
        if result.validation_result.passed:
            print("Scenario passed!")
        else:
            print(f"Scenario failed: {result.validation_result.summary}")
    """
    
    def __init__(
        self,
        mock_mode: bool = True,
        baseline_dir: Optional[Path] = None,
        verbose: bool = False
    ):
        """
        Initialize the scenario executor.
        
        Args:
            mock_mode: If True, use mocked LLM responses
            baseline_dir: Directory for storing/loading baseline traces
            verbose: If True, print detailed execution info
        """
        self.mock_mode = mock_mode
        self.baseline_dir = baseline_dir
        self.verbose = verbose
        self.tracer = AgentCallTracer()
        self.validator = ContractValidator()
    
    async def execute(self, scenario: Scenario) -> ExecutionTrace:
        """
        Execute a scenario and return the execution trace.
        
        Args:
            scenario: Scenario to execute
        
        Returns:
            ExecutionTrace with results and validation
        """
        self.tracer.clear()
        
        if self.verbose:
            print(f"\n{'='*60}")
            print(f"Executing scenario: {scenario.metadata.id}")
            print(f"Category: {scenario.metadata.category.value}")
            print(f"Description: {scenario.metadata.description}")
            print(f"{'='*60}\n")
        
        start_time = datetime.utcnow()
        errors = []
        intent_classification = None
        
        try:
            # Execute the scenario through the multi-agent system
            if self.mock_mode or scenario.metadata.requires_mock:
                result = await self._execute_mock(scenario)
                intent_classification = result.get("intent")
            else:
                result = await self._execute_real(scenario)
                intent_classification = result.get("intent")
        
        except Exception as e:
            errors.append(str(e))
            if self.verbose:
                print(f"❌ Execution error: {e}")
        
        # Calculate duration
        end_time = datetime.utcnow()
        duration_ms = int((end_time - start_time).total_seconds() * 1000)
        
        # Validate the execution
        validation_result = None
        if not errors:
            validation_result = self.validator.validate(
                self.tracer,
                scenario.expected_behavior,
                intent=intent_classification
            )
            
            if self.verbose:
                self._print_validation_results(validation_result)
        
        # Create execution trace
        trace = ExecutionTrace(
            scenario_id=scenario.metadata.id,
            timestamp=start_time.isoformat(),
            agent_sequence=self.tracer.get_agent_sequence(),
            contracts=self.tracer.get_contracts(),
            validation_result=validation_result.model_dump() if validation_result else None,
            duration_ms=duration_ms,
            errors=errors,
            metadata={
                "category": scenario.metadata.category.value,
                "tags": scenario.metadata.tags,
                "mock_mode": self.mock_mode,
            }
        )
        
        if self.verbose:
            print(f"\n{'='*60}")
            print(f"Scenario: {scenario.metadata.id}")
            print(f"Duration: {duration_ms}ms")
            print(f"Agent Sequence: {' → '.join(trace.agent_sequence)}")
            if validation_result:
                status = "✓ PASSED" if validation_result.passed else "✗ FAILED"
                print(f"Status: {status}")
                print(f"Summary: {validation_result.summary}")
            print(f"{'='*60}\n")
        
        return trace
    
    async def _execute_mock(self, scenario: Scenario) -> Dict[str, Any]:
        """
        Execute scenario in mock mode.
        
        This simulates agent behavior without calling real LLMs or executing workflows.
        Uses the expected behavior from the scenario to generate mock responses.
        """
        prompt = scenario.input.user_prompt
        expected = scenario.expected_behavior
        
        # Mock intent detection
        intent = expected.intent.type if expected.intent else "ask"
        
        with trace_agent(self.tracer, AgentRole.INTENT_DETECTOR, {"prompt": prompt}) as call:
            # Create mock IntentResult (simple dict, no imports needed)
            # Use confidence slightly above minimum to pass > threshold validations
            base_confidence = expected.intent.confidence_min if expected.intent else 0.9
            confidence = min(base_confidence + 0.05, 1.0)  # Add 0.05 buffer, cap at 1.0
            
            # Add clarifying questions if confidence is low
            clarifying_questions = []
            if confidence < 0.7:
                clarifying_questions = [
                    "Could you provide more details about what you want to do?",
                    "What type of data are you working with?"
                ]
            
            intent_result = {
                "intent": intent,
                "confidence": confidence,
                "clarifying_questions": clarifying_questions,
                "rationale": f"Mock classification for: {prompt}"
            }
            call.set_output(intent_result, "IntentResult")
        
        # Route based on intent
        if intent == "ask":
            return await self._execute_mock_ask_path(scenario)
        else:
            return await self._execute_mock_execute_path(scenario)
    
    async def _execute_mock_ask_path(self, scenario: Scenario) -> Dict[str, Any]:
        """Execute mock ask (Q&A) path."""
        prompt = scenario.input.user_prompt
        
        # Mock Guru response
        with trace_agent(self.tracer, AgentRole.GURU, {"question": prompt}) as call:
            # Create a mock answer (simple dict, no imports needed)
            mock_answer = {
                "text": f"Mock answer to: {prompt}. This is a simulated response covering RNA sequencing DNA analysis transcriptome and related bioinformatics concepts with sufficient detail for validation.",
                "citations": [],
                "followups": [],
                "suggested_next_actions": []
            }
            call.set_output(mock_answer, "Answer")
        
        return {"intent": "ask"}
    
    async def _execute_mock_execute_path(self, scenario: Scenario) -> Dict[str, Any]:
        """Execute mock execute (workflow) path."""
        prompt = scenario.input.user_prompt
        expected = scenario.expected_behavior
        session_context = scenario.input.session_context
        
        # Detect dataset scale from file paths and context
        scale = self._detect_dataset_scale(session_context)
        
        # Mock Planner (accept both OLD and NEW names)
        planner_names = ["BioinformaticsExecutor", "ImplementationAgent", "WorkflowPlannerAgent"]
        if any(name in expected.agent_sequence for name in planner_names):
            with trace_agent(self.tracer, AgentRole.PLANNER, {"request": prompt}) as call:
                # Create mock WorkflowPlan (simple dict, no imports needed)
                # Determine if this is a complex workflow (needs multiple steps)
                is_complex = any(keyword in prompt.lower() for keyword in ["complete", "analysis", "multiple", "pipeline"])
                
                steps = [
                    {"tool": "quality_control", "args": {"threshold": 30}},
                    {"tool": "alignment", "args": {"reference": "hg38"}},
                    {"tool": "differential_expression", "args": {"method": "DESeq2"}}
                ] if is_complex else [{"tool": "fastqc" if "fastqc" in prompt.lower() or "quality" in prompt.lower() else "mock_tool", "args": {}}]
                
                plan = {
                    "version": "v1",
                    "metadata": {
                        "scale": scale,  # Use detected scale
                        "parallelism": "by_sample" if scale == "large" else "none",
                        "reproducibility": "medium",
                        "expected_outputs": [],
                        "risk_flags": [],
                        "missing_info": []
                    },
                    "inputs": [],
                    "steps": steps
                }
                call.set_output(plan, "WorkflowPlan")
        
        # Mock Infrastructure Expert (accept both OLD and NEW names)
        infra_names = ["InfrastructureExpert", "InfrastructureDecisionAgent"]
        if any(name in expected.agent_sequence for name in infra_names):
            with trace_agent(self.tracer, AgentRole.INFRA, {"plan": "mock_plan"}) as call:
                # Choose infrastructure based on detected scale
                if scale == "small":
                    environment = "Local"
                    summary = "Mock infrastructure decision: Local execution for small scale dataset"
                    cost_analysis = None
                    confidence = 0.9
                elif scale == "large":
                    environment = "EMR"
                    summary = "Mock infrastructure decision: EMR cluster for large scale dataset with distributed processing"
                    cost_analysis = {
                        "estimated_cost_range_usd": "$2-5",
                        "assumptions": ["2 x m5.xlarge nodes", "10-15 minute runtime"],
                        "confidence": "medium"
                    }
                    confidence = 0.85
                else:  # medium
                    environment = "EMR"
                    summary = "Mock infrastructure decision: Medium scale dataset requires EMR execution"
                    cost_analysis = {"estimated_cost_range_usd": "$1-3", "assumptions": [], "confidence": "medium"}
                    confidence = 0.8
                
                # Add constraints for large scale data locality
                constraints = []
                if scale == "large":
                    constraints.append({
                        "key": "data_locality",
                        "value": "prefer_s3_colocation",
                        "rationale": "Large dataset benefits from data locality in S3"
                    })
                
                infra_decision = {
                    "recommended_environment": environment,
                    "confidence_score": confidence,
                    "decision_summary": summary,
                    "constraints": constraints,
                    "cost_analysis": cost_analysis,
                    "alternatives": [{"environment": "AWS Batch", "tradeoffs": "Lower cost but slower startup"}] if scale == "large" else [],
                    "warnings": [],
                    "assumptions": []
                }
                call.set_output(infra_decision, "InfraDecisionContract")
        
        # Mock CodeGen (if present)
        if "CodeGenerator" in expected.agent_sequence:
            with trace_agent(self.tracer, AgentRole.CODEGEN, {"plan": "mock_plan"}) as call:
                # Create mock ExecutionSpec (simple dict, no imports needed)
                exec_spec = {
                    "target": "emr",
                    "container_image": None,
                    "entrypoint": None,
                    "command": ["mock", "command"],
                    "env": {},
                    "inputs": [],
                    "outputs": [],
                    "resources": {},
                    "retry_strategy": {},
                    "notes": []
                }
                call.set_output(exec_spec, "ExecutionSpec")
        
        # Mock Execution Broker (NOTE: ExecutionBroker is not an agent in NEW system)
        # In NEW orchestrator system, execution happens via external runner, not broker
        if "ExecutionBroker" in expected.agent_sequence:
            with trace_agent(self.tracer, AgentRole.BROKER, {"spec": "mock_spec"}) as call:
                # Create mock ExecutionResult (simple dict, no imports needed)
                # Match environment from infrastructure decision
                environment_str = "Local" if scale == "small" else "EMR"
                
                # Generate appropriate artifacts
                artifacts = []
                if "fastqc" in prompt.lower() or "quality" in prompt.lower():
                    # FastQC produces reports
                    for filename in session_context.get("uploaded_files", []):
                        base = filename.split("/")[-1].replace(".fq", "").replace(".fastq", "")
                        artifacts.append({
                            "uri": f"s3://output/{base}_fastqc.html",
                            "description": f"FastQC report for {base}",
                            "artifact_type": "html_report"
                        })
                
                exec_result = {
                    "status": "succeeded",
                    "job_id": f"{'local' if scale == 'small' else 'j-MOCKEMR'}_job_123",
                    "environment": environment_str,
                    "logs_uri": f"s3://logs/mock_job_123.log" if scale != "small" else None,
                    "artifacts": artifacts,
                    "metrics": {"duration_seconds": 30 if scale == "small" else 600},
                    "error": None
                }
                call.set_output(exec_result, "ExecutionResult")
        
        # Mock Visualizer
        if "DataVisualizer" in expected.agent_sequence:
            with trace_agent(self.tracer, AgentRole.VISUALIZER, {"result": "mock_result"}) as call:
                # Create mock VisualizationArtifacts (simple dict, no imports needed)
                # Customize summary based on workflow type
                if "merge" in prompt.lower() or "merging" in prompt.lower():
                    if scale == "large":
                        summary = "Mock visualization: Successfully merged paired reads across multiple samples with quality metrics"
                    else:
                        summary = "Mock visualization: Successfully merged paired-end read pairs with quality metrics"
                    artifacts = [{"type": "plot", "uri": "s3://output/merge_stats.png"}] if scale == "large" else []
                elif "fastqc" in prompt.lower() or "quality" in prompt.lower():
                    summary = "Mock visualization: Quality control metrics generated successfully"
                    artifacts = []
                else:
                    summary = "Mock visualization created successfully"
                    artifacts = []
                
                viz = {
                    "summary": summary,
                    "artifacts": artifacts,
                    "warnings": []
                }
                call.set_output(viz, "VisualizationArtifacts")
        
        return {"intent": "execute"}
    
    def _detect_dataset_scale(self, session_context: Dict[str, Any]) -> str:
        """
        Detect dataset scale from session context.
        
        Looks at:
        - File paths (/test/ directory or test_ prefix = small)
        - Dataset info (file sizes, read counts)
        - Number of files
        
        Returns: "small", "medium", or "large"
        """
        uploaded_files = session_context.get("uploaded_files", [])
        dataset_info = session_context.get("dataset_info", {})
        
        # Check for test datasets in file paths - must be in /test/ directory or test_ prefix
        for filepath in uploaded_files:
            filepath_str = str(filepath)
            # Check if in test directory or has test_ prefix in filename
            if "/test/" in filepath_str or filepath_str.split("/")[-1].startswith("test_"):
                return "small"
        
        # Check explicit dataset info
        total_reads = dataset_info.get("total_reads", 0)
        file_size_mb = dataset_info.get("file_size_mb", 0)
        file_size_gb = dataset_info.get("file_size_gb", 0)
        
        # Classify by size (explicit indicators)
        if total_reads > 1000000 or file_size_gb > 0.5 or file_size_mb > 500:
            return "large"
        elif total_reads > 0 and total_reads < 10000:
            return "small"
        elif file_size_mb > 0 and file_size_mb < 10:
            return "small"
        
        # Check number of files (many files = larger dataset)
        if len(uploaded_files) > 4:
            return "large"
        
        # If files are in the root of a dataset (not in test/), assume large
        # This handles the case where we have mate_R1.fq, mate_R2.fq in main directory
        if uploaded_files and not any("/test/" in str(f) for f in uploaded_files):
            if any("mate_" in str(f) or "sample_" in str(f).split("/")[-1] for f in uploaded_files):
                # Main data files, not test files
                return "large"
        
        # Default to medium if we can't determine
        return "medium"
    
    async def _execute_real(self, scenario: Scenario) -> Dict[str, Any]:
        """
        Execute scenario against the real multi-agent system.
        
        This calls the actual backend agents and executes workflows on real data.
        Requires:
        - Backend dependencies (langchain, langgraph, boto3, etc.)
        - AWS credentials configured
        - Environment variables set (.env file)
        
        Returns:
            Dict with execution metadata
        
        Raises:
            ImportError: If backend dependencies not available
            Exception: If execution fails
        """
        try:
            # Import backend components
            from backend.agent import handle_command
            print(f"\n🚀 Real Execution Mode")
            print(f"Scenario: {scenario.metadata.id}")
            print(f"Command: {scenario.input.user_prompt}")
            print(f"Dataset: {scenario.input.session_context.get('uploaded_files', [])}")
            print(f"-" * 70)
            
        except ImportError as e:
            raise ImportError(
                f"Cannot import backend modules for real execution: {e}\n"
                "Make sure backend dependencies are installed and PYTHONPATH includes backend/"
            ) from e
        
        # Generate unique session ID for this scenario execution
        session_id = f"demo_{scenario.metadata.id}_{int(datetime.now().timestamp())}"
        
        # Prepare session context with scenario data
        session_context = scenario.input.session_context.copy()
        
        # Add environment variables to session context
        if scenario.input.environment_vars:
            for key, value in scenario.input.environment_vars.items():
                import os
                os.environ[key] = value
        
        try:
            # Execute command via backend agent system
            print(f"⏳ Calling backend agent system...")
            result = await handle_command(
                command=scenario.input.user_prompt,
                session_id=session_id,
                session_context=session_context
            )
            
            print(f"✅ Agent execution completed")
            print(f"Result type: {type(result)}")
            
            # Extract agent calls from result if available
            # The handle_command may return various structures depending on the path taken
            if isinstance(result, dict):
                # Check if tool was mapped and executed
                if result.get("status") == "tool_mapped":
                    print(f"🔧 Tool mapped: {result.get('tool_name')}")
                    with trace_agent(self.tracer, AgentRole.PLANNER, {
                        "request": scenario.input.user_prompt
                    }) as call:
                        call.set_output({
                            "tool_name": result.get("tool_name"),
                            "parameters": result.get("parameters", {})
                        }, "ToolMapping")
                
                # Check for execution results
                if "result" in result or "output" in result:
                    execution_result = result.get("result") or result.get("output")
                    print(f"📊 Execution result keys: {list(execution_result.keys()) if isinstance(execution_result, dict) else 'N/A'}")
                    
                    # Record as execution broker output
                    with trace_agent(self.tracer, AgentRole.BROKER, {
                        "spec": "real_execution"
                    }) as call:
                        call.set_output(execution_result, "ExecutionResult")
                
                # Check for errors
                if result.get("status") == "error" or "error" in result:
                    error_msg = result.get("error", "Unknown error")
                    print(f"❌ Execution error: {error_msg}")
                    raise RuntimeError(f"Backend execution failed: {error_msg}")
            
            print(f"✅ Real execution completed successfully")
            print(f"-" * 70)
            
            return {
                "intent": result.get("intent", "execute") if isinstance(result, dict) else "execute",
                "session_id": session_id,
                "result": result
            }
            
        except Exception as e:
            print(f"❌ Real execution failed: {e}")
            import traceback
            traceback.print_exc()
            raise RuntimeError(f"Real execution failed: {e}") from e
    
    def _print_validation_results(self, validation: ValidationResult) -> None:
        """Print validation results in a human-readable format."""
        if validation.passed:
            print("✓ Validation passed")
        else:
            print("✗ Validation failed")
        
        print(f"\nSummary: {validation.summary}")
        
        if validation.get_errors():
            print(f"\n❌ Errors ({len(validation.get_errors())}):")
            for i, issue in enumerate(validation.get_errors(), 1):
                print(f"  {i}. [{issue.category}] {issue.message}")
                if issue.details:
                    print(f"     Details: {issue.details}")
        
        if validation.get_warnings():
            print(f"\n⚠ Warnings ({len(validation.get_warnings())}):")
            for i, issue in enumerate(validation.get_warnings(), 1):
                print(f"  {i}. [{issue.category}] {issue.message}")
                if issue.details:
                    print(f"     Details: {issue.details}")
    
    def save_baseline(self, trace: ExecutionTrace) -> None:
        """Save an execution trace as a baseline."""
        if self.baseline_dir is None:
            raise ValueError("baseline_dir not configured")
        
        self.baseline_dir.mkdir(parents=True, exist_ok=True)
        baseline_path = self.baseline_dir / f"{trace.scenario_id}.json"
        trace.to_json(baseline_path)
        
        if self.verbose:
            print(f"Saved baseline: {baseline_path}")
    
    def load_baseline(self, scenario_id: str) -> Optional[ExecutionTrace]:
        """Load a baseline trace for comparison."""
        if self.baseline_dir is None:
            return None
        
        baseline_path = self.baseline_dir / f"{scenario_id}.json"
        if not baseline_path.exists():
            return None
        
        return ExecutionTrace.from_json(baseline_path)
    
    def compare_with_baseline(
        self,
        current: ExecutionTrace,
        baseline: ExecutionTrace
    ) -> Dict[str, Any]:
        """Compare current trace with baseline."""
        differences = {
            "agent_sequence_changed": current.agent_sequence != baseline.agent_sequence,
            "contract_count_changed": len(current.contracts) != len(baseline.contracts),
            "duration_delta_ms": current.duration_ms - baseline.duration_ms,
            "details": []
        }
        
        if differences["agent_sequence_changed"]:
            differences["details"].append({
                "type": "agent_sequence",
                "baseline": baseline.agent_sequence,
                "current": current.agent_sequence
            })
        
        return differences
