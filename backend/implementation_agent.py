"""
Implementation Agent - Plans HOW to execute a workflow (but does not execute).

This agent is the second agent in the multi-agent pipeline:
1. InfrastructureDecisionAgent decides WHERE to execute (Local, EC2, EMR, etc.)
2. ImplementationAgent decides HOW to execute (container, commands, retries)
3. External runner actually executes

Input:
- WorkflowPlan (what operation to perform)
- InfraDecision (where to execute)

Output:
- ExecutionToolSpec (how to execute)

Key principle: Separation of planning vs execution.
This agent only plans, it does not execute anything.
"""

import logging
import hashlib
import json
import os
from typing import Optional

from backend.contracts.workflow_plan import WorkflowPlan
from backend.contracts.infra_decision import InfraDecision
from backend.contracts.execution_spec import (
    ExecutionToolSpec,
    ContainerSpec,
    CommandSpec,
    RetryPolicy,
    ResourceRequirements,
    OutputSpec
)

logger = logging.getLogger(__name__)


# Load system prompt from file
from pathlib import Path

IMPLEMENTATION_AGENT_SYSTEM_PROMPT_FILE = Path(__file__).parent.parent / "agents" / "implementation-agent.md"
if IMPLEMENTATION_AGENT_SYSTEM_PROMPT_FILE.exists():
    with open(IMPLEMENTATION_AGENT_SYSTEM_PROMPT_FILE, "r") as f:
        IMPLEMENTATION_AGENT_SYSTEM_PROMPT = f.read()
else:
    # Fallback to inline prompt if file not found
    IMPLEMENTATION_AGENT_SYSTEM_PROMPT = """You are an Implementation Agent that plans HOW to execute bioinformatics workflows.

Your role:
- Consume WorkflowPlan (WHAT to do) and InfraDecision (WHERE to execute)
- Produce ExecutionToolSpec (HOW to execute: container, commands, retries)
- You do NOT execute - you only plan

Key responsibilities:
1. **Container Selection**: Choose appropriate container image (biocontainers, custom, or native)
2. **Command Construction**: Build correct shell commands with proper flags and file paths
3. **Resource Estimation**: Estimate CPU, memory, disk requirements
4. **Retry Strategy**: Define retry policy based on failure modes
5. **Output Specification**: Specify expected outputs for validation

Output schema (JSON):
{json_schema}

Best practices:
- Use biocontainers when available (e.g., biocontainers/fastqc:0.11.9)
- For EMR: Commands should work with S3 URIs directly
- For Local/EC2: Include data download/upload steps if needed
- Set realistic timeouts (most bioinformatics tools: 5-60min)
- Define clear success criteria (exit_code==0, output_file_exists)
- Estimate resources conservatively (better to over-provision than fail)
- Include warnings for uncertainty or assumptions

Confidence scoring:
- 0.9-1.0: Well-known tool with documented container
- 0.7-0.9: Standard tool but custom container or uncommon flags
- 0.5-0.7: Tool requires significant customization
- <0.5: High uncertainty, suggest human review

Think step-by-step:
1. Identify the tool/operation from WorkflowPlan
2. Check if container is available (biocontainers, quay.io, docker.io)
3. Build commands with correct input/output paths
4. Estimate resources based on file sizes and tool requirements
5. Define retry policy (most tools: max_retries=2, retry_on=['exit_code', 'timeout'])
6. List expected outputs with formats
7. Calculate confidence based on certainty of each decision

Return ONLY valid JSON matching the schema. No markdown, no explanations outside JSON.
"""

# Store prompt file path for reference
IMPLEMENTATION_AGENT_PROMPT_FILE = IMPLEMENTATION_AGENT_SYSTEM_PROMPT_FILE


def _get_llm():
    """Get LLM instance for Implementation Agent."""
    model_name = os.getenv("LLM_MODEL", "claude-sonnet-4-20250514")
    
    try:
        from langchain.chat_models import init_chat_model
        llm = init_chat_model(model_name, temperature=0.0)
        return llm
    except Exception:
        # Fallback to DeepSeek if init_chat_model fails
        try:
            from langchain_deepseek import ChatDeepSeek
            return ChatDeepSeek(
                model="deepseek-chat",
                temperature=0.0,
                timeout=60.0
            )
        except Exception:
            raise RuntimeError(
                f"Failed to initialize LLM (model={model_name}). "
                "Check your LangChain installation and API keys."
            )


def _hash_contract(obj) -> str:
    """Compute hash of Pydantic model for traceability."""
    json_str = obj.model_dump_json()
    return hashlib.sha256(json_str.encode()).hexdigest()[:16]


async def plan_implementation(
    workflow_plan: WorkflowPlan,
    infra_decision: InfraDecision,
    request_id: Optional[str] = None
) -> ExecutionToolSpec:
    """
    Plan HOW to execute a workflow step.
    
    This is the main entry point for the Implementation Agent.
    
    Args:
        workflow_plan: What operation to perform
        infra_decision: Where to execute (infrastructure choice)
        request_id: Optional request ID for tracing
    
    Returns:
        ExecutionToolSpec with container, commands, retries, resources
    
    Raises:
        Exception: If planning fails after retries
    """
    logger.info(f"ImplementationAgent: Planning execution for {workflow_plan.description}")
    
    # Get JSON schema for LLM
    json_schema = json.dumps(ExecutionToolSpec.model_json_schema(), indent=2)
    
    # Build user prompt
    user_prompt = f"""Plan the execution for this workflow:

**Workflow Plan:**
```json
{workflow_plan.model_dump_json(indent=2)}
```

**Infrastructure Decision:**
```json
{infra_decision.model_dump_json(indent=2)}
```

Target infrastructure: {infra_decision.infrastructure}

Generate an ExecutionToolSpec that describes HOW to execute this workflow on {infra_decision.infrastructure}.

Return ONLY valid JSON matching the schema (no markdown formatting).
"""
    
    # Try LLM first, fall back to heuristic if fails
    max_retries = 2
    validation_error = None
    
    for attempt in range(max_retries + 1):
        try:
            llm = _get_llm()
            
            # Include validation error in prompt if retrying
            if validation_error:
                user_prompt += f"\n\n**Previous attempt failed validation:**\n{validation_error}\n\nPlease fix the issues and return valid JSON."
            
            messages = [
                {"role": "system", "content": IMPLEMENTATION_AGENT_SYSTEM_PROMPT.format(json_schema=json_schema)},
                {"role": "user", "content": user_prompt}
            ]
            
            # Run LLM call in thread pool to avoid blocking async event loop
            import asyncio
            if hasattr(asyncio, "to_thread"):
                response = await asyncio.to_thread(llm.invoke, messages)
            else:
                loop = asyncio.get_event_loop()
                response = await loop.run_in_executor(None, llm.invoke, messages)
            
            response_content = response.content.strip()
            
            # Parse JSON (handle markdown code blocks)
            json_str = response_content
            if "```json" in json_str:
                json_str = json_str.split("```json")[1].split("```")[0].strip()
            elif "```" in json_str:
                json_str = json_str.split("```")[1].split("```")[0].strip()
            
            result_dict = json.loads(json_str)
            
            # Validate with Pydantic
            execution_spec = ExecutionToolSpec(**result_dict)
            
            # Add traceability hashes
            execution_spec.request_id = request_id
            execution_spec.workflow_plan_hash = _hash_contract(workflow_plan)
            execution_spec.infra_decision_hash = _hash_contract(infra_decision)
            
            logger.info(f"ImplementationAgent: Generated execution plan with confidence {execution_spec.confidence_score}")
            return execution_spec
        
        except Exception as e:
            error_type = type(e).__name__
            validation_error = f"{error_type}: {str(e)[:200]}"
            
            if attempt < max_retries:
                logger.warning(f"ImplementationAgent: Attempt {attempt+1} failed: {validation_error}, retrying...")
                continue
            else:
                logger.error(f"ImplementationAgent: All attempts failed, using heuristic fallback")
                return _heuristic_implementation_plan(workflow_plan, infra_decision, request_id)
    
    # Should not reach here
    return _heuristic_implementation_plan(workflow_plan, infra_decision, request_id)


def _heuristic_implementation_plan(
    workflow_plan: WorkflowPlan,
    infra_decision: InfraDecision,
    request_id: Optional[str] = None
) -> ExecutionToolSpec:
    """
    Heuristic fallback for implementation planning.
    
    Used when LLM fails to produce valid ExecutionToolSpec.
    Returns a conservative, generic execution plan.
    """
    logger.warning("ImplementationAgent: Using heuristic fallback (LLM failed)")
    
    # Extract operation name from description or operations list
    if workflow_plan.operations:
        operation_name = workflow_plan.operations[0].operation_name.lower()
    else:
        operation_name = workflow_plan.description.split()[0].lower() if workflow_plan.description else "unknown_operation"
    
    # Generic container for bioinformatics (if containerized environment)
    container_spec = None
    if infra_decision.infrastructure in ["Batch", "Lambda", "EMR"]:
        container_spec = ContainerSpec(
            image="biocontainers/biocontainers:latest",
            image_type="docker",
            pull_policy="IfNotPresent",
            env_vars={},
            mount_paths=["/data"],
            working_dir="/data"
        )
    
    # Build generic command
    input_uris = [inp.uri for inp in workflow_plan.data_inputs]
    command_str = f"# Placeholder command for {operation_name}\necho 'Execution planned but requires manual implementation'"
    
    commands = [
        CommandSpec(
            name=f"run_{operation_name}",
            command=command_str,
            inputs=input_uris,
            outputs=[],
            success_criteria="exit_code==0",
            timeout_minutes=30.0
        )
    ]
    
    # Conservative resource requirements
    resources = ResourceRequirements(
        min_cpu_cores=2.0,
        min_memory_gb=4.0,
        min_disk_gb=20.0,
        gpu_required=False
    )
    
    # Standard retry policy
    retry_policy = RetryPolicy(
        max_retries=2,
        retry_on=["all"],
        backoff_multiplier=2.0,
        initial_delay_seconds=10.0
    )
    
    return ExecutionToolSpec(
        tool_name=operation_name,
        infrastructure=infra_decision.infrastructure,
        container_spec=container_spec,
        commands=commands,
        retry_policy=retry_policy,
        resource_requirements=resources,
        expected_outputs=[],
        confidence_score=0.3,
        reasoning="Heuristic fallback: LLM failed to generate execution plan. This is a placeholder requiring manual review.",
        warnings=[
            "⚠️ This execution plan was generated by heuristic fallback",
            "⚠️ Manual review required before execution",
            "⚠️ Commands are placeholders and will not execute successfully"
        ],
        estimated_runtime_minutes=None,
        request_id=request_id,
        workflow_plan_hash=_hash_contract(workflow_plan),
        infra_decision_hash=_hash_contract(infra_decision)
    )


# For backward compatibility / testing
async def decide_implementation(*args, **kwargs):
    """Alias for plan_implementation."""
    return await plan_implementation(*args, **kwargs)
