"""
Infrastructure Decision Agent v2 - With Pydantic Validation and Repair Mechanism

Key improvements over v1:
- Uses Pydantic contracts for strict validation
- Repair mechanism: re-prompts LLM on validation errors (max 2 retries)
- Confidence scoring (0-1) on all decisions
- Cost ranges instead of exact values
- Uncertainty tracking (warnings for unknown inputs)

This module maintains backward compatibility with v1 by keeping InputAsset/OutputAsset
dataclasses, but internally uses Pydantic contracts for validation.
"""

import os
import json
import logging
from pathlib import Path
from typing import Dict, Any, Optional, List
from dataclasses import dataclass
from pydantic import ValidationError

from backend.contracts.infra_decision import (
    InfraDecision,
    FileAnalysis,
    ComputationalRequirements,
    CostAnalysis,
    InfraAlternative,
)

logger = logging.getLogger(__name__)

# Load the infrastructure-decision-agent prompt
PROJECT_ROOT = Path(__file__).resolve().parent.parent
INFRASTRUCTURE_DECISION_PROMPT_PATH = PROJECT_ROOT / "agents" / "infrastructure-decision-agent.md"

try:
    # Load infrastructure decision agent prompt
    if INFRASTRUCTURE_DECISION_PROMPT_PATH.exists():
        INFRASTRUCTURE_DECISION_SYSTEM_PROMPT = INFRASTRUCTURE_DECISION_PROMPT_PATH.read_text(encoding='utf-8')
        if not INFRASTRUCTURE_DECISION_SYSTEM_PROMPT or not INFRASTRUCTURE_DECISION_SYSTEM_PROMPT.strip():
            raise ValueError(f"Infrastructure decision prompt file is empty: {INFRASTRUCTURE_DECISION_PROMPT_PATH}")
        logger.info(f"Loaded infrastructure decision prompt from {INFRASTRUCTURE_DECISION_PROMPT_PATH} ({len(INFRASTRUCTURE_DECISION_SYSTEM_PROMPT)} chars)")
    else:
        raise FileNotFoundError(f"Infrastructure decision prompt file not found at: {INFRASTRUCTURE_DECISION_PROMPT_PATH}")
except Exception as e:
    logger.warning(f"Could not load infrastructure-decision-agent prompt: {e}")
    logger.warning("Using fallback prompt")
    INFRASTRUCTURE_DECISION_SYSTEM_PROMPT = """You are an Infrastructure Decision Agent. Analyze bioinformatics operations and recommend execution infrastructure (Local, EC2, EMR, Batch, Lambda) based on file sizes, locations, and computational requirements. Return JSON with infrastructure choice and reasoning."""


# Keep backward-compatible dataclasses for InputAsset/OutputAsset
# (Used by other modules like execution_broker.py and tool_generator_agent.py)
@dataclass
class InputAsset:
    """Represents an input file/asset for infrastructure decision."""
    uri: str
    size_bytes: Optional[int] = None
    source: str = "unknown"


@dataclass
class OutputAsset:
    """Represents an output file/asset for infrastructure decision."""
    uri: str
    estimated_size_bytes: Optional[int] = None


def _get_llm():
    """
    Lazily initialize and return the LLM for infrastructure decision.
    
    Uses the same pattern as backend.agent._get_llm() for consistency.
    Handles HELIX_MOCK_MODE and supports both OpenAI and DeepSeek.
    """
    if os.getenv("HELIX_MOCK_MODE") == "1":
        raise RuntimeError("LLM is disabled in HELIX_MOCK_MODE")

    openai_key = os.getenv("OPENAI_API_KEY", "").strip()
    deepseek_key = os.getenv("DEEPSEEK_API_KEY", "").strip()

    openai_enabled = openai_key and openai_key not in ["", "disabled", "your_openai_api_key_here", "none"]
    deepseek_enabled = deepseek_key and deepseek_key not in ["", "disabled", "your_deepseek_api_key_here", "none"]

    if openai_enabled:
        from langchain.chat_models import init_chat_model
        openai_model = os.getenv("HELIX_INFRASTRUCTURE_OPENAI_MODEL", "openai:gpt-5.5").strip()
        if ":" not in openai_model:
            openai_model = f"openai:{openai_model}"
        return init_chat_model(openai_model, temperature=0)

    if deepseek_enabled:
        from langchain_deepseek import ChatDeepSeek
        return ChatDeepSeek(
            model="deepseek-chat",
            temperature=0,
            max_tokens=None,
            timeout=30.0,
            max_retries=1,
        )

    raise ValueError(
        "No API keys found for infrastructure decision. "
        "Please set either OPENAI_API_KEY or DEEPSEEK_API_KEY "
        "in your environment variables."
    )


def _analyze_files(inputs: List[InputAsset], outputs: List[OutputAsset]) -> FileAnalysis:
    """
    Analyze file characteristics for infrastructure decision.
    
    Returns a FileAnalysis Pydantic model (v2) instead of plain dict (v1).
    
    Uses size_bytes and source fields from InputAsset objects, which are now
    populated from WorkflowPlan.data_inputs (size_bytes, location_type).
    """
    total_size_bytes = 0
    file_count = 0
    unknown_sizes = 0
    locations = set()
    largest_file_bytes = 0
    
    for inp in inputs:
        file_count += 1
        if inp.size_bytes is not None and inp.size_bytes >= 0:
            total_size_bytes += inp.size_bytes
            largest_file_bytes = max(largest_file_bytes, inp.size_bytes)
        else:
            unknown_sizes += 1
        
        # Use location from InputAsset.source (populated from workflow_plan.data_inputs[].location_type)
        # If source is "unknown", infer from URI as fallback
        if inp.source and inp.source != "unknown":
            locations.add(inp.source)
        elif inp.uri.startswith("s3://"):
            locations.add("S3")
        elif inp.uri.startswith("/") or inp.uri.startswith("./") or inp.uri.startswith("../"):
            locations.add("Local")
        else:
            locations.add("Unknown")
    
    total_size_mb = total_size_bytes / (1024 * 1024) if total_size_bytes > 0 else 0.0
    largest_file_mb = largest_file_bytes / (1024 * 1024) if largest_file_bytes > 0 else 0.0
    
    # Return Pydantic model
    return FileAnalysis(
        total_size_bytes=total_size_bytes,
        total_size_mb=round(total_size_mb, 2),
        file_count=file_count,
        unknown_sizes=unknown_sizes,
        locations=list(locations),
        largest_file_mb=round(largest_file_mb, 2)
    )


def _build_decision_prompt(
    command: str,
    file_analysis: FileAnalysis,
    inputs: List[InputAsset],
    outputs: List[OutputAsset],
    session_context: Optional[Dict[str, Any]] = None,
    validation_error: Optional[str] = None
) -> str:
    """
    Build the prompt for infrastructure decision.
    
    Args:
        command: User command/request
        file_analysis: FileAnalysis Pydantic model
        inputs: List of input assets
        outputs: List of output assets
        session_context: Optional session context
        validation_error: If retry, include validation error from previous attempt
    
    Returns:
        Formatted prompt string
    """
    
    # Build input file details
    input_details = []
    for inp in inputs:
        size_str = f"{inp.size_bytes / (1024*1024):.2f} MB" if inp.size_bytes else "size unknown"
        input_details.append(f"  - {inp.uri} ({size_str}, source: {inp.source})")
    
    output_details = []
    for out in outputs:
        size_str = f"{out.estimated_size_bytes / (1024*1024):.2f} MB" if out.estimated_size_bytes else "size unknown"
        output_details.append(f"  - {out.uri} ({size_str})")
    
    # Get JSON schema for InfraDecision (Pydantic V2)
    schema_dict = InfraDecision.model_json_schema()
    schema_json = json.dumps(schema_dict, indent=2)
    
    prompt = f"""Analyze this bioinformatics operation and determine the optimal execution infrastructure.

**Operation**: {command}

**Input Files**:
{chr(10).join(input_details) if input_details else "  - No input files specified"}

**Output Files**:
{chr(10).join(output_details) if output_details else "  - No output files specified"}

**File Analysis Summary**:
- Total size: {file_analysis.total_size_mb} MB ({file_analysis.total_size_bytes:,} bytes)
- File count: {file_analysis.file_count}
- Unknown sizes: {file_analysis.unknown_sizes}
- Locations: {', '.join(file_analysis.locations)}
- Largest file: {file_analysis.largest_file_mb} MB

**Environment Configuration**:
- HELIX_USE_EC2: {os.getenv('HELIX_USE_EC2', 'false')}
- HELIX_ASYNC_BYTES_THRESHOLD: {os.getenv('HELIX_ASYNC_BYTES_THRESHOLD', '104857600')} bytes ({int(os.getenv('HELIX_ASYNC_BYTES_THRESHOLD', '104857600')) / (1024*1024):.0f} MB)

Based on the infrastructure decision matrix and decision rules, determine the optimal infrastructure and return a JSON object that conforms to this schema:

{schema_json}

**CRITICAL REQUIREMENTS**:

1. **confidence_score** (0-1): REQUIRED. Reflects your certainty in this recommendation.
   - 0.9-1.0: High confidence (known sizes, clear thresholds, well-understood operation)
   - 0.7-0.89: Medium confidence (some unknowns, but reasonable estimate)
   - 0.5-0.69: Low confidence (multiple unknowns, unclear requirements)
   - <0.5: Very low confidence (recommend human review)

2. **decision_summary**: REQUIRED. 1-2 sentence high-level justification.

3. **cost_analysis.estimated_cost_range_usd**: REQUIRED. Provide a RANGE [min, max], not exact value.
   - Include explicit assumptions (region, instance type, runtime)
   - Set cost_confidence to reflect uncertainty

4. **warnings**: REQUIRED if:
   - File sizes are unknown (unknown_sizes > 0)
   - Confidence < 0.7
   - Operation may require human review

5. **alternatives**: At least 1-2 alternative options with tradeoffs.

"""

    # If this is a retry after validation error, include the error
    if validation_error:
        prompt += f"""
**VALIDATION ERROR FROM PREVIOUS ATTEMPT**:
Your previous response failed validation with this error:
{validation_error}

Please correct the error and return valid JSON that matches the schema above.
"""
    
    prompt += """
Return ONLY valid JSON that matches the schema exactly. No additional text before or after the JSON.
"""
    
    return prompt


async def decide_infrastructure(
    command: str,
    inputs: Optional[List[InputAsset]] = None,
    outputs: Optional[List[OutputAsset]] = None,
    session_context: Optional[Dict[str, Any]] = None,
    operation_type: Optional[str] = None,
    workflow_plan: Optional[Any] = None,  # WorkflowPlan Pydantic model
    request_id: Optional[str] = None
) -> InfraDecision:
    """
    Decide optimal infrastructure for executing a bioinformatics operation.
    
    This is the main entry point for the Infrastructure Decision Agent v2.
    
    Key improvements over v1:
    - Returns Pydantic InfraDecision model (v1 returned dataclass)
    - Automatic validation with repair mechanism (max 2 retries)
    - Confidence scoring on all decisions
    - Cost ranges instead of exact values
    
    Args:
        command: The user's command/request
        inputs: List of input files/assets with sizes (DEPRECATED: use workflow_plan)
        outputs: List of output files/assets (optional)
        session_context: Session context (optional)
        operation_type: Type of operation (e.g., "read_merging") - optional
        workflow_plan: WorkflowPlan Pydantic model with rich dataset metadata (PREFERRED)
        request_id: Request ID for tracing (optional)
    
    Returns:
        InfraDecision Pydantic model with validated infrastructure recommendation
    
    Raises:
        ValueError: If LLM fails to produce valid output after max retries
    
    Note:
        If workflow_plan is provided, it takes precedence over inputs/outputs.
        workflow_plan.data_inputs contains size_bytes, location_type, and other metadata
        that enable more accurate infrastructure decisions.
    """
    # Extract dataset info from workflow_plan if provided (PREFERRED)
    if workflow_plan is not None:
        # Convert WorkflowPlan.data_inputs to InputAsset objects with full metadata
        inputs = []
        for data_input in workflow_plan.data_inputs:
            input_asset = InputAsset(
                uri=data_input.uri,
                size_bytes=data_input.size_bytes,
                source=data_input.location_type if data_input.location_type != "Unknown" else "unknown"
            )
            inputs.append(input_asset)
        
        # Also extract command from workflow_plan if not provided
        if not command or command == "":
            command = workflow_plan.description
        
        logger.info(
            f"[{request_id}] Using workflow_plan data: {len(inputs)} inputs "
            f"({sum(1 for i in inputs if i.size_bytes is not None)} with known sizes)"
        )
    else:
        inputs = inputs or []
        outputs = outputs or []
    
    # Analyze files (returns Pydantic FileAnalysis)
    file_analysis = _analyze_files(inputs, outputs or [])
    
    # Try LLM-based decision with repair mechanism
    max_retries = 2
    validation_error = None
    
    for attempt in range(max_retries + 1):
        try:
            llm = _get_llm()
            
            # Build prompt (includes validation error if retry)
            user_prompt = _build_decision_prompt(
                command, file_analysis, inputs, outputs, session_context, validation_error
            )
            
            messages = [
                {"role": "system", "content": INFRASTRUCTURE_DECISION_SYSTEM_PROMPT},
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
            
            # Parse JSON response
            json_start = response_content.find("{")
            json_end = response_content.rfind("}") + 1
            
            if json_start >= 0 and json_end > json_start:
                json_str = response_content[json_start:json_end]
                result_dict = json.loads(json_str)
            else:
                result_dict = json.loads(response_content)
            
            # Validate with Pydantic (this is the key improvement!)
            infra_decision = InfraDecision(**result_dict)
            
            # Success! Return validated decision
            logger.info(
                f"✅ Infrastructure Decision (attempt {attempt+1}): "
                f"{infra_decision.infrastructure} "
                f"(confidence: {infra_decision.confidence_score:.2f})"
            )
            return infra_decision
            
        except ValidationError as e:
            # Validation failed - prepare error message for retry
            validation_error = _format_validation_error(e)
            logger.warning(
                f"⚠️  Validation error on attempt {attempt+1}/{max_retries+1}: {validation_error}"
            )
            
            if attempt < max_retries:
                logger.info(f"🔄 Retrying with validation error feedback (attempt {attempt+2})...")
                continue
            else:
                # Max retries exceeded - log and fall back to heuristic
                logger.error(
                    f"❌ Failed to get valid infrastructure decision after {max_retries+1} attempts. "
                    "Falling back to heuristic decision."
                )
                return _heuristic_infrastructure_decision(command, file_analysis, inputs, outputs)
        
        except json.JSONDecodeError as e:
            logger.warning(f"⚠️  JSON decode error on attempt {attempt+1}: {e}")
            validation_error = f"Invalid JSON: {str(e)}"
            
            if attempt < max_retries:
                logger.info(f"🔄 Retrying with JSON error feedback (attempt {attempt+2})...")
                continue
            else:
                logger.error(f"❌ Failed to parse JSON after {max_retries+1} attempts. Falling back to heuristic.")
                return _heuristic_infrastructure_decision(command, file_analysis, inputs, outputs)
        
        except Exception as e:
            logger.warning(f"❌ LLM infrastructure decision failed on attempt {attempt+1}: {e}")
            if attempt < max_retries:
                validation_error = f"Unexpected error: {str(e)}"
                logger.info(f"🔄 Retrying (attempt {attempt+2})...")
                continue
            else:
                logger.error(f"❌ All attempts failed. Falling back to heuristic decision.")
                return _heuristic_infrastructure_decision(command, file_analysis, inputs, outputs)
    
    # Should never reach here, but fallback just in case
    return _heuristic_infrastructure_decision(command, file_analysis, inputs, outputs)


def _format_validation_error(error: ValidationError) -> str:
    """
    Format Pydantic ValidationError into a concise error message for LLM feedback.
    
    Args:
        error: Pydantic ValidationError
    
    Returns:
        Formatted error string
    """
    errors = []
    for err in error.errors():
        loc = ".".join(str(l) for l in err['loc'])
        msg = err['msg']
        errors.append(f"  - {loc}: {msg}")
    
    return f"Pydantic validation failed:\n" + "\n".join(errors)


def _heuristic_infrastructure_decision(
    command: str,
    file_analysis: FileAnalysis,
    inputs: List[InputAsset],
    outputs: List[OutputAsset]
) -> InfraDecision:
    """
    Fallback heuristic infrastructure decision when LLM is unavailable or fails validation.
    
    Uses simple rules based on file sizes and locations, but now returns Pydantic InfraDecision.
    
    Note: This is a fallback - confidence will be low (<0.5) to signal heuristic use.
    """
    total_size_bytes = file_analysis.total_size_bytes
    total_size_mb = file_analysis.total_size_mb
    locations = file_analysis.locations
    unknown_sizes = file_analysis.unknown_sizes
    
    # Default threshold (100MB)
    emr_threshold = int(os.getenv('HELIX_ASYNC_BYTES_THRESHOLD', 100 * 1024 * 1024))
    use_ec2 = os.getenv('HELIX_USE_EC2', 'false').lower() == 'true'
    
    infrastructure = "Local"
    reasoning = ""
    warnings: list[str] = []
    # For many inputs (known size + location), the heuristic is deterministic and reliable.
    confidence = 0.75
    
    # Check if S3 files are large
    if "S3" in locations:
        if total_size_bytes > emr_threshold:
            infrastructure = "EMR"
            reasoning = f"Files on S3 are large ({total_size_mb:.2f} MB > {emr_threshold / (1024*1024):.0f} MB). EMR recommended to avoid large data transfers and process where data lives."
            confidence = 0.9 if unknown_sizes == 0 else 0.6
        elif use_ec2:
            infrastructure = "EC2"
            reasoning = f"Files on S3 are small ({total_size_mb:.2f} MB < {emr_threshold / (1024*1024):.0f} MB). EC2 recommended for fast execution with pre-installed tools."
            confidence = 0.85 if unknown_sizes == 0 else 0.55
        else:
            infrastructure = "Local"
            reasoning = f"Files on S3 are small ({total_size_mb:.2f} MB < {emr_threshold / (1024*1024):.0f} MB). Local execution recommended after downloading files."
            confidence = 0.8 if unknown_sizes == 0 else 0.5
    elif "Local" in locations:
        if total_size_bytes > 10 * 1024 * 1024 * 1024:  # >10GB
            infrastructure = "EMR"
            reasoning = f"Local files are very large ({total_size_mb:.2f} MB > 10GB). Recommend uploading to S3 and using EMR for distributed processing."
            warnings.append("Large local files should be uploaded to S3 before EMR execution")
            confidence = 0.75
        elif use_ec2 and total_size_bytes > 0:
            infrastructure = "EC2"
            reasoning = f"Local files ({total_size_mb:.2f} MB) can be processed on EC2 with pre-installed tools."
            confidence = 0.75
        else:
            infrastructure = "Local"
            reasoning = f"Local files are small ({total_size_mb:.2f} MB). Local execution is fastest."
            confidence = 0.9
    
    if unknown_sizes > 0:
        warnings.append(f"{unknown_sizes} file(s) have unknown sizes - decision based on known files only")
        confidence *= 0.8  # Reduce confidence when sizes are unknown
    
    if infrastructure == "EMR" and not use_ec2:
        warnings.append("EMR routing for tool-generator-agent may not be fully implemented")
    
    # Build cost analysis (rough estimates for heuristic)
    if infrastructure == "EMR":
        cost_range = (2.0, 10.0)
    elif infrastructure == "EC2":
        cost_range = (0.5, 2.0)
    else:
        cost_range = (0.0, 0.1)
    
    decision_summary = f"Heuristic decision: {infrastructure} based on file size ({total_size_mb:.2f} MB) and location ({', '.join(locations)})"

    # Provide at least one reasonable alternative (useful for UI + tests).
    alternatives: list[InfraAlternative] = []
    if infrastructure == "EMR":
        alternatives.append(
            InfraAlternative(
                infrastructure="EC2",
                reasoning="Could run on a single EC2 instance after downloading inputs (simpler operationally).",
                tradeoffs="Lower setup overhead, but may be slower and requires data transfer / more manual tuning.",
                confidence=0.4,
            )
        )
    elif infrastructure == "EC2":
        alternatives.append(
            InfraAlternative(
                infrastructure="Local",
                reasoning="Could run locally for development or small inputs.",
                tradeoffs="Fast iteration, but limited resources and may not scale for larger datasets.",
                confidence=0.3,
            )
        )
    elif infrastructure == "Local":
        alternatives.append(
            InfraAlternative(
                infrastructure="EC2",
                reasoning="Could run on EC2 if local resources are constrained.",
                tradeoffs="More compute options, but introduces AWS setup and cost.",
                confidence=0.3,
            )
        )
    
    return InfraDecision(
        infrastructure=infrastructure,
        confidence_score=round(confidence, 2),
        decision_summary=decision_summary,
        reasoning=reasoning,
        file_analysis=file_analysis,
        computational_requirements=ComputationalRequirements(
            estimated_cpu_cores=1,
            estimated_memory_gb=max(2.0, total_size_mb / 500),
            estimated_runtime_minutes=5.0,
            parallelizable=False
        ),
        cost_analysis=CostAnalysis(
            estimated_cost_range_usd=cost_range,
            cost_assumptions="Rough heuristic estimate, not based on detailed analysis",
            cost_confidence=0.3,  # Low confidence in heuristic costs
            data_transfer_cost_usd=0.0
        ),
        alternatives=alternatives,
        warnings=warnings,
        inputs_analyzed=file_analysis.file_count
    )


def infrastructure_decision_to_dict(decision: InfraDecision) -> Dict[str, Any]:
    """
    Convert InfraDecision Pydantic model to dictionary.
    
    Use .model_dump() for Pydantic V2 models (v1 used .dict()).
    """
    return decision.model_dump()


# For backward compatibility, keep the old dataclass type name
# But map it to the new Pydantic model
InfrastructureDecision = InfraDecision
