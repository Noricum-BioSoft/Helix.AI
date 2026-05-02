"""
Tool Generator Agent - Dynamically generates and executes bioinformatics tools
when no pre-existing tool exists in the toolbox.

This agent follows the workflow defined in agents/tool-generator-agent.md:
1. Task Analysis
2. Tool Research
3. Infrastructure Decision
4. Tool Function Generation
5. Code Execution
"""

import os
import sys
import uuid
import asyncio
import logging
import tempfile
import subprocess
import time
import shutil
from pathlib import Path
from typing import Dict, Any, Optional, List
import json

logger = logging.getLogger(__name__)

# Load the tool-generator-agent prompt
PROJECT_ROOT = Path(__file__).resolve().parent.parent
TOOL_GENERATOR_PROMPT_PATH = PROJECT_ROOT / "agents" / "tool-generator-agent.md"


def _discover_inputs_from_args(arguments: Dict[str, Any], session_context: Optional[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
    """
    Discover input files from arguments and session context.
    Returns a list of dicts with 'uri' and 'size_bytes' keys.
    
    This is a simplified version of ExecutionBroker._discover_inputs that doesn't require
    instantiating ExecutionBroker.
    """
    from backend.execution_broker import ExecutionBroker
    
    # Create a temporary broker instance just to use its _discover_inputs method
    # We pass a dummy tool executor since we only need input discovery
    async def dummy_executor(tool_name: str, args: Dict[str, Any]) -> Dict[str, Any]:
        return {"status": "error", "message": "dummy executor"}
    
    broker = ExecutionBroker(tool_executor=dummy_executor)
    input_assets = broker._discover_inputs(arguments or {}, session_context or {})
    
    # Convert InputAsset objects to dicts for easier serialization
    return [
        {
            "uri": asset.uri,
            "size_bytes": asset.size_bytes,
            "source": asset.source
        }
        for asset in input_assets
    ]


def _discover_outputs_from_args(arguments: Dict[str, Any], command: Optional[str] = None) -> List[Dict[str, Any]]:
    """
    Discover output files/paths from arguments and command text.
    Returns a list of dicts with 'uri' and optional metadata.
    
    Looks for:
    - Explicit output_path or output arguments
    - S3 paths in command text following patterns like "output ... on s3://..." or "save to s3://..."
    """
    import re
    outputs = []
    
    # 1. Check explicit arguments
    output_path = arguments.get("output_path") or arguments.get("output")
    if output_path and isinstance(output_path, str):
        if output_path.startswith("s3://") or output_path.startswith("/"):
            outputs.append({
                "uri": output_path,
                "source": "args"
            })
    
    # 2. Extract from command text if provided
    if command:
        # Pattern 1: "output ... on s3://..." or "output: s3://..."
        output_patterns = [
            r'(?:output|save|write|upload)[^\n]*?(s3://[^\s]+)',
            r'\bon\s+(s3://[^\s]+)',
            r'output\s*:\s*(s3://[^\s]+|/[^\s]+)',
            r'(?:qc\s*reports?|fastqc)\s*[^\n:]*:\s*(s3://[^\s]+|/[^\s]+)',
            r'(?:alignments?|aligned)\s*[^\n:]*:\s*(s3://[^\s]+|/[^\s]+)',
            r'(?:counts?|gene\s*counts)\s*[^\n:]*:\s*(s3://[^\s]+|/[^\s]+)',
        ]
        
        for pattern in output_patterns:
            matches = re.finditer(pattern, command, re.IGNORECASE)
            for match in matches:
                uri = match.group(1).strip()
                # Avoid duplicates
                if not any(o.get("uri") == uri for o in outputs):
                    outputs.append({
                        "uri": uri,
                        "source": "command_text"
                    })
    
    return outputs

try:
    TOOL_GENERATOR_SYSTEM_PROMPT = TOOL_GENERATOR_PROMPT_PATH.read_text()
except Exception as e:
    logger.warning(f"Could not load tool-generator-agent.md: {e}")
    TOOL_GENERATOR_SYSTEM_PROMPT = """You are an advanced coding assistant specialized in bioinformatics workflows.
Your core capability is to dynamically research, design, and implement solutions for bioinformatics tasks."""


def _get_llm():
    """Initialize and return the LLM for tool generation."""
    try:
        openai_key = os.getenv("OPENAI_API_KEY", "").strip()
        deepseek_key = os.getenv("DEEPSEEK_API_KEY", "").strip()
        # Allow overriding the model without code changes.
        # LangChain init_chat_model expects provider-prefixed names like "openai:<model>".
        # Default to a chat-capable model. Some code-focused models are *not* exposed via chat-completions
        # and will 404 when called through the Chat Completions endpoint.
        openai_model = os.getenv("HELIX_TOOLGEN_OPENAI_MODEL", "openai:gpt-5.5").strip()
        
        openai_enabled = openai_key and openai_key not in ["", "disabled", "your_openai_api_key_here", "none"]
        deepseek_enabled = deepseek_key and deepseek_key not in ["", "disabled", "your_deepseek_api_key_here", "none"]
        
        if openai_enabled:
            # Local import to avoid importing optional SSL/cert deps during test collection.
            from langchain.chat_models import init_chat_model
            # Be forgiving: allow HELIX_TOOLGEN_OPENAI_MODEL="gpt-5.5" as well (add prefix if missing).
            if ":" not in openai_model:
                openai_model = f"openai:{openai_model}"
            return init_chat_model(openai_model, temperature=0)
        elif deepseek_enabled:
            # Local import to avoid importing optional SSL/cert deps during test collection.
            from langchain_deepseek import ChatDeepSeek
            return ChatDeepSeek(
                model="deepseek-chat",
                temperature=0,
                max_tokens=None,
                timeout=120.0,  # Longer timeout for code generation
                max_retries=1,
            )
        else:
            raise ValueError("No API keys found for tool generation")
    except Exception as e:
        logger.error(f"Error initializing LLM for tool generation: {e}")
        raise


async def generate_and_execute_tool(
    command: str,
    user_request: str,
    session_id: Optional[str] = None,
    inputs: Optional[List[Any]] = None,
    outputs: Optional[List[Any]] = None,
    session_context: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Generate and execute a bioinformatics tool for the given command.
    
    Args:
        command: The user's bioinformatics command/request
        user_request: The original user request (for context)
        session_id: Optional session ID for tracking
        inputs: Optional list of InputAsset objects with discovered input files and their sizes
        outputs: Optional list of output paths/URIs where results should be written
        
    Returns:
        Dictionary with execution results
    """
    logger.info(f"🔧 Tool Generator Agent: Generating tool for command: {command}")

    # Phase 4 safety: do not generate tools for pure Q&A intent.
    try:
        from backend.intent_classifier import classify_intent

        intent = classify_intent(user_request or command)
        if intent.intent != "execute":
            return {
                "status": "skipped",
                "tool_generated": False,
                "result": {},
                "explanation": (
                    "Tool generation skipped because the request looks like Q&A intent. "
                    "Rephrase as an execution request (e.g. 'run ...', 'analyze ...')."
                ),
                "intent": intent.intent,
                "intent_reason": intent.reason,
            }
    except Exception:
        # If classifier fails for any reason, fall through to previous behavior.
        pass
    
    try:
        llm = _get_llm()
        
        # Build input metadata section if inputs are provided
        input_metadata_section = ""
        if inputs:
            total_size = 0
            known_sizes = []
            unknown_sizes = []
            
            for inp in inputs:
                # Handle both InputAsset objects and dicts
                if hasattr(inp, 'uri'):
                    uri = inp.uri
                    size_bytes = inp.size_bytes
                elif isinstance(inp, dict):
                    uri = inp.get('uri', '')
                    size_bytes = inp.get('size_bytes')
                else:
                    continue
                
                if size_bytes is not None:
                    total_size += size_bytes
                    size_mb = size_bytes / (1024 * 1024)
                    known_sizes.append(f"  - {uri}: {size_mb:.2f} MB ({size_bytes:,} bytes)")
                else:
                    unknown_sizes.append(f"  - {uri}: size unknown")
            
            if known_sizes or unknown_sizes:
                input_metadata_section = "\n\n**Input File Information:**\n"
                if known_sizes:
                    input_metadata_section += "Files with known sizes:\n" + "\n".join(known_sizes) + "\n"
                    total_mb = total_size / (1024 * 1024)
                    input_metadata_section += f"\nTotal size of known files: {total_mb:.2f} MB ({total_size:,} bytes)\n"
                if unknown_sizes:
                    input_metadata_section += "\nFiles with unknown sizes:\n" + "\n".join(unknown_sizes) + "\n"
                    input_metadata_section += "(Size unknown due to permissions or file not found)\n"
                
                # Add infrastructure recommendation based on file sizes
                if total_size > 0:
                    if total_size > 100 * 1024 * 1024:  # >100MB
                        input_metadata_section += "\n**Infrastructure Recommendation**: Based on file sizes (>100MB), AWS EMR is recommended to avoid large data transfers.\n"
                    elif total_size > 10 * 1024 * 1024:  # >10MB
                        input_metadata_section += "\n**Infrastructure Recommendation**: Files are medium-sized (10-100MB). Consider AWS Batch or local execution depending on compute requirements.\n"
                    else:
                        input_metadata_section += "\n**Infrastructure Recommendation**: Files are small (<10MB). Local execution or AWS Batch are suitable.\n"
        
        # Build output metadata section if outputs are provided
        output_metadata_section = ""
        if outputs:
            output_paths = []
            for out in outputs:
                # Handle both dicts and strings
                if isinstance(out, dict):
                    uri = out.get('uri', '')
                elif isinstance(out, str):
                    uri = out
                else:
                    continue
                
                if uri:
                    output_paths.append(f"  - {uri}")
            
            if output_paths:
                output_metadata_section = "\n\n**Output File Information:**\n"
                output_metadata_section += "Expected output locations:\n" + "\n".join(output_paths) + "\n"
                output_metadata_section += "\n**IMPORTANT**: Ensure your generated code writes the results to these exact output paths. Verify the output location matches the user's request.\n"
        
        # Get infrastructure decision from Infrastructure Decision Agent
        infrastructure_decision = None
        infrastructure_section = ""
        try:
            from backend.infrastructure_decision_agent import (
                decide_infrastructure,
                InputAsset,
                OutputAsset,
                infrastructure_decision_to_dict
            )
            
            # Convert inputs/outputs to InputAsset/OutputAsset objects
            input_assets = []
            for inp in inputs or []:
                if hasattr(inp, 'uri'):
                    input_assets.append(InputAsset(uri=inp.uri, size_bytes=getattr(inp, 'size_bytes', None), source=getattr(inp, 'source', 'unknown')))
                elif isinstance(inp, dict):
                    input_assets.append(InputAsset(uri=inp.get('uri', ''), size_bytes=inp.get('size_bytes'), source=inp.get('source', 'unknown')))
            
            output_assets = []
            for out in outputs or []:
                if isinstance(out, dict):
                    output_assets.append(OutputAsset(uri=out.get('uri', ''), estimated_size_bytes=out.get('estimated_size_bytes')))
                elif isinstance(out, str):
                    output_assets.append(OutputAsset(uri=out))
            
            logger.info(f"🔧 Calling Infrastructure Decision Agent for: {command}")
            infrastructure_decision = await decide_infrastructure(
                command=command,
                inputs=input_assets,
                outputs=output_assets,
                session_context=session_context
            )
            
            infra_dict = infrastructure_decision_to_dict(infrastructure_decision)
            infrastructure_section = f"""

**Infrastructure Decision (from Infrastructure Decision Agent):**
- **Infrastructure**: {infrastructure_decision.infrastructure}
- **Reasoning**: {infrastructure_decision.reasoning}
- **File Analysis**: Total size: {infrastructure_decision.file_analysis.total_size_mb:.2f} MB, Locations: {', '.join(infrastructure_decision.file_analysis.locations)}
- **Warnings**: {', '.join(infrastructure_decision.warnings) if infrastructure_decision.warnings else 'None'}

**IMPORTANT**: You MUST generate code that executes on **{infrastructure_decision.infrastructure}**. Do NOT make infrastructure decisions - use the provided decision.
"""
            logger.info(f"✅ Infrastructure Decision: {infrastructure_decision.infrastructure} - {infrastructure_decision.reasoning}")
        except Exception as e:
            logger.warning(f"⚠️  Infrastructure Decision Agent failed: {e}, will proceed without infrastructure decision")
            infrastructure_section = "\n\n**Note**: Infrastructure decision not available. Use file size information to make a reasonable infrastructure choice."
        
        # Build the prompt for tool generation
        user_prompt = f"""Generate and execute a solution for this bioinformatics operation:

User Command: {command}

Original Request: {user_request}
{input_metadata_section}{output_metadata_section}{infrastructure_section}
Please follow the workflow:
1. Analyze the task and identify the biological operation
2. Research appropriate bioinformatics tools
3. Generate complete, executable Python code for the specified infrastructure ({infrastructure_decision.infrastructure if infrastructure_decision else 'determine based on file sizes'})
4. Execute the code and return results

Your response should include:
- Tool selection and justification
- Complete Python implementation (for {infrastructure_decision.infrastructure if infrastructure_decision else 'the chosen infrastructure'})
- Execution results

CRITICAL REQUIREMENTS:
1. **PRIORITIZE established bioinformatics tools**: Use well-established, industry-standard tools FIRST:
   - For read merging: Prefer BBMerge (BBTools), FLASH, or PEAR - these are industry-standard tools
   - For sequence operations: Prefer established tools (samtools, bcftools, BWA, etc.)
   - Only use Python implementations as a LAST RESORT when tools cannot be installed/accessed
2. **Tool installation**: If tools are not available:
   - First check if conda is available: `shutil.which("conda")`
   - If conda is available, attempt installation: `conda install -c bioconda bbtools`
   - If conda is NOT available, skip installation and use Python fallback immediately
   - Do not attempt installation if conda is not available - it will always fail
3. **ALWAYS check tool availability**: Before using any external tool, check if it exists using `shutil.which()`, attempt installation if missing, provide Python fallback only as last resort
4. **For read merging with S3 files**: If BBMerge/FLASH/PEAR are not available and cannot be installed, you MUST use the built-in Python implementation. **IT IS FORBIDDEN TO RETURN EARLY WITHOUT MERGING**. Here is the EXACT code pattern you MUST follow:
   ```python
   import boto3
   import logging
   import os
   import shutil
   import subprocess
   
   logging.basicConfig(level=logging.INFO)
   logger = logging.getLogger(__name__)
   
   def check_tool_available(tool_name: str) -> bool:
       return shutil.which(tool_name) is not None
   
   def install_tool_via_conda(tool_name: str, conda_package: str = None) -> bool:
       if not shutil.which("conda"):
           return False
       package = conda_package or tool_name
       try:
           result = subprocess.run(
               ["conda", "install", "-c", "bioconda", package, "-y"],
               capture_output=True, text=True, timeout=300
           )
           return result.returncode == 0
       except Exception:
           return False
   
   # Main execution - this is the CRITICAL part
   # Extract input and output paths from the provided arguments/outputs
   # Example: r1_path = "s3://bucket/path/R1.fq"  # Extract from inputs or command
   # Example: r2_path = "s3://bucket/path/R2.fq"  # Extract from inputs or command
   # Example: output_path = "s3://bucket/path/merged.fq"  # Extract from outputs or command
   
   # IMPORTANT: Extract output_path from the outputs list if provided
   # If outputs is a list, get the first output URI
   # If outputs is a single string/dict, extract the URI
   # Example:
   #   if outputs and len(outputs) > 0:
   #       output_path = outputs[0].get('uri') if isinstance(outputs[0], dict) else outputs[0]
   #   else:
   #       output_path = "s3://bucket/path/merged.fq"  # Fallback or extract from command
   
   # Check if BBMerge is available
   if not check_tool_available("bbmerge.sh"):
       logger.info("BBMerge not found, attempting installation via conda...")
       if not install_tool_via_conda("bbtools"):
           logger.warning("Could not install BBTools via conda, using Python implementation")
           # CRITICAL: You MUST call merge_reads_from_s3 here - DO NOT return!
           try:
               from tools.read_merging import merge_reads_from_s3
               result = merge_reads_from_s3(
                   r1_path=r1_path,
                   r2_path=r2_path,
                   output_path=output_path,
                   min_overlap=12
               )
               if result.get("status") == "success":
                   print(f"SUCCESS: Merged reads written to {{output_path}}")
               else:
                   raise Exception(result.get("error", "Merge failed"))
           except ImportError:
               raise RuntimeError("merge_reads_from_s3 not available and inline implementation not provided")
       else:
           # BBMerge was installed, use it (implementation here)
           pass
   else:
       # BBMerge is available, use it (implementation here)
       pass
   ```
   
   **CRITICAL RULES**:
   - NEVER write `return` after "falling back to Python implementation" - you MUST call merge_reads_from_s3
   - The function merge_reads_from_s3 MUST be called when BBMerge is unavailable
   - You MUST extract output_path from the outputs list or command before using it
   - You MUST print "SUCCESS: Merged reads written to {{output_path}}" when merge completes (use the actual output_path variable)
   - If you return early without merging, the operation will be detected as a failure
   The `merge_reads_from_s3` function handles:
   - S3 file downloads
   - FASTQ parsing with quality scores
   - Sequence merging with overlap detection
   - Quality score merging (taking max in overlap regions)
   - FASTQ output writing
   - S3 upload
   
   **CRITICAL RULE**: The merge function MUST complete the merge operation. It is FORBIDDEN to return early without performing the merge. If BBMerge is unavailable, you MUST call merge_reads_from_s3 or implement the merge inline. Returning without merging is a fatal error that will cause the task to fail.
   
   **VALIDATION**: Your code MUST print "SUCCESS: Merged reads written to {{output_path}}" when the merge completes (where output_path is the actual variable containing the output path). If this message is not printed, the merge did not execute.
4. **BioPython usage**: If using BioPython, ALWAYS create new SeqRecord objects - never modify existing ones:
   ```python
   # CORRECT: Create new SeqRecord
   record = SeqRecord(Seq(sequence), id="read1", description="")
   
   # WRONG: Don't modify existing record.seq if it has annotations
   # record.seq = Seq(new_sequence)  # This will fail if record has letter annotations
   ```
5. **Error handling**: If an external tool is not found, gracefully fall back to a Python implementation
6. **Tool availability check pattern**:
   ```python
   import shutil
   if shutil.which("tool_name") is None:
       # Use Python implementation instead
       logger.warning("External tool not available, using Python fallback")
   ```

IMPORTANT: Generate ONLY executable Python code that can be run directly. Include all necessary imports, error handling, and AWS credential handling if needed. DO NOT assume external bioinformatics tools are installed - always check availability and provide Python fallbacks.

**AWS Credentials and S3 Access:**
- When executing on EC2/ECS, AWS credentials are automatically available via IAM instance/task role
- Simply use `boto3.client('s3')` without explicit credentials - boto3 automatically uses IAM role credentials
- **IMPORTANT for S3 operations:**
  1. S3 buckets can be in different regions - always detect bucket region before operations
  2. Use `s3.get_bucket_location()` to detect the bucket's region, then create a client for that region
  3. Handle permission errors gracefully - 403 Forbidden means the IAM role lacks permissions
  4. Example for S3 operations:
  ```python
  import boto3
  import os
  from botocore.exceptions import ClientError
  
  # Default region from environment
  default_region = os.getenv('AWS_REGION', os.getenv('AWS_DEFAULT_REGION', 'us-east-1'))
  
  # Create S3 client (S3 works with any region, but use default for metadata operations)
  s3_meta = boto3.client('s3', region_name=default_region)
  
  # For bucket operations, detect the bucket's actual region
  def get_bucket_region(bucket_name: str) -> str:
      try:
          response = s3_meta.head_bucket(Bucket=bucket_name)
          # HeadBucket doesn't return region in response, use get_bucket_location instead
          location = s3_meta.get_bucket_location(Bucket=bucket_name)
          # get_bucket_location returns None for us-east-1, 'EU' for Europe, etc.
          region = location.get('LocationConstraint') or 'us-east-1'
          return region
      except ClientError as e:
          error_code = e.response.get('Error', {{}}).get('Code', '')
          if error_code == '403':
              raise PermissionError(f"Access denied to bucket {{bucket_name}}. Check IAM role permissions.")
          raise
  
  # Use the bucket's region for operations (better performance and avoids issues)
  bucket_region = get_bucket_region('my-bucket')
  s3 = boto3.client('s3', region_name=bucket_region)
  # No credentials needed - uses IAM role automatically
  ```
"""
        
        # Get LLM response
        messages = [
            {"role": "system", "content": TOOL_GENERATOR_SYSTEM_PROMPT},
            {"role": "user", "content": user_prompt}
        ]
        
        logger.info("🔧 Tool Generator Agent: Requesting tool generation from LLM...")
        response = await asyncio.to_thread(llm.invoke, messages)
        
        # Extract content from response
        if hasattr(response, 'content'):
            content = response.content
        elif isinstance(response, dict):
            content = response.get('content', str(response))
        else:
            content = str(response)
        
        logger.info(f"🔧 Tool Generator Agent: Received response ({len(content)} chars)")
        
        # Try to extract Python code from the response
        code = _extract_python_code(content)
        
        if not code:
            # If no code block found, try to use the entire response as code
            # (the LLM might have generated code without markdown formatting)
            logger.warning("⚠️  No code block found in response, attempting to use full response")
            code = content
        
        # Check if we should route to EMR based on file sizes
        # Calculate total file size from inputs for execution routing
        total_size = 0
        if inputs:
            for inp in inputs:
                if hasattr(inp, 'size_bytes') and isinstance(inp.size_bytes, int):
                    total_size += inp.size_bytes
                elif isinstance(inp, dict) and isinstance(inp.get('size_bytes'), int):
                    total_size += inp['size_bytes']
        
        # Self-healing execution with retry on failure
        max_retries = 2
        execution_result = None
        
        for attempt in range(max_retries + 1):
            if attempt > 0:
                logger.info(f"🔄 Retry attempt {attempt}/{max_retries} - Regenerating code with error feedback...")
                
                # Build error feedback for LLM
                error_feedback = f"""
**PREVIOUS ATTEMPT FAILED** (Attempt {attempt}):

**Error Details:**
- Exit Code: {execution_result.get('exit_code', 'unknown')}
- Error Message: {execution_result.get('error', 'unknown')}

**stderr Output:**
```
{execution_result.get('stderr_full', execution_result.get('stderr', 'N/A'))[:2000]}
```

**stdout Output:**
```
{execution_result.get('stdout', 'N/A')[:1000]}
```

**Common Issues to Fix:**
1. If "STAR not found" or "tool not available": DO NOT try to install tools - they are not in the Docker image. Use pure Python implementation with BioPython/pysam instead.
2. If "Missing output" error: Check that you're reading OUTPUT_PATH/OUTPUT_DIR environment variables correctly
3. If "Missing input" error: Check that you're reading INPUT_R1/INPUT_R2/INPUT_FILES environment variables correctly
4. If "conda not available": Skip conda installation entirely - generate pure Python fallback immediately

**CRITICAL**: Please generate a FIXED version of the code that addresses the above error. Focus on the root cause.
"""
                
                # Regenerate with error feedback
                retry_prompt = user_prompt + error_feedback
                response = await llm.ainvoke([
                    {"role": "system", "content": TOOL_GENERATOR_SYSTEM_PROMPT},
                    {"role": "user", "content": retry_prompt}
                ])
                
                # Extract code from retry response
                retry_content = response.content if hasattr(response, 'content') else str(response)
                code = _extract_python_code(retry_content)
                if not code:
                    logger.error("❌ Failed to extract code from retry response")
                    # Try using full response as code
                    code = retry_content
                
                logger.info(f"✅ Regenerated code ({len(code)} chars) for retry attempt {attempt}")
            
            # Execute code (common for both first attempt and retries)
            if infrastructure_decision:
                infra = infrastructure_decision.infrastructure.upper()
                if attempt == 0:
                    logger.info(f"🔧 Tool Generator Agent: Executing on {infra} (from Infrastructure Decision Agent)")
                
                if infra == "EMR":
                    logger.warning(f"⚠️  Infrastructure Decision Agent recommended EMR, but EMR routing for tool-generator-agent is not yet fully implemented.")
                    logger.warning(f"⚠️  This should be routed through execution_broker. For now, executing on EC2 (may be slow for large files).")
                    execution_result = await _execute_generated_code(
                        code,
                        command,
                        session_id,
                        total_file_size_bytes=total_size,
                        inputs=inputs,
                        outputs=outputs
                    )
                elif infra == "EC2":
                    if attempt == 0:
                        logger.info("🔧 Executing generated code on EC2...")
                    execution_result = await _execute_generated_code(
                        code,
                        command,
                        session_id,
                        total_file_size_bytes=total_size,
                        inputs=inputs,
                        outputs=outputs
                    )
                elif infra == "LOCAL":
                    if attempt == 0:
                        logger.info("🔧 Executing generated code locally...")
                    execution_result = await _execute_generated_code(
                        code,
                        command,
                        session_id,
                        total_file_size_bytes=total_size,
                        inputs=inputs,
                        outputs=outputs
                    )
                else:
                    # Batch, Lambda, or unknown - default to local/EC2
                    if attempt == 0:
                        logger.info(f"🔧 Infrastructure {infra} not yet fully supported, executing on EC2/Local...")
                    execution_result = await _execute_generated_code(
                        code,
                        command,
                        session_id,
                        total_file_size_bytes=total_size,
                        inputs=inputs,
                        outputs=outputs
                    )
            else:
                # Fallback: use simple threshold-based routing if infrastructure decision unavailable
                emr_threshold = int(os.getenv('HELIX_ASYNC_BYTES_THRESHOLD', 100 * 1024 * 1024))
                if attempt == 0 and total_size > emr_threshold:
                    logger.warning(f"⚠️  Files are large ({total_size / (1024*1024):.2f} MB > {emr_threshold / (1024*1024):.2f} MB)")
                    logger.warning(f"⚠️  EMR routing for tool-generator-agent is not yet implemented. Executing on EC2 (may be slow).")
                if attempt == 0:
                    logger.info("🔧 Tool Generator Agent: Executing generated code...")
                execution_result = await _execute_generated_code(
                    code,
                    command,
                    session_id,
                    total_file_size_bytes=total_size,
                    inputs=inputs,
                    outputs=outputs
                )
            
            # Check if execution succeeded
            if isinstance(execution_result, dict):
                exec_status = execution_result.get("status", "error")
                if exec_status == "success":
                    if attempt > 0:
                        logger.info(f"✅ Code fixed successfully after {attempt} retry(ies)!")
                    break  # Success - exit retry loop
                else:
                    # Execution failed
                    if attempt < max_retries:
                        logger.warning(f"⚠️  Execution failed (attempt {attempt + 1}/{max_retries + 1})")
                        # Continue to retry
                    else:
                        logger.error(f"❌ All {max_retries + 1} execution attempts failed")
                        # Exit retry loop with final failure
            else:
                # Non-dict result (unexpected) - treat as success
                break

        # If execution failed, propagate failure status so callers don't treat it as success.
        exec_status = "success"
        if isinstance(execution_result, dict):
            exec_status = execution_result.get("status", "success")
        overall_status = "success" if exec_status == "success" else "error"

        # Extract error information from execution_result for top-level error reporting
        error_message = None
        if overall_status == "error" and isinstance(execution_result, dict):
            error_message = execution_result.get("error")
            if not error_message and execution_result.get("stderr_full"):
                error_message = f"Execution failed: {execution_result.get('stderr_full', '')}"
            elif not error_message and execution_result.get("stderr"):
                error_message = f"Execution failed: {execution_result.get('stderr', '')}"
            elif not error_message:
                error_message = "Tool execution failed (see execution_result for details)"

        result = {
            "status": overall_status,
            "tool_generated": True,
            "explanation": _extract_explanation(content),
            "code_preview": code[:500] + "..." if len(code) > 500 else code,
            "execution_result": execution_result,
            "full_response": content
        }
        
        # Add error field if execution failed
        if error_message:
            result["error"] = error_message

        # Clean up temp FASTA files created from previous_plan_steps
        for inp in (inputs or []):
            if isinstance(inp, dict) and inp.get("source") == "previous_plan_step":
                uri = inp.get("uri")
                if uri and os.path.isfile(uri):
                    try:
                        os.unlink(uri)
                        logger.debug("Removed temp alignment file from previous_plan_steps: %s", uri)
                    except Exception:
                        pass

        return result
        
    except Exception as e:
        # Clean up temp FASTA from previous_plan_steps on error path too
        for inp in (inputs or []):
            if isinstance(inp, dict) and inp.get("source") == "previous_plan_step":
                uri = inp.get("uri")
                if uri and os.path.isfile(uri):
                    try:
                        os.unlink(uri)
                    except Exception:
                        pass
        logger.error(f"❌ Tool Generator Agent error: {e}", exc_info=True)
        return {
            "status": "error",
            "tool_generated": False,
            "error": str(e),
            "message": f"Failed to generate and execute tool: {str(e)}"
        }


def _extract_python_code(content: str) -> Optional[str]:
    """Extract Python code from markdown code blocks."""
    import re
    
    # Try to find Python code blocks
    patterns = [
        r'```python\n(.*?)```',
        r'```py\n(.*?)```',
        r'```\n(.*?)```',
        r'```python(.*?)```',
    ]
    
    for pattern in patterns:
        matches = re.findall(pattern, content, re.DOTALL)
        if matches:
            # Return the longest match (most likely the main code)
            return max(matches, key=len).strip()
    
    return None


def _extract_explanation(content: str) -> str:
    """Extract explanation/justification from the LLM response."""
    # Try to extract text before code blocks
    import re
    
    # Find first code block
    code_block_match = re.search(r'```', content)
    if code_block_match:
        explanation = content[:code_block_match.start()].strip()
        # Clean up common prefixes
        explanation = re.sub(r'^(I\'ll|I will|Here|This|The solution)', '', explanation, flags=re.IGNORECASE)
        return explanation.strip()
    
    # If no code block, return first paragraph
    lines = content.split('\n')
    explanation_lines = []
    for line in lines:
        if line.strip() and not line.strip().startswith('```'):
            explanation_lines.append(line)
            if len(explanation_lines) >= 3:  # First few lines
                break
    
    return '\n'.join(explanation_lines).strip()


async def _execute_generated_code(
    code: str,
    original_command: str,
    session_id: Optional[str] = None,
    total_file_size_bytes: int = 0,
    inputs: Optional[List[Any]] = None,
    outputs: Optional[List[Any]] = None
) -> Dict[str, Any]:
    """
    Execute the generated Python code in a safe environment.
    
    This function:
    1. Checks if EC2 execution is enabled and available
    2. If EC2 is available, executes on EC2 instance with pre-installed tools
    3. Otherwise, executes locally
    4. Captures output and errors
    5. Cleans up temporary files
    """
    # Check if we should skip EC2 for large files (should use EMR instead)
    emr_threshold = int(os.getenv('HELIX_ASYNC_BYTES_THRESHOLD', 100 * 1024 * 1024))
    if total_file_size_bytes > emr_threshold:
        logger.warning(f"⚠️  File size ({total_file_size_bytes / (1024*1024):.2f} MB) exceeds EMR threshold ({emr_threshold / (1024*1024):.2f} MB)")
        logger.warning(f"⚠️  This should be executed on EMR, not EC2. Skipping EC2 execution.")
        return {
            "status": "error",
            "error": f"Files are too large ({total_file_size_bytes / (1024*1024):.2f} MB) for EC2 execution. This operation should be routed to EMR, but EMR routing for tool-generator-agent is not yet implemented. Please use a pre-existing tool that supports EMR routing.",
            "stderr": f"File size {total_file_size_bytes / (1024*1024):.2f} MB exceeds EMR threshold {emr_threshold / (1024*1024):.2f} MB. EMR routing required but not implemented for tool-generator-agent.",
            "stdout": ""
        }
    
    # Build execution environment variables for generated code
    exec_env_vars: Dict[str, str] = {}
    input_uris: List[str] = []
    output_uris: List[str] = []

    if inputs:
        for inp in inputs:
            if hasattr(inp, "uri"):
                uri = inp.uri
            elif isinstance(inp, dict):
                uri = inp.get("uri")
            else:
                uri = None
            if uri:
                input_uris.append(str(uri))

    if outputs:
        for out in outputs:
            if hasattr(out, "uri"):
                uri = out.uri
            elif isinstance(out, dict):
                uri = out.get("uri")
            else:
                uri = None
            if uri:
                output_uris.append(str(uri))

    if input_uris:
        exec_env_vars["INPUT_FILES"] = ",".join(input_uris)
        exec_env_vars["INPUT_R1"] = input_uris[0]
        if len(input_uris) > 1:
            exec_env_vars["INPUT_R2"] = input_uris[1]

    if output_uris:
        exec_env_vars["OUTPUT_PATH"] = output_uris[0]
        exec_env_vars["OUTPUT_URI"] = output_uris[0]
        if output_uris[0].startswith("s3://"):
            exec_env_vars["OUTPUT_DIR"] = output_uris[0].rsplit("/", 1)[0]
        else:
            exec_env_vars["OUTPUT_DIR"] = str(Path(output_uris[0]).parent) if "/" in output_uris[0] else output_uris[0]

    exec_env_vars["HELIX_ORIGINAL_COMMAND"] = original_command or ""

    # Check if Docker sandbox execution is enabled (NEW!)
    # Docker sandbox has all bioinformatics tools pre-installed
    use_sandbox = os.getenv('HELIX_USE_SANDBOX', 'true').lower() == 'true'
    if use_sandbox and os.getenv("HELIX_MOCK_MODE") != "1":
        logger.info("🐳 Attempting to execute generated code in Docker sandbox (all bioinfo tools available)")
        try:
            # Get Docker image name
            biotools_image = os.getenv('HELIX_BIOTOOLS_IMAGE', 'helix-biotools:latest')

            # Check if Docker is available without throwing noisy tracebacks.
            docker_path = shutil.which("docker")
            if not docker_path:
                logger.warning("⚠️  Docker binary not found, falling back to local/EC2 execution")
            else:
                docker_check = subprocess.run(
                    [docker_path, "--version"],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                if docker_check.returncode != 0:
                    logger.warning("⚠️  Docker command unavailable, falling back to local/EC2 execution")
                    docker_path = None
                if not docker_path:
                    # Fall through to local/EC2 execution path below.
                    pass
                else:
                # Check if image exists
                    image_check = subprocess.run(
                        [docker_path, "images", "-q", biotools_image],
                        capture_output=True,
                        text=True,
                        timeout=5
                    )
                
                    if not image_check.stdout.strip():
                        logger.warning(f"⚠️  Docker image {biotools_image} not found, falling back to local/EC2 execution")
                        logger.warning(f"💡 Build it with: docker build -f backend/Dockerfile.biotools -t {biotools_image} .")
                    else:
                        # Save code to persistent location for debugging
                        debug_dir = PROJECT_ROOT / "sessions" / "generated_code"
                        debug_dir.mkdir(parents=True, exist_ok=True)
                        
                        # Create meaningful filename
                        import re
                        tool_name = re.sub(r'[^a-zA-Z0-9_-]', '_', original_command[:50])
                        timestamp = time.strftime("%Y%m%d_%H%M%S")
                        debug_filename = f"{tool_name}_{timestamp}.py"
                        debug_path = debug_dir / debug_filename
                        
                        # Save to both debug location and temp file for execution
                        debug_path.write_text(code)
                        
                        temp_file = tempfile.NamedTemporaryFile(
                            mode='w',
                            suffix='.py',
                            delete=False,
                            dir='/tmp'
                        )
                        temp_file.write(code)
                        temp_file.flush()
                        temp_file_path = temp_file.name
                        temp_file.close()
                        
                        logger.info(f"🐳 Executing in Docker sandbox: {biotools_image}")
                        logger.info(f"📝 Code saved to: {temp_file_path}")
                        logger.info(f"💾 Debug copy saved to: {debug_path}")
                        
                        try:
                            env_args: List[str] = []
                            for key, value in exec_env_vars.items():
                                env_args.extend(["-e", f"{key}={value}"])
                            # Forward AWS credentials/config to Docker for S3 access
                            aws_env_keys = [
                                "AWS_ACCESS_KEY_ID",
                                "AWS_SECRET_ACCESS_KEY",
                                "AWS_SESSION_TOKEN",
                                "AWS_DEFAULT_REGION",
                                "AWS_REGION",
                                "AWS_PROFILE"
                            ]
                            for key in aws_env_keys:
                                if os.getenv(key):
                                    env_args.extend(["-e", f"{key}={os.getenv(key)}"])

                            aws_credentials_dir = os.path.expanduser("~/.aws")
                            aws_volume_args: List[str] = []
                            if os.path.isdir(aws_credentials_dir):
                                aws_volume_args = ["-v", f"{aws_credentials_dir}:/root/.aws:ro"]

                            # Execute code in Docker sandbox
                            result = subprocess.run(
                                [
                                    docker_path, "run", "--rm",
                                    *env_args,
                                    *aws_volume_args,
                                    "-v", f"{temp_file_path}:/code/script.py:ro",
                                    "-w", "/sandbox/work",
                                    "--memory", "4g",
                                    "--cpus", "2",
                                    "--network", "host",  # For S3 access
                                    biotools_image,
                                    "python", "/code/script.py"
                                ],
                                capture_output=True,
                                text=True,
                                timeout=300
                            )

                            # Cleanup temp file
                            try:
                                os.unlink(temp_file_path)
                            except:
                                pass

                            logger.info(f"🐳 Docker execution completed with exit code: {result.returncode}")

                            if result.returncode == 0:
                                logger.info("✅ Code executed successfully in Docker sandbox")
                                return {
                                    "status": "success",
                                    "stdout": result.stdout,
                                    "stderr": result.stderr,
                                    "stderr_full": result.stderr,
                                    "exit_code": result.returncode
                                }

                            logger.warning(f"⚠️  Docker execution failed with exit code {result.returncode}")
                            if result.stderr:
                                # Log first 500 chars as warning
                                logger.warning(f"stderr (preview): {result.stderr[:500]}")
                                # Log full stderr as error for debugging
                                logger.error(f"stderr (FULL):\n{result.stderr}")
                            else:
                                logger.warning("stderr: N/A")

                            if result.stdout:
                                logger.info(f"stdout:\n{result.stdout}")

                            return {
                                "status": "error",
                                "error": f"Code execution failed with return code {result.returncode}",
                                "stdout": result.stdout,
                                "stderr": result.stderr,
                                "stderr_full": result.stderr,
                                "exit_code": result.returncode
                            }

                        except subprocess.TimeoutExpired:
                            logger.error("❌ Docker execution timed out (300s)")
                            try:
                                os.unlink(temp_file_path)
                            except:
                                pass
                            return {
                                "status": "error",
                                "error": "Execution timed out after 300 seconds",
                                "stdout": "",
                                "stderr": "Timeout"
                            }
                        except Exception as e:
                            logger.error(f"❌ Docker execution exception: {e}", exc_info=True)
                            try:
                                os.unlink(temp_file_path)
                            except:
                                pass
                            # Fall through to EC2/local execution
                        
        except FileNotFoundError as e:
            logger.warning(f"⚠️  Sandbox setup skipped: {e}. Falling back to local/EC2 execution")
        except Exception as e:
            logger.error(f"❌ Sandbox setup failed: {e}")
            # Fall through to EC2/local execution
    
    # Check if EC2 execution is enabled.
    # In unit tests we default to mock mode; never attempt real EC2/SSH there.
    if os.getenv("HELIX_MOCK_MODE") == "1":
        use_ec2 = False
    else:
        use_ec2 = os.getenv('HELIX_USE_EC2', 'false').lower() == 'true'
    
    if use_ec2:
        try:
            # Prefer package import. The backend is normally run as `python -m backend.main`,
            # so `backend.ec2_executor` is the correct module path.
            try:
                from backend.ec2_executor import get_ec2_executor
            except Exception:
                # Fallback for legacy execution where backend/ is on PYTHONPATH directly.
                from ec2_executor import get_ec2_executor
            executor = get_ec2_executor()
            
            logger.info("🔧 Using EC2 instance for code execution (with pre-installed bioinformatics tools)")
            
            # Get or create EC2 instance
            instance_id = executor.get_or_create_instance()
            if not instance_id:
                logger.warning("⚠️  Could not get EC2 instance, falling back to local execution")
                use_ec2 = False
            else:
                # Execute on EC2
                logger.info(f"🔧 Executing code on EC2 instance: {instance_id}")
                try:
                    result = executor.execute_code_on_instance(code, instance_id, timeout=300)
                    
                    if result.get("status") == "success":
                        logger.info("✅ Code executed successfully on EC2 instance")
                    else:
                        error_msg = result.get('error', 'Unknown error')
                        error_details = result.get('details', '')
                        stderr = result.get('stderr', '')
                        stdout = result.get('stdout', '')
                        
                        # Build comprehensive error message
                        error_parts = [error_msg]
                        if error_details:
                            error_parts.append(error_details)
                        if stderr:
                            error_parts.append(f"stderr: {stderr[:500]}")  # Limit length
                        
                        full_error = " - ".join(error_parts)
                        logger.warning(f"⚠️  EC2 execution failed: {full_error}")
                        
                        # Log full result for debugging (at debug level to avoid noise)
                        logger.debug(f"EC2 execution result: {result}")
                    
                    return result
                except Exception as e:
                    logger.error(f"❌ Exception during EC2 execution: {e}", exc_info=True)
                    return {
                        "status": "error",
                        "error": f"Exception during EC2 execution: {str(e)}",
                        "error_type": type(e).__name__
                    }
                
        except ImportError as e:
            # This can happen if paramiko is missing, but can also be triggered by *any* import-time
            # failure inside backend/ec2_executor.py (boto3, botocore, etc). Log the real reason.
            logger.warning(
                f"⚠️  EC2 executor not available (import failed: {type(e).__name__}: {e}). "
                "Falling back to local execution."
            )
            use_ec2 = False
        except Exception as e:
            logger.error(f"❌ EC2 execution error: {e}", exc_info=True)
            logger.warning("⚠️  Falling back to local execution")
            use_ec2 = False
    
    # Local execution (fallback or default)
    if not use_ec2:
        logger.info("🔧 Executing generated code locally")
    
    temp_file = None
    try:
        # Create temporary Python file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as f:
            f.write(code)
            temp_file = f.name
        
        logger.info(f"🔧 Executing generated code from: {temp_file}")
        
        # Set up environment
        env = os.environ.copy()
        # Add project paths
        project_root = str(PROJECT_ROOT)
        tools_path = str(PROJECT_ROOT / "tools")
        backend_path = str(PROJECT_ROOT / "backend")
        
        if 'PYTHONPATH' in env:
            env['PYTHONPATH'] = f"{project_root}:{tools_path}:{backend_path}:{env['PYTHONPATH']}"
        else:
            env['PYTHONPATH'] = f"{project_root}:{tools_path}:{backend_path}"
        
        # Add information about available tools to help generated code
        # Check for common bioinformatics tools
        import shutil
        available_tools = []
        common_tools = ['bbmerge.sh', 'bbmerge', 'samtools', 'bcftools', 'bwa', 'bowtie2', 'fastqc', 'Rscript']
        for tool in common_tools:
            if shutil.which(tool):
                available_tools.append(tool)
        
        # Set environment variable with available tools info (comma-separated)
        if available_tools:
            env['HELIX_AVAILABLE_TOOLS'] = ','.join(available_tools)
        else:
            env['HELIX_AVAILABLE_TOOLS'] = ''
        
        # Add execution context variables for generated tools
        if exec_env_vars:
            env.update(exec_env_vars)

        # Indicate that Python-based solutions should be preferred
        env['HELIX_PREFER_PYTHON'] = '1'

        # Create an isolated per-execution working directory so generated scripts
        # write their output files there instead of polluting the repo root.
        run_tmp = Path(tempfile.gettempdir()) / "helix-tool-runs" / uuid.uuid4().hex[:16]
        run_tmp.mkdir(parents=True, exist_ok=True)

        # Expose the run directory so generated scripts can reference it explicitly.
        env.setdefault("HELIX_OUTPUT_DIR", str(run_tmp))
        env.setdefault("OUTPUT_DIR", str(run_tmp))

        # Execute the code
        result = subprocess.run(
            [sys.executable, temp_file],
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout
            env=env,
            cwd=str(run_tmp)
        )
        
        stdout = result.stdout
        stderr = result.stderr
        
        if result.returncode == 0:
            return {
                "status": "success",
                "stdout": stdout,
                "stderr": stderr,
                "stderr_full": stderr,
                "returncode": result.returncode
            }
        else:
            return {
                "status": "error",
                "stdout": stdout,
                "stderr": stderr,
                "stderr_full": stderr,
                "returncode": result.returncode,
                "error": f"Code execution failed with return code {result.returncode}"
            }
            
    except subprocess.TimeoutExpired:
        return {
            "status": "error",
            "error": "Code execution timed out after 5 minutes"
        }
    except Exception as e:
        logger.error(f"Error executing generated code: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }
    finally:
        # Clean up temporary file
        if temp_file and os.path.exists(temp_file):
            try:
                os.unlink(temp_file)
            except Exception as e:
                logger.warning(f"Could not delete temporary file {temp_file}: {e}")


def generate_tool_sync(
    command: str,
    user_request: str,
    session_id: Optional[str] = None
) -> Dict[str, Any]:
    """
    Synchronous wrapper for generate_and_execute_tool.
    Useful for integration with existing synchronous code.
    """
    return asyncio.run(generate_and_execute_tool(command, user_request, session_id))

