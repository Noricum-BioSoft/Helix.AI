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
import asyncio
import logging
import tempfile
import subprocess
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
        openai_model = os.getenv("HELIX_TOOLGEN_OPENAI_MODEL", "openai:gpt-4o").strip()
        
        openai_enabled = openai_key and openai_key not in ["", "disabled", "your_openai_api_key_here", "none"]
        deepseek_enabled = deepseek_key and deepseek_key not in ["", "disabled", "your_deepseek_api_key_here", "none"]
        
        if openai_enabled:
            # Local import to avoid importing optional SSL/cert deps during test collection.
            from langchain.chat_models import init_chat_model
            # Be forgiving: allow HELIX_TOOLGEN_OPENAI_MODEL="gpt-4o" as well.
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
    outputs: Optional[List[Any]] = None
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
        
        # Build the prompt for tool generation
        user_prompt = f"""Generate and execute a solution for this bioinformatics operation:

User Command: {command}

Original Request: {user_request}
{input_metadata_section}{output_metadata_section}
Please follow the workflow:
1. Analyze the task and identify the biological operation
2. Research appropriate bioinformatics tools
3. Decide on infrastructure:
   - **EC2 instance with pre-installed tools** (if HELIX_USE_EC2 is enabled): Best for operations requiring established tools
   - **AWS EMR**: For large S3-hosted files (>100MB), distributed processing
   - **AWS Batch**: For containerized tools, medium-sized jobs
   - **Local execution**: Fallback when EC2/cloud not available
   
   **IMPORTANT**: Use the file size information provided above to make an informed infrastructure decision. If file sizes are unknown, assume files may be large and prefer EMR for S3-hosted files.
4. Generate complete, executable Python code
5. Execute the code and return results

Your response should include:
- Tool selection and justification
- Infrastructure choice and reasoning
- Complete Python implementation
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
        # Calculate total file size from inputs
        total_size = 0
        if inputs:
            for inp in inputs:
                if hasattr(inp, 'size_bytes') and isinstance(inp.size_bytes, int):
                    total_size += inp.size_bytes
                elif isinstance(inp, dict) and isinstance(inp.get('size_bytes'), int):
                    total_size += inp['size_bytes']
        
        # Get EMR threshold (default 100MB, matching execution_broker)
        emr_threshold = int(os.getenv('HELIX_ASYNC_BYTES_THRESHOLD', 100 * 1024 * 1024))
        
        # Route to EMR if files are large
        if total_size > emr_threshold:
            logger.warning(f"⚠️  Tool Generator Agent: Files are large ({total_size / (1024*1024):.2f} MB > {emr_threshold / (1024*1024):.2f} MB)")
            logger.warning(f"⚠️  EMR routing for tool-generator-agent is not yet implemented. This should be routed through execution_broker.")
            logger.warning(f"⚠️  For now, executing on EC2 (may be slow for large files).")
            logger.info("🔧 Tool Generator Agent: Executing generated code on EC2 (EMR routing pending implementation)...")
            # Pass total_size to _execute_generated_code so it can skip EC2 if needed
            execution_result = await _execute_generated_code(code, command, session_id, total_file_size_bytes=total_size)
        else:
            # Execute the generated code locally or on EC2
            logger.info("🔧 Tool Generator Agent: Executing generated code...")
            execution_result = await _execute_generated_code(code, command, session_id, total_file_size_bytes=total_size)

        # If execution failed, propagate failure status so callers don't treat it as success.
        exec_status = "success"
        if isinstance(execution_result, dict):
            exec_status = execution_result.get("status", "success")
        overall_status = "success" if exec_status == "success" else "error"

        # Extract error information from execution_result for top-level error reporting
        error_message = None
        if overall_status == "error" and isinstance(execution_result, dict):
            error_message = execution_result.get("error")
            if not error_message and execution_result.get("stderr"):
                error_message = f"Execution failed: {execution_result.get('stderr', '')[:500]}"
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

        return result
        
    except Exception as e:
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
    total_file_size_bytes: int = 0
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
    
    # Check if EC2 execution is enabled.
    # In unit tests we default to mock mode; never attempt real EC2/SSH there.
    if os.getenv("HELIX_MOCK_MODE") == "1":
        use_ec2 = False
    else:
        use_ec2 = os.getenv('HELIX_USE_EC2', 'false').lower() == 'true'
    
    if use_ec2:
        try:
            # Prefer package import. The backend is normally run as `python -m backend.main_with_mcp`,
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
        
        # Indicate that Python-based solutions should be preferred
        env['HELIX_PREFER_PYTHON'] = '1'
        
        # Execute the code
        result = subprocess.run(
            [sys.executable, temp_file],
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout
            env=env,
            cwd=str(PROJECT_ROOT)
        )
        
        stdout = result.stdout
        stderr = result.stderr
        
        if result.returncode == 0:
            return {
                "status": "success",
                "stdout": stdout,
                "stderr": stderr,
                "returncode": result.returncode
            }
        else:
            return {
                "status": "error",
                "stdout": stdout,
                "stderr": stderr,
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

