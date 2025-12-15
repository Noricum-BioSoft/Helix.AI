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
from typing import Dict, Any, Optional
import json

logger = logging.getLogger(__name__)

# Load the tool-generator-agent prompt
PROJECT_ROOT = Path(__file__).resolve().parent.parent
TOOL_GENERATOR_PROMPT_PATH = PROJECT_ROOT / "agents" / "tool-generator-agent.md"

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
        
        openai_enabled = openai_key and openai_key not in ["", "disabled", "your_openai_api_key_here", "none"]
        deepseek_enabled = deepseek_key and deepseek_key not in ["", "disabled", "your_deepseek_api_key_here", "none"]
        
        if openai_enabled:
            # Local import to avoid importing optional SSL/cert deps during test collection.
            from langchain.chat_models import init_chat_model
            return init_chat_model("openai:gpt-4o", temperature=0)
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
    session_id: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate and execute a bioinformatics tool for the given command.
    
    Args:
        command: The user's bioinformatics command/request
        user_request: The original user request (for context)
        session_id: Optional session ID for tracking
        
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
        
        # Build the prompt for tool generation
        user_prompt = f"""Generate and execute a solution for this bioinformatics operation:

User Command: {command}

Original Request: {user_request}

Please follow the workflow:
1. Analyze the task and identify the biological operation
2. Research appropriate bioinformatics tools
3. Decide on infrastructure:
   - **EC2 instance with pre-installed tools** (if HELIX_USE_EC2 is enabled): Best for operations requiring established tools
   - **AWS EMR**: For large S3-hosted files (>100MB), distributed processing
   - **AWS Batch**: For containerized tools, medium-sized jobs
   - **Local execution**: Fallback when EC2/cloud not available
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
4. **For read merging**: If BBMerge/FLASH/PEAR are not available and cannot be installed, use simple Python string manipulation for overlap detection - NO BioPython SeqRecord manipulation needed:
   ```python
   def reverse_complement(seq: str) -> str:
       complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
       return seq.translate(complement)[::-1]
   
   def merge_pair(forward: str, reverse: str, min_overlap: int) -> str:
       rc_reverse = reverse_complement(reverse)
       for overlap in range(min(len(forward), len(rc_reverse)), min_overlap - 1, -1):
           if forward.endswith(rc_reverse[:overlap]):
               return forward + rc_reverse[overlap:]
       return forward + rc_reverse
   ```
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

**AWS Credentials on EC2:**
- When executing on EC2 instances, AWS credentials are automatically available via IAM instance role
- Simply use `boto3.client('s3')` without explicit credentials - boto3 will automatically use IAM role credentials
- The AWS region is set via environment variables (AWS_REGION, AWS_DEFAULT_REGION)
- Example:
  ```python
  import boto3
  import os
  
  # Get region from environment (set automatically on EC2)
  region = os.getenv('AWS_REGION', 'us-east-1')
  s3 = boto3.client('s3', region_name=region)
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
        
        # Execute the generated code
        logger.info("🔧 Tool Generator Agent: Executing generated code...")
        execution_result = await _execute_generated_code(code, command, session_id)

        # If execution failed, propagate failure status so callers don't treat it as success.
        exec_status = "success"
        if isinstance(execution_result, dict):
            exec_status = execution_result.get("status", "success")
        overall_status = "success" if exec_status == "success" else "error"

        return {
            "status": overall_status,
            "tool_generated": True,
            "explanation": _extract_explanation(content),
            "code_preview": code[:500] + "..." if len(code) > 500 else code,
            "execution_result": execution_result,
            "full_response": content
        }
        
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
    session_id: Optional[str] = None
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
    # Check if EC2 execution is enabled
    use_ec2 = os.getenv('HELIX_USE_EC2', 'false').lower() == 'true'
    
    if use_ec2:
        try:
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
                
        except ImportError:
            logger.warning("⚠️  EC2 executor not available (paramiko not installed?), falling back to local execution")
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

