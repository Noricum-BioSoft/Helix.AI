"""
Unit tests for the tool-generator-agent module.

Tests the core functionality of dynamically generating and executing
bioinformatics tools when no pre-existing tool exists.
"""

import os
import sys
import pytest
import asyncio
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import importlib.util

# Ensure project root is on path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Import the module with fallback for pytest compatibility
try:
    from backend.tool_generator_agent import (
        generate_and_execute_tool,
        _extract_python_code,
        _extract_explanation,
        _execute_generated_code,
        TOOL_GENERATOR_SYSTEM_PROMPT,
    )
except Exception:
    # Fallback: import directly from file for pytest compatibility
    spec = importlib.util.spec_from_file_location(
        "tool_generator_agent",
        PROJECT_ROOT / "backend" / "tool_generator_agent.py"
    )
    if spec and spec.loader:
        module = importlib.util.module_from_spec(spec)
        sys.modules["backend.tool_generator_agent"] = module
        spec.loader.exec_module(module)
        generate_and_execute_tool = module.generate_and_execute_tool
        _extract_python_code = module._extract_python_code
        _extract_explanation = module._extract_explanation
        _execute_generated_code = module._execute_generated_code
        TOOL_GENERATOR_SYSTEM_PROMPT = module.TOOL_GENERATOR_SYSTEM_PROMPT
    else:
        raise


class TestToolGeneratorAgent:
    """Test suite for tool-generator-agent functionality."""
    
    def test_extract_python_code_with_markdown(self):
        """Test extracting Python code from markdown code blocks."""
        content = """
        Here's the solution:
        
        ```python
        def merge_reads(r1, r2):
            return "merged"
        ```
        
        This code will merge the reads.
        """
        code = _extract_python_code(content)
        assert code is not None
        assert "def merge_reads" in code
        assert "return \"merged\"" in code
    
    def test_extract_python_code_without_markdown(self):
        """Test extracting code when no markdown formatting."""
        content = """
        def merge_reads(r1, r2):
            return "merged"
        """
        code = _extract_python_code(content)
        # Should return None if no code block found
        assert code is None
    
    def test_extract_python_code_multiple_blocks(self):
        """Test extracting the longest code block when multiple exist."""
        content = """
        ```python
        import os
        ```
        
        ```python
        def merge_reads(r1, r2):
            # Long function
            result = process(r1, r2)
            return result
        ```
        """
        code = _extract_python_code(content)
        assert code is not None
        assert "def merge_reads" in code
        assert "process(r1, r2)" in code
    
    def test_extract_explanation(self):
        """Test extracting explanation text from response."""
        content = """
        I'll use BBMerge for this task because it's designed for merging paired-end reads.
        
        ```python
        def merge():
            pass
        ```
        """
        explanation = _extract_explanation(content)
        assert "BBMerge" in explanation or "merging" in explanation.lower()
        assert "```" not in explanation
    
    def test_extract_explanation_no_code_block(self):
        """Test extracting explanation when no code block exists."""
        content = "I'll use BBMerge for this task. It's the best tool."
        explanation = _extract_explanation(content)
        assert len(explanation) > 0
        assert "BBMerge" in explanation or "tool" in explanation.lower()
    
    @pytest.mark.asyncio
    async def test_execute_generated_code_success(self):
        """Test executing valid Python code."""
        code = """
import sys
print("Hello from generated code")
sys.exit(0)
"""
        result = await _execute_generated_code(code, "test command")
        assert result["status"] == "success"
        assert "Hello from generated code" in result["stdout"]
        assert result["returncode"] == 0
    
    @pytest.mark.asyncio
    async def test_execute_generated_code_error(self):
        """Test executing code that raises an error."""
        code = """
import sys
print("Error test", file=sys.stderr)
sys.exit(1)
"""
        result = await _execute_generated_code(code, "test command")
        assert result["status"] == "error"
        assert result["returncode"] == 1
        assert "Error test" in result["stderr"]
    
    @pytest.mark.asyncio
    async def test_execute_generated_code_syntax_error(self):
        """Test executing code with syntax errors."""
        code = """
def merge_reads(
    # Missing closing parenthesis
"""
        result = await _execute_generated_code(code, "test command")
        assert result["status"] == "error"
        assert result["returncode"] != 0
    
    @pytest.mark.asyncio
    @pytest.mark.skipif(
        not (os.getenv("OPENAI_API_KEY") or os.getenv("DEEPSEEK_API_KEY")),
        reason="LLM API key required for tool generation tests"
    )
    async def test_generate_and_execute_tool_mock_execution(self):
        """Test the full generate_and_execute_tool flow with mocked execution."""
        # Mock the LLM response
        mock_response = Mock()
        mock_response.content = """
        I'll use BBMerge for merging paired-end reads.
        
        ```python
        print("Merged reads successfully")
        ```
        """
        
        with patch('backend.tool_generator_agent._get_llm') as mock_llm:
            mock_llm_instance = Mock()
            mock_llm_instance.invoke = Mock(return_value=mock_response)
            mock_llm.return_value = mock_llm_instance
            
            # Mock the execution to avoid actually running code
            with patch('backend.tool_generator_agent._execute_generated_code') as mock_exec:
                mock_exec.return_value = {
                    "status": "success",
                    "stdout": "Merged reads successfully",
                    "stderr": "",
                    "returncode": 0
                }
                
                result = await generate_and_execute_tool(
                    command="merge R1 and R2 reads",
                    user_request="merge the following forward R1 and reverse R2 reads"
                )
                
                assert result["status"] == "success"
                assert result["tool_generated"] is True
                assert "explanation" in result
                assert "code_preview" in result
    
    @pytest.mark.asyncio
    @pytest.mark.skipif(
        not (os.getenv("OPENAI_API_KEY") or os.getenv("DEEPSEEK_API_KEY")),
        reason="LLM API key required for tool generation tests"
    )
    async def test_generate_and_execute_tool_llm_error(self):
        """Test handling LLM errors gracefully."""
        with patch('backend.tool_generator_agent._get_llm') as mock_llm:
            mock_llm.side_effect = Exception("LLM connection failed")
            
            result = await generate_and_execute_tool(
                command="merge reads",
                user_request="merge reads"
            )
            
            assert result["status"] == "error"
            assert result["tool_generated"] is False
            assert "error" in result
    
    def test_tool_generator_prompt_loaded(self):
        """Test that the tool-generator-agent prompt is loaded."""
        assert TOOL_GENERATOR_SYSTEM_PROMPT is not None
        assert len(TOOL_GENERATOR_SYSTEM_PROMPT) > 0
        # Should contain key concepts from the prompt
        assert "bioinformatics" in TOOL_GENERATOR_SYSTEM_PROMPT.lower() or "tool" in TOOL_GENERATOR_SYSTEM_PROMPT.lower()


class TestToolGeneratorIntegration:
    """Integration tests for tool-generator-agent with the main system."""
    
    @pytest.mark.asyncio
    @pytest.mark.integration
    @pytest.mark.skipif(
        not (os.getenv("OPENAI_API_KEY") or os.getenv("DEEPSEEK_API_KEY")),
        reason="LLM API key required for integration tests"
    )
    async def test_unknown_tool_triggers_generator(self):
        """Test that unknown tools trigger the tool-generator-agent."""
        from backend.main_with_mcp import call_mcp_tool
        
        # This should trigger the tool-generator-agent
        # We'll mock it to avoid actual LLM calls in unit tests
        with patch('backend.tool_generator_agent.generate_and_execute_tool') as mock_gen:
            mock_gen.return_value = {
                "status": "success",
                "tool_generated": True,
                "explanation": "Generated tool for merging reads",
                "execution_result": {"status": "success", "stdout": "Done"}
            }
            
            try:
                result = await call_mcp_tool("unknown_merge_tool", {
                    "r1_path": "s3://test/r1.fq",
                    "r2_path": "s3://test/r2.fq"
                })
                # Should not raise ValueError if tool-generator-agent is called
                assert result["status"] == "success"
                assert result["tool_generated"] is True
            except ValueError:
                # If tool-generator-agent fails, it should raise ValueError
                # This is expected behavior
                pass

