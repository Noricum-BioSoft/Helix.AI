"""
Proof that read_merging refactoring works correctly.

This test demonstrates that:
1. Agent tool mapping correctly identifies read_merging
2. Execution directly calls merge_reads_from_s3 (no tool-generator-agent)
3. The flow is simplified and efficient
"""

import pytest
from unittest.mock import Mock, AsyncMock, patch, MagicMock
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

from backend.main_with_mcp import call_mcp_tool
from backend.agent_tools import read_merging
from backend.execution_broker import ExecutionBroker, ExecutionRequest


class TestReadMergingRefactor:
    """Test that read_merging refactoring eliminated redundancy."""
    
    @pytest.mark.asyncio
    async def test_read_merging_s3_paths_direct_execution(self):
        """Test that S3 paths directly call merge_reads_from_s3, not tool-generator-agent."""
        arguments = {
            "forward_reads": "s3://bucket/test_R1.fq",
            "reverse_reads": "s3://bucket/test_R2.fq",
            "output": "s3://bucket/merged.fq",
            "min_overlap": 12,
            "_from_broker": True  # Simulate call from execution broker
        }
        
        # Mock the read_merging module import inside call_mcp_tool
        with patch('read_merging.merge_reads_from_s3') as mock_merge_from_s3:
            mock_merge_from_s3.return_value = {
                "status": "success",
                "text": "Read merging completed successfully.",
                "summary": {"total_pairs": 100},
                "output_path": "s3://bucket/merged.fq"
            }
            
            # Ensure tool-generator-agent is NOT called
            with patch('backend.tool_generator_agent.generate_and_execute_tool') as mock_tool_gen:
                result = await call_mcp_tool("read_merging", arguments)
                
                # Verify merge_reads_from_s3 was called directly
                mock_merge_from_s3.assert_called_once_with(
                    r1_path="s3://bucket/test_R1.fq",
                    r2_path="s3://bucket/test_R2.fq",
                    output_path="s3://bucket/merged.fq",
                    min_overlap=12
                )
                
                # Verify tool-generator-agent was NOT called
                mock_tool_gen.assert_not_called()
                
                # Verify result structure
                assert result["status"] == "success"
                assert "summary" in result
    
    @pytest.mark.asyncio
    async def test_read_merging_fastq_content_direct_execution(self):
        """Test that FASTQ content directly calls run_read_merging_raw."""
        arguments = {
            "forward_reads": "@read1\nATCG\n+\nIIII",
            "reverse_reads": "@read1\nCGAT\n+\nIIII",
            "min_overlap": 12,
            "_from_broker": True
        }
        
        with patch('read_merging.run_read_merging_raw') as mock_run_raw:
            mock_run_raw.return_value = {
                "text": "Read merging completed successfully.",
                "merged_sequences": ">merged_1\nATCG",
                "summary": {"total_pairs": 1}
            }
            
            with patch('read_merging.merge_reads_from_s3') as mock_merge_from_s3:
                result = await call_mcp_tool("read_merging", arguments)
                
                # Verify run_read_merging_raw was called directly
                mock_run_raw.assert_called_once_with(
                    "@read1\nATCG\n+\nIIII",
                    "@read1\nCGAT\n+\nIIII",
                    12
                )
                
                # Verify merge_reads_from_s3 was NOT called
                mock_merge_from_s3.assert_not_called()
    
    @pytest.mark.asyncio
    async def test_read_merging_output_path_inference(self):
        """Test that output path is inferred from R1 path if not provided."""
        arguments = {
            "forward_reads": "s3://bucket/data_R1.fq",
            "reverse_reads": "s3://bucket/data_R2.fq",
            "min_overlap": 12,
            "_from_broker": True
        }
        
        with patch('read_merging.merge_reads_from_s3') as mock_merge_from_s3:
            mock_merge_from_s3.return_value = {
                "status": "success",
                "text": "Read merging completed successfully.",
                "summary": {"total_pairs": 100},
                "output_path": "s3://bucket/data_R1_merged.fq"
            }
            
            result = await call_mcp_tool("read_merging", arguments)
            
            # Verify merge_reads_from_s3 was called with inferred output path
            call_args = mock_merge_from_s3.call_args
            assert call_args[1]["output_path"] == "s3://bucket/data_R1_merged.fq"
    
    def test_read_merging_tool_wrapper_documentation(self):
        """Test that the agent tool wrapper exists for mapping purposes."""
        # The tool wrapper should exist and have correct signature
        # Note: read_merging is a @tool decorated function, so we need to check the underlying function
        import inspect
        from langchain_core.tools import tool
        
        # Get the actual function from the tool wrapper
        if hasattr(read_merging, 'func'):
            actual_func = read_merging.func
        elif hasattr(read_merging, '__wrapped__'):
            actual_func = read_merging.__wrapped__
        else:
            actual_func = read_merging
        
        sig = inspect.signature(actual_func)
        params = list(sig.parameters.keys())
        
        assert "forward_reads" in params
        assert "reverse_reads" in params
        assert "min_overlap" in params
        assert "output" in params
        
        # Verify it has documentation
        assert actual_func.__doc__ is not None
        assert "Merge paired-end reads" in actual_func.__doc__
        assert "S3 path" in actual_func.__doc__
    
    def test_no_circular_routing(self):
        """Test that there's no circular routing back to execution broker."""
        # This test verifies the refactored code structure
        import ast
        import inspect
        
        # Read the call_mcp_tool function source
        source = inspect.getsource(call_mcp_tool)
        
        # Parse and analyze the AST
        tree = ast.parse(source)
        
        # Check that for read_merging, we don't call broker.execute_tool recursively
        # We should find direct calls to merge_reads_from_s3 or run_read_merging_raw
        
        class BrokerCallChecker(ast.NodeVisitor):
            def __init__(self):
                self.has_broker_call_in_read_merging = False
                self.in_read_merging_branch = False
                
            def visit_If(self, node):
                # Check if we're in the read_merging elif branch
                if isinstance(node.test, ast.Compare):
                    if any(isinstance(comp, ast.Eq) for comp in node.test.ops):
                        if any(isinstance(left, ast.Name) and left.id == "tool_name" 
                               for left in ast.walk(node.test.left)):
                            if any(isinstance(s, ast.Str) and s.s == "read_merging" 
                                   for s in ast.walk(node.test)):
                                self.in_read_merging_branch = True
                                self.generic_visit(node)
                                self.in_read_merging_branch = False
                                return
                
                self.generic_visit(node)
            
            def visit_Call(self, node):
                if self.in_read_merging_branch:
                    # Check for broker.execute_tool calls
                    if isinstance(node.func, ast.Attribute):
                        if node.func.attr == "execute_tool":
                            self.has_broker_call_in_read_merging = True
                self.generic_visit(node)
        
        checker = BrokerCallChecker()
        checker.visit(tree)
        
        # After refactoring, we should NOT have broker.execute_tool calls in read_merging branch
        assert not checker.has_broker_call_in_read_merging, \
            "Found circular routing: read_merging still calls broker.execute_tool"


if __name__ == "__main__":
    # Run the tests
    pytest.main([__file__, "-v"])

