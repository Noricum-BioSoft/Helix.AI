"""
Test that stdout/stderr logs from EC2 execution are properly captured and included in responses.
"""

import pytest
from backend.main_with_mcp import _extract_execution_logs, build_standard_response


def test_extract_logs_from_direct_ec2_result():
    """Test extracting logs from direct EC2 execution result."""
    result = {
        "status": "success",
        "stdout": "Processing 100 reads...\nMerging complete!",
        "stderr": "Warning: Quality score below threshold",
        "returncode": 0
    }
    
    logs = _extract_execution_logs(result)
    
    assert len(logs) == 2
    assert logs[0]["type"] == "stdout"
    assert "Processing 100 reads" in logs[0]["content"]
    assert logs[1]["type"] == "stderr"
    assert "Warning" in logs[1]["content"]


def test_extract_logs_from_tool_generator_result():
    """Test extracting logs from tool-generator-agent result."""
    result = {
        "status": "success",
        "tool_generated": True,
        "execution_result": {
            "status": "success",
            "stdout": "BBMerge installed successfully\nMerging reads...",
            "stderr": "",
            "returncode": 0
        }
    }
    
    logs = _extract_execution_logs(result)
    
    assert len(logs) == 1  # Only stdout (stderr is empty)
    assert logs[0]["type"] == "stdout"
    assert "BBMerge installed" in logs[0]["content"]


def test_extract_logs_from_nested_result():
    """Test extracting logs from broker-wrapped result."""
    result = {
        "status": "success",
        "result": {
            "execution_result": {
                "status": "success",
                "stdout": "Read merging output",
                "stderr": "Some warnings here",
                "returncode": 0
            }
        }
    }
    
    logs = _extract_execution_logs(result)
    
    assert len(logs) == 2
    assert any(log["type"] == "stdout" and "Read merging" in log["content"] for log in logs)
    assert any(log["type"] == "stderr" and "warnings" in log["content"] for log in logs)


def test_extract_logs_handles_empty_strings():
    """Test that empty strings are not added as logs."""
    result = {
        "status": "success",
        "stdout": "Some output",
        "stderr": "",  # Empty stderr should not be added
        "returncode": 0
    }
    
    logs = _extract_execution_logs(result)
    
    assert len(logs) == 1
    assert logs[0]["type"] == "stdout"


def test_extract_logs_handles_no_logs():
    """Test that results without logs return empty list."""
    result = {
        "status": "success",
        "text": "Operation completed",
        "summary": {"total": 100}
    }
    
    logs = _extract_execution_logs(result)
    
    assert logs == []


def test_build_standard_response_includes_logs():
    """Test that build_standard_response includes logs in the response."""
    result = {
        "status": "success",
        "text": "Read merging completed",
        "execution_result": {
            "stdout": "Processed 1000 reads",
            "stderr": "Warning: low quality"
        }
    }
    
    response = build_standard_response(
        prompt="Merge reads",
        tool="read_merging",
        result=result,
        session_id="test-session",
        mcp_route="/execute",
        success=True
    )
    
    assert "logs" in response
    assert len(response["logs"]) == 2
    assert any(log["type"] == "stdout" for log in response["logs"])
    assert any(log["type"] == "stderr" for log in response["logs"])


def test_build_standard_response_logs_have_timestamps():
    """Test that logs include timestamps."""
    result = {
        "status": "success",
        "stdout": "Test output"
    }
    
    response = build_standard_response(
        prompt="Test",
        tool="test_tool",
        result=result,
        session_id="test-session",
        mcp_route="/execute",
        success=True
    )
    
    assert len(response["logs"]) == 1
    assert "timestamp" in response["logs"][0]
    assert response["logs"][0]["timestamp"]  # Should not be empty


if __name__ == "__main__":
    pytest.main([__file__, "-v"])




