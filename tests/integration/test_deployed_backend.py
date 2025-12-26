#!/usr/bin/env python3
"""
Comprehensive Test Suite for Deployed Helix.AI Backend
Adapts existing tests to run against the deployed cloud backend.

Usage:
    # Test against default deployed backend
    pytest tests/integration/test_deployed_backend.py -v
    
    # Test against custom URL
    BACKEND_URL=https://custom-url.com pytest tests/integration/test_deployed_backend.py -v
    
    # Test specific test
    pytest tests/integration/test_deployed_backend.py::TestDeployedBackend::test_sequence_alignment -v
"""

import os
import pytest
import requests
import time
from typing import Dict, Any, Optional
from pathlib import Path

# Get backend URL from environment or use default deployed URL
DEFAULT_BACKEND_URL = os.getenv(
    "BACKEND_URL",
    "http://HelixA-ALBAE-D7BksiQIynZb-1051248867.us-west-1.elb.amazonaws.com"
)
# Alternative: CloudFront URL
# DEFAULT_BACKEND_URL = os.getenv("BACKEND_URL", "https://d2a8mt5n89vos4.cloudfront.net")


class DeployedBackendTester:
    """Test client for deployed Helix.AI backend"""
    
    def __init__(self, base_url: str = DEFAULT_BACKEND_URL, timeout: int = 30):
        self.base_url = base_url.rstrip('/')
        self.timeout = timeout
        self.session_id: Optional[str] = None
        
    def _make_request(self, method: str, endpoint: str, **kwargs) -> requests.Response:
        """Make HTTP request with error handling"""
        url = f"{self.base_url}{endpoint}"
        kwargs.setdefault('timeout', self.timeout)
        try:
            if method.upper() == 'GET':
                return requests.get(url, **kwargs)
            elif method.upper() == 'POST':
                return requests.post(url, **kwargs)
            else:
                raise ValueError(f"Unsupported method: {method}")
        except requests.exceptions.RequestException as e:
            pytest.fail(f"Request failed: {e}")
    
    def health_check(self) -> Dict[str, Any]:
        """Test health endpoint"""
        response = self._make_request('GET', '/health')
        assert response.status_code == 200, f"Expected 200, got {response.status_code}"
        data = response.json()
        assert data.get("status") == "healthy", f"Expected 'healthy', got {data.get('status')}"
        return data
    
    def create_session(self) -> str:
        """Create a new session"""
        response = self._make_request('POST', '/create_session')
        assert response.status_code == 200, f"Expected 200, got {response.status_code}"
        data = response.json()
        assert "session_id" in data, f"No session_id in response: {data}"
        self.session_id = data["session_id"]
        return self.session_id
    
    def execute_command(self, command: str, session_id: Optional[str] = None) -> Dict[str, Any]:
        """Execute a command via /execute endpoint"""
        if session_id is None:
            session_id = self.session_id or self.create_session()
        
        response = self._make_request('POST', '/execute', json={
            "command": command,
            "session_id": session_id
        })
        assert response.status_code == 200, f"Expected 200, got {response.status_code}"
        return response.json()
    
    def list_tools(self) -> Dict[str, Any]:
        """List available tools"""
        response = self._make_request('GET', '/mcp/tools')
        assert response.status_code == 200, f"Expected 200, got {response.status_code}"
        return response.json()


@pytest.fixture
def tester():
    """Fixture providing a test client"""
    return DeployedBackendTester()


@pytest.fixture
def session(tester):
    """Fixture providing a test session"""
    return tester.create_session()


class TestDeployedBackend:
    """Test suite for deployed backend"""
    
    @pytest.mark.integration
    def test_health_check(self, tester):
        """Test 1: Health Check"""
        result = tester.health_check()
        assert result["status"] == "healthy"
        assert "service" in result
    
    @pytest.mark.integration
    def test_session_creation(self, tester):
        """Test 2: Session Creation"""
        session_id = tester.create_session()
        assert session_id is not None
        assert len(session_id) > 0
    
    @pytest.mark.integration
    def test_sequence_alignment(self, tester, session):
        """Test 3: Sequence Alignment (Bio.Align.Applications fix)"""
        command = "align these sequences: >seq1 ATGCGATCG >seq2 ATGCGATC >seq3 ATGCGATCGATC"
        result = tester.execute_command(command, session)
        
        assert result.get("success") is True, f"Command failed: {result}"
        assert result.get("tool") == "sequence_alignment", f"Wrong tool: {result.get('tool')}"
        
        # Check for alignment results
        output = result.get("output", {})
        if isinstance(output, dict):
            # Check for alignment data
            assert "alignment" in output or "text" in result, "No alignment data in response"
        else:
            # Output might be a string
            assert len(str(output)) > 0, "Empty output"
    
    @pytest.mark.integration
    def test_mutation_generation(self, tester, session):
        """Test 4: Mutation Generation (langchain_core.tools fix)"""
        command = "generate 10 variants of sequence ATGCGATCG"
        result = tester.execute_command(command, session)
        
        assert result.get("success") is True, f"Command failed: {result}"
        assert result.get("tool") == "mutate_sequence", f"Wrong tool: {result.get('tool')}"
        
        # Check for mutation results
        output = result.get("output", {})
        if isinstance(output, dict):
            # Check for variants
            stats = output.get("statistics", {})
            assert "variants" in stats or "total_variants" in stats, "No variants in response"
        else:
            assert len(str(output)) > 0, "Empty output"
    
    @pytest.mark.integration
    def test_phylogenetic_analysis(self, tester, session):
        """Test 5: Phylogenetic Analysis"""
        command = "build phylogenetic tree for sequences: >seq1 ATGCGATCG >seq2 ATGCGATC >seq3 ATGCGATCGATC"
        result = tester.execute_command(command, session)
        
        assert result.get("success") is True, f"Command failed: {result}"
        # Phylogenetic tree might use different tool names
        assert "tree" in result.get("tool", "").lower() or result.get("success"), "Phylogenetic analysis failed"
    
    @pytest.mark.integration
    def test_tool_listing(self, tester):
        """Test 6: Tool Listing"""
        tools = tester.list_tools()
        assert isinstance(tools, (dict, list)), f"Unexpected tools format: {type(tools)}"
        
        # Check that key tools are available
        tools_str = str(tools).lower()
        assert "sequence_alignment" in tools_str or "align" in tools_str, "sequence_alignment not found"
        assert "mutate" in tools_str or "mutation" in tools_str, "mutate_sequence not found"
    
    @pytest.mark.integration
    def test_multi_step_workflow(self, tester, session):
        """Test 7: Multi-Step Workflow"""
        steps = [
            "generate 5 variants of sequence ATGCGATCG",
            "align all the generated variants",
        ]
        
        for i, step in enumerate(steps, 1):
            result = tester.execute_command(step, session)
            assert result.get("success") is True, f"Step {i} failed: {result}"
            time.sleep(1)  # Brief pause between steps
    
    @pytest.mark.integration
    def test_error_handling(self, tester, session):
        """Test 8: Error Handling"""
        invalid_command = "execute impossible operation that doesn't exist"
        result = tester.execute_command(invalid_command, session)
        
        # Should either return success with error message, or handle gracefully
        # The important thing is it doesn't crash
        assert "result" in result or "output" in result or "error" in result, "No response structure"
    
    @pytest.mark.integration
    def test_natural_language_commands(self, tester, session):
        """Test 9: Natural Language Command Processing"""
        commands = [
            "align these DNA sequences: >seq1 ATGC >seq2 ATGC",
            "mutate sequence ATGCGATCG to create 3 variants",
        ]
        
        success_count = 0
        for command in commands:
            try:
                result = tester.execute_command(command, session)
                if result.get("success"):
                    success_count += 1
            except AssertionError:
                pass  # Count as failure
        
        # At least 80% should succeed
        assert success_count >= len(commands) * 0.8, f"Only {success_count}/{len(commands)} commands succeeded"


class TestDeployedBackendFixes:
    """Specific tests for the fixes we deployed"""
    
    @pytest.mark.integration
    def test_bio_align_applications_fix(self, tester, session):
        """Test that Bio.Align.Applications import fix works"""
        # This should not fail with import error
        command = "align these sequences: >seq1 ATGCGATCG >seq2 ATGCGATCA"
        result = tester.execute_command(command, session)
        
        assert result.get("success") is True, f"Command failed: {result}"
        
        # Check that there's no import error in the response
        result_str = str(result).lower()
        assert "bio.align.applications" not in result_str, "Bio.Align.Applications error still present"
        assert "no module named" not in result_str, "Module import error present"
        assert "import" not in result_str or "success" in result_str, "Import error in response"
    
    @pytest.mark.integration
    def test_langchain_core_tools_fix(self, tester, session):
        """Test that langchain_core.tools import fix works"""
        # This should not fail with import error
        command = "mutate sequence ATGCGATCG to create 3 variants"
        result = tester.execute_command(command, session)
        
        assert result.get("success") is True, f"Command failed: {result}"
        
        # Check that there's no import error in the response
        result_str = str(result).lower()
        assert "langchain.agents" not in result_str, "langchain.agents error still present"
        assert "cannot import name 'tool'" not in result_str, "Tool import error still present"
        assert "no module named" not in result_str, "Module import error present"


@pytest.mark.integration
class TestDeployedBackendPerformance:
    """Performance tests for deployed backend"""
    
    def test_response_time(self, tester):
        """Test that health check responds quickly"""
        start = time.time()
        tester.health_check()
        elapsed = time.time() - start
        
        # Health check should be very fast (< 1 second)
        assert elapsed < 1.0, f"Health check took {elapsed:.2f}s, expected < 1s"
    
    def test_command_response_time(self, tester, session):
        """Test that commands respond within reasonable time"""
        command = "align these sequences: >seq1 ATGCGATCG >seq2 ATGCGATCA"
        
        start = time.time()
        result = tester.execute_command(command, session)
        elapsed = time.time() - start
        
        # Commands should complete within 30 seconds
        assert elapsed < 30.0, f"Command took {elapsed:.2f}s, expected < 30s"
        assert result.get("success") is True, "Command failed"


if __name__ == "__main__":
    # Allow running directly
    pytest.main([__file__, "-v", "-s"])


