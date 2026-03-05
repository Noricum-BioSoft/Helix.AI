"""
Test that all tools with medium/large files route to async/EMR execution via infrastructure decision agent.
"""

import pytest
from unittest.mock import Mock, AsyncMock, patch, MagicMock
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

from backend.execution_broker import ExecutionBroker, ExecutionRequest, InputAsset


class TestInfrastructureBasedAsyncRouting:
    """Test that infrastructure decision agent routes medium/large files to async."""
    
    @pytest.mark.asyncio
    async def test_medium_files_route_to_async_via_infrastructure_agent(self):
        """Test that medium files (100MB-10GB) route to async when infrastructure agent recommends EMR."""
        broker = ExecutionBroker(tool_executor=AsyncMock())
        
        # Mock infrastructure decision agent to recommend EMR for medium files
        mock_infra_decision = MagicMock()
        mock_infra_decision.infrastructure = "EMR"
        mock_infra_decision.reasoning = "Files on S3 are medium-sized, EMR recommended"
        
        inputs = [
            InputAsset(uri="s3://bucket/file.fq", size_bytes=500 * 1024 * 1024, source="args"),  # 500MB (medium)
        ]
        
        with patch('backend.infrastructure_decision_agent.decide_infrastructure', return_value=mock_infra_decision):
            decision = await broker._evaluate_routing_policy(
                tool_name="any_tool",
                estimated_bytes=500 * 1024 * 1024,  # 500MB (medium size)
                unknown_inputs=0,
                inputs=inputs,
                command="Process medium files",
                session_context={}
            )
        
        # Should route to async because infrastructure agent recommended EMR
        assert decision.mode == "async"
        assert "emr" in decision.reason.lower()
        print(f"✅ PASS: Medium files route to async via infrastructure agent (reason: {decision.reason})")
    
    @pytest.mark.asyncio
    async def test_large_files_route_to_async_via_infrastructure_agent(self):
        """Test that large files (>10GB) route to async when infrastructure agent recommends EMR."""
        broker = ExecutionBroker(tool_executor=AsyncMock())
        
        # Mock infrastructure decision agent to recommend EMR for large files
        mock_infra_decision = MagicMock()
        mock_infra_decision.infrastructure = "EMR"
        mock_infra_decision.reasoning = "Files are very large, EMR recommended for distributed processing"
        
        inputs = [
            InputAsset(uri="s3://bucket/file.bam", size_bytes=15 * 1024 * 1024 * 1024, source="args"),  # 15GB (large)
        ]
        
        with patch('backend.infrastructure_decision_agent.decide_infrastructure', return_value=mock_infra_decision):
            decision = await broker._evaluate_routing_policy(
                tool_name="any_tool",
                estimated_bytes=15 * 1024 * 1024 * 1024,  # 15GB (large size)
                unknown_inputs=0,
                inputs=inputs,
                command="Process large files",
                session_context={}
            )
        
        # Should route to async because infrastructure agent recommended EMR
        assert decision.mode == "async"
        assert "emr" in decision.reason.lower()
        print(f"✅ PASS: Large files route to async via infrastructure agent (reason: {decision.reason})")
    
    @pytest.mark.asyncio
    async def test_small_files_route_to_sync_via_infrastructure_agent(self):
        """Test that small files (<100MB) route to sync when infrastructure agent recommends Local/EC2."""
        broker = ExecutionBroker(tool_executor=AsyncMock())
        
        # Mock infrastructure decision agent to recommend Local for small files
        mock_infra_decision = MagicMock()
        mock_infra_decision.infrastructure = "LOCAL"
        mock_infra_decision.reasoning = "Small files, local execution is fastest"
        
        inputs = [
            InputAsset(uri="/local/path/file.fq", size_bytes=50 * 1024 * 1024, source="args"),  # 50MB (small)
        ]
        
        with patch('backend.infrastructure_decision_agent.decide_infrastructure', return_value=mock_infra_decision):
            decision = await broker._evaluate_routing_policy(
                tool_name="any_tool",
                estimated_bytes=50 * 1024 * 1024,  # 50MB (small, below threshold)
                unknown_inputs=0,
                inputs=inputs,
                command="Process small files",
                session_context={}
            )
        
        # Should route to sync because infrastructure agent recommended Local
        assert decision.mode == "sync"
        assert "local" in decision.reason.lower()
        print(f"✅ PASS: Small files route to sync via infrastructure agent (reason: {decision.reason})")
    
    @pytest.mark.asyncio
    async def test_fallback_threshold_routing_for_medium_files(self):
        """Test that fallback threshold routing works when infrastructure agent is unavailable."""
        broker = ExecutionBroker(tool_executor=AsyncMock())
        
        inputs = [
            InputAsset(uri="s3://bucket/file.fq", size_bytes=200 * 1024 * 1024, source="args"),  # 200MB (medium)
        ]
        
        # Simulate infrastructure agent failure - should fall back to threshold
        with patch('backend.infrastructure_decision_agent.decide_infrastructure', side_effect=Exception("Agent unavailable")):
            decision = await broker._evaluate_routing_policy(
                tool_name="any_tool",
                estimated_bytes=200 * 1024 * 1024,  # 200MB (exceeds 100MB threshold)
                unknown_inputs=0,
                inputs=inputs,
                command="Process files",
                session_context={}
            )
        
        # Should route to async via fallback threshold logic
        assert decision.mode == "async"
        assert "threshold" in decision.reason.lower()
        print(f"✅ PASS: Fallback threshold routing works for medium files (reason: {decision.reason})")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

