"""
Tests for Phase 2 Read-Only Tools:
- FileMetadataInspector
- EnvironmentCapabilityCatalog
- CostHeuristicTable

These tools provide factual grounding to the Infrastructure Decision Agent.
"""

import pytest
import json
import tempfile
from pathlib import Path

from backend.tools.file_metadata import FileMetadataInspector, FileMetadata
from backend.tools.env_catalog import EnvironmentCapabilityCatalog, EnvironmentCapability
from backend.tools.cost_heuristics import CostHeuristicTable, CostHeuristic


class TestFileMetadataInspector:
    """Test FileMetadataInspector tool."""
    
    def test_mock_mode_with_catalog(self):
        """Inspector in mock mode should use catalog."""
        catalog = {
            "s3://bucket/file1.fq": 1000000,
            "/local/file2.fq": 500000
        }
        
        inspector = FileMetadataInspector(mock_mode=True, mock_catalog=catalog)
        
        # S3 file from catalog
        metadata = inspector.inspect_file("s3://bucket/file1.fq")
        assert metadata.size_bytes == 1000000
        assert metadata.size_confidence == 1.0
        assert metadata.source == "mock_catalog"
        assert metadata.location_type == "S3"
        assert metadata.accessible is True
        
        # Local file from catalog
        metadata = inspector.inspect_file("/local/file2.fq")
        assert metadata.size_bytes == 500000
        assert metadata.location_type == "Local"
    
    def test_mock_mode_file_not_in_catalog(self):
        """File not in catalog should return inaccessible."""
        inspector = FileMetadataInspector(mock_mode=True, mock_catalog={})
        
        metadata = inspector.inspect_file("s3://bucket/missing.fq")
        assert metadata.accessible is False
        assert metadata.error == "File not in mock catalog"
        assert metadata.size_bytes is None
    
    def test_inspect_multiple_files(self):
        """Inspector should handle batch of files."""
        catalog = {
            "s3://bucket/file1.fq": 1000000,
            "s3://bucket/file2.fq": 2000000,
            "s3://bucket/file3.fq": 3000000
        }
        
        inspector = FileMetadataInspector(mock_mode=True, mock_catalog=catalog)
        
        uris = ["s3://bucket/file1.fq", "s3://bucket/file2.fq", "s3://bucket/file3.fq"]
        results = inspector.inspect_files(uris)
        
        assert len(results) == 3
        assert results[0].size_bytes == 1000000
        assert results[1].size_bytes == 2000000
        assert results[2].size_bytes == 3000000
        assert all(r.accessible for r in results)
    
    def test_unknown_uri_format(self):
        """Unknown URI format should return error (non-mock mode)."""
        # In non-mock mode, unknown URIs should be rejected with proper error
        inspector = FileMetadataInspector(mock_mode=False)
        
        metadata = inspector.inspect_file("http://example.com/file.fq")
        assert metadata.accessible is False
        assert "Unknown URI format" in metadata.error
        assert metadata.location_type == "Unknown"
    
    def test_load_mock_catalog_from_file(self):
        """Inspector should load catalog from JSON file."""
        # Create temporary catalog file
        catalog_data = {
            "s3://bucket/test1.fq": 100000,
            "s3://bucket/test2.fq": 200000
        }
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(catalog_data, f)
            temp_path = f.name
        
        try:
            inspector = FileMetadataInspector.load_mock_catalog_from_file(temp_path)
            
            metadata = inspector.inspect_file("s3://bucket/test1.fq")
            assert metadata.size_bytes == 100000
            assert metadata.source == "mock_catalog"
        finally:
            Path(temp_path).unlink()
    
    def test_file_metadata_pydantic_validation(self):
        """FileMetadata should validate properly."""
        # Valid metadata
        metadata = FileMetadata(
            uri="s3://bucket/file.fq",
            size_bytes=1000000,
            size_confidence=1.0,
            source="head_object",
            location_type="S3",
            accessible=True
        )
        
        assert metadata.uri == "s3://bucket/file.fq"
        assert metadata.size_bytes == 1000000
        
        # Invalid confidence (out of range)
        with pytest.raises(Exception):  # Pydantic ValidationError
            FileMetadata(
                uri="s3://bucket/file.fq",
                size_confidence=1.5,  # Invalid: >1.0
                source="head_object",
                location_type="S3"
            )


class TestEnvironmentCapabilityCatalog:
    """Test EnvironmentCapabilityCatalog tool."""
    
    def test_load_default_catalog(self):
        """Catalog should load defaults if no config file."""
        catalog = EnvironmentCapabilityCatalog(config_path="/nonexistent/path.yaml")
        
        # Should have 5 default environments
        assert len(catalog.environments) == 5
        assert "Local" in catalog.environments
        assert "EC2" in catalog.environments
        assert "EMR" in catalog.environments
        assert "Batch" in catalog.environments
        assert "Lambda" in catalog.environments
    
    def test_get_environment(self):
        """Get specific environment capabilities."""
        catalog = EnvironmentCapabilityCatalog(config_path="/nonexistent/path.yaml")
        
        local = catalog.get_environment("Local")
        assert local is not None
        assert local.name == "Local"
        assert local.cost_class == "Free"
        assert local.startup_time_seconds == 0.0
        
        emr = catalog.get_environment("EMR")
        assert emr is not None
        assert emr.name == "EMR"
        assert emr.supports_distributed is True
        assert emr.s3_native is True
    
    def test_get_available_environments(self):
        """Get only available environments."""
        catalog = EnvironmentCapabilityCatalog(config_path="/nonexistent/path.yaml")
        
        available = catalog.get_available_environments()
        
        # Local and EMR should be available by default
        names = [env.name for env in available]
        assert "Local" in names
        assert "EMR" in names
        
        # EC2, Batch, Lambda not available by default
        # (EC2 requires HELIX_USE_EC2=true, others not implemented)
    
    def test_filter_environments_by_use_case(self):
        """Filter environments by requirements."""
        catalog = EnvironmentCapabilityCatalog(config_path="/nonexistent/path.yaml")
        
        # Need distributed support
        distributed = catalog.get_environments_for_use_case(supports_distributed=True)
        assert len(distributed) >= 1
        assert all(env.supports_distributed for env in distributed)
        
        # Need S3 native
        s3_native = catalog.get_environments_for_use_case(s3_native=True)
        assert len(s3_native) >= 1
        assert all(env.s3_native for env in s3_native)
        
        # Runtime limit
        short_runtime = catalog.get_environments_for_use_case(max_runtime_minutes=10.0)
        # Should exclude Lambda (15min max)
        for env in short_runtime:
            if env.max_runtime_minutes is not None:
                assert env.max_runtime_minutes >= 10.0
    
    def test_environment_capability_pydantic_validation(self):
        """EnvironmentCapability should validate properly."""
        # Valid capability
        cap = EnvironmentCapability(
            name="Local",
            startup_time_seconds=0.0,
            startup_overhead="None",
            supports_distributed=False,
            supports_gpu=False,
            pre_installed_tools=[],
            container_support=False,
            s3_native=False,
            data_transfer_required=True,
            reproducibility="High",
            debuggability="High",
            failure_blast_radius="Low",
            cost_class="Free",
            cost_model="No cloud costs",
            description="Local execution"
        )
        
        assert cap.name == "Local"
        assert cap.cost_class == "Free"
        
        # Invalid name (not in Literal)
        with pytest.raises(Exception):  # Pydantic ValidationError
            EnvironmentCapability(
                name="Unknown",  # Invalid
                startup_time_seconds=0.0,
                startup_overhead="None",
                supports_distributed=False,
                supports_gpu=False,
                pre_installed_tools=[],
                container_support=False,
                s3_native=False,
                data_transfer_required=True,
                reproducibility="High",
                debuggability="High",
                failure_blast_radius="Low",
                cost_class="Free",
                cost_model="Test",
                description="Test"
            )


class TestCostHeuristicTable:
    """Test CostHeuristicTable tool."""
    
    def test_load_default_table(self):
        """Table should load defaults if no config file."""
        table = CostHeuristicTable(config_path="/nonexistent/path.yaml")
        
        # Should have multiple heuristics
        assert len(table.heuristics) > 0
    
    def test_get_cost_heuristic(self):
        """Get cost heuristic for environment + operation."""
        table = CostHeuristicTable(config_path="/nonexistent/path.yaml")
        
        # Local small (should be free)
        local_small = table.get_cost_heuristic("Local", "Small")
        assert local_small is not None
        assert local_small.cost_range_usd == (0.0, 0.0)
        assert local_small.cost_class == "Free"
        assert local_small.cost_confidence == 1.0
        
        # EC2 medium
        ec2_medium = table.get_cost_heuristic("EC2", "Medium")
        assert ec2_medium is not None
        assert ec2_medium.cost_class == "Low"
        assert ec2_medium.cost_range_usd[0] < ec2_medium.cost_range_usd[1]
        
        # EMR large
        emr_large = table.get_cost_heuristic("EMR", "Large")
        assert emr_large is not None
        assert emr_large.cost_class in ["Medium", "High"]
    
    def test_estimate_cost_for_files(self):
        """Estimate cost based on file sizes."""
        table = CostHeuristicTable(config_path="/nonexistent/path.yaml")
        
        # Small files (<100MB) on Local
        heuristic = table.estimate_cost_for_files("Local", total_size_mb=50)
        assert heuristic is not None
        assert heuristic.operation_class == "Small"
        
        # Medium files (100-10GB) on EC2
        heuristic = table.estimate_cost_for_files("EC2", total_size_mb=500)
        assert heuristic is not None
        assert heuristic.operation_class == "Medium"
        
        # Large files (>10GB) on EMR
        heuristic = table.estimate_cost_for_files("EMR", total_size_mb=15000)
        assert heuristic is not None
        assert heuristic.operation_class in ["Large", "Distributed"]
    
    def test_compare_costs(self):
        """Compare costs across environments."""
        table = CostHeuristicTable(config_path="/nonexistent/path.yaml")
        
        environments = ["Local", "EC2", "EMR"]
        comparison = table.compare_costs(environments, "Medium")
        
        # Should return sorted by cost (low to high)
        assert len(comparison) >= 1
        
        # Local should be cheapest for Medium operations
        if len(comparison) > 1:
            # Check sorting (average cost increasing)
            costs = [sum(h.cost_range_usd) / 2 for _, h in comparison]
            assert costs == sorted(costs)
    
    def test_cost_heuristic_pydantic_validation(self):
        """CostHeuristic should validate properly."""
        # Valid heuristic
        heuristic = CostHeuristic(
            environment="Local",
            operation_class="Small",
            cost_range_usd=(0.0, 0.0),
            cost_class="Free",
            cost_confidence=1.0,
            assumptions="No cloud costs",
            description="Test heuristic"
        )
        
        assert heuristic.environment == "Local"
        assert heuristic.cost_range_usd == (0.0, 0.0)
        
        # Invalid cost_confidence (out of range)
        with pytest.raises(Exception):  # Pydantic ValidationError
            CostHeuristic(
                environment="Local",
                operation_class="Small",
                cost_range_usd=(0.0, 0.0),
                cost_class="Free",
                cost_confidence=1.5,  # Invalid: >1.0
                assumptions="Test",
                description="Test"
            )
        
        # Invalid cost_range (negative)
        with pytest.raises(Exception):  # Pydantic ValidationError or custom validation
            CostHeuristic(
                environment="Local",
                operation_class="Small",
                cost_range_usd=(-1.0, 0.0),  # Invalid: negative
                cost_class="Free",
                cost_confidence=1.0,
                assumptions="Test",
                description="Test"
            )


class TestToolsIntegration:
    """Integration tests using all three tools together."""
    
    def test_infrastructure_decision_with_tools(self):
        """Use all three tools to inform infrastructure decision."""
        # Load all tools
        file_inspector = FileMetadataInspector(
            mock_mode=True,
            mock_catalog={
                "s3://bucket/file_R1.fq": 250000000,
                "s3://bucket/file_R2.fq": 240000000
            }
        )
        
        env_catalog = EnvironmentCapabilityCatalog(config_path="/nonexistent/path.yaml")
        cost_table = CostHeuristicTable(config_path="/nonexistent/path.yaml")
        
        # Inspect files
        files = file_inspector.inspect_files([
            "s3://bucket/file_R1.fq",
            "s3://bucket/file_R2.fq"
        ])
        
        total_size_mb = sum(f.size_bytes for f in files if f.size_bytes) / (1024 * 1024)
        assert total_size_mb == pytest.approx(467.3, rel=0.1)
        
        # Total > 100MB → consider EMR or EC2
        # Get available environments
        available = env_catalog.get_available_environments()
        assert len(available) >= 2
        
        # EMR is best for S3 files >100MB
        emr = env_catalog.get_environment("EMR")
        assert emr is not None
        assert emr.s3_native is True
        assert emr.supports_distributed is True
        
        # Get cost estimate for EMR Medium operation
        emr_cost = cost_table.estimate_cost_for_files("EMR", total_size_mb, file_count=2)
        assert emr_cost is not None
        assert emr_cost.cost_range_usd[0] > 0
        assert emr_cost.cost_class in ["Low", "Medium"]
        
        # Decision: EMR recommended for these files
        # - S3-native (no transfer cost)
        # - Supports distributed processing
        # - Cost is medium but acceptable for file size
    
    def test_small_local_files_decision(self):
        """Small local files should prefer Local execution."""
        file_inspector = FileMetadataInspector(
            mock_mode=True,
            mock_catalog={"/local/data/test.fq": 5000000}
        )
        
        env_catalog = EnvironmentCapabilityCatalog(config_path="/nonexistent/path.yaml")
        cost_table = CostHeuristicTable(config_path="/nonexistent/path.yaml")
        
        # Inspect file
        files = file_inspector.inspect_files(["/local/data/test.fq"])
        total_size_mb = sum(f.size_bytes for f in files if f.size_bytes) / (1024 * 1024)
        assert total_size_mb == pytest.approx(4.77, rel=0.1)
        
        # Get Local capabilities
        local = env_catalog.get_environment("Local")
        assert local is not None
        assert local.cost_class == "Free"
        assert local.startup_time_seconds == 0.0
        
        # Get cost estimate for Local small operation
        local_cost = cost_table.estimate_cost_for_files("Local", total_size_mb, file_count=1)
        assert local_cost is not None
        assert local_cost.cost_range_usd == (0.0, 0.0)
        assert local_cost.cost_class == "Free"
        
        # Decision: Local execution recommended
        # - No cloud costs
        # - Fast (no startup overhead)
        # - File is small (<100MB)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
