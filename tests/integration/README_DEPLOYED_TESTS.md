# Deployed Backend Test Suite

This directory contains tests adapted from the existing test suite to run against the **deployed cloud backend** instead of a local instance.

## Overview

The `test_deployed_backend.py` file adapts existing backend tests to work with the deployed backend via HTTP API calls. This allows you to:

- ✅ Test the actual deployed backend in AWS
- ✅ Verify fixes are working in production
- ✅ Run integration tests against the live system
- ✅ Monitor backend health and functionality

## Quick Start

### Run All Deployed Tests

```bash
# Using the test runner script
./tests/integration/run_deployed_tests.sh

# Or directly with pytest
pytest tests/integration/test_deployed_backend.py -v
```

### Test Against Custom Backend URL

```bash
# Set environment variable
export BACKEND_URL="https://your-custom-backend.com"
pytest tests/integration/test_deployed_backend.py -v

# Or inline
BACKEND_URL="https://your-custom-backend.com" pytest tests/integration/test_deployed_backend.py -v
```

### Run Specific Tests

```bash
# Test only health check
pytest tests/integration/test_deployed_backend.py::TestDeployedBackend::test_health_check -v

# Test only the fixes
pytest tests/integration/test_deployed_backend.py::TestDeployedBackendFixes -v

# Test only performance
pytest tests/integration/test_deployed_backend.py::TestDeployedBackendPerformance -v
```

## Test Categories

### 1. Core Functionality Tests (`TestDeployedBackend`)

Tests the main backend functionalities:

- ✅ `test_health_check` - Backend health endpoint
- ✅ `test_session_creation` - Session management
- ✅ `test_sequence_alignment` - Sequence alignment tool
- ✅ `test_mutation_generation` - Mutation generation tool
- ✅ `test_phylogenetic_analysis` - Phylogenetic tree building
- ✅ `test_tool_listing` - Available tools endpoint
- ✅ `test_multi_step_workflow` - Multi-step workflows
- ✅ `test_error_handling` - Error handling
- ✅ `test_natural_language_commands` - Natural language processing

### 2. Fix Verification Tests (`TestDeployedBackendFixes`)

Specifically tests the fixes we deployed:

- ✅ `test_bio_align_applications_fix` - Verifies Bio.Align.Applications import fix
- ✅ `test_langchain_core_tools_fix` - Verifies langchain_core.tools import fix

### 3. Performance Tests (`TestDeployedBackendPerformance`)

Tests backend performance:

- ✅ `test_response_time` - Health check response time
- ✅ `test_command_response_time` - Command execution time

## Default Backend URLs

The tests default to the deployed backend:

- **ALB URL**: `http://HelixA-ALBAE-D7BksiQIynZb-1051248867.us-west-1.elb.amazonaws.com`
- **CloudFront URL**: `https://d2a8mt5n89vos4.cloudfront.net` (can be used by setting `BACKEND_URL`)

## Test Structure

### Test Client Class

The `DeployedBackendTester` class provides a clean interface for testing:

```python
tester = DeployedBackendTester(base_url="https://your-backend.com")
tester.health_check()  # Test health
session_id = tester.create_session()  # Create session
result = tester.execute_command("align sequences", session_id)  # Execute command
```

### Pytest Fixtures

- `tester` - Provides a test client instance
- `session` - Provides a test session ID

### Example Test

```python
def test_sequence_alignment(tester, session):
    """Test sequence alignment"""
    command = "align these sequences: >seq1 ATGCGATCG >seq2 ATGCGATCA"
    result = tester.execute_command(command, session)
    
    assert result.get("success") is True
    assert result.get("tool") == "sequence_alignment"
```

## Comparison with Local Tests

| Feature | Local Tests | Deployed Tests |
|---------|------------|----------------|
| **Target** | Local backend (`localhost:8001`) | Deployed cloud backend |
| **Method** | Direct imports or TestClient | HTTP API calls |
| **Speed** | Fast (no network) | Slower (network latency) |
| **Environment** | Test environment | Production environment |
| **Use Case** | Development, CI/CD | Deployment verification |

## Running in CI/CD

You can integrate these tests into your CI/CD pipeline:

```yaml
# Example GitHub Actions
- name: Test Deployed Backend
  env:
    BACKEND_URL: ${{ secrets.BACKEND_URL }}
  run: |
    pytest tests/integration/test_deployed_backend.py -v
```

## Troubleshooting

### Connection Errors

If you get connection errors:

1. **Check backend URL**: Verify the backend is accessible
   ```bash
   curl $BACKEND_URL/health
   ```

2. **Check network**: Ensure you can reach the backend
   ```bash
   ping $(echo $BACKEND_URL | sed 's|https\?://||' | sed 's|/.*||')
   ```

3. **Check timeout**: Increase timeout if backend is slow
   ```python
   tester = DeployedBackendTester(timeout=60)  # 60 seconds
   ```

### Test Failures

If tests fail:

1. **Check backend logs**: Look at CloudWatch logs for errors
2. **Verify deployment**: Ensure latest code is deployed
3. **Check API changes**: Verify API hasn't changed
4. **Review test output**: Check detailed error messages

## Adding New Tests

To add new tests for deployed backend:

1. **Add to existing test class**:
   ```python
   @pytest.mark.integration
   def test_new_feature(self, tester, session):
       result = tester.execute_command("your command", session)
       assert result.get("success") is True
   ```

2. **Or create new test class**:
   ```python
   class TestNewFeature:
       @pytest.mark.integration
       def test_something(self, tester):
           # Your test
   ```

## Related Files

- `test_core_functionalities.py` - Original local test suite
- `test_natural_language_mapping.py` - Natural language command tests
- `test_deployed_backend.py` - This file (deployed backend tests)

## Notes

- These tests make actual HTTP requests to the deployed backend
- Tests may be slower due to network latency
- Tests use real sessions and may create data in the backend
- Consider using test-specific session IDs or cleanup after tests


