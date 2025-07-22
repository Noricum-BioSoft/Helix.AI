# Helix.AI Test Suite

This directory contains all tests for the Helix.AI bioinformatics platform.

## Test Structure

### `/backend/`
Backend-specific tests including:
- **test_agent.py**: Agent functionality tests
- **test_command_handling.py**: Command routing and execution tests
- **test_enhanced_mcp.py**: Enhanced MCP server tests
- **test_history.py**: Session and history management tests
- **test_imports.py**: Import validation tests
- **test_mcp_integration.py**: MCP integration tests
- **test_mcp_startup.py**: MCP server startup tests
- **test_mutation.py**: Sequence mutation tests
- **run_agent.py**: Agent execution runner

### `/frontend/`
Frontend-specific tests including:
- **test-integration.js**: Frontend integration tests

### `/integration/`
End-to-end integration tests including:
- **test-microservices.py**: Microservices integration tests

## Running Tests

### Backend Tests
```bash
cd tests/backend
python -m pytest test_*.py -v
```

### Frontend Tests
```bash
cd tests/frontend
npm test
```

### Integration Tests
```bash
cd tests/integration
python test-microservices.py
```

## Test Data

Test data files are located in:
- `/data/samples/`: Sample sequence files
- `/data/phylogenetic/`: Phylogenetic tree datasets

## Environment Setup

Copy the environment file for backend tests:
```bash
cp tests/backend/.env.example tests/backend/.env
```

## Test Coverage

Tests cover:
- ✅ Natural language command parsing
- ✅ Bioinformatics tool integration
- ✅ MCP server functionality
- ✅ Session management
- ✅ Plasmid visualization
- ✅ Phylogenetic tree generation
- ✅ Sequence clustering
- ✅ Frontend-backend integration 