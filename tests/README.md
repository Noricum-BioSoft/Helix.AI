# 🧪 Helix.AI Test Suite

This directory contains the comprehensive test suite for the Helix.AI bioinformatics platform.

## 📁 Test Structure

```
tests/
├── backend/              # Backend-specific unit tests
│   ├── test_agent.py    # LangChain agent tests
│   ├── test_command_handling.py  # Command handling tests
│   ├── test_history.py  # History manager tests
│   ├── test_imports.py  # Import validation tests
│   ├── test_mutation.py # Mutation tool tests
│   └── test_session_history.py  # Session history unit tests
├── frontend/            # Frontend-specific tests
│   ├── test_components/ # React component tests
│   ├── test_services/   # API service tests
│   └── test_utils/      # Utility function tests
├── integration/         # End-to-end integration tests
│   ├── test_core_functionalities.py  # Core functionality tests via API
│   ├── test_natural_language_mapping.py  # Natural language command mapping tests
│   └── test-microservices.py  # Microservices integration tests
├── workflows/          # End-to-end workflow tests
│   ├── test_rnaseq_workflow.py              # RNA-seq preprocessing workflow
│   ├── test_alignment_phylogenetic_workflow.py  # Alignment & phylogeny workflow
│   ├── test_variant_analysis_workflow.py    # Variant analysis workflow
│   └── test_workflow_permutations.py       # Comprehensive workflow permutations
├── data/               # Test data files
│   ├── sequences/      # Test sequence files
│   ├── alignments/     # Test alignment files
│   └── trees/          # Test phylogenetic trees
└── README.md           # This documentation
```

## 🚀 Running Tests

### Prerequisites

```bash
# Install test dependencies
pip install pytest pytest-asyncio pytest-cov
npm install --prefix frontend
```

### Backend Tests

```bash
# Run all backend tests
python -m pytest tests/backend/ -v

# Run specific test file
python -m pytest tests/backend/test_command_router.py -v

# Run with coverage
python -m pytest tests/backend/ --cov=backend --cov-report=html
```

### Workflow Tests (End-to-End)

```bash
# Run all workflow tests
python -m pytest tests/workflows/ -v

# Run specific workflow test
python3 tests/workflows/test_rnaseq_workflow.py
python3 tests/workflows/test_alignment_phylogenetic_workflow.py
python3 tests/workflows/test_variant_analysis_workflow.py

# Run comprehensive permutation test
python3 tests/workflows/test_workflow_permutations.py
```

**Note**: Workflow tests require the backend server to be running on `http://localhost:8001`

### Frontend Tests

```bash
# Run frontend tests
cd frontend && npm test

# Run with coverage
cd frontend && npm test -- --coverage

# Run specific test
cd frontend && npm test -- --testNamePattern="PlasmidVisualizer"
```

### Integration Tests

```bash
# Run all integration tests
python -m pytest tests/integration/ -v

# Run specific integration test
python -m pytest tests/integration/test_core_functionalities.py -v
python -m pytest tests/integration/test_natural_language_mapping.py -v

# Or run directly
python3 tests/integration/test_core_functionalities.py
python3 tests/integration/test_natural_language_mapping.py
```

**Note**: Integration tests require the backend server to be running on `http://localhost:8001`

### Complete Test Suite

```bash
# Run all tests
python -m pytest tests/ -v
cd frontend && npm test
```

## 🧬 Test Categories

### 1. **Unit Tests** (`tests/backend/`, `tests/frontend/`)

**Purpose**: Test individual functions and components in isolation

**Examples**:
- Command parsing logic
- Sequence alignment algorithms
- React component rendering
- API service functions

**Running**:
```bash
# Backend unit tests
python -m pytest tests/backend/test_command_router.py -v

# Frontend unit tests
cd frontend && npm test -- --testNamePattern="CommandParser"
```

### 2. **Integration Tests** (`tests/integration/`)

**Purpose**: Test complete workflows and system interactions

**Examples**:
- Complete phylogenetic analysis workflow
- File upload → analysis → visualization pipeline
- Session management across multiple commands
- Natural language command processing

**Running**:
```bash
python -m pytest tests/integration/test_workflows/ -v
```

### 3. **API Tests** (`tests/integration/test_api/`)

**Purpose**: Test REST API endpoints and MCP integration

**Examples**:
- Endpoint response validation
- Session creation and management
- Natural language command processing
- File upload and processing

**Running**:
```bash
python -m pytest tests/integration/test_api/ -v
```

### 4. **UI Tests** (`tests/frontend/`, `tests/integration/test_ui/`)

**Purpose**: Test user interface components and interactions

**Examples**:
- Component rendering
- User interactions (click, type, drag-drop)
- State management
- API integration

**Running**:
```bash
cd frontend && npm test
```

## 📊 Test Coverage

### Current Coverage Targets

- **Backend**: >90% line coverage
- **Frontend**: >80% line coverage
- **Integration**: >95% workflow coverage

### Coverage Reports

```bash
# Generate coverage reports
python -m pytest tests/ --cov=backend --cov=frontend --cov-report=html
cd frontend && npm test -- --coverage --coverageReporters=html
```

Reports will be generated in:
- `htmlcov/` (backend coverage)
- `frontend/coverage/` (frontend coverage)

## 🧪 Test Data

### Test Sequences

Located in `tests/data/sequences/`:
- `test_sequences.fasta`: Basic test sequences
- `phylogenetic_test.fasta`: Sequences for tree building
- `alignment_test.fasta`: Sequences for alignment testing

### Test Alignments

Located in `tests/data/alignments/`:
- `test_alignment.fasta`: Pre-aligned sequences
- `clustal_output.aln`: ClustalW output format

### Test Trees

Located in `tests/data/trees/`:
- `test_tree.newick`: Newick format tree
- `complex_tree.newick`: Multi-branch tree

## 🔧 Writing Tests

### Backend Test Example

```python
# tests/backend/test_command_router.py
import pytest
from backend.command_router import CommandRouter

def test_plasmid_command_routing():
    router = CommandRouter()
    command = "insert each of the sequences into pUC19"
    
    tool, params = router.route_command(command)
    
    assert tool == "plasmid_visualization"
    assert "vector_name" in params
    assert params["vector_name"] == "pUC19"
```

### Frontend Test Example

```typescript
// tests/frontend/test_components/PlasmidVisualizer.test.tsx
import { render, screen } from '@testing-library/react';
import { PlasmidDataVisualizer } from '../../src/components/PlasmidVisualizer';

test('renders plasmid visualization with correct data', () => {
  const testData = {
    name: 'Test Plasmid',
    sequence: 'ATGCGATCG',
    features: []
  };
  
  render(<PlasmidDataVisualizer data={testData} />);
  
  expect(screen.getByText('Test Plasmid')).toBeInTheDocument();
});
```

### Integration Test Example

```python
# tests/integration/test_workflows/test_phylogenetic_workflow.py
import pytest
from backend.main_with_mcp import app
from fastapi.testclient import TestClient

def test_complete_phylogenetic_workflow():
    client = TestClient(app)
    
    # Step 1: Create session
    response = client.post("/create_session")
    session_id = response.json()["session_id"]
    
    # Step 2: Upload sequences
    sequences = ">seq1 ATGCGATCG\n>seq2 ATGCGATC"
    response = client.post("/mcp/phylogenetic-tree", 
                          json={"aligned_sequences": sequences, "session_id": session_id})
    assert response.status_code == 200
    
    # Step 3: Select representatives
    response = client.post("/mcp/clustering-analysis",
                          json={"num_clusters": 5, "session_id": session_id})
    assert response.status_code == 200
```

## 🐛 Debugging Tests

### Common Issues

1. **Import Errors**:
   ```bash
   # Add project root to Python path
   export PYTHONPATH="${PYTHONPATH}:$(pwd)"
   ```

2. **Frontend Test Failures**:
   ```bash
   # Clear node modules and reinstall
   cd frontend && rm -rf node_modules package-lock.json
   npm install
   ```

3. **API Connection Issues**:
   ```bash
   # Ensure backend is running
   python backend/main_with_mcp.py &
   ```

### Debug Mode

```bash
# Run tests with debug output
python -m pytest tests/ -v -s --tb=long

# Run specific test with debug
python -m pytest tests/backend/test_command_router.py::test_specific_function -v -s
```

## 📈 Performance Testing

### Load Testing

```bash
# Test API performance
python -m pytest tests/integration/test_performance.py -v
```

### Memory Testing

```bash
# Test memory usage
python -m pytest tests/integration/test_memory.py -v
```

## 🔄 Continuous Integration

### GitHub Actions

Tests are automatically run on:
- Pull requests
- Push to main branch
- Release tags

### Local CI

```bash
# Run full CI locally
./scripts/run_ci.sh
```

## 📚 Test Documentation

### Adding New Tests

1. **Unit Tests**: Add to appropriate `tests/backend/` or `tests/frontend/` directory
2. **Integration Tests**: Add to `tests/integration/` with descriptive names
3. **Test Data**: Add required test files to `tests/data/`
4. **Documentation**: Update this README with new test information

### Test Naming Convention

- **Backend**: `test_<module>_<function>.py`
- **Frontend**: `test_<Component>.test.tsx`
- **Integration**: `test_<workflow>_workflow.py`

### Test Data Management

- Keep test data minimal but representative
- Use realistic but synthetic data
- Document data sources and formats
- Version control test data files

## 🎯 Test Quality Checklist

- [ ] All new features have corresponding tests
- [ ] Tests cover both success and failure cases
- [ ] Integration tests cover complete workflows
- [ ] Test data is properly documented
- [ ] Coverage targets are met
- [ ] Tests run in reasonable time (< 5 minutes)
- [ ] Tests are independent and repeatable
- [ ] Error messages are clear and helpful 