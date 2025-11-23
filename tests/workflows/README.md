# End-to-End Workflow Tests

This directory contains comprehensive end-to-end workflow tests that validate complete bioinformatics analysis pipelines in Helix.AI.

## Test Scripts

### Individual Workflow Tests

1. **`test_rnaseq_workflow.py`**
   - Tests the RNA-seq preprocessing workflow (Workflow A from NGS101)
   - Steps: Quality trimming → Adapter removal → Merging → Quality assessment
   - Uses paired-end FASTQ files from `data/rnaseq_demo/`

2. **`test_alignment_phylogenetic_workflow.py`**
   - Tests sequence alignment and phylogenetic analysis (Workflow B from NGS101)
   - Steps: Align sequences → Build tree → Select representatives
   - Generates test RNA sequences programmatically

3. **`test_variant_analysis_workflow.py`**
   - Tests variant analysis workflow (Workflow C from NGS101)
   - Steps: Generate variants → Select diverse → Visualize in plasmid
   - Generates test DNA sequence programmatically

### Comprehensive Permutation Test

4. **`test_workflow_permutations.py`**
   - Tests 13 different biologically meaningful workflow patterns
   - Validates command chaining and data passing between steps
   - Provides success rate and error reporting

## Running the Tests

### Prerequisites

1. **Backend Server**: Must be running on `http://localhost:8001`
   ```bash
   cd backend
   uvicorn main_with_mcp:app --reload
   ```

2. **Test Data**: FASTQ files should be in `data/rnaseq_demo/`
   - `sample_R1_realistic.fastq`
   - `sample_R2_realistic.fastq`

### Run Individual Tests

```bash
# From project root
python3 tests/workflows/test_rnaseq_workflow.py
python3 tests/workflows/test_alignment_phylogenetic_workflow.py
python3 tests/workflows/test_variant_analysis_workflow.py
```

### Run All Workflow Tests

```bash
# Run comprehensive permutation test
python3 tests/workflows/test_workflow_permutations.py
```

### Run All Tests with pytest

```bash
# Run all workflow tests
pytest tests/workflows/ -v

# Run specific test
pytest tests/workflows/test_rnaseq_workflow.py -v
```

## Test Structure

Each test script follows a similar structure:

1. **Setup**: Create session, load/generate test data
2. **Execution**: Run workflow steps sequentially
3. **Validation**: Check results and data passing
4. **Summary**: Display session history and statistics

## Expected Results

All tests should:
- ✅ Complete all workflow steps successfully
- ✅ Show proper session history (correct number of entries)
- ✅ Display detailed results for each step
- ✅ Demonstrate data passing between steps

## Troubleshooting

### Backend Not Running
```
✗ Request failed: Connection refused
```
**Solution**: Start the backend server first

### Missing Test Data
```
✗ Failed to load FASTQ file
```
**Solution**: Ensure FASTQ files exist in `data/rnaseq_demo/`

### Session History Empty
```
⚠ Session history is empty
```
**Solution**: Check backend logs for session management issues

### Data Not Passing Between Steps
```
✗ Error: No previous results found
```
**Solution**: Verify that previous steps completed successfully and saved results

## Adding New Workflow Tests

To add a new workflow test:

1. Create a new test file: `test_<workflow_name>_workflow.py`
2. Follow the structure of existing tests
3. Define workflow steps with natural language commands
4. Include validation and summary reporting
5. Update this README

Example structure:
```python
#!/usr/bin/env python3
"""
Test script for <Workflow Name> workflow.
"""

import requests
import json
from pathlib import Path

BACKEND_URL = "http://localhost:8001"

def main():
    # Create session
    # Execute workflow steps
    # Validate results
    # Display summary

if __name__ == "__main__":
    main()
```

## Related Documentation

- [NGS101 Demonstration Strategy](../../docs/NGS101_DEMONSTRATION_STRATEGY.md)
- [Workflow Test Scripts Documentation](../../docs/demos/WORKFLOW_TEST_SCRIPTS.md)
- [Workflow Permutations Test](../../docs/demos/WORKFLOW_PERMUTATIONS_TEST.md)

---

**Last Updated:** 2025-01-27

