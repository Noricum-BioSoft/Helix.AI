# Bioinformatics MCP Server

This directory contains the Model Context Protocol (MCP) server implementation for bioinformatics tasks, integrated with your existing Helix.AI backend.

## Overview

The MCP server provides a standardized interface for bioinformatics operations, including:
- Sequence alignment
- Sequence mutation analysis
- Sequence data analysis and visualization
- Alignment visualization

## Files

- `mcp_server.py` - Main MCP server implementation
- `main_with_mcp.py` - FastAPI server with MCP integration
- `mcp_config.json` - MCP server configuration
- `start_mcp_server.py` - Startup script for the MCP server

## Installation

1. Install the required dependencies:
```bash
pip install -r requirements.txt
```

2. Ensure you have the bioinformatics tools installed:
```bash
# For Clustal Omega (sequence alignment)
brew install clustalo  # macOS
# or
sudo apt-get install clustalo  # Ubuntu/Debian
```

## Usage

### Starting the MCP Server

```bash
cd backend
python start_mcp_server.py
```

### Using the FastAPI Server with MCP

```bash
cd backend
python main_with_mcp.py
```

The server will be available at `http://localhost:8001`

## API Endpoints

### General Command Execution
```http
POST /execute
Content-Type: application/json

{
  "command": "align the given sequences"
}
```

### MCP-Specific Endpoints

#### Sequence Alignment
```http
POST /mcp/sequence-alignment
Content-Type: application/json

{
  "sequences": ">seq1\nACTGTTGAC\n>seq2\nACTGCATCC\n>seq3\nACTGCAATGAC",
  "algorithm": "clustal"
}
```

#### Sequence Mutation
```http
POST /mcp/mutate-sequence
Content-Type: application/json

{
  "sequence": "ACTGTTGAC",
  "num_variants": 96,
  "mutation_rate": 0.1
}
```

#### Sequence Data Analysis
```http
POST /mcp/analyze-sequence-data
Content-Type: application/json

{
  "data": ">seq1\nACTGTTGAC\n>seq2\nACTGCATCC",
  "analysis_type": "alignment"
}
```

#### Alignment Visualization
```http
POST /mcp/visualize-alignment
Content-Type: application/json

{
  "alignment_file": "path/to/alignment.fasta",
  "output_format": "png"
}
```

#### List Available Tools
```http
GET /mcp/tools
```

## MCP Tools

### 1. sequence_alignment
Performs multiple sequence alignment on DNA/RNA sequences.

**Parameters:**
- `sequences` (string): FASTA format sequences
- `algorithm` (string): Alignment algorithm (clustal, muscle, mafft)

**Example:**
```python
result = await call_mcp_tool("sequence_alignment", {
    "sequences": ">seq1\nACTGTTGAC\n>seq2\nACTGCATCC",
    "algorithm": "clustal"
})
```

### 2. mutate_sequence
Generates mutations and variants of a given sequence.

**Parameters:**
- `sequence` (string): Input DNA/RNA sequence
- `num_variants` (integer): Number of variants to generate (default: 96)
- `mutation_rate` (float): Mutation rate (default: 0.1)

**Example:**
```python
result = await call_mcp_tool("mutate_sequence", {
    "sequence": "ACTGTTGAC",
    "num_variants": 96,
    "mutation_rate": 0.1
})
```

### 3. analyze_sequence_data
Analyzes sequence data and generates visualizations.

**Parameters:**
- `data` (string): Sequence data in CSV or FASTA format
- `analysis_type` (string): Type of analysis (alignment, phylogeny, composition)

**Example:**
```python
result = await call_mcp_tool("analyze_sequence_data", {
    "data": "path/to/sequences.csv",
    "analysis_type": "alignment"
})
```

### 4. visualize_alignment
Creates visual representation of sequence alignments.

**Parameters:**
- `alignment_file` (string): Path to alignment file (FASTA format)
- `output_format` (string): Output image format (png, svg, pdf)

**Example:**
```python
result = await call_mcp_tool("visualize_alignment", {
    "alignment_file": "path/to/alignment.fasta",
    "output_format": "png"
})
```

## Integration with Existing Tools

The MCP server integrates with your existing bioinformatics tools:

- `tools/alignment.py` - Sequence alignment functionality
- `tools/bio.py` - Advanced bioinformatics operations
- `tools/mutations.py` - Sequence mutation analysis
- `tools/data_science.py` - Data analysis capabilities

## Configuration

The MCP server can be configured through `mcp_config.json`:

```json
{
  "mcpServers": {
    "bioinformatics": {
      "command": "python",
      "args": ["mcp_server.py"],
      "env": {
        "PYTHONPATH": "."
      }
    }
  }
}
```

## Development

### Adding New Tools

1. Create the tool function in the appropriate `tools/` file
2. Add the tool to the MCP server in `mcp_server.py`
3. Update the API endpoints in `main_with_mcp.py`
4. Add documentation to this README

### Testing

```bash
# Test the MCP server
python -m pytest tests/test_mcp.py

# Test the API endpoints
python -m pytest tests/test_api.py
```

## Troubleshooting

### Common Issues

1. **Import Errors**: Ensure all dependencies are installed and the Python path is set correctly
2. **Tool Not Found**: Check that the tool is properly registered in the MCP server
3. **Alignment Failures**: Verify that Clustal Omega is installed and accessible

### Logging

The MCP server uses Python's logging module. Set the log level to DEBUG for detailed output:

```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

## Contributing

When adding new bioinformatics tools:

1. Follow the existing pattern in `tools/` directory
2. Add proper error handling and validation
3. Include comprehensive documentation
4. Add tests for new functionality
5. Update this README with new tool information

## License

This MCP server implementation is part of the Helix.AI project and follows the same license terms. 