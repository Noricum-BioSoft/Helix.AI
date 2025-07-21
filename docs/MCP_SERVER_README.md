# MCP Server Setup and Troubleshooting

This document explains how the Model Context Protocol (MCP) server works in Helix.AI and how to troubleshoot common issues.

## Overview

The MCP server provides bioinformatics tools through the Model Context Protocol, allowing AI models to interact with sequence analysis, mutation generation, and data visualization tools.

## Architecture

```
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│   Frontend      │    │   FastAPI       │    │   MCP Server    │
│   (React)       │◄──►│   Server        │◄──►│   (Bioinformatics│
│                 │    │   (Port 8001)   │    │    Tools)       │
└─────────────────┘    └─────────────────┘    └─────────────────┘
```

## Files

- `simple_mcp_server.py` - Simplified MCP server with better error handling
- `mcp_server.py` - Original MCP server implementation
- `start_mcp_server.py` - Startup script for the original MCP server
- `diagnose_mcp.py` - Diagnostic script to identify issues
- `test_mcp_startup.py` - Test script for MCP server functionality

## Startup Process

The `start-app.sh` script follows this sequence:

1. **Prerequisites Check** - Python, Node.js, npm
2. **Dependency Installation** - Backend and frontend dependencies
3. **MCP Server Start** - Tries multiple MCP server implementations:
   - `simple_mcp_server.py` (preferred)
   - `start_mcp_server.py` (fallback)
   - `mcp_server.py` (fallback)
4. **FastAPI Server Start** - `main_with_mcp.py`
5. **Frontend Server Start** - React development server

## Available Tools

The MCP server provides these bioinformatics tools:

### 1. Sequence Alignment
- **Tool**: `sequence_alignment`
- **Description**: Perform multiple sequence alignment
- **Parameters**:
  - `sequences`: FASTA format sequences
  - `algorithm`: clustal, muscle, or mafft

### 2. Sequence Mutation
- **Tool**: `mutate_sequence`
- **Description**: Generate sequence variants
- **Parameters**:
  - `sequence`: Input DNA/RNA sequence
  - `num_variants`: Number of variants (default: 96)
  - `mutation_rate`: Mutation rate (default: 0.1)

### 3. Sequence Analysis
- **Tool**: `analyze_sequence_data`
- **Description**: Analyze sequence data and generate visualizations
- **Parameters**:
  - `data`: CSV or FASTA format data
  - `analysis_type`: alignment, phylogeny, or composition

### 4. Alignment Visualization
- **Tool**: `visualize_alignment`
- **Description**: Create visual representations of alignments
- **Parameters**:
  - `alignment_file`: Path to FASTA file
  - `output_format`: png, svg, or pdf

## Troubleshooting

### Common Issues

#### 1. MCP Server Won't Start

**Symptoms**: 
- "MCP server failed to start" error
- No MCP server process running

**Solutions**:
1. Run diagnostics: `python diagnose_mcp.py`
2. Check MCP installation: `pip install mcp>=1.0.0`
3. Verify dependencies: `pip install -r requirements.txt`
4. Check Python version: `python --version` (requires 3.8+)

#### 2. Import Errors

**Symptoms**:
- "Failed to import MCP dependencies" error
- ModuleNotFoundError for bioinformatics tools

**Solutions**:
1. Install MCP: `pip install mcp>=1.0.0`
2. Install bioinformatics tools: `pip install biopython pandas numpy matplotlib seaborn`
3. Check tools directory exists: `ls backend/tools/`

#### 3. Port Conflicts

**Symptoms**:
- "Port 8001 is already in use" error
- FastAPI server won't start

**Solutions**:
1. Find and kill existing process: `lsof -ti:8001 | xargs kill -9`
2. Use different port: Modify `main_with_mcp.py` port number
3. Restart terminal/computer

#### 4. Frontend Can't Connect

**Symptoms**:
- Frontend shows "Server not connected"
- API calls fail

**Solutions**:
1. Verify FastAPI server is running: `curl http://localhost:8001/health`
2. Check CORS settings in `main_with_mcp.py`
3. Verify frontend is using correct API URL

### Diagnostic Commands

```bash
# Check if MCP server is running
ps aux | grep mcp_server

# Test MCP server directly
cd backend
python simple_mcp_server.py

# Run diagnostics
python diagnose_mcp.py

# Check ports
lsof -i :8001
lsof -i :5173

# Test FastAPI server
curl http://localhost:8001/health
```

### Manual Startup

If the automated startup fails, you can start components manually:

```bash
# Terminal 1: Start MCP server
cd backend
python simple_mcp_server.py

# Terminal 2: Start FastAPI server
cd backend
python main_with_mcp.py

# Terminal 3: Start frontend
cd frontend
npm run dev
```

## Configuration

### MCP Server Configuration

The MCP server can be configured through environment variables:

```bash
# Set log level
export MCP_LOG_LEVEL=INFO

# Set server name
export MCP_SERVER_NAME=bioinformatics-mcp-server

# Set server version
export MCP_SERVER_VERSION=1.0.0
```

### FastAPI Configuration

The FastAPI server configuration is in `main_with_mcp.py`:

```python
# Server settings
HOST = "0.0.0.0"
PORT = 8001
RELOAD = True

# CORS settings
ALLOW_ORIGINS = ["*"]
ALLOW_METHODS = ["*"]
ALLOW_HEADERS = ["*"]
```

## Development

### Adding New Tools

To add a new bioinformatics tool:

1. **Create the tool function** in `tools/` directory
2. **Add tool definition** in MCP server:
   ```python
   Tool(
       name="new_tool",
       description="Description of the tool",
       inputSchema={
           "type": "object",
           "properties": {
               "param1": {"type": "string"}
           },
           "required": ["param1"]
       }
   )
   ```
3. **Add handler function**:
   ```python
   async def handle_new_tool(arguments: Dict[str, Any]) -> Dict[str, Any]:
       # Tool implementation
       return {"result": "success"}
   ```
4. **Update tool routing** in `handle_call_tool()`

### Testing

Test the MCP server:

```bash
# Run unit tests
cd backend
python -m pytest tests/

# Test specific tool
python -c "
from tools.alignment import run_alignment
result = run_alignment('>seq1\nATCG\n>seq2\nGCTA')
print(result)
"
```

## Logs

MCP server logs are available at different levels:

- **INFO**: General operation messages
- **DEBUG**: Detailed debugging information
- **ERROR**: Error messages and stack traces

To enable debug logging:

```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

## Support

If you encounter issues:

1. **Run diagnostics**: `python diagnose_mcp.py`
2. **Check logs**: Look for error messages in terminal output
3. **Verify dependencies**: Ensure all required packages are installed
4. **Test manually**: Try starting components individually
5. **Check documentation**: Review this README and related files

## Related Files

- `start-app.sh` - Main startup script
- `main_with_mcp.py` - FastAPI server with MCP integration
- `frontend/src/mcpApi.ts` - Frontend MCP API client
- `frontend/src/commandParser.ts` - Command routing logic 