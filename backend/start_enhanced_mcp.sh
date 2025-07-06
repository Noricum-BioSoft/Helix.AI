#!/bin/bash

# Enhanced Bioinformatics MCP Server Startup Script

echo "🚀 Starting Enhanced Bioinformatics MCP Server..."
echo "================================================"

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is not installed"
    exit 1
fi

# Check if required files exist
if [ ! -f "mcp_server_enhanced.py" ]; then
    echo "❌ Enhanced MCP server file not found"
    exit 1
fi

# Create logs directory if it doesn't exist
mkdir -p logs

# Set environment variables
export PYTHONPATH="."
export MCP_LOG_LEVEL="INFO"
export MCP_VALIDATION_STRICT="true"

echo "✅ Environment configured"
echo "📝 Logging to logs/mcp_server.log"

# Start the enhanced MCP server
echo "🔧 Starting enhanced MCP server..."
python3 mcp_server_enhanced.py 2>&1 | tee logs/mcp_server.log

echo "🛑 Enhanced MCP server stopped" 