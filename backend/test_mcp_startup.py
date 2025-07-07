#!/usr/bin/env python3
"""
Test script to verify MCP server startup
"""

import sys
import asyncio
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def test_imports():
    """Test if all required imports are available."""
    try:
        logger.info("Testing imports...")
        
        # Test MCP imports
        from mcp.server import Server
        from mcp.server.stdio import stdio_server
        from mcp.types import (
            CallToolRequest,
            CallToolResult,
            ListToolsRequest,
            ListToolsResult,
            Tool,
            TextContent,
            ImageContent,
            EmbeddedResource,
        )
        logger.info("✅ MCP imports successful")
        
        # Test tool imports
        sys.path.append(str(Path(__file__).parent.parent / "tools"))
        from tools.alignment import run_alignment
        from tools.bio import align_and_visualize_fasta
        from tools.mutations import run_mutation
        from tools.data_science import analyze_basic_stats
        logger.info("✅ Tool imports successful")
        
        return True
        
    except ImportError as e:
        logger.error(f"❌ Import error: {e}")
        return False
    except Exception as e:
        logger.error(f"❌ Unexpected error: {e}")
        return False

def test_server_creation():
    """Test if the MCP server can be created."""
    try:
        logger.info("Testing server creation...")
        
        from mcp.server import Server
        server = Server("bioinformatics-mcp-server")
        logger.info("✅ Server creation successful")
        
        return True
        
    except Exception as e:
        logger.error(f"❌ Server creation error: {e}")
        return False

async def test_async_functions():
    """Test if async functions work properly."""
    try:
        logger.info("Testing async functions...")
        
        # Import the server module
        sys.path.append(str(Path(__file__).parent))
        from mcp_server import handle_list_tools, handle_call_tool
        
        # Test list tools
        result = await handle_list_tools()
        logger.info(f"✅ List tools successful: {len(result.tools)} tools found")
        
        # Test call tool with dummy data
        dummy_result = await handle_call_tool("sequence_alignment", {"sequences": "ATCG"})
        logger.info("✅ Call tool successful")
        
        return True
        
    except Exception as e:
        logger.error(f"❌ Async function error: {e}")
        return False

async def main():
    """Main test function."""
    logger.info("🧪 Starting MCP server tests...")
    
    # Test imports
    if not test_imports():
        logger.error("❌ Import tests failed")
        return False
    
    # Test server creation
    if not test_server_creation():
        logger.error("❌ Server creation tests failed")
        return False
    
    # Test async functions
    if not await test_async_functions():
        logger.error("❌ Async function tests failed")
        return False
    
    logger.info("✅ All tests passed! MCP server should start properly.")
    return True

if __name__ == "__main__":
    success = asyncio.run(main())
    sys.exit(0 if success else 1) 