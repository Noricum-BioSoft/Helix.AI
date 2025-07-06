#!/usr/bin/env python3
"""
Startup script for the Bioinformatics MCP Server
"""

import asyncio
import sys
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def main():
    """Main entry point for starting the MCP server."""
    try:
        # Import and run the MCP server
        from mcp_server import main as mcp_main
        
        logger.info("Starting Bioinformatics MCP Server...")
        logger.info("Available tools:")
        logger.info("- sequence_alignment: Perform multiple sequence alignment")
        logger.info("- mutate_sequence: Generate sequence mutations")
        logger.info("- analyze_sequence_data: Analyze sequence data")
        logger.info("- visualize_alignment: Create alignment visualizations")
        
        # Run the MCP server
        asyncio.run(mcp_main())
        
    except KeyboardInterrupt:
        logger.info("MCP Server stopped by user")
    except Exception as e:
        logger.error(f"Error starting MCP server: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 