#!/usr/bin/env python3
"""
MCP Server for Bioinformatics Tasks
Provides tools for sequence alignment, mutation analysis, and other bioinformatics operations.
"""

import asyncio
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from mcp.server import Server
from mcp.server.models import InitializationOptions
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

tools_path = str((Path(__file__).resolve().parent.parent / "tools").resolve())
sys.path.insert(0, tools_path)

from tools.alignment import run_alignment
from tools.bio import align_and_visualize_fasta
from tools.mutations import run_mutation_raw
from tools.data_science import analyze_basic_stats
from tools.dna_vendor_research import run_dna_vendor_research_raw

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize the MCP server
server = Server("bioinformatics-mcp-server")

@server.list_tools()
async def handle_list_tools() -> ListToolsResult:
    """List all available bioinformatics tools."""
    return ListToolsResult(
        tools=[
            Tool(
                name="sequence_alignment",
                description="Perform multiple sequence alignment on a set of DNA/RNA sequences",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "sequences": {
                            "type": "string",
                            "description": "FASTA format sequences or sequence data"
                        },
                        "algorithm": {
                            "type": "string",
                            "enum": ["clustal", "muscle", "mafft"],
                            "default": "clustal",
                            "description": "Alignment algorithm to use"
                        }
                    },
                    "required": ["sequences"]
                }
            ),
            Tool(
                name="mutate_sequence",
                description="Generate mutations and variants of a given sequence",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "sequence": {
                            "type": "string",
                            "description": "The input DNA/RNA sequence to mutate"
                        },
                        "num_variants": {
                            "type": "integer",
                            "default": 96,
                            "description": "Number of variants to generate"
                        },
                        "mutation_rate": {
                            "type": "number",
                            "default": 0.1,
                            "description": "Mutation rate (0.0 to 1.0)"
                        }
                    },
                    "required": ["sequence"]
                }
            ),
            Tool(
                name="analyze_sequence_data",
                description="Analyze sequence data and generate visualizations",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "data": {
                            "type": "string",
                            "description": "Sequence data in CSV or FASTA format"
                        },
                        "analysis_type": {
                            "type": "string",
                            "enum": ["alignment", "phylogeny", "composition"],
                            "default": "alignment",
                            "description": "Type of analysis to perform"
                        }
                    },
                    "required": ["data"]
                }
            ),
            Tool(
                name="visualize_alignment",
                description="Create visual representation of sequence alignments",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "alignment_file": {
                            "type": "string",
                            "description": "Path to alignment file (FASTA format)"
                        },
                        "output_format": {
                            "type": "string",
                            "enum": ["png", "svg", "pdf"],
                            "default": "png",
                            "description": "Output image format"
                        }
                    },
                    "required": ["alignment_file"]
                }
            ),
            Tool(
                name="dna_vendor_research",
                description="Research DNA synthesis vendors and testing options for experimental validation",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "command": {
                            "type": "string",
                            "description": "Research command (e.g., 'order sequences', 'find testing options')"
                        },
                        "sequence_length": {
                            "type": "integer",
                            "description": "Length of sequences to synthesize (optional)"
                        },
                        "quantity": {
                            "type": "string",
                            "enum": ["standard", "large", "custom"],
                            "default": "standard",
                            "description": "Quantity needed for synthesis"
                        }
                    },
                    "required": ["command"]
                }
            )
        ]
    )

@server.call_tool()
async def handle_call_tool(name: str, arguments: Dict[str, Any]) -> CallToolResult:
    """Handle tool calls for bioinformatics operations."""
    
    try:
        if name == "sequence_alignment":
            result = await handle_sequence_alignment(arguments)
        elif name == "mutate_sequence":
            result = await handle_mutate_sequence(arguments)
        elif name == "analyze_sequence_data":
            result = await handle_analyze_sequence_data(arguments)
        elif name == "visualize_alignment":
            result = await handle_visualize_alignment(arguments)
        elif name == "dna_vendor_research":
            result = await handle_dna_vendor_research(arguments)
        else:
            raise ValueError(f"Unknown tool: {name}")
        
        return CallToolResult(
            content=[
                TextContent(
                    type="text",
                    text=json.dumps(result, indent=2)
                )
            ]
        )
    
    except Exception as e:
        logger.error(f"Error in tool call {name}: {e}")
        return CallToolResult(
            content=[
                TextContent(
                    type="text",
                    text=f"Error: {str(e)}"
                )
            ]
        )

async def handle_sequence_alignment(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Handle sequence alignment requests."""
    sequences = arguments.get("sequences", "")
    algorithm = arguments.get("algorithm", "clustal")
    
    # Use the existing alignment tool
    result = run_alignment(sequences)
    
    return {
        "status": "success",
        "tool": "sequence_alignment",
        "algorithm": algorithm,
        "result": result
    }

async def handle_mutate_sequence(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Handle sequence mutation requests."""
    sequence = arguments.get("sequence", "")
    num_variants = arguments.get("num_variants", 96)
    mutation_rate = arguments.get("mutation_rate", 0.1)
    
    # Use the existing mutation tool
    result = run_mutation_raw(sequence, num_variants)
    
    return {
        "status": "success",
        "tool": "mutate_sequence",
        "input_sequence": sequence,
        "num_variants": num_variants,
        "mutation_rate": mutation_rate,
        "result": result
    }

async def handle_analyze_sequence_data(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Handle sequence data analysis requests."""
    data = arguments.get("data", "")
    analysis_type = arguments.get("analysis_type", "alignment")
    
    # Convert data to DataFrame if needed
    import pandas as pd
    if isinstance(data, str):
        # Assume CSV format for now
        try:
            df = pd.read_csv(data)
        except:
            # Try to parse as FASTA
            df = parse_fasta_to_dataframe(data)
    else:
        df = data
    
    if analysis_type == "alignment":
        result = align_and_visualize_fasta(df)
    else:
        result = analyze_basic_stats(df)
    
    return {
        "status": "success",
        "tool": "analyze_sequence_data",
        "analysis_type": analysis_type,
        "result": str(result)
    }

async def handle_visualize_alignment(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Handle alignment visualization requests."""
    alignment_file = arguments.get("alignment_file", "")
    output_format = arguments.get("output_format", "png")
    
    # Use the existing bio tool for visualization
    result = align_and_visualize_fasta(None)  # Will use default file
    
    return {
        "status": "success",
        "tool": "visualize_alignment",
        "output_format": output_format,
        "result": str(result)
    }

async def handle_dna_vendor_research(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Handle DNA vendor research requests."""
    command = arguments.get("command", "")
    sequence_length = arguments.get("sequence_length", None)
    quantity = arguments.get("quantity", "standard")
    
    # Use the DNA vendor research tool
    result = run_dna_vendor_research_raw(command, sequence_length, quantity)
    
    return {
        "status": "success",
        "tool": "dna_vendor_research",
        "command": command,
        "sequence_length": sequence_length,
        "quantity": quantity,
        "result": result
    }

def parse_fasta_to_dataframe(fasta_content: str) -> pd.DataFrame:
    """Parse FASTA content to DataFrame."""
    import pandas as pd
    
    sequences = []
    names = []
    current_seq = ""
    current_name = ""
    
    for line in fasta_content.split('\n'):
        if line.startswith('>'):
            if current_name and current_seq:
                names.append(current_name)
                sequences.append(current_seq)
            current_name = line[1:].strip()
            current_seq = ""
        else:
            current_seq += line.strip()
    
    if current_name and current_seq:
        names.append(current_name)
        sequences.append(current_seq)
    
    return pd.DataFrame({
        'name': names,
        'sequence': sequences
    })

async def main():
    """Main entry point for the MCP server."""
    async with stdio_server() as (read_stream, write_stream):
        await server.run(
            read_stream,
            write_stream,
            InitializationOptions(
                server_name="bioinformatics-mcp-server",
                server_version="1.0.0",
                capabilities=server.get_capabilities(
                    notification_options=None,
                    experimental_capabilities=None,
                ),
            ),
        )

if __name__ == "__main__":
    asyncio.run(main()) 