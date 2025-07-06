#!/usr/bin/env python3
"""
Enhanced MCP Server for Bioinformatics Tasks
Provides tools for sequence alignment, mutation analysis, and other bioinformatics operations.
"""

import asyncio
import json
import logging
import sys
import re
import base64
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence
from datetime import datetime

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

# Add the tools directory to the path
sys.path.append(str(Path(__file__).parent / "tools"))

from tools.alignment import run_alignment
from tools.bio import align_and_visualize_fasta
from tools.mutations import mutate_sequence
from tools.data_science import analyze_data

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Initialize the MCP server
server = Server("enhanced-bioinformatics-mcp-server")

class ValidationError(Exception):
    """Custom exception for validation errors."""
    pass

def validate_sequence(sequence: str) -> bool:
    """Validate DNA/RNA sequence."""
    if not sequence:
        return False
    
    # Remove whitespace and convert to uppercase
    sequence = re.sub(r'\s+', '', sequence.upper())
    
    # Check for valid DNA/RNA characters
    valid_chars = set('ATCGN-')
    return all(char in valid_chars for char in sequence)

def validate_fasta_format(fasta_content: str) -> bool:
    """Validate FASTA format content."""
    if not fasta_content:
        return False
    
    lines = fasta_content.strip().split('\n')
    if not lines:
        return False
    
    # Check if first line starts with '>'
    if not lines[0].startswith('>'):
        return False
    
    # Check for at least one sequence
    has_sequence = False
    for line in lines[1:]:
        if line.startswith('>'):
            has_sequence = False
        elif line.strip():
            has_sequence = True
            if not validate_sequence(line.strip()):
                return False
    
    return has_sequence

def parse_fasta_to_dataframe(fasta_content: str):
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
                            "minimum": 1,
                            "maximum": 1000,
                            "description": "Number of variants to generate"
                        },
                        "mutation_rate": {
                            "type": "number",
                            "default": 0.1,
                            "minimum": 0.0,
                            "maximum": 1.0,
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
                name="sequence_validation",
                description="Validate DNA/RNA sequences and FASTA format",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "sequence": {
                            "type": "string",
                            "description": "DNA/RNA sequence to validate"
                        },
                        "fasta_content": {
                            "type": "string",
                            "description": "FASTA format content to validate"
                        }
                    }
                }
            ),
            Tool(
                name="sequence_statistics",
                description="Calculate statistics for DNA/RNA sequences",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "sequence": {
                            "type": "string",
                            "description": "DNA/RNA sequence to analyze"
                        },
                        "include_composition": {
                            "type": "boolean",
                            "default": True,
                            "description": "Include nucleotide composition analysis"
                        }
                    },
                    "required": ["sequence"]
                }
            ),
            Tool(
                name="reverse_complement",
                description="Generate reverse complement of DNA sequence",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "sequence": {
                            "type": "string",
                            "description": "DNA sequence to reverse complement"
                        }
                    },
                    "required": ["sequence"]
                }
            )
        ]
    )

@server.call_tool()
async def handle_call_tool(name: str, arguments: Dict[str, Any]) -> CallToolResult:
    """Handle tool calls for bioinformatics operations."""
    
    start_time = datetime.now()
    
    try:
        # Validate input arguments based on tool type
        if name == "sequence_alignment":
            result = await handle_sequence_alignment(arguments)
        elif name == "mutate_sequence":
            result = await handle_mutate_sequence(arguments)
        elif name == "analyze_sequence_data":
            result = await handle_analyze_sequence_data(arguments)
        elif name == "visualize_alignment":
            result = await handle_visualize_alignment(arguments)
        elif name == "sequence_validation":
            result = await handle_sequence_validation(arguments)
        elif name == "sequence_statistics":
            result = await handle_sequence_statistics(arguments)
        elif name == "reverse_complement":
            result = await handle_reverse_complement(arguments)
        else:
            raise ValueError(f"Unknown tool: {name}")
        
        # Add timing information
        execution_time = (datetime.now() - start_time).total_seconds()
        result["execution_time"] = execution_time
        result["timestamp"] = datetime.now().isoformat()
        
        return CallToolResult(
            content=[
                TextContent(
                    type="text",
                    text=json.dumps(result, indent=2)
                )
            ]
        )
    
    except ValidationError as e:
        logger.error(f"Validation error in tool call {name}: {e}")
        return CallToolResult(
            content=[
                TextContent(
                    type="text",
                    text=json.dumps({
                        "error": "validation_error",
                        "message": str(e),
                        "tool": name,
                        "timestamp": datetime.now().isoformat()
                    }, indent=2)
                )
            ]
        )
    except Exception as e:
        logger.error(f"Error in tool call {name}: {e}")
        return CallToolResult(
            content=[
                TextContent(
                    type="text",
                    text=json.dumps({
                        "error": "execution_error",
                        "message": str(e),
                        "tool": name,
                        "timestamp": datetime.now().isoformat()
                    }, indent=2)
                )
            ]
        )

async def handle_sequence_alignment(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Handle sequence alignment requests."""
    sequences = arguments.get("sequences", "")
    algorithm = arguments.get("algorithm", "clustal")
    
    # Validate input
    if not sequences:
        raise ValidationError("Sequences parameter is required")
    
    if not validate_fasta_format(sequences):
        raise ValidationError("Invalid FASTA format")
    
    # Use the existing alignment tool
    result = run_alignment(sequences)
    
    return {
        "status": "success",
        "tool": "sequence_alignment",
        "algorithm": algorithm,
        "input_sequences_count": len([line for line in sequences.split('\n') if line.startswith('>')]),
        "result": result
    }

async def handle_mutate_sequence(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Handle sequence mutation requests."""
    sequence = arguments.get("sequence", "")
    num_variants = arguments.get("num_variants", 96)
    mutation_rate = arguments.get("mutation_rate", 0.1)
    
    # Validate input
    if not sequence:
        raise ValidationError("Sequence parameter is required")
    
    if not validate_sequence(sequence):
        raise ValidationError("Invalid DNA/RNA sequence")
    
    if num_variants < 1 or num_variants > 1000:
        raise ValidationError("num_variants must be between 1 and 1000")
    
    if mutation_rate < 0.0 or mutation_rate > 1.0:
        raise ValidationError("mutation_rate must be between 0.0 and 1.0")
    
    # Use the existing mutation tool
    result = mutate_sequence(sequence, num_variants)
    
    return {
        "status": "success",
        "tool": "mutate_sequence",
        "input_sequence": sequence,
        "input_sequence_length": len(sequence),
        "num_variants": num_variants,
        "mutation_rate": mutation_rate,
        "result": result
    }

async def handle_analyze_sequence_data(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Handle sequence data analysis requests."""
    data = arguments.get("data", "")
    analysis_type = arguments.get("analysis_type", "alignment")
    
    # Validate input
    if not data:
        raise ValidationError("Data parameter is required")
    
    # Convert data to DataFrame if needed
    import pandas as pd
    if isinstance(data, str):
        # Assume CSV format for now
        try:
            df = pd.read_csv(data)
        except:
            # Try to parse as FASTA
            if validate_fasta_format(data):
                df = parse_fasta_to_dataframe(data)
            else:
                raise ValidationError("Invalid data format. Expected CSV or FASTA format.")
    else:
        df = data
    
    if analysis_type == "alignment":
        result = align_and_visualize_fasta(df)
    else:
        result = analyze_data(df)
    
    return {
        "status": "success",
        "tool": "analyze_sequence_data",
        "analysis_type": analysis_type,
        "input_data_rows": len(df),
        "result": str(result)
    }

async def handle_visualize_alignment(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Handle alignment visualization requests."""
    alignment_file = arguments.get("alignment_file", "")
    output_format = arguments.get("output_format", "png")
    
    # Validate input
    if not alignment_file:
        raise ValidationError("alignment_file parameter is required")
    
    if output_format not in ["png", "svg", "pdf"]:
        raise ValidationError("output_format must be one of: png, svg, pdf")
    
    # Use the existing bio tool for visualization
    result = align_and_visualize_fasta(None)  # Will use default file
    
    return {
        "status": "success",
        "tool": "visualize_alignment",
        "alignment_file": alignment_file,
        "output_format": output_format,
        "result": str(result)
    }

async def handle_sequence_validation(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Handle sequence validation requests."""
    sequence = arguments.get("sequence", "")
    fasta_content = arguments.get("fasta_content", "")
    
    result = {
        "status": "success",
        "tool": "sequence_validation",
        "validations": {}
    }
    
    if sequence:
        is_valid = validate_sequence(sequence)
        result["validations"]["sequence"] = {
            "is_valid": is_valid,
            "length": len(sequence) if sequence else 0,
            "sequence": sequence[:50] + "..." if len(sequence) > 50 else sequence
        }
    
    if fasta_content:
        is_valid = validate_fasta_format(fasta_content)
        sequence_count = len([line for line in fasta_content.split('\n') if line.startswith('>')])
        result["validations"]["fasta"] = {
            "is_valid": is_valid,
            "sequence_count": sequence_count,
            "content_length": len(fasta_content)
        }
    
    return result

async def handle_sequence_statistics(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Handle sequence statistics requests."""
    sequence = arguments.get("sequence", "")
    include_composition = arguments.get("include_composition", True)
    
    # Validate input
    if not sequence:
        raise ValidationError("Sequence parameter is required")
    
    if not validate_sequence(sequence):
        raise ValidationError("Invalid DNA/RNA sequence")
    
    # Clean sequence
    clean_sequence = re.sub(r'\s+', '', sequence.upper())
    
    # Calculate basic statistics
    stats = {
        "length": len(clean_sequence),
        "gc_content": 0.0,
        "at_content": 0.0,
        "n_content": 0.0,
        "gap_content": 0.0
    }
    
    if include_composition:
        # Calculate nucleotide composition
        gc_count = clean_sequence.count('G') + clean_sequence.count('C')
        at_count = clean_sequence.count('A') + clean_sequence.count('T')
        n_count = clean_sequence.count('N')
        gap_count = clean_sequence.count('-')
        
        total = len(clean_sequence)
        if total > 0:
            stats["gc_content"] = gc_count / total
            stats["at_content"] = at_count / total
            stats["n_content"] = n_count / total
            stats["gap_content"] = gap_count / total
        
        # Add detailed composition
        stats["composition"] = {
            "A": clean_sequence.count('A'),
            "T": clean_sequence.count('T'),
            "C": clean_sequence.count('C'),
            "G": clean_sequence.count('G'),
            "N": clean_sequence.count('N'),
            "Gap": clean_sequence.count('-')
        }
    
    return {
        "status": "success",
        "tool": "sequence_statistics",
        "sequence_length": len(clean_sequence),
        "statistics": stats
    }

async def handle_reverse_complement(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Handle reverse complement requests."""
    sequence = arguments.get("sequence", "")
    
    # Validate input
    if not sequence:
        raise ValidationError("Sequence parameter is required")
    
    if not validate_sequence(sequence):
        raise ValidationError("Invalid DNA/RNA sequence")
    
    # Clean sequence
    clean_sequence = re.sub(r'\s+', '', sequence.upper())
    
    # Create reverse complement
    complement_map = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '-': '-'
    }
    
    complement = ''.join(complement_map.get(base, base) for base in clean_sequence)
    reverse_complement = complement[::-1]
    
    return {
        "status": "success",
        "tool": "reverse_complement",
        "original_sequence": clean_sequence,
        "reverse_complement": reverse_complement,
        "sequence_length": len(clean_sequence)
    }

async def main():
    """Main entry point for the MCP server."""
    async with stdio_server() as (read_stream, write_stream):
        await server.run(
            read_stream,
            write_stream,
            InitializationOptions(
                server_name="enhanced-bioinformatics-mcp-server",
                server_version="1.1.0",
                capabilities=server.get_capabilities(
                    notification_options=None,
                    experimental_capabilities=None,
                ),
            ),
        )

if __name__ == "__main__":
    asyncio.run(main()) 