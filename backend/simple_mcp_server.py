#!/usr/bin/env python3
"""
Simplified MCP Server for Bioinformatics Tasks
This version implements a basic MCP server without the official SDK.
"""

import asyncio
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Add the tools directory to the path (tools are in the root directory)
tools_path = str((Path(__file__).resolve().parent.parent / "tools").resolve())
sys.path.insert(0, tools_path)
logger.info(f"🔧 Added tools path: {tools_path}")

try:
    # Import tools using direct module imports
    import alignment
    import bio
    import mutations
    import data_science
    
    # Get the functions from the modules
    run_alignment = alignment.run_alignment
    align_and_visualize_fasta = bio.align_and_visualize_fasta
    # Replace all uses of run_mutation with run_mutation_raw for direct calls
    # import mutations and use mutations.run_mutation_raw(...)
    analyze_basic_stats = data_science.analyze_basic_stats
    
    logger.info("✅ Bioinformatics tools imported successfully")
except ImportError as e:
    logger.error(f"❌ Failed to import bioinformatics tools: {e}")
    logger.error(f"Current sys.path: {sys.path}")
    sys.exit(1)

class SimpleMCPServer:
    """Simple MCP server implementation without the official SDK."""
    
    def __init__(self):
        self.server_name = "bioinformatics-mcp-server"
        self.server_version = "1.0.0"
        self.tools = self._get_tools()
    
    def _get_tools(self) -> List[Dict[str, Any]]:
        """Define available tools."""
        return [
            {
                "name": "sequence_alignment",
                "description": "Perform multiple sequence alignment on a set of DNA/RNA sequences",
                "inputSchema": {
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
            },
            {
                "name": "mutate_sequence",
                "description": "Generate mutations and variants of a given sequence",
                "inputSchema": {
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
            },
            {
                "name": "analyze_sequence_data",
                "description": "Analyze sequence data and generate visualizations",
                "inputSchema": {
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
            },
            {
                "name": "select_variants",
                "description": "Select variants from previous mutation results based on criteria",
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "session_id": {
                            "type": "string",
                            "description": "Session ID to get previous results from"
                        },
                        "selection_criteria": {
                            "type": "string",
                            "enum": ["diversity", "random", "length", "custom"],
                            "default": "diversity",
                            "description": "Criteria for selecting variants"
                        },
                        "num_variants": {
                            "type": "integer",
                            "default": 10,
                            "description": "Number of variants to select"
                        },
                        "custom_filters": {
                            "type": "object",
                            "description": "Custom filtering criteria (optional)"
                        }
                    },
                    "required": ["session_id"]
                }
            },
            {
                "name": "parse_command",
                "description": "Parse natural language commands into structured operations",
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "command": {
                            "type": "string",
                            "description": "Natural language command to parse"
                        },
                        "session_id": {
                            "type": "string",
                            "description": "Optional session ID for context"
                        }
                    },
                    "required": ["command"]
                }
            },
            {
                "name": "execute_command",
                "description": "Execute parsed commands using appropriate tools",
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "parsed_command": {
                            "type": "object",
                            "description": "Parsed command object from parse_command"
                        }
                    },
                    "required": ["parsed_command"]
                }
            },
            {
                "name": "handle_natural_command",
                "description": "Handle natural language commands by parsing and executing them",
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "command": {
                            "type": "string",
                            "description": "Natural language command to handle"
                        },
                        "session_id": {
                            "type": "string",
                            "description": "Optional session ID for context"
                        }
                    },
                    "required": ["command"]
                }
            },
            {
                "name": "visualize_alignment",
                "description": "Create visual representation of sequence alignments",
                "inputSchema": {
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
            }
        ]
    
    def convert_sequences_to_fasta(self, sequences_input: str) -> str:
        """Convert various sequence input formats to FASTA format."""
        print(f"🔍 convert_sequences_to_fasta input: '{sequences_input}'")
        print(f"🔍 Input type: {type(sequences_input)}")
        print(f"🔍 Input length: {len(sequences_input)}")
        print(f"🔍 Input repr: {repr(sequences_input)}")
        
        if not sequences_input or not sequences_input.strip():
            return ""
        
        # If it's already in FASTA format (contains '>'), return as is
        if '>' in sequences_input:
            print(f"🔍 Input already in FASTA format, returning as-is")
            return sequences_input
        
        # Check if this looks like a command string (e.g., "align sequences ACTGTTGAC ACTGCATCC")
        # Extract only the sequences part
        import re
        
        # Try to extract sequences from command format
        # Pattern: "align sequences" followed by sequences
        command_pattern = r'align\s+sequences\s+(.+)'
        match = re.match(command_pattern, sequences_input.strip(), re.IGNORECASE)
        
        if match:
            # Extract just the sequences part
            sequences_part = match.group(1).strip()
            print(f"🔍 Extracted sequences part: '{sequences_part}'")
            sequences_input = sequences_part
        
        # Split by spaces, commas, or newlines
        sequences = re.split(r'[,\s\n]+', sequences_input.strip())
        sequences = [seq.strip() for seq in sequences if seq.strip()]
        
        print(f"🔍 After splitting: {sequences}")
        
        # Convert to FASTA format
        fasta_output = ""
        for i, seq in enumerate(sequences, 1):
            fasta_output += f">sequence_{i}\n{seq}\n"
        
        print(f"🔍 Final FASTA output: '{fasta_output}'")
        return fasta_output

    async def handle_sequence_alignment(self, arguments: Dict[str, Any]) -> Dict[str, Any]:
        """Handle sequence alignment requests."""
        logger.info(f"🔍 Received arguments: {arguments}")
        logger.info(f"🔍 Arguments type: {type(arguments)}")
        logger.info(f"🔍 Arguments keys: {list(arguments.keys()) if isinstance(arguments, dict) else 'Not a dict'}")
        
        sequences = arguments.get("sequences", "")
        algorithm = arguments.get("algorithm", "clustal")
        
        logger.info(f"🔍 Extracted sequences: '{sequences}' (type: {type(sequences)}, length: {len(sequences) if sequences else 0})")
        logger.info(f"🔍 Extracted algorithm: '{algorithm}' (type: {type(algorithm)})")
        
        # Convert to FASTA format if needed
        fasta_sequences = self.convert_sequences_to_fasta(sequences)
        logger.info(f"🔍 Converted to FASTA: '{fasta_sequences}'")
        
        logger.info(f"Aligning sequences using {algorithm} algorithm")
        
        try:
            # Use the existing alignment tool
            result = run_alignment(fasta_sequences)
            logger.info(f"✅ Alignment tool executed successfully")
        except Exception as e:
            logger.error(f"❌ Error in alignment tool: {e}")
            result = {
                "text": f"Error performing alignment: {str(e)}",
                "alignment": []
            }
        
        return {
            "status": "success",
            "tool": "sequence_alignment",
            "algorithm": algorithm,
            "result": result
        }

    async def handle_mutate_sequence(self, arguments: Dict[str, Any]) -> Dict[str, Any]:
        """Handle sequence mutation requests."""
        logger.info(f"🔍 Received arguments: {arguments}")
        logger.info(f"🔍 Arguments type: {type(arguments)}")
        logger.info(f"🔍 Arguments keys: {list(arguments.keys()) if isinstance(arguments, dict) else 'Not a dict'}")
        
        sequence = arguments.get("sequence", "")
        num_variants = arguments.get("num_variants", 96)
        mutation_rate = arguments.get("mutation_rate", 0.1)
        
        logger.info(f"🔍 Extracted sequence: '{sequence}' (type: {type(sequence)}, length: {len(sequence) if sequence else 0})")
        logger.info(f"🔍 Extracted num_variants: {num_variants} (type: {type(num_variants)})")
        logger.info(f"🔍 Extracted mutation_rate: {mutation_rate} (type: {type(mutation_rate)})")
        
        # Ensure proper type conversion
        try:
            num_variants = int(num_variants) if num_variants is not None else 96
        except (ValueError, TypeError):
            num_variants = 96
            logger.warning(f"Invalid num_variants value, using default: {num_variants}")
        
        try:
            mutation_rate = float(mutation_rate) if mutation_rate is not None else 0.1
        except (ValueError, TypeError):
            mutation_rate = 0.1
            logger.warning(f"Invalid mutation_rate value, using default: {mutation_rate}")
        
        logger.info(f"🔍 Final values - sequence: '{sequence}', num_variants: {num_variants}, mutation_rate: {mutation_rate}")
        logger.info(f"Generating {num_variants} variants with mutation rate {mutation_rate}")
        
        try:
            # Use the existing mutation tool
            result = mutations.run_mutation_raw(sequence, num_variants)
            logger.info(f"✅ Mutation tool executed successfully")
        except Exception as e:
            logger.error(f"❌ Error in mutation tool: {e}")
            result = {
                "text": f"Error generating mutations: {str(e)}",
                "plot": {
                    "data": [{"x": [1], "y": [1], "type": "bar"}],
                    "layout": {"title": "Error"}
                }
            }
        
        return {
            "status": "success",
            "tool": "mutate_sequence",
            "input_sequence": sequence,
            "num_variants": num_variants,
            "mutation_rate": mutation_rate,
            "result": result
        }

    async def handle_analyze_sequence_data(self, arguments: Dict[str, Any]) -> Dict[str, Any]:
        """Handle sequence data analysis requests."""
        data = arguments.get("data", "")
        analysis_type = arguments.get("analysis_type", "alignment")
        
        logger.info(f"Analyzing sequence data with type: {analysis_type}")
        
        # Convert data to DataFrame if needed
        import pandas as pd
        if isinstance(data, str):
            # Assume CSV format for now
            try:
                df = pd.read_csv(data)
            except:
                # Try to parse as FASTA
                df = self.parse_fasta_to_dataframe(data)
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

    async def handle_visualize_alignment(self, arguments: Dict[str, Any]) -> Dict[str, Any]:
        """Handle alignment visualization requests."""
        alignment_file = arguments.get("alignment_file", "")
        output_format = arguments.get("output_format", "png")
        
        logger.info(f"Creating visualization in {output_format} format")
        
        # Use the existing bio tool for visualization
        result = align_and_visualize_fasta(None)  # Will use default file
        
        return {
            "status": "success",
            "tool": "visualize_alignment",
            "output_format": output_format,
            "result": str(result)
        }

    async def handle_select_variants(self, arguments: Dict[str, Any]) -> Dict[str, Any]:
        """Handle variant selection requests."""
        session_id = arguments.get("session_id", "")
        selection_criteria = arguments.get("selection_criteria", "diversity")
        num_variants = arguments.get("num_variants", 10)
        custom_filters = arguments.get("custom_filters", None)
        
        logger.info(f"Selecting variants from session {session_id} using {selection_criteria} criteria")
        
        try:
            # Import variant selection tool
            import variant_selection
            
            result = variant_selection.run_variant_selection_raw(
                session_id, selection_criteria, num_variants, custom_filters
            )
            
            return {
                "status": "success",
                "tool": "select_variants",
                "session_id": session_id,
                "selection_criteria": selection_criteria,
                "result": result
            }
        except Exception as e:
            logger.error(f"❌ Error in variant selection: {e}")
            return {
                "status": "error",
                "tool": "select_variants",
                "error": str(e)
            }

    async def handle_parse_command(self, arguments: Dict[str, Any]) -> Dict[str, Any]:
        """Handle command parsing requests."""
        command = arguments.get("command", "")
        session_id = arguments.get("session_id", None)
        
        logger.info(f"Parsing command: '{command}'")
        
        try:
            # Import command parser tool
            import command_parser
            
            result = command_parser.parse_command_raw(command, session_id)
            
            return {
                "status": "success",
                "tool": "parse_command",
                "parsed_command": result,
                "original_command": command
            }
        except Exception as e:
            logger.error(f"❌ Error in command parsing: {e}")
            return {
                "status": "error",
                "tool": "parse_command",
                "error": str(e)
            }

    async def handle_execute_command(self, arguments: Dict[str, Any]) -> Dict[str, Any]:
        """Handle command execution requests."""
        parsed_command = arguments.get("parsed_command", {})
        
        logger.info(f"Executing parsed command: {parsed_command}")
        
        try:
            # Import command executor tool
            import command_executor
            
            result = command_executor.execute_command_raw(parsed_command)
            
            return {
                "status": "success",
                "tool": "execute_command",
                "execution_result": result
            }
        except Exception as e:
            logger.error(f"❌ Error in command execution: {e}")
            return {
                "status": "error",
                "tool": "execute_command",
                "error": str(e)
            }

    async def handle_natural_command(self, arguments: Dict[str, Any]) -> Dict[str, Any]:
        """Handle natural language command requests."""
        command = arguments.get("command", "")
        session_id = arguments.get("session_id", None)
        
        logger.info(f"Handling natural command: '{command}'")
        
        try:
            # Import command handler tool
            import command_handler
            
            result = command_handler.handle_command_raw(command, session_id)
            
            return {
                "status": "success",
                "tool": "handle_natural_command",
                "command_result": result,
                "original_command": command
            }
        except Exception as e:
            logger.error(f"❌ Error in natural command handling: {e}")
            return {
                "status": "error",
                "tool": "handle_natural_command",
                "error": str(e)
            }

    def parse_fasta_to_dataframe(self, fasta_content: str):
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

    async def handle_call_tool(self, name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
        """Handle tool calls for bioinformatics operations."""
        logger.info(f"🔧 Calling tool: {name} with arguments: {arguments}")
        
        try:
            if name == "sequence_alignment":
                result = await self.handle_sequence_alignment(arguments)
            elif name == "mutate_sequence":
                result = await self.handle_mutate_sequence(arguments)
            elif name == "analyze_sequence_data":
                result = await self.handle_analyze_sequence_data(arguments)
            elif name == "select_variants":
                result = await self.handle_select_variants(arguments)
            elif name == "parse_command":
                result = await self.handle_parse_command(arguments)
            elif name == "execute_command":
                result = await self.handle_execute_command(arguments)
            elif name == "handle_natural_command":
                result = await self.handle_natural_command(arguments)
            elif name == "visualize_alignment":
                result = await self.handle_visualize_alignment(arguments)
            else:
                raise ValueError(f"Unknown tool: {name}")
            
            logger.info(f"✅ Tool {name} completed successfully")
            return result
        
        except Exception as e:
            logger.error(f"❌ Error in tool call {name}: {e}")
            return {
                "status": "error",
                "tool": name,
                "error": str(e)
            }

    async def handle_list_tools(self) -> Dict[str, Any]:
        """List all available bioinformatics tools."""
        logger.info("📋 Listing available tools...")
        logger.info(f"✅ Listed {len(self.tools)} tools")
        return {
            "status": "success",
            "tools": self.tools
        }

    async def handle_initialize(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Handle initialization request."""
        logger.info("🚀 Initializing MCP server...")
        return {
            "protocolVersion": "2024-11-05",
            "capabilities": {
                "tools": {}
            },
            "serverInfo": {
                "name": self.server_name,
                "version": self.server_version
            }
        }

    async def process_message(self, message: Dict[str, Any]) -> Dict[str, Any]:
        """Process incoming MCP messages."""
        method = message.get("method")
        params = message.get("params", {})
        request_id = message.get("id")
        
        logger.info(f"📨 Processing message: {method}")
        
        try:
            if method == "initialize":
                result = await self.handle_initialize(params)
            elif method == "tools/list":
                result = await self.handle_list_tools()
            elif method == "tools/call":
                name = params.get("name")
                arguments = params.get("arguments", {})
                result = await self.handle_call_tool(name, arguments)
            else:
                result = {"error": f"Unknown method: {method}"}
            
            return {
                "jsonrpc": "2.0",
                "id": request_id,
                "result": result
            }
        
        except Exception as e:
            logger.error(f"❌ Error processing message: {e}")
            return {
                "jsonrpc": "2.0",
                "id": request_id,
                "error": {
                    "code": -32603,
                    "message": str(e)
                }
            }

    async def run(self):
        """Run the MCP server using stdio."""
        logger.info("🚀 Starting Bioinformatics MCP Server...")
        logger.info("Available tools:")
        logger.info("- sequence_alignment: Perform multiple sequence alignment")
        logger.info("- mutate_sequence: Generate sequence mutations")
        logger.info("- analyze_sequence_data: Analyze sequence data")
        logger.info("- visualize_alignment: Create alignment visualizations")
        
        try:
            # Read from stdin, write to stdout
            while True:
                try:
                    # Read a line from stdin
                    line = await asyncio.get_event_loop().run_in_executor(None, sys.stdin.readline)
                    if not line:
                        break
                    
                    # Parse the JSON message
                    message = json.loads(line.strip())
                    
                    # Process the message
                    response = await self.process_message(message)
                    
                    # Send the response
                    response_line = json.dumps(response) + "\n"
                    await asyncio.get_event_loop().run_in_executor(None, sys.stdout.write, response_line)
                    await asyncio.get_event_loop().run_in_executor(None, sys.stdout.flush)
                    
                except json.JSONDecodeError as e:
                    logger.error(f"❌ Invalid JSON: {e}")
                except Exception as e:
                    logger.error(f"❌ Error processing message: {e}")
                    
        except KeyboardInterrupt:
            logger.info("🛑 MCP Server stopped by user")
        except Exception as e:
            logger.error(f"❌ Error running MCP server: {e}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            raise

async def main():
    """Main entry point for the MCP server."""
    server = SimpleMCPServer()
    await server.run()

if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        logger.info("🛑 MCP Server stopped by user")
    except Exception as e:
        logger.error(f"❌ Fatal error: {e}")
        import traceback
        logger.error(f"Traceback: {traceback.format_exc()}")
        sys.exit(1) 