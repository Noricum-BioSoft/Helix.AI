import logging
import sys
from pathlib import Path
from typing import Dict, Any, Optional
import asyncio

logger = logging.getLogger(__name__)

class CommandExecutor:
    """Executes parsed commands using the appropriate tools."""
    
    def __init__(self):
        # Add tools directory to path
        tools_path = str((Path(__file__).parent).resolve())
        sys.path.insert(0, tools_path)
        
        # Add backend directory to path for history manager
        backend_path = str((Path(__file__).parent.parent / "backend").resolve())
        sys.path.insert(0, backend_path)
        
        try:
            from history_manager import history_manager
            self.history_manager = history_manager
        except ImportError:
            logger.warning("History manager not available")
            self.history_manager = None
    
    async def execute_command(self, parsed_command: Dict[str, Any]) -> Dict[str, Any]:
        """Execute a parsed command using the appropriate tool."""
        action = parsed_command.get("action")
        tool = parsed_command.get("tool")
        parameters = parsed_command.get("parameters", {})
        session_id = parsed_command.get("session_id")
        original_command = parsed_command.get("original_command", "")
        
        logger.info(f"Executing command: {action} using tool: {tool}")
        logger.info(f"Parameters: {parameters}")
        
        try:
            if tool == "select_variants":
                result = await self._execute_select_variants(parameters, session_id)
            elif tool == "mutate_sequence":
                result = await self._execute_mutate_sequence(parameters, session_id)
            elif tool == "sequence_alignment":
                result = await self._execute_sequence_alignment(parameters, session_id)
            elif tool == "analyze_sequence_data":
                result = await self._execute_analyze_data(parameters, session_id)
            elif tool == "dna_vendor_research":
                result = await self._execute_dna_vendor_research(parameters, session_id)
            elif tool == "read_trimming":
                result = await self._execute_read_trimming(parameters, session_id)
            elif tool == "read_merging":
                result = await self._execute_read_merging(parameters, session_id)
            elif tool == "quality_assessment":
                result = await self._execute_quality_assessment(parameters, session_id)
            else:
                result = {
                    "status": "error",
                    "message": f"Unknown tool: {tool}",
                    "action": action,
                    "tool": tool
                }
            
            # Track the operation in history if session_id is provided
            if session_id and self.history_manager:
                self.history_manager.add_history_entry(
                    session_id,
                    original_command,
                    tool,
                    result,
                    parameters
                )
            
            return {
                "status": "success",
                "action": action,
                "tool": tool,
                "result": result,
                "session_id": session_id,
                "confidence": parsed_command.get("confidence", 0.0)
            }
            
        except Exception as e:
            logger.error(f"Error executing command: {e}")
            return {
                "status": "error",
                "action": action,
                "tool": tool,
                "error": str(e),
                "session_id": session_id
            }
    
    async def _execute_select_variants(self, parameters: Dict[str, Any], session_id: Optional[str]) -> Dict[str, Any]:
        """Execute variant selection command."""
        try:
            import variant_selection
            
            num_variants = parameters.get("num_variants", 10)
            selection_criteria = parameters.get("selection_criteria", "random")
            custom_filters = parameters.get("custom_filters")
            
            if not session_id:
                return {
                    "status": "error",
                    "message": "Session ID required for variant selection"
                }
            
            result = variant_selection.run_variant_selection_raw(
                session_id=session_id,
                selection_criteria=selection_criteria,
                num_variants=num_variants,
                custom_filters=custom_filters
            )
            
            return result
            
        except Exception as e:
            logger.error(f"Error in variant selection: {e}")
            return {
                "status": "error",
                "message": f"Error selecting variants: {str(e)}"
            }
    
    async def _execute_mutate_sequence(self, parameters: Dict[str, Any], session_id: Optional[str]) -> Dict[str, Any]:
        """Execute sequence mutation command."""
        try:
            import mutations
            
            sequence = parameters.get("sequence", "")
            num_variants = parameters.get("num_variants", 96)
            mutation_rate = parameters.get("mutation_rate", 0.1)
            
            if not sequence:
                return {
                    "status": "error",
                    "message": "Sequence required for mutation"
                }
            
            result = mutations.run_mutation_raw(sequence, num_variants)
            
            return result
            
        except Exception as e:
            logger.error(f"Error in sequence mutation: {e}")
            return {
                "status": "error",
                "message": f"Error mutating sequence: {str(e)}"
            }
    
    async def _execute_sequence_alignment(self, parameters: Dict[str, Any], session_id: Optional[str]) -> Dict[str, Any]:
        """Execute sequence alignment command."""
        try:
            import alignment
            
            sequences = parameters.get("sequences", "")
            algorithm = parameters.get("algorithm", "clustal")
            
            if not sequences:
                return {
                    "status": "error",
                    "message": "Sequences required for alignment"
                }
            
            result = alignment.run_alignment(sequences)
            
            return result
            
        except Exception as e:
            logger.error(f"Error in sequence alignment: {e}")
            return {
                "status": "error",
                "message": f"Error aligning sequences: {str(e)}"
            }
    
    async def _execute_analyze_data(self, parameters: Dict[str, Any], session_id: Optional[str]) -> Dict[str, Any]:
        """Execute data analysis command."""
        try:
            import data_science
            
            data = parameters.get("data", "")
            analysis_type = parameters.get("analysis_type", "basic")
            
            if not data:
                return {
                    "status": "error",
                    "message": "Data required for analysis"
                }
            
            result = data_science.analyze_basic_stats(data)
            
            return {
                "status": "success",
                "result": str(result)
            }
            
        except Exception as e:
            logger.error(f"Error in data analysis: {e}")
            return {
                "status": "error",
                "message": f"Error analyzing data: {str(e)}"
            }

    async def _execute_dna_vendor_research(self, parameters: Dict[str, Any], session_id: Optional[str]) -> Dict[str, Any]:
        """Execute DNA vendor research command."""
        try:
            import dna_vendor_research
            
            command = parameters.get("command", "")
            sequence_length = parameters.get("sequence_length")
            quantity = parameters.get("quantity", "standard")
            
            if not command:
                return {
                    "status": "error",
                    "message": "Command required for DNA vendor research"
                }
            
            # Convert sequence_length to int if it's not None
            if sequence_length is not None:
                try:
                    sequence_length = int(sequence_length)
                except (ValueError, TypeError):
                    sequence_length = None
            
            result = dna_vendor_research.run_dna_vendor_research_raw(command, sequence_length, quantity)
            
            return result
            
        except Exception as e:
            logger.error(f"Error in DNA vendor research: {e}")
            return {
                "status": "error",
                "message": f"Error researching DNA vendors: {str(e)}"
            }
    
    async def _execute_read_trimming(self, parameters: Dict[str, Any], session_id: Optional[str]) -> Dict[str, Any]:
        """Execute read trimming command."""
        try:
            import read_trimming
            
            file_type = parameters.get("file_type", "single")
            quality_threshold = parameters.get("quality_threshold", 20)
            adapter = parameters.get("adapter")
            
            # Get initial values from parameters
            forward_reads = parameters.get("forward_reads", "")
            reverse_reads = parameters.get("reverse_reads", "")
            reads = parameters.get("reads", "")
            
            print(f"🔧 [EXECUTOR] Initial values - reads: {len(reads)} chars, forward: {len(forward_reads)} chars, reverse: {len(reverse_reads)} chars, adapter: {adapter}")
            
            # If reads are missing but we have an adapter and session_id, try to get trimmed reads from history
            # Check if reads are actually empty or don't look like FASTQ (FASTQ starts with "@")
            reads_empty = not reads or len(reads.strip()) == 0 or not reads.strip().startswith("@")
            forward_empty = not forward_reads or len(forward_reads.strip()) == 0 or not forward_reads.strip().startswith("@")
            reverse_empty = not reverse_reads or len(reverse_reads.strip()) == 0 or not reverse_reads.strip().startswith("@")
            
            print(f"🔧 [EXECUTOR] Empty checks - reads_empty: {reads_empty}, forward_empty: {forward_empty}, reverse_empty: {reverse_empty}")
            
            # If adapter is specified and reads are empty/don't look like FASTQ, retrieve from history
            if adapter and (reads_empty and forward_empty and reverse_empty) and session_id and self.history_manager:
                print(f"🔧 [EXECUTOR] Reads missing but adapter specified, checking session history for trimmed reads...")
                trimmed_result = self.history_manager.get_latest_result(session_id, "read_trimming")
                
                if trimmed_result and isinstance(trimmed_result, dict):
                    result_data = trimmed_result
                    # Check for nested result structure (from command handler)
                    if isinstance(result_data, dict) and "result" in result_data:
                        result_data = result_data["result"]
                    
                    print(f"🔧 [EXECUTOR] Found trimming result, file_type: {result_data.get('file_type')}")
                    
                    # Check if we have paired-end trimmed results
                    if result_data.get("file_type") == "paired_end":
                        forward_data = result_data.get("forward_reads", {})
                        reverse_data = result_data.get("reverse_reads", {})
                        
                        if isinstance(forward_data, dict) and "trimmed_reads" in forward_data:
                            forward_reads = forward_data["trimmed_reads"]
                            print(f"🔧 [EXECUTOR] Retrieved trimmed forward reads from history: {len(forward_reads)} chars")
                            parameters["forward_reads"] = forward_reads
                        
                        if isinstance(reverse_data, dict) and "trimmed_reads" in reverse_data:
                            reverse_reads = reverse_data["trimmed_reads"]
                            print(f"🔧 [EXECUTOR] Retrieved trimmed reverse reads from history: {len(reverse_reads)} chars")
                            parameters["reverse_reads"] = reverse_reads
                        
                        if forward_reads and reverse_reads:
                            file_type = "paired_end"
                            parameters["file_type"] = "paired_end"
                    else:
                        # Single file result
                        if "trimmed_reads" in result_data:
                            reads = result_data["trimmed_reads"]
                            print(f"🔧 [EXECUTOR] Retrieved trimmed reads from history: {len(reads)} chars")
                            parameters["reads"] = reads
                            file_type = "single"
                            parameters["file_type"] = "single"
            
            # Update variables from parameters if they were set above
            file_type = parameters.get("file_type", file_type)
            forward_reads = parameters.get("forward_reads", forward_reads)
            reverse_reads = parameters.get("reverse_reads", reverse_reads)
            reads = parameters.get("reads", reads)
            
            # Handle paired-end reads
            if file_type == "paired_end":
                
                print(f"🔧 [EXECUTOR] Paired-end mode detected")
                print(f"🔧 [EXECUTOR] Forward reads length: {len(forward_reads)}")
                print(f"🔧 [EXECUTOR] Reverse reads length: {len(reverse_reads)}")
                
                if not forward_reads or not reverse_reads:
                    return {
                        "status": "error",
                        "message": f"Both forward and reverse reads required for paired-end trimming. Forward: {len(forward_reads)} chars, Reverse: {len(reverse_reads)} chars"
                    }
                
                print(f"🔧 [EXECUTOR] Trimming paired-end reads - forward: {len(forward_reads)} chars, reverse: {len(reverse_reads)} chars")
                
                # Process forward reads
                forward_result = read_trimming.run_read_trimming_raw(
                    reads=forward_reads,
                    adapter=adapter,
                    quality_threshold=quality_threshold
                )
                
                # Process reverse reads
                reverse_result = read_trimming.run_read_trimming_raw(
                    reads=reverse_reads,
                    adapter=adapter,
                    quality_threshold=quality_threshold
                )
                
                # Combine results
                combined_summary = {
                    "forward": {
                        "total_reads": forward_result["summary"]["total_reads"],
                        "total_bases": forward_result["summary"]["total_bases"],
                        "trimmed_bases": forward_result["summary"]["trimmed_bases"],
                    },
                    "reverse": {
                        "total_reads": reverse_result["summary"]["total_reads"],
                        "total_bases": reverse_result["summary"]["total_bases"],
                        "trimmed_bases": reverse_result["summary"]["trimmed_bases"],
                    },
                    "quality_threshold": quality_threshold,
                    "adapter_removed": bool(adapter),
                }
                
                return {
                    "text": "Read trimming completed successfully for both forward and reverse reads.",
                    "forward_reads": {
                        "trimmed_reads": forward_result["trimmed_reads"],
                        "summary": forward_result["summary"]
                    },
                    "reverse_reads": {
                        "trimmed_reads": reverse_result["trimmed_reads"],
                        "summary": reverse_result["summary"]
                    },
                    "summary": combined_summary,
                    "file_type": "paired_end"
                }
            
            # Handle single file
            else:
                reads = parameters.get("reads", "")
                
                # Check if reads are empty or don't look like FASTQ
                if not reads or len(reads.strip()) == 0 or not reads.strip().startswith("@"):
                    # If we have an adapter, try to retrieve from history
                    if adapter and session_id and self.history_manager:
                        print(f"🔧 [EXECUTOR] Single file mode: reads empty, checking history for adapter removal...")
                        trimmed_result = self.history_manager.get_latest_result(session_id, "read_trimming")
                        
                        if trimmed_result and isinstance(trimmed_result, dict):
                            result_data = trimmed_result
                            # Check for nested result structure
                            if isinstance(result_data, dict) and "result" in result_data:
                                result_data = result_data["result"]
                            
                            # Check if we have paired-end results (use those for adapter removal)
                            if result_data.get("file_type") == "paired_end":
                                forward_data = result_data.get("forward_reads", {})
                                reverse_data = result_data.get("reverse_reads", {})
                                
                                if isinstance(forward_data, dict) and "trimmed_reads" in forward_data:
                                    forward_reads = forward_data["trimmed_reads"]
                                if isinstance(reverse_data, dict) and "trimmed_reads" in reverse_data:
                                    reverse_reads = reverse_data["trimmed_reads"]
                                
                                if forward_reads and reverse_reads:
                                    # Switch to paired-end mode
                                    file_type = "paired_end"
                                    parameters["forward_reads"] = forward_reads
                                    parameters["reverse_reads"] = reverse_reads
                                    parameters["file_type"] = "paired_end"
                                    # Re-run with paired-end logic
                                    return await self._execute_read_trimming(parameters, session_id)
                            elif "trimmed_reads" in result_data:
                                reads = result_data["trimmed_reads"]
                                parameters["reads"] = reads
                                print(f"🔧 [EXECUTOR] Retrieved trimmed reads from history: {len(reads)} chars")
                
                if not reads or len(reads.strip()) == 0 or not reads.strip().startswith("@"):
                    return {
                        "status": "error",
                        "message": "FASTQ reads required for trimming"
                    }
                
                print(f"🔧 [EXECUTOR] Trimming single file - reads length: {len(reads)}, quality_threshold: {quality_threshold}")
                result = read_trimming.run_read_trimming_raw(
                    reads=reads,
                    adapter=adapter,
                    quality_threshold=quality_threshold
                )
                
                print(f"🔧 [EXECUTOR] Read trimming result keys: {result.keys() if isinstance(result, dict) else 'Not a dict'}")
                print(f"🔧 [EXECUTOR] Read trimming summary: {result.get('summary', 'No summary') if isinstance(result, dict) else 'N/A'}")
                
                result["file_type"] = "single"
                return result
            
        except Exception as e:
            logger.error(f"Error in read trimming: {e}")
            import traceback
            traceback.print_exc()
            return {
                "status": "error",
                "message": f"Error trimming reads: {str(e)}"
            }
    
    async def _execute_read_merging(self, parameters: Dict[str, Any], session_id: Optional[str]) -> Dict[str, Any]:
        """Execute read merging command."""
        try:
            import read_merging
            
            forward_reads = parameters.get("forward_reads", "")
            reverse_reads = parameters.get("reverse_reads", "")
            min_overlap = parameters.get("min_overlap", 12)
            
            print(f"🔧 [EXECUTOR] Merging paired-end reads")
            print(f"🔧 [EXECUTOR] Forward reads length: {len(forward_reads)}")
            print(f"🔧 [EXECUTOR] Reverse reads length: {len(reverse_reads)}")
            print(f"🔧 [EXECUTOR] Minimum overlap: {min_overlap}")
            
            # If reads are missing, try to get them from previous trimming results
            # Check if reads are actually empty (not just falsy)
            forward_empty = not forward_reads or len(forward_reads.strip()) == 0
            reverse_empty = not reverse_reads or len(reverse_reads.strip()) == 0
            
            if (forward_empty or reverse_empty) and session_id and self.history_manager:
                print(f"🔧 [EXECUTOR] Reads missing, checking session history for trimmed reads...")
                trimmed_result = self.history_manager.get_latest_result(session_id, "read_trimming")
                
                if trimmed_result and isinstance(trimmed_result, dict):
                    print(f"🔧 [EXECUTOR] Found trimming result, keys: {list(trimmed_result.keys())}")
                    
                    # The result is stored directly (not wrapped), so check the structure
                    # For paired-end, it should have: forward_reads, reverse_reads, file_type
                    file_type = trimmed_result.get("file_type")
                    print(f"🔧 [EXECUTOR] File type in result: {file_type}")
                    
                    if file_type == "paired_end":
                        # Extract forward and reverse reads
                        forward_data = trimmed_result.get("forward_reads", {})
                        reverse_data = trimmed_result.get("reverse_reads", {})
                        
                        print(f"🔧 [EXECUTOR] Forward data type: {type(forward_data)}, Reverse data type: {type(reverse_data)}")
                        
                        # forward_reads and reverse_reads should be dicts with "trimmed_reads" key
                        if isinstance(forward_data, dict) and "trimmed_reads" in forward_data:
                            forward_reads = forward_data["trimmed_reads"]
                            print(f"🔧 [EXECUTOR] Extracted forward reads: {len(forward_reads)} chars")
                        
                        if isinstance(reverse_data, dict) and "trimmed_reads" in reverse_data:
                            reverse_reads = reverse_data["trimmed_reads"]
                            print(f"🔧 [EXECUTOR] Extracted reverse reads: {len(reverse_reads)} chars")
                    else:
                        # Single file result - can't use for merging
                        print(f"🔧 [EXECUTOR] Found trimming result but it's single file (file_type: {file_type}), not paired-end")
                else:
                    print(f"🔧 [EXECUTOR] No trimming results found in session history (result type: {type(trimmed_result)})")
            
            if not forward_reads or not reverse_reads:
                return {
                    "status": "error",
                    "message": f"Both forward and reverse reads required for merging. Forward: {len(forward_reads)} chars, Reverse: {len(reverse_reads)} chars. Please provide the trimmed reads or run trimming first."
                }
            
            result = read_merging.run_read_merging_raw(
                forward_reads=forward_reads,
                reverse_reads=reverse_reads,
                min_overlap=min_overlap
            )
            
            print(f"🔧 [EXECUTOR] Read merging result keys: {result.keys() if isinstance(result, dict) else 'Not a dict'}")
            print(f"🔧 [EXECUTOR] Read merging summary: {result.get('summary', 'No summary') if isinstance(result, dict) else 'N/A'}")
            
            return result
            
        except Exception as e:
            logger.error(f"Error in read merging: {e}")
            import traceback
            traceback.print_exc()
            return {
                "status": "error",
                "message": f"Error merging reads: {str(e)}"
            }
    
    async def _execute_quality_assessment(self, parameters: Dict[str, Any], session_id: Optional[str]) -> Dict[str, Any]:
        """Execute quality assessment command."""
        try:
            import quality_assessment
            
            sequences = parameters.get("sequences", "")
            
            # If sequences are missing, try to get merged sequences from history
            if (not sequences or len(sequences.strip()) == 0) and session_id and self.history_manager:
                print(f"🔧 [EXECUTOR] Sequences missing, checking session history for merged sequences...")
                merging_result = self.history_manager.get_latest_result(session_id, "read_merging", skip_errors=True)
                
                if merging_result and isinstance(merging_result, dict):
                    result_data = merging_result
                    # Check for nested result structure
                    if isinstance(result_data, dict) and "result" in result_data:
                        result_data = result_data["result"]
                    
                    if "merged_sequences" in result_data:
                        sequences = result_data["merged_sequences"]
                        print(f"🔧 [EXECUTOR] Retrieved merged sequences from history: {len(sequences)} chars")
            
            if not sequences or len(sequences.strip()) == 0:
                return {
                    "status": "error",
                    "message": "No sequences provided for quality assessment. Please provide sequences or run merging first."
                }
            
            result = quality_assessment.run_quality_assessment_raw(sequences)
            
            print(f"🔧 [EXECUTOR] Quality assessment result keys: {result.keys() if isinstance(result, dict) else 'Not a dict'}")
            
            return result
            
        except Exception as e:
            logger.error(f"Error in quality assessment: {e}")
            import traceback
            traceback.print_exc()
            return {
                "status": "error",
                "message": f"Error assessing quality: {str(e)}"
            }

def execute_command_raw(parsed_command: Dict[str, Any]) -> Dict[str, Any]:
    """Raw function for command execution (for direct calls)."""
    executor = CommandExecutor()
    
    # Check if we're already in an event loop
    try:
        loop = asyncio.get_running_loop()
        # We're in an event loop, create a task
        import concurrent.futures
        with concurrent.futures.ThreadPoolExecutor() as executor_pool:
            future = executor_pool.submit(asyncio.run, executor.execute_command(parsed_command))
            return future.result()
    except RuntimeError:
        # No event loop running, we can use asyncio.run
        return asyncio.run(executor.execute_command(parsed_command)) 