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