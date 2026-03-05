import logging
import sys
from pathlib import Path
from typing import Dict, Any, Optional
import asyncio

logger = logging.getLogger(__name__)

class CommandHandler:
    """Handles natural language commands by parsing and executing them."""
    
    def __init__(self):
        # Add tools directory to path
        tools_path = str((Path(__file__).parent).resolve())
        sys.path.insert(0, tools_path)
        
        # Add backend directory to path for history manager
        backend_path = str((Path(__file__).parent.parent / "backend").resolve())
        sys.path.insert(0, backend_path)
        
        try:
            from command_parser import CommandParser
            from command_executor import CommandExecutor
            from history_manager import history_manager
            
            self.parser = CommandParser()
            self.executor = CommandExecutor()
            self.history_manager = history_manager
            
        except ImportError as e:
            logger.error(f"Failed to import required modules: {e}")
            raise
    
    async def handle_command(self, command: str, session_id: Optional[str] = None) -> Dict[str, Any]:
        """Handle a natural language command by parsing and executing it."""
        logger.info(f"Handling command: '{command}'")
        
        try:
            # Step 1: Parse the command
            parsed_command = self.parser.parse_command(command, session_id)
            
            # Convert parsed command to dict format
            command_dict = {
                "action": parsed_command.action,
                "tool": parsed_command.tool,
                "parameters": parsed_command.parameters,
                "session_id": parsed_command.session_id,
                "confidence": parsed_command.confidence,
                "original_command": command
            }
            
            logger.info(f"Parsed command: {command_dict}")
            
            # Step 2: Execute the command
            result = await self.executor.execute_command(command_dict)
            
            logger.info(f"Command execution result: {result}")
            
            return result
            
        except Exception as e:
            logger.error(f"Error handling command: {e}")
            return {
                "status": "error",
                "message": f"Error handling command: {str(e)}",
                "original_command": command,
                "session_id": session_id
            }
    
    def create_session_if_needed(self, session_id: Optional[str] = None) -> str:
        """Create a session if one doesn't exist."""
        if not session_id:
            session_id = self.history_manager.create_session()
            logger.info(f"Created new session: {session_id}")
        
        return session_id
    
    def get_session_summary(self, session_id: str) -> Dict[str, Any]:
        """Get summary of a session."""
        return self.history_manager.get_session_summary(session_id)
    
    def get_latest_result(self, session_id: str, tool: str) -> Optional[Dict[str, Any]]:
        """Get the latest result for a specific tool in a session."""
        return self.history_manager.get_latest_result(session_id, tool)

def handle_command_raw(command: str, session_id: Optional[str] = None) -> Dict[str, Any]:
    """Raw function for command handling (for direct calls)."""
    handler = CommandHandler()
    
    # Check if we're already in an event loop
    try:
        loop = asyncio.get_running_loop()
        # We're in an event loop, create a task
        import concurrent.futures
        with concurrent.futures.ThreadPoolExecutor() as executor_pool:
            future = executor_pool.submit(asyncio.run, handler.handle_command(command, session_id))
            return future.result()
    except RuntimeError:
        # No event loop running, we can use asyncio.run
        return asyncio.run(handler.handle_command(command, session_id))

# Test function for the specific command
def test_specific_command():
    """Test the specific command: 'from the sequence variants, pick 10 sequences randomly and output them.'"""
    command = "from the sequence variants, pick 10 sequences randomly and output them."
    
    print(f"🧪 Testing specific command: '{command}'")
    print("=" * 60)
    
    try:
        # Create a test session with some mutation results
        handler = CommandHandler()
        session_id = handler.create_session_if_needed()
        
        # Add some test mutation results to the session
        test_mutation_result = {
            "status": "success",
            "variants": [
                {"name": "variant_1", "sequence": "ACTGTTGAC"},
                {"name": "variant_2", "sequence": "ACTGTTGGC"},
                {"name": "variant_3", "sequence": "ACTGTTGAT"},
                {"name": "variant_4", "sequence": "ACTGTTGCC"},
                {"name": "variant_5", "sequence": "ACTGTTGAG"},
                {"name": "variant_6", "sequence": "ACTGTTGTC"},
                {"name": "variant_7", "sequence": "ACTGTTGAA"},
                {"name": "variant_8", "sequence": "ACTGTTGCT"},
                {"name": "variant_9", "sequence": "ACTGTTGCA"},
                {"name": "variant_10", "sequence": "ACTGTTGCC"},
                {"name": "variant_11", "sequence": "ACTGTTGTA"},
                {"name": "variant_12", "sequence": "ACTGTTGGA"},
                {"name": "variant_13", "sequence": "ACTGTTGCT"},
                {"name": "variant_14", "sequence": "ACTGTTGAG"},
                {"name": "variant_15", "sequence": "ACTGTTGTC"},
            ]
        }
        
        handler.history_manager.add_history_entry(
            session_id,
            "Generate test variants",
            "mutate_sequence",
            test_mutation_result
        )
        
        print(f"✅ Created test session with {len(test_mutation_result['variants'])} variants")
        
        # Handle the specific command
        result = asyncio.run(handler.handle_command(command, session_id))
        
        print(f"\n📋 Command Analysis:")
        print(f"   Original command: '{command}'")
        print(f"   Parsed action: {result.get('action', 'unknown')}")
        print(f"   Tool used: {result.get('tool', 'unknown')}")
        print(f"   Session ID: {result.get('session_id', 'none')}")
        print(f"   Status: {result.get('status', 'unknown')}")
        
        if result.get('status') == 'success':
            result_data = result.get('result', {})
            print(f"\n✅ Command executed successfully!")
            print(f"   Selected variants: {result_data.get('num_variants_selected', 0)}")
            print(f"   Selection criteria: {result_data.get('selection_criteria', 'unknown')}")
            
            if result_data.get('selected_variants'):
                print(f"\n📝 Selected sequences:")
                for i, variant in enumerate(result_data['selected_variants'][:5], 1):
                    print(f"   {i}. {variant['name']}: {variant['sequence']}")
                if len(result_data['selected_variants']) > 5:
                    print(f"   ... and {len(result_data['selected_variants']) - 5} more")
        else:
            print(f"\n❌ Command failed: {result.get('message', 'Unknown error')}")
        
        return result
        
    except Exception as e:
        print(f"❌ Error testing command: {e}")
        import traceback
        traceback.print_exc()
        return {"status": "error", "message": str(e)}

if __name__ == "__main__":
    test_specific_command() 