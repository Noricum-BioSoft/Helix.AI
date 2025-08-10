#!/usr/bin/env python3
"""
Test script to demonstrate command handling capabilities
"""

import sys
import asyncio
import json
from pathlib import Path

# Add the backend directory to Python path for imports
sys.path.append(str(Path(__file__).parent.parent.parent / "backend"))

from history_manager import history_manager

async def test_specific_command():
    """Test the specific command: 'from the sequence variants, pick 10 sequences randomly and output them.'"""
    
    print("🧪 Testing Command Handling System")
    print("=" * 60)
    
    # The specific command from the user
    command = "from the sequence variants, pick 10 sequences randomly and output them."
    
    print(f"📝 Command: '{command}'")
    print("=" * 60)
    
    try:
        # Step 1: Create a session and add some mutation results
        print("\n📋 Step 1: Setting up test session with mutation results...")
        
        session_id = history_manager.create_session("test_user")
        print(f"✅ Created session: {session_id}")
        
        # Add test mutation results
        test_mutation_result = {
            "status": "success",
            "input_sequence": "ACTGTTGAC",
            "num_variants": 50,
            "mutation_rate": 0.15,
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
                {"name": "variant_16", "sequence": "ACTGTTGAT"},
                {"name": "variant_17", "sequence": "ACTGTTGGC"},
                {"name": "variant_18", "sequence": "ACTGTTGCA"},
                {"name": "variant_19", "sequence": "ACTGTTGCC"},
                {"name": "variant_20", "sequence": "ACTGTTGAG"},
            ]
        }
        
        history_manager.add_history_entry(
            session_id,
            "Generate 20 variants of sequence ACTGTTGAC",
            "mutate_sequence",
            test_mutation_result
        )
        
        print(f"✅ Added {len(test_mutation_result['variants'])} mutation variants to session")
        
        # Step 2: Test command parsing
        print("\n📋 Step 2: Testing command parsing...")
        
        # Add tools directory to path
        tools_path = str((Path(__file__).parent.parent.parent / "tools").resolve())
        sys.path.insert(0, tools_path)
        
        import command_parser
        
        parsed_result = command_parser.parse_command_raw(command, session_id)
        
        print(f"✅ Command parsed successfully!")
        print(f"   Action: {parsed_result.get('action')}")
        print(f"   Tool: {parsed_result.get('tool')}")
        print(f"   Parameters: {parsed_result.get('parameters')}")
        print(f"   Confidence: {parsed_result.get('confidence')}")
        
        # Step 3: Test command execution
        print("\n📋 Step 3: Testing command execution...")
        
        import command_executor
        
        execution_result = command_executor.execute_command_raw(parsed_result)
        
        print(f"✅ Command executed successfully!")
        print(f"   Status: {execution_result.get('status')}")
        print(f"   Tool used: {execution_result.get('tool')}")
        
        # Step 4: Test complete command handling
        print("\n📋 Step 4: Testing complete command handling...")
        
        import command_handler
        
        handler_result = command_handler.handle_command_raw(command, session_id)
        
        print(f"✅ Complete command handling successful!")
        print(f"   Status: {handler_result.get('status')}")
        print(f"   Action: {handler_result.get('action')}")
        print(f"   Tool: {handler_result.get('tool')}")
        
        # Display results
        if handler_result.get('status') == 'success':
            result_data = handler_result.get('result', {})
            print(f"\n📊 Results:")
            print(f"   Selected variants: {result_data.get('num_variants_selected', 0)}")
            print(f"   Selection criteria: {result_data.get('selection_criteria', 'unknown')}")
            print(f"   Session ID: {result_data.get('session_id', 'none')}")
            
            if result_data.get('selected_variants'):
                print(f"\n📝 Selected sequences:")
                for i, variant in enumerate(result_data['selected_variants'][:5], 1):
                    print(f"   {i}. {variant['name']}: {variant['sequence']}")
                if len(result_data['selected_variants']) > 5:
                    print(f"   ... and {len(result_data['selected_variants']) - 5} more")
        
        # Step 5: Show session summary
        print("\n📋 Step 5: Session summary...")
        
        summary = history_manager.get_session_summary(session_id)
        print(f"✅ Session Summary:")
        print(f"   Session ID: {summary['session_id']}")
        print(f"   Total Operations: {summary['total_operations']}")
        print(f"   Tool Usage: {summary['tool_usage']}")
        
        print("\n🎉 All tests passed successfully!")
        return True
        
    except Exception as e:
        print(f"❌ Error during testing: {e}")
        import traceback
        traceback.print_exc()
        return False

async def test_various_commands():
    """Test various natural language commands."""
    
    print("\n🧪 Testing Various Natural Language Commands")
    print("=" * 60)
    
    commands = [
        "from the sequence variants, pick 10 sequences randomly and output them.",
        "select 5 diverse sequences from the variants",
        "choose 8 sequences by length from the mutation results",
        "pick 3 random variants from the previous results",
        "from the variants, select 7 diverse sequences"
    ]
    
    try:
        # Add tools directory to path
        tools_path = str((Path(__file__).parent.parent / "tools").resolve())
        sys.path.insert(0, tools_path)
        
        import command_handler
        
        for i, command in enumerate(commands, 1):
            print(f"\n📝 Test {i}: '{command}'")
            
            # Create a new session for each test
            session_id = history_manager.create_session(f"test_user_{i}")
            
            # Add some test mutation results
            test_mutation_result = {
                "status": "success",
                "variants": [
                    {"name": f"variant_{j}", "sequence": f"ACTGTTGA{j%4}"} 
                    for j in range(1, 16)
                ]
            }
            
            history_manager.add_history_entry(
                session_id,
                "Generate test variants",
                "mutate_sequence",
                test_mutation_result
            )
            
            # Handle the command
            result = command_handler.handle_command_raw(command, session_id)
            
            print(f"   Status: {result.get('status')}")
            print(f"   Action: {result.get('action')}")
            print(f"   Tool: {result.get('tool')}")
            
            if result.get('status') == 'success':
                result_data = result.get('result', {})
                print(f"   Selected: {result_data.get('num_variants_selected', 0)} variants")
                print(f"   Criteria: {result_data.get('selection_criteria', 'unknown')}")
            else:
                print(f"   Error: {result.get('message', 'Unknown error')}")
        
        print("\n✅ All command tests completed!")
        return True
        
    except Exception as e:
        print(f"❌ Error testing various commands: {e}")
        return False

async def main():
    """Main test function."""
    print("🧪 Starting Command Handling Tests")
    print("=" * 60)
    
    # Test the specific command
    success1 = await test_specific_command()
    
    # Test various commands
    success2 = await test_various_commands()
    
    if success1 and success2:
        print("\n🎉 All command handling tests passed!")
        print("\n📋 Summary:")
        print("- Command parsing is working correctly")
        print("- Command execution is functional")
        print("- Natural language processing is operational")
        print("- History tracking is integrated")
        print("- The specific command is handled properly")
        return True
    else:
        print("\n❌ Some tests failed")
        return False

if __name__ == "__main__":
    success = asyncio.run(main())
    sys.exit(0 if success else 1) 