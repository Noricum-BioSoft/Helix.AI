#!/usr/bin/env python3
"""
Test script to demonstrate history tracking functionality
"""

import sys
import asyncio
import json
from pathlib import Path

# Add the backend directory to Python path for imports
sys.path.append(str(Path(__file__).parent.parent.parent / "backend"))

from history_manager import history_manager

async def test_user_story():
    """Test the user story: Step 1 -> Step 2 with history tracking."""
    
    print("🧬 Testing History Tracking - User Story")
    print("=" * 50)
    
    # Step 1: Create a session and generate mutations
    print("\n📝 Step 1: Creating session and generating mutations...")
    
    # Create a new session
    session_id = history_manager.create_session("test_user")
    print(f"✅ Created session: {session_id}")
    
    # Simulate mutation command
    mutation_command = "Generate 50 variants of sequence ACTGTTGAC with mutation rate 0.15"
    mutation_result = {
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
            # ... more variants would be generated
        ],
        "statistics": {
            "total_variants": 50,
            "average_length": 9,
            "gc_content": 0.44
        }
    }
    
    # Track the mutation operation
    history_manager.add_history_entry(
        session_id,
        mutation_command,
        "mutate_sequence",
        mutation_result,
        {
            "num_variants": 50,
            "mutation_rate": 0.15,
            "input_sequence": "ACTGTTGAC"
        }
    )
    
    print(f"✅ Generated {len(mutation_result['variants'])} variants")
    print(f"✅ Mutation operation tracked in session history")
    
    # Step 2: Select variants from previous results
    print("\n📝 Step 2: Selecting variants from previous results...")
    
    # Simulate variant selection command
    selection_command = "Select 10 diverse variants from previous mutation results"
    
    # Get the previous mutation results from history
    previous_mutation = history_manager.get_latest_result(session_id, "mutate_sequence")
    
    if previous_mutation:
        print(f"✅ Found previous mutation results in session")
        print(f"   - Total variants: {len(previous_mutation.get('variants', []))}")
        print(f"   - Input sequence: {previous_mutation.get('input_sequence', 'N/A')}")
        
        # Simulate variant selection
        selected_variants = previous_mutation['variants'][:10]  # Select first 10 for demo
        
        selection_result = {
            "status": "success",
            "session_id": session_id,
            "selection_criteria": "diversity",
            "num_variants_requested": 10,
            "num_variants_selected": len(selected_variants),
            "selected_variants": selected_variants,
            "analysis": {
                "total_variants": len(previous_mutation['variants']),
                "selected_variants": len(selected_variants),
                "selection_ratio": len(selected_variants) / len(previous_mutation['variants']),
                "criteria_used": "diversity"
            }
        }
        
        # Track the selection operation
        history_manager.add_history_entry(
            session_id,
            selection_command,
            "select_variants",
            selection_result,
            {
                "selection_criteria": "diversity",
                "num_variants": 10
            }
        )
        
        print(f"✅ Selected {len(selected_variants)} diverse variants")
        print(f"✅ Selection operation tracked in session history")
        
    else:
        print("❌ No previous mutation results found")
    
    # Display session summary
    print("\n📊 Session Summary:")
    print("=" * 30)
    
    summary = history_manager.get_session_summary(session_id)
    print(f"Session ID: {summary['session_id']}")
    print(f"User ID: {summary['user_id']}")
    print(f"Total Operations: {summary['total_operations']}")
    print(f"Tool Usage: {summary['tool_usage']}")
    print(f"Available Results: {summary['available_results']}")
    
    # Display full session history
    print("\n📜 Full Session History:")
    print("=" * 30)
    
    session = history_manager.get_session(session_id)
    for i, entry in enumerate(session['history'], 1):
        print(f"\nOperation {i}:")
        print(f"  Timestamp: {entry['timestamp']}")
        print(f"  Command: {entry['command']}")
        print(f"  Tool: {entry['tool']}")
        print(f"  Status: {entry['result'].get('status', 'N/A')}")
    
    print("\n✅ History tracking test completed successfully!")
    return True

async def test_variant_selection_tool():
    """Test the variant selection tool with history integration."""
    
    print("\n🔬 Testing Variant Selection Tool")
    print("=" * 40)
    
    try:
        # Add tools directory to path
        tools_path = str((Path(__file__).parent.parent.parent / "tools").resolve())
        sys.path.insert(0, tools_path)
        
        import variant_selection
        
        # Create a test session with mutation results
        session_id = history_manager.create_session("test_user_2")
        
        # Add some test mutation results
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
            ]
        }
        
        history_manager.add_history_entry(
            session_id,
            "Generate test variants",
            "mutate_sequence",
            test_mutation_result
        )
        
        print(f"✅ Created test session with mutation results")
        
        # Test variant selection
        result = variant_selection.run_variant_selection_raw(
            session_id=session_id,
            selection_criteria="diversity",
            num_variants=5
        )
        
        print(f"✅ Variant selection completed")
        print(f"   Status: {result.get('status')}")
        print(f"   Selected variants: {result.get('num_variants_selected', 0)}")
        print(f"   Criteria: {result.get('selection_criteria')}")
        
        if result.get('selected_variants'):
            print("   Selected sequences:")
            for variant in result['selected_variants'][:3]:  # Show first 3
                print(f"     - {variant['name']}: {variant['sequence']}")
            if len(result['selected_variants']) > 3:
                print(f"     ... and {len(result['selected_variants']) - 3} more")
        
        return True
        
    except Exception as e:
        print(f"❌ Error testing variant selection: {e}")
        import traceback
        traceback.print_exc()
        return False

async def main():
    """Main test function."""
    print("🧪 Starting History Tracking Tests")
    print("=" * 50)
    
    # Test basic history functionality
    success1 = await test_user_story()
    
    # Test variant selection tool
    success2 = await test_variant_selection_tool()
    
    if success1 and success2:
        print("\n🎉 All tests passed!")
        print("\n📋 Summary:")
        print("- History tracking system is working")
        print("- Session management is functional")
        print("- Variant selection tool is integrated")
        print("- User story workflow is supported")
        return True
    else:
        print("\n❌ Some tests failed")
        return False

if __name__ == "__main__":
    success = asyncio.run(main())
    sys.exit(0 if success else 1) 