#!/usr/bin/env python3
"""
Test script for mutation tool
"""

import sys
from pathlib import Path

# Add tools to path
tools_path = str((Path(__file__).resolve().parent.parent / "tools").resolve())
sys.path.insert(0, tools_path)

def test_mutation():
    """Test the mutation tool directly."""
    try:
        import mutations
        print("✅ Successfully imported mutations module")
        
        # Test the mutation function
        result = mutations.run_mutation_raw("ACTGTTGAC", 10)
        print(f"✅ Mutation test successful: {result}")
        
        return True
    except Exception as e:
        print(f"❌ Mutation test failed: {e}")
        return False

if __name__ == "__main__":
    print("🧪 Testing mutation tool...")
    success = test_mutation()
    if success:
        print("✅ All tests passed!")
    else:
        print("❌ Tests failed!")
        sys.exit(1) 