#!/usr/bin/env python3
"""
Simple test to verify tool imports work correctly
"""

import sys
from pathlib import Path

tools_path = str((Path(__file__).resolve().parent.parent / "tools").resolve())
sys.path.insert(0, tools_path)
print(f"🔧 Added tools path: {tools_path}")

def test_imports():
    """Test if all tool imports work."""
    print("🧪 Testing tool imports...")
    
    try:
        import alignment
        import bio
        import mutations
        import data_science
        print("✅ All tool modules imported successfully")
    except ImportError as e:
        print(f"❌ Tool module import failed: {e}")
        print(f"Current sys.path: {sys.path}")
        return False
    
    try:
        # Get the functions from the modules
        run_alignment = alignment.run_alignment
        align_and_visualize_fasta = bio.align_and_visualize_fasta
        run_mutation = mutations.run_mutation
        analyze_basic_stats = data_science.analyze_basic_stats
        print("✅ All tool functions imported successfully")
    except AttributeError as e:
        print(f"❌ Tool function import failed: {e}")
        return False
    
    print("🎉 All tool imports successful!")
    return True

def test_function_calls():
    """Test if the functions can be called."""
    print("\n🧪 Testing function calls...")
    
    try:
        import alignment
        result = alignment.run_alignment(">seq1\nATCG\n>seq2\nGCTA")
        print("✅ run_alignment call successful")
    except Exception as e:
        print(f"❌ run_alignment call failed: {e}")
        return False
    
    try:
        import mutations
        result = mutations.run_mutation_raw("ATCG", 5)
        print("✅ run_mutation call successful")
    except Exception as e:
        print(f"❌ run_mutation call failed: {e}")
        return False
    
    print("🎉 All function calls successful!")
    return True

if __name__ == "__main__":
    print("🔍 Testing MCP Server Tool Imports")
    print("=" * 40)
    
    if test_imports() and test_function_calls():
        print("\n✅ All tests passed! MCP server should work correctly.")
        sys.exit(0)
    else:
        print("\n❌ Some tests failed. Check the errors above.")
        sys.exit(1) 