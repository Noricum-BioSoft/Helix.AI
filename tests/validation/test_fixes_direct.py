#!/usr/bin/env python3
"""
Direct test of the fixes by importing and testing the tools.
This tests the actual code that was deployed.
"""

import sys
from pathlib import Path

# Add tools directory to path (same as deployed environment)
tools_path = str((Path(__file__).resolve().parent.parent.parent / "tools").resolve())
sys.path.insert(0, tools_path)

# Add backend directory to path
backend_path = str((Path(__file__).resolve().parent.parent.parent / "backend").resolve())
sys.path.insert(0, backend_path)

def test_sequence_alignment_import():
    """Test that sequence_alignment can be imported without Bio.Align.Applications error."""
    print("1️⃣ Testing sequence_alignment import (Bio.Align.Applications fix)...")
    try:
        import alignment
        print("   ✅ alignment module imported successfully")
        
        # Check if the tool decorator is available
        if hasattr(alignment, 'run_alignment_tool'):
            print("   ✅ run_alignment_tool function found")
        else:
            print("   ⚠️  run_alignment_tool function not found")
        
        # Test that we can call the function
        test_sequences = ">seq1\nATGCGATCG\n>seq2\nATGCGATCA"
        try:
            result = alignment.run_alignment(test_sequences)
            if result and "text" in result:
                print("   ✅ run_alignment executed successfully")
                print(f"   ✅ Result preview: {result.get('text', '')[:100]}...")
                return True
            else:
                print("   ⚠️  run_alignment returned unexpected result")
                return False
        except Exception as e:
            error_msg = str(e)
            if "Bio.Align.Applications" in error_msg or "No module named" in error_msg:
                print(f"   ❌ Import error still present: {error_msg}")
                return False
            else:
                print(f"   ⚠️  Execution error (not import-related): {error_msg}")
                return False
                
    except ImportError as e:
        error_msg = str(e)
        if "Bio.Align.Applications" in error_msg:
            print(f"   ❌ Bio.Align.Applications import error: {error_msg}")
            return False
        else:
            print(f"   ❌ Import error: {error_msg}")
            return False
    except Exception as e:
        print(f"   ❌ Unexpected error: {e}")
        return False

def test_mutate_sequence_import():
    """Test that mutate_sequence can be imported without langchain.agents error."""
    print("\n2️⃣ Testing mutate_sequence import (langchain_core.tools fix)...")
    try:
        import mutations
        print("   ✅ mutations module imported successfully")
        
        # Check if the tool decorator is available
        if hasattr(mutations, 'run_mutation'):
            print("   ✅ run_mutation function found (with @tool decorator)")
        else:
            print("   ⚠️  run_mutation function not found")
        
        # Test that we can call the function
        try:
            result = mutations.run_mutation_raw("ATGCGATCG", num_variants=3)
            if result and "variants" in result.get("statistics", {}):
                print("   ✅ run_mutation_raw executed successfully")
                variants = result["statistics"]["variants"]
                print(f"   ✅ Generated {len(variants)} variants")
                return True
            else:
                print("   ⚠️  run_mutation_raw returned unexpected result")
                return False
        except Exception as e:
            error_msg = str(e)
            if "langchain.agents" in error_msg or "cannot import name 'tool'" in error_msg:
                print(f"   ❌ Import error still present: {error_msg}")
                return False
            else:
                print(f"   ⚠️  Execution error (not import-related): {error_msg}")
                return False
                
    except ImportError as e:
        error_msg = str(e)
        if "langchain.agents" in error_msg or "cannot import name 'tool'" in error_msg:
            print(f"   ❌ langchain.agents import error: {error_msg}")
            return False
        else:
            print(f"   ❌ Import error: {error_msg}")
            return False
    except Exception as e:
        print(f"   ❌ Unexpected error: {e}")
        return False

def test_other_tools_import():
    """Test that other tools can still be imported."""
    print("\n3️⃣ Testing other tool imports...")
    tools_to_test = [
        ("phylogenetic_tree", "run_phylogenetic_tree"),
        ("alignment", "run_alignment_tool"),
    ]
    
    all_passed = True
    for tool_module, tool_func in tools_to_test:
        try:
            module = __import__(tool_module)
            if hasattr(module, tool_func):
                print(f"   ✅ {tool_module}.{tool_func} imported successfully")
            else:
                print(f"   ⚠️  {tool_module}.{tool_func} not found")
                all_passed = False
        except ImportError as e:
            error_msg = str(e)
            if "langchain.agents" in error_msg or "cannot import name 'tool'" in error_msg:
                print(f"   ❌ {tool_module} has import error: {error_msg}")
                all_passed = False
            else:
                print(f"   ⚠️  {tool_module} import error: {error_msg}")
    
    return all_passed

def main():
    print("🧪 Testing Deployed Fixes (Direct Import Test)")
    print("=" * 60)
    print("This test directly imports and tests the fixed modules.")
    print("=" * 60)
    
    results = []
    
    # Test 1: sequence_alignment
    results.append(("sequence_alignment", test_sequence_alignment_import()))
    
    # Test 2: mutate_sequence
    results.append(("mutate_sequence", test_mutate_sequence_import()))
    
    # Test 3: Other tools
    results.append(("other_tools", test_other_tools_import()))
    
    # Summary
    print("\n" + "=" * 60)
    print("📊 Test Summary")
    print("=" * 60)
    
    all_passed = True
    for test_name, passed in results:
        status = "✅ PASSED" if passed else "❌ FAILED"
        print(f"{test_name}: {status}")
        if not passed:
            all_passed = False
    
    print("=" * 60)
    if all_passed:
        print("🎉 All tests passed! Your fixes are working correctly!")
        return 0
    else:
        print("⚠️  Some tests failed. Please check the errors above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())



