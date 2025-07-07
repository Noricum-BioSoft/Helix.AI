#!/usr/bin/env python3
"""
Comprehensive test script for the Enhanced Bioinformatics MCP Server
"""

import asyncio
import json
import sys
from pathlib import Path

# Add the current directory to Python path
sys.path.append(str(Path(__file__).parent))

async def test_enhanced_mcp_server():
    """Test the enhanced MCP server functionality."""
    
    print("🧪 Testing Enhanced Bioinformatics MCP Server")
    print("=" * 60)
    
    # Test 1: Import and basic functionality
    print("\n1. Testing server import and initialization...")
    try:
        from mcp_server_enhanced import (
            validate_sequence, 
            validate_fasta_format, 
            parse_fasta_to_dataframe
        )
        print("✅ Enhanced MCP server modules imported successfully")
    except Exception as e:
        print(f"❌ Import failed: {e}")
        return
    
    # Test 2: Sequence validation
    print("\n2. Testing sequence validation...")
    test_sequences = [
        "ACTGTTGAC",  # Valid
        "actgttgac",  # Valid (case insensitive)
        "ACTG N TTGAC",  # Valid (with spaces)
        "ACTGTTGAC123",  # Invalid
        "",  # Invalid (empty)
        "ACTGTTGACX",  # Invalid character
    ]
    
    for seq in test_sequences:
        is_valid = validate_sequence(seq)
        status = "✅" if is_valid else "❌"
        print(f"   {status} '{seq}' -> {is_valid}")
    
    # Test 3: FASTA validation
    print("\n3. Testing FASTA format validation...")
    valid_fasta = """>seq1
ACTGTTGAC
>seq2
ACTGCATCC"""
    
    invalid_fasta = """seq1
ACTGTTGAC
seq2
ACTGCATCC"""
    
    print(f"   ✅ Valid FASTA: {validate_fasta_format(valid_fasta)}")
    print(f"   ❌ Invalid FASTA: {validate_fasta_format(invalid_fasta)}")
    
    # Test 4: FASTA parsing
    print("\n4. Testing FASTA parsing...")
    try:
        df = parse_fasta_to_dataframe(valid_fasta)
        print(f"   ✅ Parsed {len(df)} sequences")
        print(f"   Sequence names: {list(df['name'])}")
    except Exception as e:
        print(f"   ❌ FASTA parsing failed: {e}")
    
    # Test 5: Enhanced tools functionality
    print("\n5. Testing enhanced tools...")
    
    # Test sequence statistics
    print("   Testing sequence statistics...")
    try:
        from mcp_server_enhanced import handle_sequence_statistics
        result = await handle_sequence_statistics({
            "sequence": "ACTGTTGAC",
            "include_composition": True
        })
        print(f"   ✅ Statistics calculated: {result['statistics']['length']} bp")
        print(f"   GC content: {result['statistics']['gc_content']:.2%}")
    except Exception as e:
        print(f"   ❌ Statistics failed: {e}")
    
    # Test reverse complement
    print("   Testing reverse complement...")
    try:
        from mcp_server_enhanced import handle_reverse_complement
        result = await handle_reverse_complement({
            "sequence": "ACTGTTGAC"
        })
        print(f"   ✅ Original: {result['original_sequence']}")
        print(f"   Reverse complement: {result['reverse_complement']}")
    except Exception as e:
        print(f"   ❌ Reverse complement failed: {e}")
    
    # Test sequence validation tool
    print("   Testing sequence validation tool...")
    try:
        from mcp_server_enhanced import handle_sequence_validation
        result = await handle_sequence_validation({
            "sequence": "ACTGTTGAC",
            "fasta_content": valid_fasta
        })
        print(f"   ✅ Validation results: {result['validations']}")
    except Exception as e:
        print(f"   ❌ Validation tool failed: {e}")
    
    # Test 6: Error handling
    print("\n6. Testing error handling...")
    
    # Test invalid sequence
    try:
        result = await handle_sequence_statistics({
            "sequence": "ACTGTTGAC123"  # Invalid
        })
        print("   ❌ Should have failed validation")
    except Exception as e:
        print(f"   ✅ Properly caught validation error: {type(e).__name__}")
    
    # Test missing required parameter
    try:
        result = await handle_sequence_statistics({})
        print("   ❌ Should have failed with missing parameter")
    except Exception as e:
        print(f"   ✅ Properly caught missing parameter error: {type(e).__name__}")
    
    # Test 7: Performance and timing
    print("\n7. Testing performance features...")
    
    import time
    start_time = time.time()
    
    try:
        from mcp_server_enhanced import handle_sequence_alignment
        result = await handle_sequence_alignment({
            "sequences": valid_fasta,
            "algorithm": "clustal"
        })
        execution_time = time.time() - start_time
        print(f"   ✅ Alignment completed in {execution_time:.2f}s")
        print(f"   Execution time recorded: {result.get('execution_time', 'N/A')}")
    except Exception as e:
        print(f"   ❌ Performance test failed: {e}")
    
    print("\n" + "=" * 60)
    print("🎉 Enhanced MCP Server Test Completed!")

async def test_integration_with_existing_tools():
    """Test integration with existing bioinformatics tools."""
    
    print("\n🔗 Testing Integration with Existing Tools")
    print("-" * 50)
    
    # Test integration with alignment tool
    print("\n1. Testing alignment tool integration...")
    try:
        from tools.alignment import run_alignment
        from mcp_server_enhanced import validate_fasta_format
        
        test_sequences = """>seq1
ACTGTTGAC
>seq2
ACTGCATCC"""
        
        if validate_fasta_format(test_sequences):
            result = run_alignment(test_sequences)
            print(f"   ✅ Alignment integration successful: {result['text']}")
        else:
            print("   ❌ FASTA validation failed")
    except Exception as e:
        print(f"   ❌ Alignment integration failed: {e}")
    
    # Test integration with mutation tool
    print("\n2. Testing mutation tool integration...")
    try:
        from tools.mutations import run_mutation_raw
        from mcp_server_enhanced import validate_sequence
        
        test_sequence = "ACTGTTGAC"
        
        if validate_sequence(test_sequence):
            result = run_mutation_raw(test_sequence, 5)
            print(f"   ✅ Mutation integration successful: {result['text']}")
        else:
            print("   ❌ Sequence validation failed")
    except Exception as e:
        print(f"   ❌ Mutation integration failed: {e}")
    
    # Test integration with bio tool
    print("\n3. Testing bio tool integration...")
    try:
        from tools.bio import align_and_visualize_fasta
        import pandas as pd
        
        test_data = pd.DataFrame({
            'name': ['seq1', 'seq2'],
            'sequence': ['ACTGTTGAC', 'ACTGCATCC']
        })
        
        result = align_and_visualize_fasta(test_data)
        print(f"   ✅ Bio tool integration successful: {result}")
    except Exception as e:
        print(f"   ❌ Bio tool integration failed: {e}")

async def main():
    """Main test function."""
    print("Enhanced Bioinformatics MCP Server Test Suite")
    print("=" * 60)
    
    await test_enhanced_mcp_server()
    await test_integration_with_existing_tools()
    
    print("\n" + "=" * 60)
    print("✅ All tests completed!")

if __name__ == "__main__":
    asyncio.run(main()) 