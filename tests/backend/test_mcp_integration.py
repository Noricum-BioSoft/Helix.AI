#!/usr/bin/env python3
"""
Test script for MCP server integration
"""

import asyncio
import json
import sys
from pathlib import Path

# Add the current directory to Python path
sys.path.append(str(Path(__file__).parent))

async def test_mcp_tools():
    """Test the MCP tool integration."""
    
    print("Testing MCP Server Integration...")
    print("=" * 50)
    
    # Test sequence alignment
    print("\n1. Testing sequence alignment...")
    try:
        from tools.alignment import run_alignment
        result = run_alignment(">seq1\nACTGTTGAC\n>seq2\nACTGCATCC")
        print(f"✓ Sequence alignment successful: {result}")
    except Exception as e:
        print(f"✗ Sequence alignment failed: {e}")
    
    # Test mutation generation
    print("\n2. Testing sequence mutation...")
    try:
        from tools.mutations import mutate_sequence
        result = mutate_sequence("ACTGTTGAC", 5)
        print(f"✓ Sequence mutation successful: {result}")
    except Exception as e:
        print(f"✗ Sequence mutation failed: {e}")
    
    # Test bioinformatics analysis
    print("\n3. Testing bioinformatics analysis...")
    try:
        from tools.bio import align_and_visualize_fasta
        import pandas as pd
        
        # Create test data
        test_data = pd.DataFrame({
            'name': ['seq1', 'seq2', 'seq3'],
            'sequence': ['ACTGTTGAC', 'ACTGCATCC', 'ACTGCAATGAC']
        })
        
        result = align_and_visualize_fasta(test_data)
        print(f"✓ Bioinformatics analysis successful: {result}")
    except Exception as e:
        print(f"✗ Bioinformatics analysis failed: {e}")
    
    # Test API endpoints (simulated)
    print("\n4. Testing API endpoint simulation...")
    try:
        # Simulate API call
        async def simulate_api_call():
            return {
                "success": True,
                "result": {
                    "tool": "sequence_alignment",
                    "status": "completed",
                    "data": "Test alignment result"
                }
            }
        
        result = await simulate_api_call()
        print(f"✓ API endpoint simulation successful: {result}")
    except Exception as e:
        print(f"✗ API endpoint simulation failed: {e}")
    
    print("\n" + "=" * 50)
    print("MCP Integration Test Complete!")

def test_mcp_server_import():
    """Test that the MCP server can be imported."""
    print("\nTesting MCP server import...")
    try:
        import mcp_server
        print("✓ MCP server module imported successfully")
        return True
    except ImportError as e:
        print(f"✗ MCP server import failed: {e}")
        return False
    except Exception as e:
        print(f"✗ Unexpected error importing MCP server: {e}")
        return False

def test_dependencies():
    """Test that all required dependencies are available."""
    print("\nTesting dependencies...")
    
    dependencies = [
        "mcp",
        "fastapi", 
        "uvicorn",
        "pandas",
        "biopython",
        "pymsaviz"
    ]
    
    missing_deps = []
    
    for dep in dependencies:
        try:
            __import__(dep)
            print(f"✓ {dep} available")
        except ImportError:
            print(f"✗ {dep} missing")
            missing_deps.append(dep)
    
    if missing_deps:
        print(f"\nMissing dependencies: {missing_deps}")
        print("Install them with: pip install " + " ".join(missing_deps))
        return False
    else:
        print("\n✓ All dependencies available")
        return True

async def main():
    """Main test function."""
    print("Helix.AI MCP Server Integration Test")
    print("=" * 50)
    
    # Test dependencies first
    deps_ok = test_dependencies()
    if not deps_ok:
        print("\nPlease install missing dependencies before running tests.")
        return
    
    # Test MCP server import
    mcp_ok = test_mcp_server_import()
    if not mcp_ok:
        print("\nMCP server import failed. Check the implementation.")
        return
    
    # Test MCP tools
    await test_mcp_tools()
    
    print("\n" + "=" * 50)
    print("All tests completed!")

if __name__ == "__main__":
    asyncio.run(main()) 