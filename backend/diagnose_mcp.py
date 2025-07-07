#!/usr/bin/env python3
"""
Diagnostic script for MCP server issues
"""

import sys
import subprocess
import importlib
from pathlib import Path

def check_python_version():
    """Check Python version."""
    print(f"🐍 Python version: {sys.version}")
    if sys.version_info < (3, 8):
        print("❌ Python 3.8+ required")
        return False
    print("✅ Python version OK")
    return True

def check_mcp_installation():
    """Check if MCP is installed."""
    try:
        import mcp
        print("✅ MCP package found")
        
        # Check specific MCP modules
        try:
            from mcp.server import Server
            print("✅ mcp.server.Server available")
        except ImportError as e:
            print(f"❌ mcp.server.Server not available: {e}")
            return False
            
        try:
            from mcp.server.stdio import stdio_server
            print("✅ mcp.server.stdio available")
        except ImportError as e:
            print(f"❌ mcp.server.stdio not available: {e}")
            return False
            
        try:
            from mcp.types import Tool, CallToolResult
            print("✅ mcp.types available")
        except ImportError as e:
            print(f"❌ mcp.types not available: {e}")
            return False
            
        return True
        
    except ImportError:
        print("❌ MCP package not installed")
        print("   Install with: pip install mcp>=1.0.0")
        return False

def check_bioinformatics_tools():
    """Check if bioinformatics tools are available."""
    tools_path = str((Path(__file__).resolve().parent.parent / "tools").resolve())
    sys.path.insert(0, tools_path)
    
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
        # Test function access
        alignment.run_alignment
        bio.align_and_visualize_fasta
        mutations.run_mutation
        data_science.analyze_basic_stats
        print("✅ All tool functions accessible")
    except AttributeError as e:
        print(f"❌ Tool function access failed: {e}")
        return False
    
    return True

def check_dependencies():
    """Check if all required dependencies are installed."""
    required_packages = [
        "pandas",
        "numpy", 
        "matplotlib",
        "seaborn",
        "biopython"
    ]
    
    print("📦 Checking required packages...")
    
    for package in required_packages:
        try:
            importlib.import_module(package)
            print(f"✅ {package} available")
        except ImportError:
            print(f"❌ {package} not available")
            return False
    
    return True

def test_simple_mcp_server():
    """Test if the simple MCP server can be imported."""
    try:
        print("🧪 Testing simple MCP server import...")
        
        # Add current directory to path
        sys.path.insert(0, str(Path(__file__).parent))
        
        # Try to import the simple MCP server
        import simple_mcp_server
        print("✅ simple_mcp_server import successful")
        
        # Test server creation
        from simple_mcp_server import server
        print("✅ MCP server object created successfully")
        
        return True
        
    except Exception as e:
        print(f"❌ Simple MCP server test failed: {e}")
        return False

def main():
    """Run all diagnostics."""
    print("🔍 MCP Server Diagnostics")
    print("=" * 40)
    
    checks = [
        ("Python Version", check_python_version),
        ("MCP Installation", check_mcp_installation),
        ("Bioinformatics Tools", check_bioinformatics_tools),
        ("Dependencies", check_dependencies),
        ("Simple MCP Server", test_simple_mcp_server)
    ]
    
    all_passed = True
    
    for check_name, check_func in checks:
        print(f"\n📋 {check_name}:")
        try:
            if check_func():
                print(f"✅ {check_name} passed")
            else:
                print(f"❌ {check_name} failed")
                all_passed = False
        except Exception as e:
            print(f"❌ {check_name} error: {e}")
            all_passed = False
    
    print("\n" + "=" * 40)
    if all_passed:
        print("🎉 All checks passed! MCP server should work.")
        print("\n💡 If the server still doesn't start, try:")
        print("   1. Restart the terminal")
        print("   2. Run: pip install --upgrade mcp")
        print("   3. Check if any antivirus is blocking the process")
    else:
        print("❌ Some checks failed. Please fix the issues above.")
        print("\n💡 Common solutions:")
        print("   1. Install MCP: pip install mcp>=1.0.0")
        print("   2. Install dependencies: pip install -r requirements.txt")
        print("   3. Check Python version: python --version")

if __name__ == "__main__":
    main() 