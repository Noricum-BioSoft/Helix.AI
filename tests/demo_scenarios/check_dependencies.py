#!/usr/bin/env python3
"""
Quick dependency checker for demo/eval system.

Checks if all required dependencies are installed before running real execution.
"""

import sys
from pathlib import Path

def check_dependency(module_name, package_name=None):
    """Check if a module can be imported."""
    package_name = package_name or module_name
    try:
        __import__(module_name)
        print(f"✅ {package_name}")
        return True
    except ImportError:
        print(f"❌ {package_name} - NOT INSTALLED")
        return False

def main():
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║           Dependency Check - Helix.AI Demo/Eval System              ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")
    print()
    
    # Check Python version
    print("🐍 Python Version:")
    print(f"   {sys.version}")
    print()
    
    # Check framework dependencies (always needed)
    print("📦 Framework Dependencies (required for mock mode):")
    framework_deps = {
        'yaml': 'pyyaml',
        'pydantic': 'pydantic',
    }
    
    framework_ok = True
    for module, package in framework_deps.items():
        if not check_dependency(module, package):
            framework_ok = False
    print()
    
    # Check backend dependencies (needed for real execution)
    print("🔧 Backend Dependencies (required for real execution):")
    backend_deps = {
        'anthropic': 'anthropic',
        'langchain': 'langchain',
        'langgraph': 'langgraph',
        'boto3': 'boto3',
        'dotenv': 'python-dotenv',
    }
    
    backend_ok = True
    for module, package in backend_deps.items():
        if not check_dependency(module, package):
            backend_ok = False
    print()
    
    # Summary
    print("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    print()
    
    if framework_ok and backend_ok:
        print("✅ All dependencies installed!")
        print()
        print("You can run:")
        print("  • Mock tests: pytest test_scenarios.py -v")
        print("  • Mock demos: python demo_cli.py")
        print("  • Real execution: python test_real_execution.py --scenario real_fastqc_small")
        return 0
    elif framework_ok:
        print("⚠️  Framework dependencies OK, but backend dependencies missing")
        print()
        print("You can run:")
        print("  • Mock tests: pytest test_scenarios.py -v")
        print("  • Mock demos: python demo_cli.py")
        print()
        print("To enable real execution, install backend dependencies:")
        print("  cd /Users/eoberortner/git/Helix.AI")
        print("  pip install -r requirements.txt")
        return 1
    else:
        print("❌ Missing dependencies")
        print()
        print("To install all dependencies:")
        print("  cd /Users/eoberortner/git/Helix.AI")
        print("  pip install -r requirements.txt")
        return 1

if __name__ == "__main__":
    sys.exit(main())
