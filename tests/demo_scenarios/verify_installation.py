#!/usr/bin/env python3
"""
Verification script for the demo/eval framework.

Run this to verify everything is set up correctly:
    python tests/demo_scenarios/verify_installation.py
"""

import sys
from pathlib import Path


def verify_imports():
    """Verify all framework imports work."""
    print("Checking framework imports...")
    try:
        from framework import (
            Scenario,
            ScenarioLoader,
            ScenarioExecutor,
            ScenarioReporter,
            ContractValidator,
            ScenarioCategory,
        )
        print("  ✓ All imports successful")
        return True
    except ImportError as e:
        print(f"  ✗ Import error: {e}")
        return False


def verify_scenarios():
    """Verify scenarios can be loaded."""
    print("\nChecking scenario loading...")
    try:
        from framework import ScenarioLoader
        
        scenarios_dir = Path(__file__).parent / "scenarios"
        loader = ScenarioLoader(scenarios_dir)
        scenarios = loader.load_all_scenarios()
        
        print(f"  ✓ Loaded {len(scenarios)} scenarios")
        
        # Count by category
        from framework import ScenarioCategory
        by_category = {}
        for scenario in scenarios:
            cat = scenario.metadata.category
            by_category[cat] = by_category.get(cat, 0) + 1
        
        for cat, count in sorted(by_category.items()):
            print(f"    - {cat.value}: {count}")
        
        return True
    except Exception as e:
        print(f"  ✗ Error loading scenarios: {e}")
        return False


def verify_execution():
    """Verify a simple scenario can be executed."""
    print("\nChecking scenario execution...")
    try:
        import asyncio
        from pathlib import Path
        from framework import ScenarioLoader, ScenarioExecutor
        
        async def run_test():
            loader = ScenarioLoader(Path(__file__).parent / "scenarios")
            scenarios = loader.load_all_scenarios()
            
            if not scenarios:
                print("  ✗ No scenarios found to test")
                return False
            
            # Try to execute first scenario
            scenario = scenarios[0]
            executor = ScenarioExecutor(mock_mode=True, verbose=False)
            trace = await executor.execute(scenario)
            
            print(f"  ✓ Executed test scenario: {scenario.metadata.id}")
            print(f"    - Duration: {trace.duration_ms}ms")
            print(f"    - Agents: {' → '.join(trace.agent_sequence)}")
            print(f"    - Contracts: {len(trace.contracts)}")
            
            if trace.validation_result:
                from framework.validator import ValidationResult
                validation = ValidationResult(**trace.validation_result)
                status = "passed" if validation.passed else "failed"
                print(f"    - Validation: {status}")
            
            return True
        
        return asyncio.run(run_test())
    
    except Exception as e:
        print(f"  ✗ Error executing scenario: {e}")
        import traceback
        traceback.print_exc()
        return False


def verify_contracts():
    """Verify contract imports work."""
    print("\nChecking contract imports...")
    try:
        # These might not be available in standalone mode, so it's optional
        from shared.contracts import (
            IntentResult,
            WorkflowPlan,
            InfraDecisionContract,
            ExecutionResult,
        )
        print("  ✓ Contract imports successful")
        return True
    except ImportError:
        print("  ⚠ Contract imports not available (expected in standalone mode)")
        return True  # This is OK


def main():
    """Run all verification checks."""
    print("="*70)
    print("Helix.AI Demo/Eval Framework Verification")
    print("="*70)
    
    checks = [
        ("Framework Imports", verify_imports),
        ("Scenario Loading", verify_scenarios),
        ("Scenario Execution", verify_execution),
        ("Contract Imports", verify_contracts),
    ]
    
    results = []
    for name, check_fn in checks:
        try:
            results.append(check_fn())
        except Exception as e:
            print(f"\n✗ {name} failed with unexpected error: {e}")
            results.append(False)
    
    print("\n" + "="*70)
    print("Verification Summary")
    print("="*70)
    
    passed = sum(results)
    total = len(results)
    
    print(f"\nPassed: {passed}/{total}")
    
    if all(results):
        print("\n✓ All checks passed! The framework is ready to use.")
        print("\nNext steps:")
        print("  1. Try the CLI demo: python tests/demo_scenarios/demo_cli.py")
        print("  2. Run pytest tests: pytest tests/demo_scenarios/test_scenarios.py -v")
        print("  3. Read QUICK_START.md for more information")
        return 0
    else:
        print("\n✗ Some checks failed. Please review the errors above.")
        print("\nTroubleshooting:")
        print("  - Ensure you're in the Helix.AI project root")
        print("  - Check that dependencies are installed: pip install -r requirements.txt")
        print("  - Verify Python path includes project root")
        return 1


if __name__ == "__main__":
    sys.exit(main())
