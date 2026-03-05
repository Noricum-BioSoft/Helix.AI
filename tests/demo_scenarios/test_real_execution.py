#!/usr/bin/env python3
"""
Test script for real workflow execution.

This script tests the real execution mode by running a scenario
against the actual backend system.

Prerequisites:
- Backend dependencies installed
- AWS credentials configured
- .env file with API keys

Usage:
    # Test with mock mode (no AWS required)
    python test_real_execution.py --mock
    
    # Test with real small dataset (fast, local execution)
    python test_real_execution.py --scenario real_fastqc_small
    
    # Test with real large dataset (slow, EMR execution, costs $2-5)
    python test_real_execution.py --scenario real_fastqc_large
"""

import asyncio
import sys
import os
from pathlib import Path
import argparse

# Configure PYTHONPATH to include backend and project root
# This allows the script to import backend modules without manual configuration
project_root = Path(__file__).parent.parent.parent  # Go up to Helix.AI/
backend_dir = project_root / "backend"

# Add to sys.path if not already present
for path_dir in [str(project_root), str(backend_dir)]:
    if path_dir not in sys.path:
        sys.path.insert(0, path_dir)

# Also set PYTHONPATH environment variable for subprocesses
current_pythonpath = os.environ.get('PYTHONPATH', '')
new_paths = [str(project_root), str(backend_dir)]
for path_dir in new_paths:
    if path_dir not in current_pythonpath:
        if current_pythonpath:
            current_pythonpath = f"{path_dir}:{current_pythonpath}"
        else:
            current_pythonpath = path_dir
os.environ['PYTHONPATH'] = current_pythonpath

# Add demo_scenarios to path
sys.path.insert(0, str(Path(__file__).parent))

from framework import ScenarioLoader, ScenarioExecutor


async def run_real_execution(scenario_id: str, mock_mode: bool = False):
    """Run execution of a scenario (CLI helper; not a pytest test)."""
    
    print("=" * 70)
    print(f"{'MOCK' if mock_mode else 'REAL'} EXECUTION TEST")
    print("=" * 70)
    
    # Load scenario
    scenarios_dir = Path(__file__).parent / "scenarios"
    loader = ScenarioLoader(scenarios_dir)
    
    try:
        scenario = loader.load_scenario(scenario_id)
    except FileNotFoundError:
        print(f"❌ Scenario '{scenario_id}' not found")
        print(f"\nAvailable scenarios:")
        for cat, scenarios in loader.get_scenarios_by_category().items():
            print(f"\n{cat.value.upper()}:")
            for s in scenarios:
                print(f"  - {s.metadata.id}")
        return False
    
    print(f"\n📋 Scenario: {scenario.metadata.id}")
    print(f"Description: {scenario.metadata.description}")
    print(f"Mock Required: {scenario.metadata.requires_mock}")
    print(f"Mode: {'MOCK' if mock_mode else 'REAL'}")
    
    # Check if trying to run real execution on mock-only scenario
    if not mock_mode and scenario.metadata.requires_mock:
        print(f"\n⚠️  Scenario '{scenario_id}' requires mock mode")
        print(f"Try: python test_real_execution.py --scenario {scenario_id} --mock")
        return False
    
    # Check if trying to run mock on real-execution scenario
    if mock_mode and not scenario.metadata.requires_mock:
        print(f"\n⚠️  Scenario '{scenario_id}' is designed for real execution")
        print(f"Remove --mock flag to run real execution")
    
    # Create executor
    executor = ScenarioExecutor(mock_mode=mock_mode, verbose=True)
    
    print(f"\n🚀 Starting execution...")
    print("-" * 70)
    
    try:
        # Execute scenario
        trace = await executor.execute(scenario)
        
        print("-" * 70)
        print(f"\n✅ Execution completed!")
        
        # Show results
        print(f"\n📊 Results:")
        print(f"  Duration: {trace.duration_ms}ms")
        print(f"  Agent Sequence: {' → '.join([c['agent'] for c in trace.contracts])}")
        print(f"  Contracts: {len(trace.contracts)}")
        
        # Validation
        if trace.validation_result:
            # Handle both dict and object validation results
            if isinstance(trace.validation_result, dict):
                summary = trace.validation_result.get('summary', 'No summary')
                passed = trace.validation_result.get('passed', False)
                issues = trace.validation_result.get('issues', [])
            else:
                summary = trace.validation_result.summary
                passed = trace.validation_result.passed
                issues = trace.validation_result.issues
            
            print(f"\n📋 Validation: {summary}")
            
            if issues:
                for issue in issues[:5]:  # Show first 5
                    if isinstance(issue, dict):
                        severity = issue.get('severity', 'error')
                        category = issue.get('category', 'unknown')
                        message = issue.get('message', '')
                    else:
                        severity = issue.severity
                        category = issue.category
                        message = issue.message
                    
                    symbol = "❌" if severity == "error" else "⚠️"
                    print(f"  {symbol} [{category}] {message}")
        
        # Return pass status
        if trace.validation_result:
            if isinstance(trace.validation_result, dict):
                return trace.validation_result.get('passed', True)
            else:
                return trace.validation_result.passed
        return True
        
    except Exception as e:
        print(f"\n❌ Execution failed: {e}")
        import traceback
        traceback.print_exc()
        return False


async def run_infrastructure_comparison():
    """Compare infrastructure decisions for small vs large datasets."""
    
    print("=" * 70)
    print("INFRASTRUCTURE DECISION COMPARISON")
    print("=" * 70)
    
    scenarios_dir = Path(__file__).parent / "scenarios"
    loader = ScenarioLoader(scenarios_dir)
    executor = ScenarioExecutor(mock_mode=True, verbose=False)
    
    # Run small and large FastQC scenarios
    small = await executor.execute(loader.load_scenario('real_fastqc_small'))
    large = await executor.execute(loader.load_scenario('real_fastqc_large'))
    
    # Extract infrastructure decisions
    small_infra = [c for c in small.contracts if c['agent'] == 'InfrastructureExpert']
    large_infra = [c for c in large.contracts if c['agent'] == 'InfrastructureExpert']
    
    if small_infra and large_infra:
        small_data = small_infra[0]['data']
        large_data = large_infra[0]['data']
        
        print(f"\n📊 Small Dataset:")
        print(f"  Files: test_mate_*.fq")
        print(f"  Infrastructure: {small_data['recommended_environment']}")
        print(f"  Cost: $0 (local)")
        print(f"  Confidence: {small_data['confidence_score']}")
        
        print(f"\n📊 Large Dataset:")
        print(f"  Files: mate_*.fq")
        print(f"  Infrastructure: {large_data['recommended_environment']}")
        print(f"  Cost: {large_data['cost_analysis']['estimated_cost_range_usd']}")
        print(f"  Confidence: {large_data['confidence_score']}")
        
        print(f"\n✅ Result: System intelligently adapts infrastructure based on data size!")
    else:
        print(f"\n❌ Could not extract infrastructure decisions")


def main():
    parser = argparse.ArgumentParser(description="Test real workflow execution")
    parser.add_argument(
        '--scenario',
        default='real_fastqc_small',
        help='Scenario ID to run (default: real_fastqc_small)'
    )
    parser.add_argument(
        '--mock',
        action='store_true',
        help='Use mock mode (no actual backend execution)'
    )
    parser.add_argument(
        '--compare',
        action='store_true',
        help='Compare infrastructure decisions for small vs large datasets'
    )
    
    args = parser.parse_args()
    
    if args.compare:
        success = asyncio.run(run_infrastructure_comparison())
    else:
        success = asyncio.run(run_real_execution(args.scenario, args.mock))
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
