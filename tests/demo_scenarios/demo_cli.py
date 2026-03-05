#!/usr/bin/env python3
"""
Interactive CLI demo for Helix.AI multi-agent scenarios.

This provides a human-friendly way to:
- Explore available scenarios
- Run scenarios interactively
- See detailed execution traces
- Compare with baselines
- Demonstrate system capabilities

Usage:
    # Interactive mode
    python tests/demo_scenarios/demo_cli.py
    
    # Run specific scenario
    python tests/demo_scenarios/demo_cli.py --scenario ask_basic_question
    
    # Run all scenarios in a category
    python tests/demo_scenarios/demo_cli.py --category ask
    
    # Show diff with baseline
    python tests/demo_scenarios/demo_cli.py --scenario ask_basic_question --show-diff
    
    # Run and save as baseline
    python tests/demo_scenarios/demo_cli.py --scenario ask_basic_question --save-baseline
"""

import asyncio
import sys
import os
from pathlib import Path
import argparse
from typing import Optional, List

# Configure PYTHONPATH to include backend and project root
# This allows the script to import backend modules for real execution without manual configuration
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

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from framework import (
    Scenario,
    ScenarioLoader,
    ScenarioExecutor,
    ScenarioReporter,
    ScenarioCategory,
)


class DemoCLI:
    """Interactive CLI for scenario demos."""
    
    def __init__(self):
        self.scenarios_dir = Path(__file__).parent / "scenarios"
        self.baselines_dir = Path(__file__).parent / "baselines"
        self.loader = ScenarioLoader(self.scenarios_dir)
        self.executor = ScenarioExecutor(
            mock_mode=True,
            baseline_dir=self.baselines_dir,
            verbose=True
        )
        self.reporter = ScenarioReporter(verbose=True)
    
    def list_scenarios(self) -> None:
        """List all available scenarios."""
        print("\n" + "="*70)
        print("Available Scenarios")
        print("="*70 + "\n")
        
        scenarios = self.loader.load_all_scenarios()
        
        if not scenarios:
            print("⚠ No scenarios found in", self.scenarios_dir)
            return
        
        # Group by category
        by_category = {}
        for scenario in scenarios:
            cat = scenario.metadata.category
            by_category.setdefault(cat, []).append(scenario)
        
        for category, cat_scenarios in sorted(by_category.items()):
            print(f"\n{category.value.upper()} ({len(cat_scenarios)} scenarios):")
            for scenario in cat_scenarios:
                tags = f" [{', '.join(scenario.metadata.tags)}]" if scenario.metadata.tags else ""
                print(f"  • {scenario.metadata.id}{tags}")
                print(f"    {scenario.metadata.description}")
    
    async def run_scenario(
        self,
        scenario_id: str,
        show_diff: bool = False,
        save_baseline: bool = False
    ) -> None:
        """Run a specific scenario."""
        try:
            scenario = self.loader.load_scenario(scenario_id)
        except FileNotFoundError as e:
            print(f"❌ Error: {e}")
            return
        
        print(f"\n{'='*70}")
        print(f"Running Scenario: {scenario.metadata.id}")
        print(f"{'='*70}\n")
        print(f"Description: {scenario.metadata.description}")
        print(f"Category: {scenario.metadata.category.value}")
        print(f"User Prompt: \"{scenario.input.user_prompt}\"")
        print(f"\nExecuting...\n")
        
        # Execute
        trace = await self.executor.execute(scenario)
        
        # Show results
        report = self.reporter.report_single(trace, show_contracts=True)
        print(report)
        
        # Show diff if requested
        if show_diff:
            baseline = self.executor.load_baseline(scenario_id)
            if baseline:
                diff = self.reporter.report_diff(trace, baseline)
                print(diff)
            else:
                print(f"⚠ No baseline found for {scenario_id}")
        
        # Save baseline if requested
        if save_baseline:
            self.executor.save_baseline(trace)
            print(f"✓ Saved baseline for {scenario_id}")
    
    async def run_category(self, category: ScenarioCategory) -> None:
        """Run all scenarios in a category."""
        scenarios = self.loader.load_by_category(category)
        
        if not scenarios:
            print(f"⚠ No scenarios found in category: {category.value}")
            return
        
        print(f"\n{'='*70}")
        print(f"Running {category.value.upper()} Scenarios ({len(scenarios)} total)")
        print(f"{'='*70}\n")
        
        traces = []
        for i, scenario in enumerate(scenarios, 1):
            print(f"\n[{i}/{len(scenarios)}] {scenario.metadata.id}")
            trace = await self.executor.execute(scenario)
            traces.append(trace)
        
        # Show summary
        summary = self.reporter.report_batch(traces, group_by="category")
        print(summary)
    
    async def run_all(self) -> None:
        """Run all scenarios."""
        scenarios = self.loader.load_all_scenarios()
        
        if not scenarios:
            print("⚠ No scenarios found")
            return
        
        print(f"\n{'='*70}")
        print(f"Running ALL Scenarios ({len(scenarios)} total)")
        print(f"{'='*70}\n")
        
        traces = []
        for i, scenario in enumerate(scenarios, 1):
            print(f"\n[{i}/{len(scenarios)}] {scenario.metadata.id}")
            trace = await self.executor.execute(scenario)
            traces.append(trace)
        
        # Show summary
        summary = self.reporter.report_batch(traces, group_by="category")
        print(summary)
    
    async def interactive_mode(self) -> None:
        """Run in interactive mode."""
        print("\n" + "="*70)
        print("Helix.AI Multi-Agent Demo")
        print("="*70)
        print("\nInteractive mode - explore agent behaviors\n")
        
        while True:
            print("\nOptions:")
            print("  1. List all scenarios")
            print("  2. Run a scenario")
            print("  3. Run all scenarios in a category")
            print("  4. Run all scenarios")
            print("  5. Exit")
            
            choice = input("\nEnter choice (1-5): ").strip()
            
            if choice == "1":
                self.list_scenarios()
            
            elif choice == "2":
                scenario_id = input("Enter scenario ID: ").strip()
                show_diff = input("Show diff with baseline? (y/n): ").strip().lower() == "y"
                save_baseline = input("Save as baseline? (y/n): ").strip().lower() == "y"
                await self.run_scenario(scenario_id, show_diff, save_baseline)
            
            elif choice == "3":
                print("\nCategories:")
                for cat in ScenarioCategory:
                    print(f"  • {cat.value}")
                cat_input = input("Enter category: ").strip()
                try:
                    category = ScenarioCategory(cat_input)
                    await self.run_category(category)
                except ValueError:
                    print(f"❌ Invalid category: {cat_input}")
            
            elif choice == "4":
                confirm = input("Run all scenarios? This may take a while. (y/n): ").strip().lower()
                if confirm == "y":
                    await self.run_all()
            
            elif choice == "5":
                print("\nExiting. Thanks for using Helix.AI Demo!")
                break
            
            else:
                print("❌ Invalid choice")


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Interactive demo for Helix.AI multi-agent scenarios",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode
  python demo_cli.py
  
  # Run specific scenario
  python demo_cli.py --scenario ask_basic_question
  
  # Run category
  python demo_cli.py --category ask
  
  # Show diff
  python demo_cli.py --scenario ask_basic_question --show-diff
        """
    )
    
    parser.add_argument(
        "--scenario",
        help="Run a specific scenario by ID"
    )
    
    parser.add_argument(
        "--category",
        choices=[c.value for c in ScenarioCategory],
        help="Run all scenarios in a category"
    )
    
    parser.add_argument(
        "--list",
        action="store_true",
        help="List all available scenarios"
    )
    
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all scenarios"
    )
    
    parser.add_argument(
        "--show-diff",
        action="store_true",
        help="Show diff with baseline (requires --scenario)"
    )
    
    parser.add_argument(
        "--save-baseline",
        action="store_true",
        help="Save execution as baseline (requires --scenario)"
    )
    
    args = parser.parse_args()
    
    cli = DemoCLI()
    
    # Dispatch based on arguments
    if args.list:
        cli.list_scenarios()
    
    elif args.scenario:
        asyncio.run(cli.run_scenario(
            args.scenario,
            show_diff=args.show_diff,
            save_baseline=args.save_baseline
        ))
    
    elif args.category:
        category = ScenarioCategory(args.category)
        asyncio.run(cli.run_category(category))
    
    elif args.all:
        asyncio.run(cli.run_all())
    
    else:
        # Interactive mode
        asyncio.run(cli.interactive_mode())


if __name__ == "__main__":
    main()
