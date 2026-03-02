"""
Pytest configuration for demo scenarios.

This file is automatically loaded by pytest and provides custom command-line options.
"""

import sys
from pathlib import Path
import pytest

# Add the demo_scenarios directory to Python path so 'framework' can be imported
demo_scenarios_dir = Path(__file__).parent
if str(demo_scenarios_dir) not in sys.path:
    sys.path.insert(0, str(demo_scenarios_dir))


def pytest_addoption(parser):
    """Add custom pytest options."""
    parser.addoption(
        "--show-traces",
        action="store_true",
        default=False,
        help="Show detailed execution traces"
    )
    parser.addoption(
        "--update-baseline",
        action="store_true",
        default=False,
        help="Update baseline traces for all scenarios"
    )
    parser.addoption(
        "--compare-baseline",
        action="store_true",
        default=False,
        help="Compare current execution with baseline"
    )
