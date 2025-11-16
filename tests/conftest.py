"""
Pytest configuration to ensure project packages are importable during tests.
"""

import sys
from pathlib import Path

# Make the repository root importable when running tests from within `tests/`
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


