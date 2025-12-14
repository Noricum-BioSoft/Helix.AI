"""
Pytest configuration to ensure project packages are importable during tests.
"""

import sys
from pathlib import Path
import pytest
import os

# Make the repository root importable when running tests from within `tests/`
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Ensure `tests/` itself does NOT shadow top-level packages like `backend/`.
# Pytest (and some environments) may prepend the tests directory to sys.path,
# which would make `tests/backend/` shadow the real `backend/` package.
TESTS_DIR = PROJECT_ROOT / "tests"
tests_dir_str = str(TESTS_DIR)
if tests_dir_str in sys.path:
    sys.path.remove(tests_dir_str)
    sys.path.append(tests_dir_str)

# Default to mock mode in tests unless a test explicitly opts out.
# This prevents optional LLM / cloud dependencies from being imported/executed.
os.environ.setdefault("HELIX_MOCK_MODE", "1")

# Also add tools directory explicitly for tests that import tool modules directly.
TOOLS_PATH = PROJECT_ROOT / "tools"
if str(TOOLS_PATH) not in sys.path:
    # Append so `PROJECT_ROOT` stays earlier on sys.path (so `import backend...` works).
    sys.path.append(str(TOOLS_PATH))


def pytest_configure(config):
    """Register custom pytest marks."""
    config.addinivalue_line(
        "markers", "integration: marks tests as integration tests (requires backend running)"
    )

