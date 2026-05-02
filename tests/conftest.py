"""
Pytest configuration to ensure project packages are importable during tests.
"""

import asyncio
import importlib
import sys
import warnings
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

# Redirect all tool and agent output files to tmp/test-results/ so they never
# land in the repo root or pollute the working tree.
_TEST_OUTPUT_DIR = PROJECT_ROOT / "tmp" / "test-results"
_TEST_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("HELIX_OUTPUT_DIR", str(_TEST_OUTPUT_DIR))
os.environ.setdefault("OUTPUT_DIR", str(_TEST_OUTPUT_DIR))

# Also add tools directory explicitly for tests that import tool modules directly.
TOOLS_PATH = PROJECT_ROOT / "tools"
if str(TOOLS_PATH) not in sys.path:
    # Append so `PROJECT_ROOT` stays earlier on sys.path (so `import backend...` works).
    sys.path.append(str(TOOLS_PATH))

# Keep `tools` bound to this repository's package even if another dependency
# imports a third-party module named `tools` first.
_tools_mod = sys.modules.get("tools")
if _tools_mod is not None:
    _origin = str(getattr(_tools_mod, "__file__", "") or "")
    if str(PROJECT_ROOT / "tools") not in _origin:
        sys.modules.pop("tools", None)
        sys.modules.pop("tools.ncbi_tools", None)

# Eager import makes `patch("tools.ncbi_tools....")` deterministic across suite order.
# Guard with try/except: optional heavy deps (biopython, pysam…) may be absent
# in lightweight CI environments that only run unit/benchmark tests.
importlib.import_module("tools")
try:
    importlib.import_module("tools.ncbi_tools")
except ImportError as _e:
    import warnings
    warnings.warn(
        f"tools.ncbi_tools could not be imported ({_e}). "
        "Tests that patch this module will be skipped or may fail.",
        ImportWarning,
        stacklevel=1,
    )

@pytest.fixture(autouse=True)
def _ensure_default_event_loop():
    """
    Some sync tests call asyncio.get_event_loop().run_until_complete(...).
    On Python 3.13, prior tests may leave no current loop; ensure one exists.
    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            loop = asyncio.get_event_loop_policy().get_event_loop()
        if loop.is_closed():
            raise RuntimeError("event loop was closed")
    except Exception:
        asyncio.get_event_loop_policy().set_event_loop(asyncio.new_event_loop())
    yield


def pytest_configure(config):
    """Register custom pytest marks."""
    config.addinivalue_line(
        "markers", "integration: marks tests as integration tests (requires backend running)"
    )

