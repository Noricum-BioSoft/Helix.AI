"""
Unit tests for context builder module.
Tests session-aware context generation from session state.
"""

import sys
from pathlib import Path
import pytest

# Add project root to path (must be before any backend imports)
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Import after path setup
try:
    from backend.context_builder import build_context_snippet
except ImportError:
    # Fallback: try direct import if backend is in path
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "context_builder",
        PROJECT_ROOT / "backend" / "context_builder.py"
    )
    context_builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(context_builder)
    build_context_snippet = context_builder.build_context_snippet


def test_build_context_snippet_empty():
    """Test context builder with empty session context."""
    result = build_context_snippet({})
    assert result == ""


def test_build_context_snippet_mutated_sequences():
    """Test context builder with mutated sequences."""
    session_context = {
        "mutated_sequences": [
            "ATGCGATCG",
            "ATGCGATCC",
            "ATGCGATCA"
        ]
    }
    result = build_context_snippet(session_context)
    
    assert "Context:" in result
    assert "Mutated sequences:" in result
    assert ">mutant_1" in result
    assert ">mutant_2" in result
    assert ">mutant_3" in result
    assert "ATGCGATCG" in result
    assert "ATGCGATCC" in result
    assert "ATGCGATCA" in result


def test_build_context_snippet_aligned_sequences():
    """Test context builder with aligned sequences."""
    session_context = {
        "aligned_sequences": ">seq1 ATGCGATCG\n>seq2 ATGCGATC"
    }
    result = build_context_snippet(session_context)
    
    assert "Context:" in result
    assert "Aligned sequences (FASTA, truncated):" in result
    assert ">seq1" in result
    assert ">seq2" in result
    assert "ATGCGATCG" in result


def test_build_context_snippet_selected_sequences():
    """Test context builder with selected sequences."""
    session_context = {
        "selected_sequences": [
            "ATGCGATCG",
            "ATGCGATCC"
        ]
    }
    result = build_context_snippet(session_context)
    
    assert "Context:" in result
    assert "Selected sequences:" in result
    assert ">selected_1" in result
    assert ">selected_2" in result
    assert "ATGCGATCG" in result


def test_build_context_snippet_uploaded_files():
    """Test context builder with uploaded files."""
    session_context = {
        "uploaded_files": [
            {"name": "test.fasta", "size": 1000},
            {"name": "data.csv", "size": 2000}
        ]
    }
    result = build_context_snippet(session_context)
    
    assert "Context:" in result
    assert "Uploaded files in session:" in result
    assert "test.fasta" in result
    assert "data.csv" in result


def test_build_context_snippet_all_contexts():
    """Test context builder with all context types."""
    session_context = {
        "mutated_sequences": ["ATGCGATCG"],
        "aligned_sequences": ">seq1 ATGCGATCG",
        "selected_sequences": ["ATGCGATCC"],
        "uploaded_files": [{"name": "test.fasta"}]
    }
    result = build_context_snippet(session_context)
    
    assert "Context:" in result
    assert "Mutated sequences:" in result
    assert "Aligned sequences (FASTA, truncated):" in result
    assert "Selected sequences:" in result
    assert "Uploaded files in session:" in result


def test_build_context_snippet_none_values():
    """Test context builder handles None values gracefully."""
    session_context = {
        "mutated_sequences": None,
        "aligned_sequences": None,
        "selected_sequences": None,
        "uploaded_files": None
    }
    result = build_context_snippet(session_context)
    assert result == ""


def test_build_context_snippet_empty_lists():
    """Test context builder handles empty lists."""
    session_context = {
        "mutated_sequences": [],
        "selected_sequences": [],
        "uploaded_files": []
    }
    result = build_context_snippet(session_context)
    assert result == ""

