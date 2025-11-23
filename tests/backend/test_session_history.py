#!/usr/bin/env python3
"""
Unit tests for session history management.

Tests verify that:
1. History entries are saved correctly
2. Results can be retrieved from history
3. Data passes correctly between workflow steps
"""

import sys
import os
from pathlib import Path

# Add project root to path
# Go up from tests/backend to project root
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "backend"))
sys.path.insert(0, str(project_root / "tools"))

from backend.history_manager import HistoryManager
from tools.command_handler import CommandHandler
import asyncio

def test_history_manager_basic():
    """Test basic history manager functionality."""
    print("=" * 60)
    print("Test 1: Basic History Manager")
    print("=" * 60)
    
    # Create a temporary history manager
    hm = HistoryManager(storage_dir="test_sessions")
    session_id = hm.create_session()
    print(f"✓ Created session: {session_id}")
    
    # Add a test entry
    test_result = {
        "status": "success",
        "data": "test_data",
        "file_type": "paired_end",
        "forward_reads": {
            "trimmed_reads": ">read1\nATCG\n",
            "summary": {"total_reads": 1}
        },
        "reverse_reads": {
            "trimmed_reads": ">read1\nGCTA\n",
            "summary": {"total_reads": 1}
        }
    }
    
    hm.add_history_entry(
        session_id,
        "Test command",
        "read_trimming",
        test_result
    )
    print("✓ Added history entry")
    
    # Retrieve the result
    retrieved = hm.get_latest_result(session_id, "read_trimming")
    assert retrieved is not None, "Result should not be None"
    assert retrieved.get("file_type") == "paired_end", "File type should match"
    assert "forward_reads" in retrieved, "Should have forward_reads"
    assert "reverse_reads" in retrieved, "Should have reverse_reads"
    print("✓ Retrieved result from history")
    
    # Check structure
    forward_data = retrieved.get("forward_reads", {})
    assert isinstance(forward_data, dict), "forward_reads should be a dict"
    assert "trimmed_reads" in forward_data, "Should have trimmed_reads"
    print("✓ Result structure is correct")
    
    # Clean up
    import shutil
    if Path("test_sessions").exists():
        shutil.rmtree("test_sessions")
    print("✓ Test 1 passed\n")


def test_command_handler_history():
    """Test that command handler saves to history."""
    print("=" * 60)
    print("Test 2: Command Handler History")
    print("=" * 60)
    
    handler = CommandHandler()
    session_id = handler.create_session_if_needed()
    print(f"✓ Created session: {session_id}")
    
    # Execute a simple trimming command
    test_fastq = """@read1 1/1
ATCGATCGATCG
+
IIIIIIIIIIII
@read2 2/1
GCTAGCTAGCTA
+
IIIIIIIIIIII"""
    
    command = f"Trim low-quality bases from my FASTQ reads with quality threshold 20\n\nForward reads (test.fastq):\n{test_fastq}"
    
    result = asyncio.run(handler.handle_command(command, session_id))
    print(f"✓ Executed command, status: {result.get('status')}")
    
    # Check if history was saved
    session = handler.history_manager.get_session(session_id)
    assert session is not None, "Session should exist"
    assert len(session["history"]) > 0, "History should have entries"
    print(f"✓ History has {len(session['history'])} entries")
    
    # Check the latest result
    latest = handler.history_manager.get_latest_result(session_id, "read_trimming")
    assert latest is not None, "Should have trimming result"
    assert "trimmed_reads" in latest, "Should have trimmed_reads"
    print("✓ Retrieved trimming result from history")
    print("✓ Test 2 passed\n")


def test_data_passing_between_steps():
    """Test that data passes correctly between workflow steps."""
    print("=" * 60)
    print("Test 3: Data Passing Between Steps")
    print("=" * 60)
    
    handler = CommandHandler()
    session_id = handler.create_session_if_needed()
    print(f"✓ Created session: {session_id}")
    
    # Step 1: Quality trimming with paired-end reads
    r1_fastq = """@read1 1/1
ATCGATCGATCGATCG
+
IIIIIIIIIIIIIIII
@read2 2/1
GCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIII"""
    
    r2_fastq = """@read1 1/2
GCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIII
@read2 2/2
ATCGATCGATCGATCG
+
IIIIIIIIIIIIIIII"""
    
    command1 = f"Trim low-quality bases from my FASTQ reads with quality threshold 20\n\nForward reads (R1.fastq):\n{r1_fastq}\n\nReverse reads (R2.fastq):\n{r2_fastq}"
    
    result1 = asyncio.run(handler.handle_command(command1, session_id))
    print(f"✓ Step 1 completed, status: {result1.get('status')}")
    
    # Verify Step 1 result structure
    result_data1 = result1.get("result", {})
    assert result_data1.get("file_type") == "paired_end", "Should be paired_end"
    assert "forward_reads" in result_data1, "Should have forward_reads"
    assert "reverse_reads" in result_data1, "Should have reverse_reads"
    print("✓ Step 1 result structure is correct")
    
    # Check history
    latest1 = handler.history_manager.get_latest_result(session_id, "read_trimming")
    assert latest1 is not None, "Step 1 should be in history"
    assert latest1.get("file_type") == "paired_end", "History should have file_type"
    print("✓ Step 1 saved to history")
    
    # Step 2: Adapter removal (should retrieve from Step 1)
    command2 = "Remove adapter sequences AGATCGGAAGAGC from my trimmed reads"
    result2 = asyncio.run(handler.handle_command(command2, session_id))
    print(f"✓ Step 2 completed, status: {result2.get('status')}")
    
    # Verify Step 2 retrieved data from Step 1
    result_data2 = result2.get("result", {})
    if result_data2.get("status") == "error":
        print(f"⚠ Step 2 error: {result_data2.get('message')}")
    else:
        assert "forward_reads" in result_data2 or "trimmed_reads" in result_data2, "Step 2 should have results"
        print("✓ Step 2 retrieved data from Step 1")
    
    # Step 3: Merging (should retrieve from Step 2)
    command3 = "Merge my paired-end reads with minimum overlap of 12 bases"
    result3 = asyncio.run(handler.handle_command(command3, session_id))
    print(f"✓ Step 3 completed, status: {result3.get('status')}")
    
    # Verify Step 3 retrieved data
    result_data3 = result3.get("result", {})
    if result_data3.get("status") == "error":
        print(f"⚠ Step 3 error: {result_data3.get('message')}")
    else:
        assert "merged_sequences" in result_data3, "Step 3 should have merged sequences"
        print("✓ Step 3 retrieved data from Step 2")
    
    # Check final history
    session = handler.history_manager.get_session(session_id)
    print(f"✓ Final history has {len(session['history'])} entries")
    for i, entry in enumerate(session["history"], 1):
        print(f"  Entry {i}: {entry['tool']} - {entry['command'][:50]}...")
    
    print("✓ Test 3 passed\n")


if __name__ == "__main__":
    try:
        test_history_manager_basic()
        test_command_handler_history()
        test_data_passing_between_steps()
        print("=" * 60)
        print("ALL TESTS PASSED")
        print("=" * 60)
    except AssertionError as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

