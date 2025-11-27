#!/usr/bin/env python3
"""
Test script for RNA-seq preprocessing workflow via Helix.AI backend API.

This script:
1. Loads paired-end FASTQ files
2. Executes the workflow commands from RNASEQ_COMMANDS.txt
3. Validates backend responses
4. Checks session handling
"""

import requests
import json
import sys
from pathlib import Path
from typing import Dict, Any, Optional
import time

# Configuration
BACKEND_URL = "http://localhost:8001"
TIMEOUT = 60  # seconds

# Colors for terminal output
class Colors:
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    RESET = '\033[0m'
    BOLD = '\033[1m'

def print_success(msg: str):
    print(f"{Colors.GREEN}✓ {msg}{Colors.RESET}")

def print_error(msg: str):
    print(f"{Colors.RED}✗ {msg}{Colors.RESET}")

def print_info(msg: str):
    print(f"{Colors.BLUE}ℹ {msg}{Colors.RESET}")

def print_warning(msg: str):
    print(f"{Colors.YELLOW}⚠ {msg}{Colors.RESET}")

def print_step(step_num: int, step_name: str):
    print(f"\n{Colors.BOLD}{'='*60}{Colors.RESET}")
    print(f"{Colors.BOLD}Step {step_num}: {step_name}{Colors.RESET}")
    print(f"{Colors.BOLD}{'='*60}{Colors.RESET}")

def load_fastq_file(filepath: Path) -> str:
    """Load FASTQ file content."""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        print_success(f"Loaded {filepath.name} ({len(content)} characters)")
        return content
    except Exception as e:
        print_error(f"Failed to load {filepath}: {e}")
        sys.exit(1)

def create_session() -> str:
    """Create a new session."""
    try:
        response = requests.post(
            f"{BACKEND_URL}/session/create",
            json={"user_id": None},
            timeout=TIMEOUT
        )
        response.raise_for_status()
        data = response.json()
        session_id = data.get("session_id")
        if session_id:
            print_success(f"Created session: {session_id}")
            return session_id
        else:
            raise ValueError("No session_id in response")
    except Exception as e:
        print_error(f"Failed to create session: {e}")
        sys.exit(1)

def get_session_info(session_id: str) -> Dict[str, Any]:
    """Get session information."""
    try:
        response = requests.get(
            f"{BACKEND_URL}/session/{session_id}",
            timeout=TIMEOUT
        )
        response.raise_for_status()
        return response.json()
    except Exception as e:
        print_warning(f"Failed to get session info: {e}")
        return {}

def execute_command(command: str, session_id: str, files: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
    """
    Execute a natural language command via the /execute endpoint.
    
    Args:
        command: Natural language command
        session_id: Session ID
        files: Optional dict of {filename: content} to include in command
    
    Returns:
        Response data from backend
    """
    # If files are provided, append them to the command
    if files:
        for filename, content in files.items():
            # Detect R1/R2 pattern
            if 'R1' in filename or '_1' in filename or filename.endswith('_1.fastq'):
                command = f"{command}\n\nForward reads ({filename}):\n{content}"
            elif 'R2' in filename or '_2' in filename or filename.endswith('_2.fastq'):
                command = f"{command}\n\nReverse reads ({filename}):\n{content}"
            else:
                command = f"{command}\n\nFile content ({filename}):\n{content}"
    
    print_info(f"Command: {command[:100]}...")
    
    try:
        response = requests.post(
            f"{BACKEND_URL}/execute",
            json={
                "command": command,
                "session_id": session_id
            },
            timeout=TIMEOUT
        )
        response.raise_for_status()
        data = response.json()
        
        # Validate response structure
        if "success" not in data:
            print_warning("Response missing 'success' field")
        if "result" not in data:
            print_warning("Response missing 'result' field")
        if "session_id" not in data:
            print_warning("Response missing 'session_id' field")
        
        return data
    except requests.exceptions.RequestException as e:
        print_error(f"Request failed: {e}")
        if hasattr(e, 'response') and e.response is not None:
            try:
                error_data = e.response.json()
                print_error(f"Error response: {json.dumps(error_data, indent=2)}")
            except:
                print_error(f"Error response text: {e.response.text}")
        raise
    except Exception as e:
        print_error(f"Unexpected error: {e}")
        raise

def validate_trimming_result(result: Dict[str, Any]) -> bool:
    """Validate read trimming result."""
    success = True
    
    if not result.get("success"):
        print_error("Command was not successful")
        success = False
    
    result_data = result.get("result", {})
    
    # Check for nested result structure
    if isinstance(result_data, dict) and "result" in result_data:
        result_data = result_data["result"]
    
    # Check for trimming-specific fields
    if "trimmed_reads" in result_data or "forward_reads" in result_data:
        print_success("Found trimmed reads in result")
    else:
        print_warning("No trimmed reads found in result")
        success = False
    
    if "summary" in result_data:
        summary = result_data["summary"]
        if isinstance(summary, dict):
            if "total_reads" in summary or "forward" in summary:
                print_success("Found summary statistics")
            else:
                print_warning("Summary missing expected fields")
        else:
            print_warning("Summary is not a dict")
    else:
        print_warning("No summary found in result")
    
    return success

def validate_merging_result(result: Dict[str, Any]) -> bool:
    """Validate read merging result."""
    success = True
    
    if not result.get("success"):
        print_error("Command was not successful")
        success = False
    
    result_data = result.get("result", {})
    
    # Check for nested result structure
    if isinstance(result_data, dict) and "result" in result_data:
        result_data = result_data["result"]
    
    # Check for merging-specific fields
    if "merged_sequences" in result_data:
        print_success("Found merged sequences in result")
        merged = result_data["merged_sequences"]
        if merged and len(merged) > 0:
            print_success(f"Merged sequences length: {len(merged)} characters")
        else:
            print_warning("Merged sequences is empty")
            success = False
    else:
        print_error("No merged sequences found in result")
        success = False
    
    if "summary" in result_data:
        summary = result_data["summary"]
        if isinstance(summary, dict):
            if "total_pairs" in summary:
                print_success(f"Found {summary.get('total_pairs', 0)} total pairs")
            if "merged_pairs" in summary:
                print_success(f"Found {summary.get('merged_pairs', 0)} merged pairs")
        else:
            print_warning("Summary is not a dict")
    else:
        print_warning("No summary found in result")
    
    return success

def display_trimming_sequences(result: Dict[str, Any], max_lines: int = 20):
    """Display trimmed sequences from trimming result."""
    result_data = result.get("result", {})
    
    # Check for nested result structure
    if isinstance(result_data, dict) and "result" in result_data:
        result_data = result_data["result"]
    
    print(f"\n{Colors.BOLD}Trimmed Sequences:{Colors.RESET}")
    
    # Handle paired-end results
    if "forward_reads" in result_data and "reverse_reads" in result_data:
        forward_data = result_data["forward_reads"]
        reverse_data = result_data["reverse_reads"]
        
        if isinstance(forward_data, dict) and "trimmed_reads" in forward_data:
            forward_sequences = forward_data["trimmed_reads"]
            print(f"\n{Colors.BLUE}Forward Reads (R1) - First {max_lines} lines:{Colors.RESET}")
            forward_lines = forward_sequences.split('\n')
            lines = forward_lines[:max_lines]
            print('\n'.join(lines))
            if len(forward_lines) > max_lines:
                remaining = len(forward_lines) - max_lines
                print(f"{Colors.YELLOW}... ({remaining} more lines){Colors.RESET}")
        
        if isinstance(reverse_data, dict) and "trimmed_reads" in reverse_data:
            reverse_sequences = reverse_data["trimmed_reads"]
            print(f"\n{Colors.BLUE}Reverse Reads (R2) - First {max_lines} lines:{Colors.RESET}")
            reverse_lines = reverse_sequences.split('\n')
            lines = reverse_lines[:max_lines]
            print('\n'.join(lines))
            if len(reverse_lines) > max_lines:
                remaining = len(reverse_lines) - max_lines
                print(f"{Colors.YELLOW}... ({remaining} more lines){Colors.RESET}")
    
    # Handle single file results
    elif "trimmed_reads" in result_data:
        sequences = result_data["trimmed_reads"]
        print(f"\n{Colors.BLUE}Trimmed Reads - First {max_lines} lines:{Colors.RESET}")
        seq_lines = sequences.split('\n')
        lines = seq_lines[:max_lines]
        print('\n'.join(lines))
        if len(seq_lines) > max_lines:
            remaining = len(seq_lines) - max_lines
            print(f"{Colors.YELLOW}... ({remaining} more lines){Colors.RESET}")
    else:
        print_warning("No trimmed sequences found to display")

def display_merging_sequences(result: Dict[str, Any], max_sequences: int = 5):
    """Display merged sequences from merging result."""
    result_data = result.get("result", {})
    
    # Check for nested result structure
    if isinstance(result_data, dict) and "result" in result_data:
        result_data = result_data["result"]
    
    if "merged_sequences" in result_data:
        merged_sequences = result_data["merged_sequences"]
        print(f"\n{Colors.BOLD}Merged Sequences:{Colors.RESET}")
        
        # Parse FASTA format (sequences start with >)
        sequences = merged_sequences.split('>')
        sequences = [s.strip() for s in sequences if s.strip()]
        
        print(f"\n{Colors.BLUE}Showing first {min(max_sequences, len(sequences))} merged sequences:{Colors.RESET}\n")
        
        for i, seq in enumerate(sequences[:max_sequences], 1):
            lines = seq.split('\n')
            header = lines[0] if lines else ""
            sequence = '\n'.join(lines[1:]) if len(lines) > 1 else ""
            
            print(f"{Colors.GREEN}Sequence {i}:{Colors.RESET}")
            print(f"  Header: {header}")
            if sequence:
                # Show first 80 characters of sequence
                seq_preview = sequence[:80] if len(sequence) > 80 else sequence
                print(f"  Sequence: {seq_preview}")
                if len(sequence) > 80:
                    print(f"  ... ({len(sequence) - 80} more characters)")
            print()
        
        if len(sequences) > max_sequences:
            print(f"{Colors.YELLOW}... ({len(sequences) - max_sequences} more sequences){Colors.RESET}")
        
        # Show summary statistics
        if "summary" in result_data:
            summary = result_data["summary"]
            print(f"\n{Colors.BOLD}Merging Statistics:{Colors.RESET}")
            if isinstance(summary, dict):
                print(f"  Total Pairs: {summary.get('total_pairs', 'N/A')}")
                print(f"  Merged Pairs: {summary.get('merged_pairs', 'N/A')}")
                print(f"  Average Overlap: {summary.get('average_overlap', 0):.2f} bases" if summary.get('average_overlap') else "  Average Overlap: N/A")
                print(f"  Minimum Overlap: {summary.get('min_overlap', 'N/A')} bases")
    else:
        print_warning("No merged sequences found to display")

def display_adapter_removal_sequences(result: Dict[str, Any], max_lines: int = 20):
    """Display sequences after adapter removal."""
    result_data = result.get("result", {})
    
    # Check for nested result structure
    if isinstance(result_data, dict) and "result" in result_data:
        result_data = result_data["result"]
    
    print(f"\n{Colors.BOLD}Sequences After Adapter Removal:{Colors.RESET}")
    
    # Similar structure to trimming
    if "trimmed_reads" in result_data or "forward_reads" in result_data:
        display_trimming_sequences(result, max_lines)
    else:
        print_warning("No sequences found to display")
        print_info("Result structure:")
        print(json.dumps(result_data, indent=2, default=str)[:500])

def display_quality_report(result: Dict[str, Any]):
    """Display quality assessment report."""
    result_data = result.get("result", {})
    
    # Check for nested result structure
    if isinstance(result_data, dict) and "result" in result_data:
        result_data = result_data["result"]
    
    print(f"\n{Colors.BOLD}Quality Assessment Report:{Colors.RESET}")
    
    if "text" in result_data:
        print(f"\n{Colors.BLUE}{result_data['text']}{Colors.RESET}")
    
    if "summary" in result_data or "report" in result_data:
        report_data = result_data.get("summary") or result_data.get("report", {})
        if isinstance(report_data, dict):
            print(f"\n{Colors.BLUE}Report Data:{Colors.RESET}")
            for key, value in report_data.items():
                print(f"  {key}: {value}")
        else:
            print(f"\n{Colors.BLUE}Report:{Colors.RESET}")
            print(str(report_data))
    else:
        print_warning("No quality report data found")
        print_info("Result structure:")
        print(json.dumps(result_data, indent=2, default=str)[:500])

def main():
    """Main test workflow."""
    print(f"{Colors.BOLD}{'='*60}{Colors.RESET}")
    print(f"{Colors.BOLD}RNA-seq Preprocessing Workflow Test{Colors.RESET}")
    print(f"{Colors.BOLD}{'='*60}{Colors.RESET}\n")
    
    # Check backend is running
    print_info("Checking backend connection...")
    try:
        response = requests.get(f"{BACKEND_URL}/health", timeout=5)
        response.raise_for_status()
        print_success("Backend is running")
    except Exception as e:
        print_error(f"Backend is not accessible at {BACKEND_URL}: {e}")
        print_info("Make sure the backend is running: cd backend && uvicorn main_with_mcp:app --reload")
        sys.exit(1)
    
    # Load FASTQ files
    base_dir = Path(__file__).parent.parent.parent  # Go up to project root
    r1_file = base_dir / "data" / "rnaseq_demo" / "sample_R1_realistic.fastq"
    r2_file = base_dir / "data" / "rnaseq_demo" / "sample_R2_realistic.fastq"
    
    print_info("Loading FASTQ files...")
    r1_content = load_fastq_file(r1_file)
    r2_content = load_fastq_file(r2_file)
    
    # Create session
    print_info("Creating session...")
    session_id = create_session()
    
    # Step 1: Quality Trimming
    print_step(1, "Quality Trimming")
    files = {
        "sample_R1_realistic.fastq": r1_content,
        "sample_R2_realistic.fastq": r2_content
    }
    
    command1 = "Trim low-quality bases from my FASTQ reads with quality threshold 20"
    result1 = execute_command(command1, session_id, files)
    
    print_info("Validating trimming result...")
    if validate_trimming_result(result1):
        print_success("Step 1 passed validation")
    else:
        print_error("Step 1 validation failed")
        print_info("Result structure:")
        print(json.dumps(result1, indent=2, default=str))
    
    # Display trimmed sequences
    display_trimming_sequences(result1, max_lines=30)
    
    # Check session after step 1
    session_info = get_session_info(session_id)
    if session_info:
        session_data = session_info.get("session", {})
        history = session_data.get("history", [])
        if history:
            print_success(f"Session has {len(history)} history entries")
        else:
            print_warning("Session history is empty")
    else:
        print_warning("Could not retrieve session info")
    
    # Wait a bit for processing
    time.sleep(1)
    
    # Step 2: Adapter Removal (if supported)
    print_step(2, "Adapter Removal")
    # Use trimmed reads from step 1 (don't pass files, let it retrieve from session history)
    command2 = "Remove adapter sequences AGATCGGAAGAGC from my trimmed reads"
    result2 = execute_command(command2, session_id)  # No files - should use trimmed reads from history
    
    print_info("Validating adapter removal result...")
    if result2.get("success"):
        print_success("Step 2 completed")
        # Display sequences after adapter removal
        display_adapter_removal_sequences(result2, max_lines=30)
    else:
        print_warning("Step 2 may not be fully implemented")
        print_info("Result:")
        print(json.dumps(result2, indent=2, default=str))
    
    # Step 3: Merge Paired-End Reads
    print_step(3, "Merge Paired-End Reads")
    # The merge command should automatically use trimmed reads from session history
    command3 = "Merge my paired-end reads with minimum overlap of 12 bases"
    result3 = execute_command(command3, session_id)
    
    print_info("Validating merging result...")
    if validate_merging_result(result3):
        print_success("Step 3 passed validation")
    else:
        print_error("Step 3 validation failed")
        print_info("Result structure:")
        print(json.dumps(result3, indent=2, default=str))
    
    # Display merged sequences
    display_merging_sequences(result3, max_sequences=10)
    
    # Step 4: Quality Assessment (if supported)
    print_step(4, "Quality Assessment")
    command4 = "Generate a quality report for the merged sequences"
    result4 = execute_command(command4, session_id)
    
    print_info("Validating quality report result...")
    if result4.get("success"):
        print_success("Step 4 completed")
        # Display quality report
        display_quality_report(result4)
    else:
        print_warning("Step 4 may not be fully implemented")
        print_info("Result:")
        print(json.dumps(result4, indent=2, default=str))
    
    # Session Summary
    print_step(5, "Session Summary")
    final_session_info = get_session_info(session_id)
    
    if final_session_info:
        # The API returns {"success": True, "session": {...}, "summary": {...}}
        session_data = final_session_info.get("session", {})
        summary_data = final_session_info.get("summary", {})
        history = session_data.get("history", [])
        
        print(f"\n{Colors.BOLD}Session Overview:{Colors.RESET}")
        print(f"  Session ID: {session_data.get('session_id', 'N/A')}")
        print(f"  Created: {session_data.get('created_at', 'N/A')}")
        print(f"  Updated: {session_data.get('updated_at', 'N/A')}")
        print(f"  Total Operations: {len(history)}")
        
        # Tools used summary
        tools_used = {}
        for entry in history:
            tool = entry.get("tool", "unknown")
            tools_used[tool] = tools_used.get(tool, 0) + 1
        
        print(f"\n{Colors.BOLD}Tools Executed:{Colors.RESET}")
        for tool, count in tools_used.items():
            print(f"  • {tool}: {count} time{'s' if count > 1 else ''}")
        
        # Detailed operation history
        print(f"\n{Colors.BOLD}Operation History:{Colors.RESET}")
        for i, entry in enumerate(history, 1):
            print(f"\n  {Colors.BOLD}Operation {i}:{Colors.RESET}")
            print(f"    Tool: {entry.get('tool', 'unknown')}")
            print(f"    Command: {entry.get('command', 'N/A')[:80]}...")
            print(f"    Timestamp: {entry.get('timestamp', 'N/A')}")
            
            # Show result summary
            result = entry.get("result", {})
            if isinstance(result, dict):
                # Handle nested result structure
                actual_result = result.get("result", result)
                
                if entry.get("tool") == "read_trimming":
                    if "forward_reads" in actual_result and "reverse_reads" in actual_result:
                        fwd_data = actual_result.get("forward_reads", {})
                        rev_data = actual_result.get("reverse_reads", {})
                        fwd_summary = fwd_data.get("summary", {}) if isinstance(fwd_data, dict) else {}
                        rev_summary = rev_data.get("summary", {}) if isinstance(rev_data, dict) else {}
                        print(f"    Input: Paired-end FASTQ files")
                        print(f"    Output: Trimmed reads (Forward: {fwd_summary.get('total_reads', 0)} reads, Reverse: {rev_summary.get('total_reads', 0)} reads)")
                    elif "trimmed_reads" in actual_result:
                        summary = actual_result.get("summary", {})
                        print(f"    Input: FASTQ file")
                        print(f"    Output: {summary.get('total_reads', 0)} trimmed reads, {summary.get('total_bases', 0)} bases")
                
                elif entry.get("tool") == "read_merging":
                    summary = actual_result.get("summary", {})
                    print(f"    Input: Paired-end trimmed reads")
                    print(f"    Output: {summary.get('total_pairs', 0)} read pairs processed, {summary.get('merged_pairs', 0)} merged")
                    if "merged_sequences" in actual_result:
                        merged = actual_result["merged_sequences"]
                        seq_count = merged.count(">") if isinstance(merged, str) else 0
                        print(f"             {seq_count} merged sequences generated")
                
                elif entry.get("tool") == "quality_assessment":
                    metrics = actual_result.get("metrics", {})
                    print(f"    Input: Merged sequences")
                    print(f"    Output: Quality report ({metrics.get('total_sequences', 0)} sequences, {metrics.get('total_bases', 0)} bases)")
                    print(f"             Avg length: {metrics.get('average_length', 0):.1f} bp, GC: {metrics.get('gc_content', 0):.2f}%")
                
                else:
                    # Generic result display
                    status = actual_result.get("status", "unknown")
                    if status == "error":
                        print(f"    Status: {Colors.RED}Error{Colors.RESET}")
                        print(f"    Message: {actual_result.get('message', 'N/A')}")
                    else:
                        print(f"    Status: {Colors.GREEN}Success{Colors.RESET}")
                        # Show key result fields
                        result_keys = [k for k in actual_result.keys() if k not in ["status", "text"]]
                        if result_keys:
                            print(f"    Output keys: {', '.join(result_keys[:5])}")
        
        # Session data structure
        print(f"\n{Colors.BOLD}Session Data Structure:{Colors.RESET}")
        print(f"  Session keys: {', '.join(session_data.keys())}")
        print(f"  Results stored: {len(session_data.get('results', {}))} result entries")
        metadata = session_data.get('metadata', {})
        if metadata:
            print(f"  Metadata keys: {', '.join(metadata.keys())}")
        else:
            print(f"  Metadata: empty")
        
        # Summary statistics
        if summary_data:
            print(f"\n{Colors.BOLD}Summary Statistics:{Colors.RESET}")
            print(f"  Total operations: {summary_data.get('total_operations', 0)}")
            tool_usage = summary_data.get('tool_usage', {})
            if tool_usage:
                print(f"  Tool usage breakdown:")
                for tool, count in tool_usage.items():
                    print(f"    • {tool}: {count}")
            available_results = summary_data.get('available_results', [])
            if available_results:
                print(f"  Available results: {', '.join(available_results)}")
        
        # Show raw session data (first level only for readability)
        print(f"\n{Colors.BOLD}Raw Session Data (Structure):{Colors.RESET}")
        print(f"  Session ID: {session_data.get('session_id', 'N/A')}")
        print(f"  User ID: {session_data.get('user_id', 'N/A')}")
        print(f"  History entries: {len(history)}")
        print(f"  Results dict keys: {list(session_data.get('results', {}).keys())}")
        print(f"  Metadata: {json.dumps(metadata, indent=4, default=str)[:200]}...")
    else:
        print_error("Could not retrieve session information")
    
    # Summary
    print(f"\n{Colors.BOLD}{'='*60}{Colors.RESET}")
    print(f"{Colors.BOLD}Workflow Test Summary{Colors.RESET}")
    print(f"{Colors.BOLD}{'='*60}{Colors.RESET}")
    
    results_summary = []
    if result1.get("success"):
        results_summary.append("✓ Step 1: Quality Trimming")
    else:
        results_summary.append("✗ Step 1: Quality Trimming")
    
    if result2.get("success"):
        results_summary.append("✓ Step 2: Adapter Removal")
    else:
        results_summary.append("⚠ Step 2: Adapter Removal (may not be implemented)")
    
    if result3.get("success"):
        results_summary.append("✓ Step 3: Read Merging")
    else:
        results_summary.append("✗ Step 3: Read Merging")
    
    if result4.get("success"):
        results_summary.append("✓ Step 4: Quality Assessment")
    else:
        results_summary.append("⚠ Step 4: Quality Assessment (may not be implemented)")
    
    for summary in results_summary:
        print(f"  {summary}")
    
    print(f"\n{Colors.BOLD}Session Information:{Colors.RESET}")
    print(f"  Session ID: {session_id}")
    print(f"  Backend URL: {BACKEND_URL}")
    session_data = final_session_info.get("session", {}) if final_session_info else {}
    history = session_data.get("history", [])
    print(f"  History Entries: {len(history)}")
    
    print(f"\n{Colors.BOLD}{'='*60}{Colors.RESET}")

if __name__ == "__main__":
    main()

