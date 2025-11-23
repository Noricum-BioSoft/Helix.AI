#!/usr/bin/env python3
"""
Test script for Sequence Alignment and Phylogenetic Analysis workflow via Helix.AI backend API.

This script demonstrates Workflow B from NGS101_DEMONSTRATION_STRATEGY.md:
- Align RNA sequences
- Build phylogenetic tree
- Select representative sequences from each cluster

Commands:
1. "Align these RNA sequences"
2. "Build phylogenetic tree"
3. "Select 10 representative sequences from each cluster"
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

def execute_command(command: str, session_id: str) -> Dict[str, Any]:
    """
    Execute a natural language command via the /execute endpoint.
    
    Args:
        command: Natural language command
        session_id: Session ID
    
    Returns:
        Response data from backend
    """
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

def generate_test_sequences() -> str:
    """Generate test RNA sequences for alignment."""
    # Create diverse RNA sequences that will form clusters
    sequences = [
        "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",  # Identical to first
        "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCT",  # One mutation
        "GTGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",  # Different start
        "GTGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCT",  # Different start + mutation
        "CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",  # Different sequence
        "CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",  # Identical to previous
        "CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCA",  # One mutation
        "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",  # Another variant
        "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCT",  # Variant with mutation
    ]
    
    # Format as FASTA
    fasta_lines = []
    for i, seq in enumerate(sequences, 1):
        fasta_lines.append(f">RNA_seq_{i}")
        fasta_lines.append(seq)
    
    return "\n".join(fasta_lines)

def display_alignment_result(result: Dict[str, Any], max_lines: int = 30):
    """Display alignment result."""
    result_data = result.get("result", {})
    
    # Check for nested result structure
    if isinstance(result_data, dict) and "result" in result_data:
        result_data = result_data["result"]
    
    print(f"\n{Colors.BOLD}Alignment Result:{Colors.RESET}")
    
    if "alignment" in result_data:
        alignment = result_data["alignment"]
        if isinstance(alignment, list) and len(alignment) > 0:
            print_success(f"Aligned {len(alignment)} sequences")
            print(f"\nFirst {min(3, len(alignment))} aligned sequences:")
            for i, seq in enumerate(alignment[:3], 1):
                name = seq.get("name", f"seq_{i}")
                sequence = seq.get("sequence", "")
                print(f"  {name}: {sequence[:80]}...")
        else:
            print_warning("Alignment is empty or invalid")
    else:
        print_warning("No alignment found in result")
    
    if "statistics" in result_data:
        stats = result_data["statistics"]
        print(f"\nAlignment Statistics:")
        for key, value in stats.items():
            print(f"  {key}: {value}")

def display_tree_result(result: Dict[str, Any]):
    """Display phylogenetic tree result."""
    result_data = result.get("result", {})
    
    # Check for nested result structure
    if isinstance(result_data, dict) and "result" in result_data:
        result_data = result_data["result"]
    
    print(f"\n{Colors.BOLD}Phylogenetic Tree Result:{Colors.RESET}")
    
    # Try multiple paths to find tree_newick
    newick = None
    if "tree_newick" in result_data:
        newick = result_data["tree_newick"]
    elif "output" in result_data and isinstance(result_data["output"], dict):
        newick = result_data["output"].get("tree_newick")
    
    if newick:
        print_success("Tree Newick format generated")
        print(f"Newick: {newick[:200]}...")
    else:
        print_warning("No tree_newick found in result")
        print_info(f"Available keys in result: {list(result_data.keys())}")
    
    # Check for ETE3 visualization
    ete = None
    if "ete_visualization" in result_data:
        ete = result_data["ete_visualization"]
    elif "output" in result_data and isinstance(result_data["output"], dict):
        ete = result_data["output"].get("ete_visualization")
    
    if ete and isinstance(ete, dict):
        if "svg" in ete:
            print_success("ETE3 SVG visualization generated")
        elif "error" in ete:
            print_warning(f"ETE3 visualization error: {ete['error']}")
        else:
            print_info("ETE3 visualization data present")
    
    # Show text if available
    if "text" in result_data:
        print(f"\nResult text: {result_data['text'][:200]}...")

def display_selection_result(result: Dict[str, Any]):
    """Display sequence selection result."""
    result_data = result.get("result", {})
    
    # Check for nested result structure
    if isinstance(result_data, dict) and "result" in result_data:
        result_data = result_data["result"]
    
    print(f"\n{Colors.BOLD}Sequence Selection Result:{Colors.RESET}")
    
    # Try multiple paths to find selected sequences/representatives
    selected = None
    representatives = None
    
    # Check for selected_sequences (from sequence_selection tool)
    if "selected_sequences" in result_data:
        selected = result_data["selected_sequences"]
    elif "output" in result_data and isinstance(result_data["output"], list):
        selected = result_data["output"]
    elif "output" in result_data and isinstance(result_data["output"], dict):
        selected = result_data["output"].get("selected_sequences")
    
    # Check for representatives (from clustering_analysis tool)
    if "clustering_result" in result_data:
        clustering = result_data["clustering_result"]
        if isinstance(clustering, dict):
            representatives = clustering.get("representatives", [])
            print_success(f"Clustering analysis completed")
            print(f"  Number of clusters: {clustering.get('num_clusters', 'N/A')}")
            print(f"  Total sequences: {clustering.get('total_sequences', 'N/A')}")
    
    # Display selected sequences or representatives
    if selected and isinstance(selected, list):
        print_success(f"Selected {len(selected)} representative sequences")
        for i, seq in enumerate(selected[:5], 1):
            if isinstance(seq, dict):
                name = seq.get("name", f"seq_{i}")
                sequence = seq.get("sequence", "")
                print(f"  {i}. {name}: {sequence[:60]}...")
            elif isinstance(seq, str):
                print(f"  {i}. {seq[:60]}...")
    elif representatives and isinstance(representatives, list):
        print_success(f"Found {len(representatives)} representative sequences from clusters")
        for i, rep in enumerate(representatives[:5], 1):
            if isinstance(rep, dict):
                name = rep.get("name", f"rep_{i}")
                print(f"  {i}. {name}")
            elif isinstance(rep, str):
                print(f"  {i}. {rep}")
    else:
        print_warning("No selected_sequences or representatives found in result")
        print_info(f"Available keys in result: {list(result_data.keys())}")
        
        # Check for error message
        if "error" in result_data:
            print_error(f"Error in clustering: {result_data['error']}")
        elif "text" in result_data and "error" in result_data["text"].lower():
            print_error(f"Error message: {result_data['text']}")
        
        if "clustering_result" in result_data:
            clustering = result_data["clustering_result"]
            if isinstance(clustering, dict):
                if "error" in clustering:
                    print_error(f"Clustering error: {clustering['error']}")
                else:
                    print_info(f"Clustering result keys: {list(clustering.keys())}")
    
    # Show summary if available
    if "summary" in result_data:
        summary = result_data["summary"]
        if isinstance(summary, dict):
            print(f"\nSelection Summary:")
            for key, value in summary.items():
                print(f"  {key}: {value}")
    
    # Show text if available
    if "text" in result_data:
        print(f"\nResult text: {result_data['text'][:200]}...")

def main():
    """Main test workflow."""
    print(f"\n{Colors.BOLD}{'='*60}{Colors.RESET}")
    print(f"{Colors.BOLD}Sequence Alignment and Phylogenetic Analysis Workflow Test{Colors.RESET}")
    print(f"{Colors.BOLD}{'='*60}{Colors.RESET}")
    
    # Create session
    session_id = create_session()
    
    # Generate test sequences
    test_sequences = generate_test_sequences()
    print_success(f"Generated {len(test_sequences.split('>')) - 1} test sequences")
    
    # Step 1: Align sequences
    print_step(1, "Sequence Alignment")
    command1 = f"Align these RNA sequences:\n\n{test_sequences}"
    result1 = execute_command(command1, session_id)
    
    if result1.get("success"):
        print_success("Step 1 completed successfully")
        display_alignment_result(result1)
    else:
        print_error("Step 1 failed")
        print(json.dumps(result1, indent=2, default=str))
        return
    
    # Step 2: Build phylogenetic tree
    print_step(2, "Build Phylogenetic Tree")
    command2 = "Build phylogenetic tree"
    result2 = execute_command(command2, session_id)
    
    if result2.get("success"):
        print_success("Step 2 completed successfully")
        display_tree_result(result2)
    else:
        print_error("Step 2 failed")
        print(json.dumps(result2, indent=2, default=str))
    
    # Step 3: Select representative sequences
    print_step(3, "Select Representative Sequences")
    command3 = "Select 10 representative sequences from each cluster"
    result3 = execute_command(command3, session_id)
    
    if result3.get("success"):
        print_success("Step 3 completed successfully")
        display_selection_result(result3)
    else:
        print_warning("Step 3 may not be fully implemented")
        print_info("Result:")
        print(json.dumps(result3, indent=2, default=str))
    
    # Session Summary
    print_step(4, "Session Summary")
    final_session_info = get_session_info(session_id)
    
    if final_session_info:
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
                actual_result = result.get("result", result)
                
                if entry.get("tool") == "sequence_alignment":
                    alignment = actual_result.get("alignment", [])
                    print(f"    Input: RNA sequences")
                    print(f"    Output: {len(alignment)} aligned sequences")
                
                elif entry.get("tool") == "phylogenetic_tree":
                    newick = actual_result.get("tree_newick", "")
                    print(f"    Input: Aligned sequences")
                    print(f"    Output: Phylogenetic tree (Newick format: {len(newick)} chars)")
                
                elif entry.get("tool") == "sequence_selection":
                    selected = actual_result.get("selected_sequences", [])
                    print(f"    Input: Aligned sequences / tree")
                    print(f"    Output: {len(selected) if isinstance(selected, list) else 0} representative sequences")
                
                else:
                    status = actual_result.get("status", "unknown")
                    if status == "error":
                        print(f"    Status: {Colors.RED}Error{Colors.RESET}")
                        print(f"    Message: {actual_result.get('message', 'N/A')}")
                    else:
                        print(f"    Status: {Colors.GREEN}Success{Colors.RESET}")
        
        # Summary statistics
        if summary_data:
            print(f"\n{Colors.BOLD}Summary Statistics:{Colors.RESET}")
            print(f"  Total operations: {summary_data.get('total_operations', 0)}")
            tool_usage = summary_data.get('tool_usage', {})
            if tool_usage:
                print(f"  Tool usage breakdown:")
                for tool, count in tool_usage.items():
                    print(f"    • {tool}: {count}")
    else:
        print_error("Could not retrieve session information")
    
    # Final summary
    print(f"\n{Colors.BOLD}{'='*60}{Colors.RESET}")
    print(f"{Colors.BOLD}Workflow Test Summary{Colors.RESET}")
    print(f"{Colors.BOLD}{'='*60}{Colors.RESET}")
    
    results_summary = []
    if result1.get("success"):
        results_summary.append("✓ Step 1: Sequence Alignment")
    else:
        results_summary.append("✗ Step 1: Sequence Alignment")
    
    if result2.get("success"):
        results_summary.append("✓ Step 2: Phylogenetic Tree")
    else:
        results_summary.append("✗ Step 2: Phylogenetic Tree")
    
    if result3.get("success"):
        results_summary.append("✓ Step 3: Representative Selection")
    else:
        results_summary.append("⚠ Step 3: Representative Selection (may not be implemented)")
    
    for summary in results_summary:
        print(f"  {summary}")
    
    print(f"\n{Colors.BOLD}Session Information:{Colors.RESET}")
    print(f"  Session ID: {session_id}")
    print(f"  Backend URL: {BACKEND_URL}")
    
    print(f"\n{Colors.BOLD}{'='*60}{Colors.RESET}")

if __name__ == "__main__":
    main()

