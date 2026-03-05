#!/usr/bin/env python3
"""
Test script for Variant Analysis Workflow via Helix.AI backend API.

This script demonstrates Workflow C from NGS101_DEMONSTRATION_STRATEGY.md:
- Generate sequence variants (mutations)
- Select diverse variants for further analysis
- Visualize variants in plasmid context

Commands:
1. "Generate 50 variants"
2. "Select 10 most diverse"
3. "Insert into pUC19 plasmid"
4. "Visualize"
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

def generate_test_sequence() -> str:
    """Generate a test DNA sequence for variant generation."""
    # Create a realistic DNA sequence (100bp)
    return "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"

def display_mutation_result(result: Dict[str, Any], max_variants: int = 5):
    """Display mutation/variant generation result."""
    result_data = result.get("result", {})
    
    # Check for nested result structure
    if isinstance(result_data, dict) and "result" in result_data:
        result_data = result_data["result"]
    
    print(f"\n{Colors.BOLD}Variant Generation Result:{Colors.RESET}")
    
    # Try multiple paths to find variants
    variants = None
    
    # Path 1: Direct variants key
    if "variants" in result_data:
        variants = result_data["variants"]
    # Path 2: Variants in statistics
    elif "statistics" in result_data and isinstance(result_data["statistics"], dict):
        variants = result_data["statistics"].get("variants")
    # Path 3: Variants in output
    elif "output" in result_data:
        output = result_data["output"]
        if isinstance(output, dict):
            variants = output.get("variants")
        elif isinstance(output, list):
            variants = output
    
    if variants and isinstance(variants, list):
        print_success(f"Generated {len(variants)} variants")
        print(f"\nFirst {min(max_variants, len(variants))} variants:")
        for i, variant in enumerate(variants[:max_variants], 1):
            if isinstance(variant, str):
                print(f"  {i}. {variant[:60]}...")
            elif isinstance(variant, dict):
                seq = variant.get("sequence", variant.get("variant", ""))
                print(f"  {i}. {seq[:60]}...")
    else:
        print_warning("No variants found in result")
        print_info(f"Available keys in result: {list(result_data.keys())}")
        if "statistics" in result_data:
            print_info(f"Statistics keys: {list(result_data['statistics'].keys()) if isinstance(result_data['statistics'], dict) else 'N/A'}")
    
    # Display statistics
    if "statistics" in result_data:
        stats = result_data["statistics"]
        if isinstance(stats, dict):
            print(f"\nMutation Statistics:")
            for key, value in stats.items():
                if key != "variants":  # Don't print variants again
                    if isinstance(value, (int, float, str, bool)):
                        print(f"  {key}: {value}")
                    elif isinstance(value, list) and len(value) <= 5:
                        print(f"  {key}: {value}")
    elif "output" in result_data:
        output = result_data["output"]
        if isinstance(output, dict):
            print(f"\nOutput Statistics:")
            for key, value in output.items():
                if key != "variants" and isinstance(value, (int, float, str, bool)):
                    print(f"  {key}: {value}")

def display_selection_result(result: Dict[str, Any], max_selected: int = 5):
    """Display variant selection result."""
    result_data = result.get("result", {})
    
    # Check for nested result structure
    if isinstance(result_data, dict) and "result" in result_data:
        result_data = result_data["result"]
    
    print(f"\n{Colors.BOLD}Variant Selection Result:{Colors.RESET}")
    
    # Try multiple paths to find selected variants
    selected = None
    
    # Path 1: Direct selected_variants key
    if "selected_variants" in result_data:
        selected = result_data["selected_variants"]
    # Path 2: Selected variants in output
    elif "output" in result_data:
        output = result_data["output"]
        if isinstance(output, dict):
            selected = output.get("selected_variants")
        elif isinstance(output, list):
            selected = output
    
    if selected and isinstance(selected, list):
        print_success(f"Selected {len(selected)} diverse variants")
        print(f"\nSelected variants:")
        for i, variant in enumerate(selected[:max_selected], 1):
            if isinstance(variant, str):
                print(f"  {i}. {variant[:60]}...")
            elif isinstance(variant, dict):
                seq = variant.get("sequence", variant.get("variant", variant.get("name", "")))
                name = variant.get("name", f"variant_{i}")
                if seq:
                    print(f"  {i}. {name}: {seq[:60]}...")
                else:
                    print(f"  {i}. {name}")
    else:
        print_warning("No selected_variants found in result")
        print_info(f"Available keys in result: {list(result_data.keys())}")
        
        # Check for error
        if "status" in result_data and result_data.get("status") == "error":
            print_error(f"Error: {result_data.get('message', 'Unknown error')}")
        elif "error" in result_data:
            print_error(f"Error: {result_data['error']}")
    
    # Display summary or analysis
    if "summary" in result_data:
        summary = result_data["summary"]
        if isinstance(summary, dict):
            print(f"\nSelection Summary:")
            for key, value in summary.items():
                if isinstance(value, (int, float, str, bool)):
                    print(f"  {key}: {value}")
    elif "analysis" in result_data:
        analysis = result_data["analysis"]
        if isinstance(analysis, dict):
            print(f"\nSelection Analysis:")
            for key, value in analysis.items():
                if isinstance(value, (int, float, str, bool)):
                    print(f"  {key}: {value}")
    
    # Show selection criteria and counts
    if "selection_criteria" in result_data:
        print(f"\nSelection Criteria: {result_data['selection_criteria']}")
    if "num_variants_selected" in result_data:
        print(f"Variants Selected: {result_data['num_variants_selected']}")
    if "num_variants_requested" in result_data:
        print(f"Variants Requested: {result_data['num_variants_requested']}")

def display_plasmid_result(result: Dict[str, Any]):
    """Display plasmid visualization result."""
    result_data = result.get("result", {})
    
    # Check for nested result structure
    if isinstance(result_data, dict) and "result" in result_data:
        result_data = result_data["result"]
    
    print(f"\n{Colors.BOLD}Plasmid Visualization Result:{Colors.RESET}")
    
    if "plasmid_data" in result_data:
        plasmid = result_data["plasmid_data"]
        if isinstance(plasmid, dict):
            print_success("Plasmid data generated")
            print(f"  Name: {plasmid.get('name', 'N/A')}")
            print(f"  Size: {plasmid.get('size', 'N/A')} bp")
            print(f"  Features: {len(plasmid.get('features', []))}")
        else:
            print_warning("Plasmid data is not a dict")
    elif "output" in result_data:
        output = result_data["output"]
        if isinstance(output, dict) and "plasmid_data" in output:
            plasmid = output["plasmid_data"]
            print_success("Plasmid data generated")
            print(f"  Name: {plasmid.get('name', 'N/A')}")
            print(f"  Size: {plasmid.get('size', 'N/A')} bp")
    else:
        print_warning("No plasmid_data found in result")
    
    if "visualization_type" in result_data:
        print(f"  Visualization type: {result_data['visualization_type']}")

def main():
    """Main test workflow."""
    print(f"\n{Colors.BOLD}{'='*60}{Colors.RESET}")
    print(f"{Colors.BOLD}Variant Analysis Workflow Test{Colors.RESET}")
    print(f"{Colors.BOLD}{'='*60}{Colors.RESET}")
    
    # Create session
    session_id = create_session()
    
    # Generate test sequence
    test_sequence = generate_test_sequence()
    print_success(f"Generated test sequence ({len(test_sequence)} bp)")
    
    # Step 1: Generate variants
    print_step(1, "Generate Variants")
    command1 = f"Generate 50 variants from this sequence:\n\n{test_sequence}"
    result1 = execute_command(command1, session_id)
    
    if result1.get("success"):
        print_success("Step 1 completed successfully")
        display_mutation_result(result1)
    else:
        print_error("Step 1 failed")
        print(json.dumps(result1, indent=2, default=str))
        return
    
    # Step 2: Select diverse variants
    print_step(2, "Select Diverse Variants")
    command2 = "Select 10 most diverse variants"
    result2 = execute_command(command2, session_id)
    
    if result2.get("success"):
        print_success("Step 2 completed successfully")
        display_selection_result(result2)
    else:
        print_warning("Step 2 may not be fully implemented")
        print_info("Result:")
        print(json.dumps(result2, indent=2, default=str))
    
    # Step 3: Insert into plasmid and visualize
    print_step(3, "Plasmid Visualization")
    # Get the first selected variant for visualization
    selected_variants = None
    if result2.get("success"):
        result_data = result2.get("result", {})
        if isinstance(result_data, dict) and "result" in result_data:
            result_data = result_data["result"]
        selected_variants = result_data.get("selected_variants", [])
    
    # If no selected variants, use a sample variant from step 1
    if not selected_variants or len(selected_variants) == 0:
        # Use a variant from step 1
        result_data = result1.get("result", {})
        if isinstance(result_data, dict) and "result" in result_data:
            result_data = result_data["result"]
        
        # Try multiple paths to find variants
        variants = result_data.get("variants", [])
        if not variants and "statistics" in result_data:
            stats = result_data["statistics"]
            if isinstance(stats, dict):
                variants = stats.get("variants", [])
        
        if variants and len(variants) > 0:
            # Take the first variant
            first_variant = variants[0]
            if isinstance(first_variant, str):
                selected_variants = [first_variant]
            elif isinstance(first_variant, dict):
                seq = first_variant.get("sequence", first_variant.get("variant", ""))
                if seq:
                    selected_variants = [seq]
    
    if selected_variants and len(selected_variants) > 0:
        # Get the first variant sequence
        variant_seq = selected_variants[0]
        if isinstance(variant_seq, dict):
            variant_seq = variant_seq.get("sequence", variant_seq.get("variant", ""))
        
        command3 = f"Insert this sequence into pUC19 plasmid and visualize:\n\n{variant_seq}"
        result3 = execute_command(command3, session_id)
        
        if result3.get("success"):
            print_success("Step 3 completed successfully")
            display_plasmid_result(result3)
        else:
            print_warning("Step 3 may not be fully implemented")
            print_info("Result:")
            print(json.dumps(result3, indent=2, default=str))
    else:
        print_warning("No variants available for plasmid visualization")
    
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
                
                if entry.get("tool") == "mutate_sequence":
                    # Try multiple paths to find variants
                    variants = actual_result.get("variants", [])
                    if not variants and "statistics" in actual_result:
                        stats = actual_result["statistics"]
                        if isinstance(stats, dict):
                            variants = stats.get("variants", [])
                    if not variants and "output" in actual_result:
                        output = actual_result["output"]
                        if isinstance(output, dict):
                            variants = output.get("variants", [])
                    print(f"    Input: DNA sequence")
                    print(f"    Output: {len(variants) if isinstance(variants, list) else 0} variants generated")
                
                elif entry.get("tool") == "select_variants":
                    selected = actual_result.get("selected_variants", [])
                    print(f"    Input: Generated variants")
                    print(f"    Output: {len(selected) if isinstance(selected, list) else 0} diverse variants selected")
                
                elif entry.get("tool") == "plasmid_visualization":
                    plasmid = actual_result.get("plasmid_data", {})
                    print(f"    Input: Selected variant sequence")
                    print(f"    Output: Plasmid visualization ({plasmid.get('size', 'N/A')} bp)")
                
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
        results_summary.append("✓ Step 1: Variant Generation")
    else:
        results_summary.append("✗ Step 1: Variant Generation")
    
    if result2.get("success"):
        results_summary.append("✓ Step 2: Variant Selection")
    else:
        results_summary.append("⚠ Step 2: Variant Selection (may not be implemented)")
    
    if result3.get("success") if 'result3' in locals() else False:
        results_summary.append("✓ Step 3: Plasmid Visualization")
    else:
        results_summary.append("⚠ Step 3: Plasmid Visualization (may not be implemented)")
    
    for summary in results_summary:
        print(f"  {summary}")
    
    print(f"\n{Colors.BOLD}Session Information:{Colors.RESET}")
    print(f"  Session ID: {session_id}")
    print(f"  Backend URL: {BACKEND_URL}")
    
    print(f"\n{Colors.BOLD}{'='*60}{Colors.RESET}")

if __name__ == "__main__":
    main()

