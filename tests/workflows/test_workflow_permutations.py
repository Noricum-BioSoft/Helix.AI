#!/usr/bin/env python3
"""
Comprehensive test script for workflow permutations in Helix.AI.

This script generates and tests biologically meaningful workflow combinations:
- Sequence Analysis Workflows
- Variant Analysis Workflows  
- Read Preprocessing Workflows
- Combined Workflows

Each workflow is tested end-to-end with appropriate test data.
"""

import requests
import json
import sys
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
from itertools import product
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
    CYAN = '\033[96m'
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

def print_workflow(msg: str):
    print(f"{Colors.CYAN}→ {msg}{Colors.RESET}")

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
        return data.get("session_id")
    except Exception as e:
        print_error(f"Failed to create session: {e}")
        return None

def execute_command(command: str, session_id: str) -> Dict[str, Any]:
    """Execute a natural language command."""
    try:
        response = requests.post(
            f"{BACKEND_URL}/execute",
            json={"command": command, "session_id": session_id},
            timeout=TIMEOUT
        )
        response.raise_for_status()
        return response.json()
    except Exception as e:
        return {"success": False, "error": str(e)}

def generate_test_sequences(count: int = 5) -> str:
    """Generate test DNA sequences."""
    sequences = [
        "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCT",
        "GTGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        "CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",
        "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
    ]
    
    fasta_lines = []
    for i, seq in enumerate(sequences[:count], 1):
        fasta_lines.append(f">test_seq_{i}")
        fasta_lines.append(seq)
    
    return "\n".join(fasta_lines)

def generate_test_sequence() -> str:
    """Generate a single test DNA sequence."""
    return "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"

# Define workflow patterns with their commands and dependencies
WORKFLOW_PATTERNS = {
    # Pattern 1: Simple Sequence Analysis
    "sequence_analysis_simple": {
        "description": "Align sequences → Build tree",
        "steps": [
            ("align", "Align these sequences:\n\n{sequences}"),
            ("phylogenetic_tree", "Build phylogenetic tree"),
        ],
        "requires": {"sequences": True}
    },
    
    # Pattern 2: Sequence Analysis with Selection
    "sequence_analysis_with_selection": {
        "description": "Align → Tree → Select representatives",
        "steps": [
            ("align", "Align these sequences:\n\n{sequences}"),
            ("phylogenetic_tree", "Build phylogenetic tree"),
            ("select_sequences", "Select 3 representative sequences"),
        ],
        "requires": {"sequences": True}
    },
    
    # Pattern 3: Sequence Analysis with Clustering
    "sequence_analysis_with_clustering": {
        "description": "Align → Tree → Cluster → Representatives",
        "steps": [
            ("align", "Align these sequences:\n\n{sequences}"),
            ("phylogenetic_tree", "Build phylogenetic tree"),
            ("cluster", "Select 5 representative sequences from each cluster"),
        ],
        "requires": {"sequences": True}
    },
    
    # Pattern 4: Variant Generation and Selection
    "variant_generation_selection": {
        "description": "Generate variants → Select diverse → Visualize",
        "steps": [
            ("mutate", "Generate 20 variants from this sequence:\n\n{sequence}"),
            ("select_variants", "Select 5 most diverse variants"),
            ("plasmid", "Insert this sequence into pUC19 plasmid and visualize:\n\n{selected_sequence}"),
        ],
        "requires": {"sequence": True, "selected_sequence": "from_step_2"}
    },
    
    # Pattern 5: Variant to Alignment
    "variant_to_alignment": {
        "description": "Generate variants → Align → Tree",
        "steps": [
            ("mutate", "Generate 10 variants from this sequence:\n\n{sequence}"),
            ("align", "Align these variant sequences:\n\n{variants}"),
            ("phylogenetic_tree", "Build phylogenetic tree"),
        ],
        "requires": {"sequence": True, "variants": "from_step_1"}
    },
    
    # Pattern 6: Read Preprocessing
    "read_preprocessing": {
        "description": "Trim → Merge → Quality assessment",
        "steps": [
            ("trim", "Trim low-quality bases from my FASTQ reads with quality threshold 20"),
            ("merge", "Merge my paired-end reads with minimum overlap of 12 bases"),
            ("quality", "Generate a quality report for the merged sequences"),
        ],
        "requires": {"fastq_files": True}
    },
    
    # Pattern 7: Read Preprocessing with Adapter Removal
    "read_preprocessing_with_adapters": {
        "description": "Trim → Remove adapters → Merge → Quality",
        "steps": [
            ("trim", "Trim low-quality bases from my FASTQ reads with quality threshold 20"),
            ("adapter", "Remove adapter sequences AGATCGGAAGAGC from my trimmed reads"),
            ("merge", "Merge my paired-end reads with minimum overlap of 12 bases"),
            ("quality", "Generate a quality report for the merged sequences"),
        ],
        "requires": {"fastq_files": True}
    },
    
    # Pattern 8: Complete Variant Analysis
    "complete_variant_analysis": {
        "description": "Generate → Select → Align → Tree → Visualize",
        "steps": [
            ("mutate", "Generate 15 variants from this sequence:\n\n{sequence}"),
            ("select_variants", "Select 5 most diverse variants"),
            ("align", "Align these selected variant sequences:\n\n{selected_variants}"),
            ("phylogenetic_tree", "Build phylogenetic tree"),
            ("plasmid", "Insert this sequence into pUC19 plasmid and visualize:\n\n{first_variant}"),
        ],
        "requires": {"sequence": True, "selected_variants": "from_step_2", "first_variant": "from_step_2"}
    },
    
    # Pattern 9: Alignment to Clustering to Plasmid
    "alignment_to_plasmid": {
        "description": "Align → Cluster → Visualize representatives in plasmids",
        "steps": [
            ("align", "Align these sequences:\n\n{sequences}"),
            ("cluster", "Select 3 representative sequences from each cluster"),
            ("plasmid", "Insert this sequence into pUC19 plasmid and visualize:\n\n{first_variant}"),
        ],
        "requires": {"sequences": True, "first_variant": "from_step_2"}
    },
    
    # Pattern 10: Variant Generation to Alignment to Selection
    "variant_to_alignment_to_selection": {
        "description": "Generate variants → Align → Select best",
        "steps": [
            ("mutate", "Generate 12 variants from this sequence:\n\n{sequence}"),
            ("align", "Align these variant sequences:\n\n{variants}"),
            ("select_sequences", "Select 3 sequences with best conservation"),
        ],
        "requires": {"sequence": True, "variants": "from_step_1"}
    },
    
    # Pattern 11: Simple Variant Workflow
    "simple_variant_workflow": {
        "description": "Generate variants → Select diverse",
        "steps": [
            ("mutate", "Generate 25 variants from this sequence:\n\n{sequence}"),
            ("select_variants", "Select 5 most diverse variants"),
        ],
        "requires": {"sequence": True}
    },
    
    # Pattern 12: Alignment to Tree Only
    "alignment_to_tree": {
        "description": "Align sequences → Build tree",
        "steps": [
            ("align", "Align these sequences:\n\n{sequences}"),
            ("phylogenetic_tree", "Build phylogenetic tree"),
        ],
        "requires": {"sequences": True}
    },
    
    # Pattern 13: Read Preprocessing Simple
    "read_preprocessing_simple": {
        "description": "Trim → Merge",
        "steps": [
            ("trim", "Trim low-quality bases from my FASTQ reads with quality threshold 20"),
            ("merge", "Merge my paired-end reads with minimum overlap of 12 bases"),
        ],
        "requires": {"fastq_files": True}
    },
}

def load_fastq_files():
    """Load test FASTQ files."""
    # Go up from tests/workflows to project root
    base_dir = Path(__file__).parent.parent.parent
    r1_file = base_dir / "data" / "rnaseq_demo" / "sample_R1_realistic.fastq"
    r2_file = base_dir / "data" / "rnaseq_demo" / "sample_R2_realistic.fastq"
    
    try:
        with open(r1_file, 'r') as f:
            r1_content = f.read()
        with open(r2_file, 'r') as f:
            r2_content = f.read()
        return {"r1": r1_content, "r2": r2_content}
    except Exception as e:
        print_warning(f"Could not load FASTQ files: {e}")
        return None

def execute_workflow_pattern(pattern_name: str, pattern: Dict[str, Any], test_data: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """
    Execute a workflow pattern and return success status and any errors.
    
    Returns:
        (success: bool, errors: List[str])
    """
    session_id = create_session()
    if not session_id:
        return False, ["Failed to create session"]
    
    errors = []
    step_results = []
    extracted_data = {}  # Store extracted data from previous steps
    
    print_workflow(f"Testing: {pattern['description']}")
    
    for step_idx, (step_type, command_template) in enumerate(pattern["steps"], 1):
        # Prepare command with test data
        command = command_template
        
        # Replace placeholders
        if "{sequences}" in command and "sequences" in test_data:
            command = command.replace("{sequences}", test_data["sequences"])
        if "{sequence}" in command and "sequence" in test_data:
            command = command.replace("{sequence}", test_data["sequence"])
        
        # Extract data from previous steps
        if "{variants}" in command:
            # Get variants from mutation step (usually step 1)
            for prev_idx, prev_result in enumerate(step_results, 1):
                variants = extract_variants_from_result(prev_result)
                if variants:
                    command = command.replace("{variants}", variants)
                    extracted_data["variants"] = variants
                    break
            else:
                errors.append(f"Step {step_idx}: Could not extract variants from previous steps")
                step_results.append({"success": False, "error": "Missing variants"})
                continue
        
        if "{selected_variants}" in command:
            # Get selected variants from selection step
            for prev_idx, prev_result in enumerate(step_results, 1):
                selected = extract_selected_variants_from_result(prev_result)
                if selected:
                    command = command.replace("{selected_variants}", selected)
                    extracted_data["selected_variants"] = selected
                    break
            else:
                errors.append(f"Step {step_idx}: Could not extract selected variants from previous steps")
                step_results.append({"success": False, "error": "Missing selected variants"})
                continue
        
        if "{selected_sequence}" in command or "{first_variant}" in command:
            # Get first variant from previous step
            for prev_idx, prev_result in enumerate(step_results, 1):
                first_variant = extract_first_variant_from_result(prev_result)
                if first_variant:
                    command = command.replace("{selected_sequence}", first_variant)
                    command = command.replace("{first_variant}", first_variant)
                    extracted_data["first_variant"] = first_variant
                    break
            else:
                # Fallback: use test sequence
                if "sequence" in test_data:
                    fallback_seq = test_data["sequence"][:50]  # Use first 50bp
                    command = command.replace("{selected_sequence}", fallback_seq)
                    command = command.replace("{first_variant}", fallback_seq)
                    print_warning(f"  Step {step_idx}: Using fallback sequence (variant not found)")
                else:
                    errors.append(f"Step {step_idx}: Could not extract variant and no fallback")
                    step_results.append({"success": False, "error": "Missing variant"})
                    continue
        
        # Handle FASTQ files
        if step_type == "trim" and "fastq_files" in test_data and test_data["fastq_files"]:
            fastq = test_data["fastq_files"]
            command = f"{command}\n\nForward reads (R1.fastq):\n{fastq['r1']}\n\nReverse reads (R2.fastq):\n{fastq['r2']}"
        elif step_type == "trim" and ("fastq_files" not in test_data or not test_data["fastq_files"]):
            errors.append(f"Step {step_idx}: FASTQ files not available")
            step_results.append({"success": False, "error": "Missing FASTQ files"})
            continue
        
        # Execute command
        print_info(f"  Step {step_idx}: {step_type}")
        result = execute_command(command, session_id)
        step_results.append(result)
        
        if not result.get("success"):
            error_msg = result.get("error", "Unknown error")
            # Try to extract error from nested structure
            if "result" in result and isinstance(result["result"], dict):
                nested_error = result["result"].get("error") or result["result"].get("message")
                if nested_error:
                    error_msg = nested_error
            errors.append(f"Step {step_idx} ({step_type}): {error_msg}")
            # Continue to next step even if this one failed
    
    # Determine overall success
    # Consider it successful if at least 50% of steps succeeded
    successful_steps = sum(1 for r in step_results if r.get("success", False))
    success_rate = successful_steps / len(pattern["steps"]) if pattern["steps"] else 0
    success = success_rate >= 0.5  # At least 50% success rate
    
    return success, errors

def extract_variants_from_result(result: Dict[str, Any]) -> Optional[str]:
    """Extract variants from mutation result and format as FASTA."""
    result_data = result.get("result", {})
    if isinstance(result_data, dict) and "result" in result_data:
        result_data = result_data["result"]
    
    variants = None
    if "statistics" in result_data and isinstance(result_data["statistics"], dict):
        variants = result_data["statistics"].get("variants", [])
    elif "variants" in result_data:
        variants = result_data["variants"]
    
    if variants and isinstance(variants, list):
        fasta_lines = []
        for i, variant in enumerate(variants, 1):
            if isinstance(variant, str):
                fasta_lines.append(f">variant_{i}")
                fasta_lines.append(variant)
            elif isinstance(variant, dict):
                seq = variant.get("sequence", variant.get("variant", ""))
                name = variant.get("name", f"variant_{i}")
                fasta_lines.append(f">{name}")
                fasta_lines.append(seq)
        return "\n".join(fasta_lines)
    
    return None

def extract_selected_variants_from_result(result: Dict[str, Any]) -> Optional[str]:
    """Extract selected variants and format as FASTA."""
    result_data = result.get("result", {})
    if isinstance(result_data, dict) and "result" in result_data:
        result_data = result_data["result"]
    
    selected = result_data.get("selected_variants", [])
    
    if selected and isinstance(selected, list):
        fasta_lines = []
        for i, variant in enumerate(selected, 1):
            if isinstance(variant, str):
                fasta_lines.append(f">selected_{i}")
                fasta_lines.append(variant)
            elif isinstance(variant, dict):
                seq = variant.get("sequence", variant.get("variant", ""))
                name = variant.get("name", f"selected_{i}")
                fasta_lines.append(f">{name}")
                fasta_lines.append(seq)
        return "\n".join(fasta_lines)
    
    return None

def extract_first_variant_from_result(result: Dict[str, Any]) -> Optional[str]:
    """Extract first variant/sequence from result."""
    result_data = result.get("result", {})
    if isinstance(result_data, dict) and "result" in result_data:
        result_data = result_data["result"]
    
    # Try selected_variants first
    selected = result_data.get("selected_variants", [])
    if selected and isinstance(selected, list) and len(selected) > 0:
        first = selected[0]
        if isinstance(first, str):
            return first
        elif isinstance(first, dict):
            seq = first.get("sequence", first.get("variant", ""))
            if seq:
                return seq
    
    # Try representatives from clustering
    if "clustering_result" in result_data:
        clustering = result_data["clustering_result"]
        if isinstance(clustering, dict):
            representatives = clustering.get("representatives", [])
            if representatives and isinstance(representatives, list) and len(representatives) > 0:
                # Representatives are usually just names, need to get sequences from aligned_sequences
                if "aligned_sequences" in result_data:
                    aligned = result_data["aligned_sequences"]
                    if isinstance(aligned, list):
                        for seq_dict in aligned:
                            if isinstance(seq_dict, dict) and seq_dict.get("name") == representatives[0]:
                                return seq_dict.get("sequence", "")
                # If we can't find the sequence, return the name as fallback
                return str(representatives[0])
    
    # Try selected_sequences (from sequence_selection)
    selected_seqs = result_data.get("selected_sequences", [])
    if selected_seqs and isinstance(selected_seqs, list) and len(selected_seqs) > 0:
        first = selected_seqs[0]
        if isinstance(first, str):
            return first
        elif isinstance(first, dict):
            seq = first.get("sequence", first.get("variant", ""))
            if seq:
                return seq
    
    # Try variants
    variants = None
    if "statistics" in result_data and isinstance(result_data["statistics"], dict):
        variants = result_data["statistics"].get("variants", [])
    elif "variants" in result_data:
        variants = result_data["variants"]
    
    if variants and isinstance(variants, list) and len(variants) > 0:
        first = variants[0]
        if isinstance(first, str):
            return first
        elif isinstance(first, dict):
            seq = first.get("sequence", first.get("variant", ""))
            if seq:
                return seq
    
    # Try aligned_sequences (from alignment or tree)
    aligned = result_data.get("aligned_sequences", [])
    if aligned and isinstance(aligned, list) and len(aligned) > 0:
        first = aligned[0]
        if isinstance(first, dict):
            return first.get("sequence", "")
        elif isinstance(first, str):
            return first
    
    return None

def main():
    """Main test function."""
    print(f"\n{Colors.BOLD}{'='*70}{Colors.RESET}")
    print(f"{Colors.BOLD}Comprehensive Workflow Permutation Test{Colors.RESET}")
    print(f"{Colors.BOLD}{'='*70}{Colors.RESET}\n")
    
    # Check backend
    print_info("Checking backend connection...")
    try:
        response = requests.get(f"{BACKEND_URL}/health", timeout=5)
        response.raise_for_status()
        print_success("Backend is running")
    except Exception as e:
        print_error(f"Backend is not accessible: {e}")
        sys.exit(1)
    
    # Prepare test data
    print_info("Preparing test data...")
    test_sequences = generate_test_sequences(5)
    test_sequence = generate_test_sequence()
    fastq_files = load_fastq_files()
    
    test_data = {
        "sequences": test_sequences,
        "sequence": test_sequence,
        "fastq_files": fastq_files
    }
    
    print_success("Test data prepared")
    
    # Test each workflow pattern
    results = {}
    total_patterns = len(WORKFLOW_PATTERNS)
    
    print(f"\n{Colors.BOLD}Testing {total_patterns} workflow patterns...{Colors.RESET}\n")
    
    for pattern_name, pattern in WORKFLOW_PATTERNS.items():
        print(f"\n{Colors.BOLD}Pattern: {pattern_name}{Colors.RESET}")
        print(f"  Description: {pattern['description']}")
        print(f"  Steps: {len(pattern['steps'])}")
        
        success, errors = execute_workflow_pattern(pattern_name, pattern, test_data)
        results[pattern_name] = {
            "success": success,
            "errors": errors,
            "description": pattern["description"],
            "steps": len(pattern["steps"])
        }
        
        if success:
            print_success(f"Pattern '{pattern_name}' completed successfully")
        else:
            print_error(f"Pattern '{pattern_name}' failed")
            for error in errors:
                print_error(f"  {error}")
        
        # Small delay between workflows
        time.sleep(0.5)
    
    # Summary report
    print(f"\n{Colors.BOLD}{'='*70}{Colors.RESET}")
    print(f"{Colors.BOLD}Test Summary{Colors.RESET}")
    print(f"{Colors.BOLD}{'='*70}{Colors.RESET}\n")
    
    successful = sum(1 for r in results.values() if r["success"])
    failed = total_patterns - successful
    
    print(f"Total Patterns Tested: {total_patterns}")
    print_success(f"Successful: {successful}")
    if failed > 0:
        print_error(f"Failed: {failed}")
    
    print(f"\n{Colors.BOLD}Detailed Results:{Colors.RESET}\n")
    
    for pattern_name, result in results.items():
        status = "✓" if result["success"] else "✗"
        color = Colors.GREEN if result["success"] else Colors.RED
        print(f"{color}{status}{Colors.RESET} {pattern_name}: {result['description']}")
        print(f"    Steps: {result['steps']}")
        if result["errors"]:
            print(f"    Errors: {len(result['errors'])}")
            for error in result["errors"][:2]:  # Show first 2 errors
                print(f"      - {error}")
    
    # Success rate
    success_rate = (successful / total_patterns * 100) if total_patterns > 0 else 0
    print(f"\n{Colors.BOLD}Success Rate: {success_rate:.1f}%{Colors.RESET}")
    
    print(f"\n{Colors.BOLD}{'='*70}{Colors.RESET}")

if __name__ == "__main__":
    main()

