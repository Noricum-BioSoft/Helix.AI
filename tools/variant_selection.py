import random
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from io import StringIO
import logging

logger = logging.getLogger(__name__)

def select_variants_from_history(session_id: str, 
                               selection_criteria: str = "diversity",
                               num_variants: int = 10,
                               custom_filters: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Select variants from previous mutation results based on specified criteria.
    
    Args:
        session_id: Session ID to get previous results from
        selection_criteria: Criteria for selection ("diversity", "random", "length", "custom")
        num_variants: Number of variants to select
        custom_filters: Custom filtering criteria
    
    Returns:
        Dictionary with selected variants and analysis
    """
    try:
        # Import history manager
        import sys
        from pathlib import Path
        backend_path = Path(__file__).parent.parent / "backend"
        sys.path.insert(0, str(backend_path))
        
        from history_manager import history_manager
        
        # Get previous mutation results
        mutation_results = history_manager.get_latest_result(session_id, "mutate_sequence")
        
        if not mutation_results:
            return {
                "status": "error",
                "message": "No previous mutation results found in session",
                "session_id": session_id
            }
        
        # DEBUG: Print mutation result structure
        logger.info(f"🔧 [DEBUG] Mutation result keys: {list(mutation_results.keys()) if isinstance(mutation_results, dict) else 'Not a dict'}")
        if isinstance(mutation_results, dict) and "statistics" in mutation_results:
            logger.info(f"🔧 [DEBUG] Statistics keys: {list(mutation_results['statistics'].keys()) if isinstance(mutation_results['statistics'], dict) else 'Not a dict'}")
            if isinstance(mutation_results['statistics'], dict) and "variants" in mutation_results['statistics']:
                variants_list = mutation_results['statistics']['variants']
                logger.info(f"🔧 [DEBUG] Found {len(variants_list) if isinstance(variants_list, list) else 'N/A'} variants in statistics.variants")
        
        # Extract variants from the result
        variants = extract_variants_from_result(mutation_results)
        
        logger.info(f"🔧 [DEBUG] Extracted {len(variants)} variants from mutation result")
        
        if not variants:
            return {
                "status": "error", 
                "message": "No variants found in previous mutation results",
                "session_id": session_id
            }
        
        if len(variants) < 2:
            return {
                "status": "error",
                "message": f"At least 2 variants are required for variant selection, but only {len(variants)} were found",
                "session_id": session_id
            }
        
        # Select variants based on criteria
        selected_variants = select_variants_by_criteria(
            variants, selection_criteria, num_variants, custom_filters
        )
        
        # Analyze selection
        analysis = analyze_variant_selection(variants, selected_variants, selection_criteria)
        
        # Create visualization
        plot_data = create_selection_visualization(variants, selected_variants, selection_criteria)
        
        return {
            "status": "success",
            "session_id": session_id,
            "selection_criteria": selection_criteria,
            "num_variants_requested": num_variants,
            "num_variants_selected": len(selected_variants),
            "selected_variants": selected_variants,
            "analysis": analysis,
            "plot": plot_data
        }
        
    except Exception as e:
        logger.error(f"Error in variant selection: {e}")
        return {
            "status": "error",
            "message": f"Error selecting variants: {str(e)}",
            "session_id": session_id
        }

def extract_variants_from_result(mutation_result: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Extract variants from mutation result."""
    variants = []
    
    # Handle different result formats
    if "statistics" in mutation_result and "variants" in mutation_result["statistics"]:
        # New format: variants stored in statistics.variants
        variant_sequences = mutation_result["statistics"]["variants"]
        variants = []
        for i, sequence in enumerate(variant_sequences):
            variants.append({
                "name": f"variant_{i+1}",
                "sequence": sequence
            })
    elif "output" in mutation_result:
        output = mutation_result["output"]
        if isinstance(output, list):
            variants = output
        elif isinstance(output, str):
            # Parse FASTA format
            variants = parse_fasta_variants(output)
    elif "variants" in mutation_result:
        variants = mutation_result["variants"]
    elif "sequences" in mutation_result:
        variants = mutation_result["sequences"]
    
    # Ensure variants have proper structure
    normalized_variants = []
    for i, variant in enumerate(variants):
        if isinstance(variant, dict):
            normalized_variants.append(variant)
        elif isinstance(variant, str):
            normalized_variants.append({
                "name": f"variant_{i+1}",
                "sequence": variant
            })
        else:
            normalized_variants.append({
                "name": f"variant_{i+1}",
                "sequence": str(variant)
            })
    
    return normalized_variants

def parse_fasta_variants(fasta_content: str) -> List[Dict[str, Any]]:
    """Parse variants from FASTA content."""
    variants = []
    lines = fasta_content.strip().split('\n')
    current_name = ""
    current_sequence = ""
    
    for line in lines:
        if line.startswith('>'):
            if current_name and current_sequence:
                variants.append({
                    "name": current_name,
                    "sequence": current_sequence
                })
            current_name = line[1:].strip()
            current_sequence = ""
        else:
            current_sequence += line.strip()
    
    # Add the last variant
    if current_name and current_sequence:
        variants.append({
            "name": current_name,
            "sequence": current_sequence
        })
    
    return variants

def select_variants_by_criteria(variants: List[Dict[str, Any]], 
                               criteria: str, 
                               num_variants: int,
                               custom_filters: Optional[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
    """Select variants based on specified criteria."""
    
    if criteria == "random":
        return random.sample(variants, min(num_variants, len(variants)))
    
    elif criteria == "diversity":
        return select_diverse_variants(variants, num_variants)
    
    elif criteria == "length":
        return select_by_length(variants, num_variants, custom_filters)
    
    elif criteria == "custom":
        return select_by_custom_filters(variants, num_variants, custom_filters)
    
    else:
        # Default to random selection
        return random.sample(variants, min(num_variants, len(variants)))

def select_diverse_variants(variants: List[Dict[str, Any]], num_variants: int) -> List[Dict[str, Any]]:
    """Select diverse variants using sequence similarity."""
    if len(variants) <= num_variants:
        return variants
    
    # Calculate pairwise distances
    sequences = [v["sequence"] for v in variants]
    distances = calculate_sequence_distances(sequences)
    
    # Use hierarchical clustering to select diverse variants
    selected_indices = select_diverse_indices(distances, num_variants)
    
    return [variants[i] for i in selected_indices]

def calculate_sequence_distances(sequences: List[str]) -> np.ndarray:
    """Calculate pairwise distances between sequences."""
    n = len(sequences)
    distances = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            # Simple Hamming distance for DNA sequences
            seq1, seq2 = sequences[i], sequences[j]
            min_len = min(len(seq1), len(seq2))
            distance = sum(1 for k in range(min_len) if seq1[k] != seq2[k])
            distance += abs(len(seq1) - len(seq2))  # Penalty for length differences
            
            distances[i, j] = distance
            distances[j, i] = distance
    
    return distances

def select_diverse_indices(distances: np.ndarray, num_variants: int) -> List[int]:
    """Select diverse indices using farthest point sampling."""
    n = distances.shape[0]
    if n <= num_variants:
        return list(range(n))
    
    # Start with the first sequence
    selected = [0]
    
    # Iteratively select the sequence farthest from all selected sequences
    for _ in range(num_variants - 1):
        # Calculate minimum distance to selected sequences for each unselected sequence
        unselected = [i for i in range(n) if i not in selected]
        min_distances = []
        
        for i in unselected:
            min_dist = min(distances[i, j] for j in selected)
            min_distances.append((min_dist, i))
        
        # Select the sequence with maximum minimum distance
        selected.append(max(min_distances)[1])
    
    return selected

def select_by_length(variants: List[Dict[str, Any]], num_variants: int, 
                    filters: Optional[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
    """Select variants based on sequence length."""
    # Sort by length
    sorted_variants = sorted(variants, key=lambda x: len(x["sequence"]))
    
    if filters and "length_range" in filters:
        min_len, max_len = filters["length_range"]
        filtered_variants = [
            v for v in sorted_variants 
            if min_len <= len(v["sequence"]) <= max_len
        ]
        return filtered_variants[:num_variants]
    
    # Default: select from middle range
    start_idx = len(sorted_variants) // 4
    end_idx = start_idx + num_variants
    return sorted_variants[start_idx:end_idx]

def select_by_custom_filters(variants: List[Dict[str, Any]], num_variants: int,
                           filters: Optional[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
    """Select variants using custom filters."""
    if not filters:
        return random.sample(variants, min(num_variants, len(variants)))
    
    filtered_variants = variants.copy()
    
    # Apply filters
    if "min_length" in filters:
        filtered_variants = [
            v for v in filtered_variants 
            if len(v["sequence"]) >= filters["min_length"]
        ]
    
    if "max_length" in filters:
        filtered_variants = [
            v for v in filtered_variants 
            if len(v["sequence"]) <= filters["max_length"]
        ]
    
    if "gc_content_range" in filters:
        min_gc, max_gc = filters["gc_content_range"]
        filtered_variants = [
            v for v in filtered_variants 
            if min_gc <= calculate_gc_content(v["sequence"]) <= max_gc
        ]
    
    if "substring" in filters:
        substring = filters["substring"]
        filtered_variants = [
            v for v in filtered_variants 
            if substring.upper() in v["sequence"].upper()
        ]
    
    return filtered_variants[:num_variants]

def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content of a sequence."""
    if not sequence:
        return 0.0
    
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return gc_count / len(sequence)

def analyze_variant_selection(all_variants: List[Dict[str, Any]], 
                            selected_variants: List[Dict[str, Any]],
                            criteria: str) -> Dict[str, Any]:
    """Analyze the variant selection."""
    analysis = {
        "total_variants": len(all_variants),
        "selected_variants": len(selected_variants),
        "selection_ratio": len(selected_variants) / len(all_variants) if all_variants else 0,
        "criteria_used": criteria
    }
    
    # Length statistics
    all_lengths = [len(v["sequence"]) for v in all_variants]
    selected_lengths = [len(v["sequence"]) for v in selected_variants]
    
    analysis["length_stats"] = {
        "all_variants": {
            "mean": np.mean(all_lengths),
            "std": np.std(all_lengths),
            "min": min(all_lengths),
            "max": max(all_lengths)
        },
        "selected_variants": {
            "mean": np.mean(selected_lengths),
            "std": np.std(selected_lengths),
            "min": min(selected_lengths),
            "max": max(selected_lengths)
        }
    }
    
    # GC content statistics
    all_gc = [calculate_gc_content(v["sequence"]) for v in all_variants]
    selected_gc = [calculate_gc_content(v["sequence"]) for v in selected_variants]
    
    analysis["gc_content_stats"] = {
        "all_variants": {
            "mean": np.mean(all_gc),
            "std": np.std(all_gc),
            "min": min(all_gc),
            "max": max(all_gc)
        },
        "selected_variants": {
            "mean": np.mean(selected_gc),
            "std": np.std(selected_gc),
            "min": min(selected_gc),
            "max": max(selected_gc)
        }
    }
    
    return analysis

def create_selection_visualization(all_variants: List[Dict[str, Any]], 
                                 selected_variants: List[Dict[str, Any]],
                                 criteria: str) -> Dict[str, Any]:
    """Create visualization for variant selection."""
    try:
        # Create comparison plots
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f"Variant Selection Analysis - {criteria.capitalize()}")
        
        # Length distribution
        all_lengths = [len(v["sequence"]) for v in all_variants]
        selected_lengths = [len(v["sequence"]) for v in selected_variants]
        
        axes[0, 0].hist(all_lengths, alpha=0.7, label="All Variants", bins=20)
        axes[0, 0].hist(selected_lengths, alpha=0.7, label="Selected Variants", bins=20)
        axes[0, 0].set_xlabel("Sequence Length")
        axes[0, 0].set_ylabel("Frequency")
        axes[0, 0].set_title("Length Distribution")
        axes[0, 0].legend()
        
        # GC content distribution
        all_gc = [calculate_gc_content(v["sequence"]) for v in all_variants]
        selected_gc = [calculate_gc_content(v["sequence"]) for v in selected_variants]
        
        axes[0, 1].hist(all_gc, alpha=0.7, label="All Variants", bins=20)
        axes[0, 1].hist(selected_gc, alpha=0.7, label="Selected Variants", bins=20)
        axes[0, 1].set_xlabel("GC Content")
        axes[0, 1].set_ylabel("Frequency")
        axes[0, 1].set_title("GC Content Distribution")
        axes[0, 1].legend()
        
        # Length vs GC content scatter
        axes[1, 0].scatter(all_lengths, all_gc, alpha=0.6, label="All Variants")
        axes[1, 0].scatter(selected_lengths, selected_gc, alpha=0.8, 
                           color='red', label="Selected Variants")
        axes[1, 0].set_xlabel("Sequence Length")
        axes[1, 0].set_ylabel("GC Content")
        axes[1, 0].set_title("Length vs GC Content")
        axes[1, 0].legend()
        
        # Selection summary
        summary_text = f"""
        Total Variants: {len(all_variants)}
        Selected Variants: {len(selected_variants)}
        Selection Ratio: {len(selected_variants)/len(all_variants)*100:.1f}%
        Criteria: {criteria}
        """
        axes[1, 1].text(0.1, 0.5, summary_text, transform=axes[1, 1].transAxes,
                        fontsize=12, verticalalignment='center')
        axes[1, 1].set_title("Selection Summary")
        axes[1, 1].axis('off')
        
        plt.tight_layout()
        
        # Convert to plotly format
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        
        fig_plotly = make_subplots(
            rows=2, cols=2,
            subplot_titles=("Length Distribution", "GC Content Distribution", 
                           "Length vs GC Content", "Selection Summary")
        )
        
        # Add traces (simplified for plotly)
        fig_plotly.add_trace(
            go.Histogram(x=all_lengths, name="All Variants", opacity=0.7),
            row=1, col=1
        )
        fig_plotly.add_trace(
            go.Histogram(x=selected_lengths, name="Selected Variants", opacity=0.7),
            row=1, col=1
        )
        
        fig_plotly.add_trace(
            go.Histogram(x=all_gc, name="All Variants", opacity=0.7),
            row=1, col=2
        )
        fig_plotly.add_trace(
            go.Histogram(x=selected_gc, name="Selected Variants", opacity=0.7),
            row=1, col=2
        )
        
        fig_plotly.add_trace(
            go.Scatter(x=all_lengths, y=all_gc, mode='markers', 
                      name="All Variants", opacity=0.6),
            row=2, col=1
        )
        fig_plotly.add_trace(
            go.Scatter(x=selected_lengths, y=selected_gc, mode='markers',
                      name="Selected Variants", marker=dict(color='red'), opacity=0.8),
            row=2, col=1
        )
        
        fig_plotly.update_layout(
            title=f"Variant Selection Analysis - {criteria.capitalize()}",
            height=600
        )
        
        return {
            "data": fig_plotly.to_dict()["data"],
            "layout": fig_plotly.to_dict()["layout"]
        }
        
    except Exception as e:
        logger.error(f"Error creating visualization: {e}")
        return {
            "data": [{"x": [1], "y": [1], "type": "bar"}],
            "layout": {"title": "Visualization Error"}
        }

def run_variant_selection_raw(session_id: str, 
                            selection_criteria: str = "diversity",
                            num_variants: int = 10,
                            custom_filters: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Raw function for variant selection (for direct calls)."""
    return select_variants_from_history(session_id, selection_criteria, num_variants, custom_filters) 