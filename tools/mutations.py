# /backend/tools/mutation.py

import random
import numpy as np
from typing import List, Dict, Any

def run_mutation_raw(sequence: str, num_variants: int = 96):
    """Generate mutations of a given sequence and return variants with statistics."""
    
    if not sequence:
        return {
            "text": "Error: No sequence provided",
            "plot": {
                "data": [{"x": [], "y": [], "type": "bar"}],
                "layout": {"title": "Error - No Sequence"}
            }
        }
    
    # Store original sequence for display (before cleaning)
    original_sequence = sequence
    
    # Clean the sequence (remove spaces, convert to uppercase)
    sequence = sequence.replace(" ", "").upper()
    
    # Validate sequence contains only valid DNA/RNA characters
    valid_bases = set("ATCGU")
    if not all(base in valid_bases for base in sequence):
        return {
            "text": f"Error: Invalid sequence '{sequence}'. Only A, T, C, G, U allowed.",
            "plot": {
                "data": [{"x": [], "y": [], "type": "bar"}],
                "layout": {"title": "Error - Invalid Sequence"}
            }
        }
    
    # Generate mutations
    variants = []
    mutation_positions = []
    mutation_types = []
    
    for i in range(num_variants):
        # Create a copy of the original sequence
        variant = list(sequence)
        
        # Randomly select mutation type and position
        mutation_type = random.choice(['substitution', 'insertion', 'deletion'])
        position = random.randint(0, len(sequence) - 1)
        
        if mutation_type == 'substitution':
            # Replace with a different base
            original_base = variant[position]
            new_bases = [b for b in "ATCG" if b != original_base]
            variant[position] = random.choice(new_bases)
            mutation_positions.append(position)
            mutation_types.append('substitution')
            
        elif mutation_type == 'insertion':
            # Insert a random base
            new_base = random.choice("ATCG")
            variant.insert(position, new_base)
            mutation_positions.append(position)
            mutation_types.append('insertion')
            
        elif mutation_type == 'deletion':
            # Delete a base
            if len(variant) > 1:  # Only delete if sequence has more than 1 base
                deleted_base = variant.pop(position)
                mutation_positions.append(position)
                mutation_types.append('deletion')
        
        variants.append(''.join(variant))
    
    # Calculate statistics
    total_mutations = len(mutation_positions)
    substitution_count = mutation_types.count('substitution')
    insertion_count = mutation_types.count('insertion')
    deletion_count = mutation_types.count('deletion')
    
    # Calculate mutation rates
    mutation_rate = total_mutations / (len(sequence) * num_variants)
    
    # Create visualization data
    mutation_counts = {
        'Substitution': substitution_count,
        'Insertion': insertion_count,
        'Deletion': deletion_count
    }
    
    # Position distribution for mutations
    position_counts = {}
    for pos in mutation_positions:
        position_counts[pos] = position_counts.get(pos, 0) + 1
    
    # Generate plot data
    plot_data = [
        {
            "x": list(mutation_counts.keys()),
            "y": list(mutation_counts.values()),
            "type": "bar",
            "name": "Mutation Types"
        },
        {
            "x": list(position_counts.keys()),
            "y": list(position_counts.values()),
            "type": "scatter",
            "mode": "markers",
            "name": "Mutation Positions"
        }
    ]
    
    # Create result text - use original sequence for display, but truncate if too long
    display_sequence = original_sequence if len(original_sequence) <= 50 else original_sequence[:47] + "..."
    result_text = f"""Generated {num_variants} mutated variants from sequence '{display_sequence}'.
    
Statistics:
- Total mutations: {total_mutations}
- Substitutions: {substitution_count}
- Insertions: {insertion_count}
- Deletions: {deletion_count}
- Average mutation rate: {mutation_rate:.4f}

First 5 variants:
{chr(10).join(variants[:5])}

... and {len(variants) - 5} more variants generated."""
    
    return {
        "text": result_text,
        "plot": {
            "data": plot_data,
            "layout": {
                "title": f"Mutation Analysis - {num_variants} Variants",
                "xaxis": {"title": "Mutation Type / Position"},
                "yaxis": {"title": "Count"},
                "barmode": "group"
            }
        },
        "statistics": {
            "total_variants": num_variants,
            "total_mutations": total_mutations,
            "substitution_count": substitution_count,
            "insertion_count": insertion_count,
            "deletion_count": deletion_count,
            "mutation_rate": mutation_rate,
            "original_sequence": sequence,
            "variants": variants  # Return all variants
        }
    }

from langchain_core.tools import tool

@tool
def run_mutation(sequence: str, num_variants: int = 96):
    """Mutates a given sequence and returns the specified number of variants."""
    return run_mutation_raw(sequence, num_variants)
