import json
import re
from typing import List, Dict, Any
from datetime import datetime, timedelta

def parse_sequences_for_synthesis(sequences_data: str) -> List[Dict[str, Any]]:
    """Parse sequences for synthesis submission."""
    sequences = []
    current_name = ""
    current_sequence = ""
    
    lines = sequences_data.strip().split('\n')
    
    # Check if input is just a sequence without FASTA headers
    if len(lines) == 1 and not lines[0].startswith('>'):
        # Single sequence without header
        sequence = lines[0].strip()
        if sequence:
            sequences.append({
                "name": "sequence1",
                "sequence": sequence,
                "length": len(sequence)
            })
        return sequences
    
    # Parse FASTA format
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            # Save previous sequence if exists
            if current_name and current_sequence:
                sequences.append({
                    "name": current_name,
                    "sequence": current_sequence,
                    "length": len(current_sequence)
                })
            # Start new sequence
            current_name = line[1:].strip()
            current_sequence = ""
        else:
            # Add to current sequence
            current_sequence += line
    
    # Add the last sequence
    if current_name and current_sequence:
        sequences.append({
            "name": current_name,
            "sequence": current_sequence,
            "length": len(current_sequence)
        })
    
    return sequences

def validate_sequences_for_synthesis(sequences: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Validate sequences for synthesis requirements."""
    validation_results = {
        "valid": True,
        "errors": [],
        "warnings": [],
        "statistics": {
            "total_sequences": len(sequences),
            "total_length": 0,
            "average_length": 0,
            "max_length": 0,
            "min_length": 0
        }
    }
    
    if not sequences:
        validation_results["valid"] = False
        validation_results["errors"].append("No sequences provided")
        return validation_results
    
    lengths = []
    total_length = 0
    
    for i, seq_data in enumerate(sequences):
        sequence = seq_data["sequence"].upper()
        name = seq_data["name"]
        length = len(sequence)
        
        # Check for invalid characters
        invalid_chars = [char for char in sequence if char not in "ATCG"]
        if invalid_chars:
            validation_results["valid"] = False
            validation_results["errors"].append(f"Sequence '{name}' contains invalid characters: {set(invalid_chars)}")
        
        # Check length constraints
        if length < 20:
            validation_results["warnings"].append(f"Sequence '{name}' is very short ({length} bp) - may be difficult to synthesize")
        elif length > 10000:
            validation_results["valid"] = False
            validation_results["errors"].append(f"Sequence '{name}' is too long ({length} bp) - exceeds typical synthesis limits")
        
        # Check GC content
        gc_content = (sequence.count('G') + sequence.count('C')) / length * 100 if length > 0 else 0
        if gc_content < 20 or gc_content > 80:
            validation_results["warnings"].append(f"Sequence '{name}' has extreme GC content ({gc_content:.1f}%) - may affect synthesis efficiency")
        
        # Check for repetitive sequences
        if length > 50:
            # Check for simple repeats
            for repeat_len in [3, 4, 5]:
                for start in range(length - repeat_len + 1):
                    repeat = sequence[start:start + repeat_len]
                    if sequence.count(repeat) > 3:
                        validation_results["warnings"].append(f"Sequence '{name}' contains repetitive motif '{repeat}' - may cause synthesis issues")
                        break
        
        lengths.append(length)
        total_length += length
    
    # Calculate statistics
    if lengths:
        validation_results["statistics"]["total_length"] = total_length
        validation_results["statistics"]["average_length"] = total_length / len(lengths)
        validation_results["statistics"]["max_length"] = max(lengths)
        validation_results["statistics"]["min_length"] = min(lengths)
    
    return validation_results

def generate_synthesis_quote(sequences: List[Dict[str, Any]], 
                           vendor_preference: str = None,
                           quantity: str = "standard",
                           delivery_time: str = "standard") -> Dict[str, Any]:
    """Generate synthesis quote for sequences."""
    
    # Validate sequences first
    validation = validate_sequences_for_synthesis(sequences)
    
    if not validation["valid"]:
        return {
            "text": f"Error: Sequences failed validation:\n" + "\n".join(validation["errors"]),
            "error": "Validation failed",
            "validation_results": validation
        }
    
    # Calculate pricing based on sequence characteristics
    total_cost = 0
    cost_breakdown = []
    
    for seq_data in sequences:
        length = seq_data["length"]
        name = seq_data["name"]
        
        # Base pricing model (simplified)
        if length <= 100:
            base_price = 0.08  # $0.08 per bp for short sequences
        elif length <= 1000:
            base_price = 0.06  # $0.06 per bp for medium sequences
        else:
            base_price = 0.04  # $0.04 per bp for long sequences
        
        # Quantity discount
        if quantity == "large":
            base_price *= 0.8
        elif quantity == "custom":
            base_price *= 1.2
        
        # Rush delivery premium
        if delivery_time == "rush":
            base_price *= 1.5
        
        sequence_cost = length * base_price
        total_cost += sequence_cost
        
        cost_breakdown.append({
            "name": name,
            "length": length,
            "base_price_per_bp": base_price,
            "total_cost": sequence_cost
        })
    
    # Add setup fees
    setup_fee = 50 if len(sequences) <= 5 else 100
    total_cost += setup_fee
    
    # Calculate delivery time
    base_delivery_days = 7
    if delivery_time == "rush":
        base_delivery_days = 3
    elif delivery_time == "standard":
        base_delivery_days = 7
    elif delivery_time == "economy":
        base_delivery_days = 14
    
    # Adjust for sequence complexity
    max_length = validation["statistics"]["max_length"]
    if max_length > 5000:
        base_delivery_days += 3
    
    delivery_date = datetime.now() + timedelta(days=base_delivery_days)
    
    # Generate quote text
    quote_text = f"""DNA Synthesis Quote

Sequences: {len(sequences)}
Total length: {validation['statistics']['total_length']:,} bp
Average length: {validation['statistics']['average_length']:.1f} bp

Cost Breakdown:"""
    
    for item in cost_breakdown:
        quote_text += f"""
- {item['name']}: {item['length']} bp × ${item['base_price_per_bp']:.3f}/bp = ${item['total_cost']:.2f}"""
    
    quote_text += f"""

Setup fee: ${setup_fee:.2f}
Total cost: ${total_cost:.2f}

Estimated delivery: {delivery_date.strftime('%B %d, %Y')}
Delivery time: {base_delivery_days} business days

Vendor recommendations:"""
    
    # Vendor recommendations based on sequence characteristics
    if validation["statistics"]["max_length"] > 5000:
        quote_text += """
- Twist Bioscience (specializes in long sequences)
- GenScript (high-quality long gene synthesis)"""
    else:
        quote_text += """
- IDT (fast turnaround for short sequences)
- Eurofins Genomics (good value for standard sequences)"""
    
    if validation["warnings"]:
        quote_text += "\n\nWarnings:"
        for warning in validation["warnings"]:
            quote_text += f"\n- {warning}"
    
    return {
        "text": quote_text,
        "quote": {
            "total_cost": total_cost,
            "setup_fee": setup_fee,
            "delivery_days": base_delivery_days,
            "delivery_date": delivery_date.isoformat(),
            "cost_breakdown": cost_breakdown,
            "validation_results": validation
        }
    }

def run_synthesis_submission_raw(sequences: str, 
                                vendor_preference: str = None,
                                quantity: str = "standard",
                                delivery_time: str = "standard"):
    """Submit sequences for DNA synthesis and get quote."""
    
    if not sequences or not sequences.strip():
        return {
            "text": "Error: No sequences provided for synthesis",
            "error": "No sequences"
        }
    
    # Parse sequences
    try:
        parsed_sequences = parse_sequences_for_synthesis(sequences)
    except Exception as e:
        return {
            "text": f"Error parsing sequences: {str(e)}",
            "error": str(e)
        }
    
    if not parsed_sequences:
        return {
            "text": "Error: No valid sequences found",
            "error": "No valid sequences"
        }
    
    # Generate synthesis quote
    return generate_synthesis_quote(parsed_sequences, vendor_preference, quantity, delivery_time)

from langchain.agents import tool

@tool
def run_synthesis_submission(sequences: str, vendor_preference: str = None, quantity: str = "standard", delivery_time: str = "standard"):
    """Submit sequences for DNA synthesis and get pricing quote. Quantity options: standard, large, custom. Delivery options: rush, standard, economy."""
    return run_synthesis_submission_raw(sequences, vendor_preference, quantity, delivery_time) 