#!/usr/bin/env python3
"""
Example usage of the Bioinformatics MCP Server
"""

import asyncio
import json
import sys
from pathlib import Path

# Add the current directory to Python path
sys.path.append(str(Path(__file__).parent))

async def example_sequence_alignment():
    """Example of sequence alignment using MCP."""
    print("Example: Sequence Alignment")
    print("-" * 30)
    
    sequences = """>seq1
ACTGTTGAC
>seq2
ACTGCATCC
>seq3
ACTGCAATGAC"""
    
    try:
        from tools.alignment import run_alignment
        result = run_alignment(sequences)
        print(f"Input sequences: {sequences}")
        print(f"Alignment result: {json.dumps(result, indent=2)}")
    except Exception as e:
        print(f"Error: {e}")

async def example_mutation_analysis():
    """Example of sequence mutation analysis using MCP."""
    print("\nExample: Sequence Mutation Analysis")
    print("-" * 40)
    
    sequence = "ACTGTTGAC"
    num_variants = 5
    
    try:
        from tools.mutations import mutate_sequence
        result = mutate_sequence(sequence, num_variants)
        print(f"Original sequence: {sequence}")
        print(f"Number of variants: {num_variants}")
        print(f"Mutation result: {json.dumps(result, indent=2)}")
    except Exception as e:
        print(f"Error: {e}")

async def example_bioinformatics_analysis():
    """Example of bioinformatics analysis using MCP."""
    print("\nExample: Bioinformatics Analysis")
    print("-" * 35)
    
    try:
        import pandas as pd
        from tools.bio import align_and_visualize_fasta
        
        # Create test data
        test_data = pd.DataFrame({
            'name': ['seq1', 'seq2', 'seq3'],
            'sequence': ['ACTGTTGAC', 'ACTGCATCC', 'ACTGCAATGAC']
        })
        
        print(f"Input data: {test_data}")
        result = align_and_visualize_fasta(test_data)
        print(f"Analysis result: {result}")
    except Exception as e:
        print(f"Error: {e}")

async def main():
    """Run all examples."""
    print("Bioinformatics MCP Server Examples")
    print("=" * 40)
    
    await example_sequence_alignment()
    await example_mutation_analysis()
    await example_bioinformatics_analysis()
    
    print("\n" + "=" * 40)
    print("Examples completed!")

if __name__ == "__main__":
    asyncio.run(main()) 