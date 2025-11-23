#!/usr/bin/env python3
"""
Generate realistic FASTQ files with quality degradation for RNA-seq demo.

This script creates paired-end FASTQ files with:
- Quality scores that degrade toward read ends (realistic sequencing behavior)
- Some reads with adapter contamination
- Overlapping regions for successful read merging
"""

import random
import string

def generate_quality_scores(length, start_quality=40, end_quality=15):
    """Generate quality scores that degrade from start to end."""
    scores = []
    for i in range(length):
        # Linear degradation from start to end
        quality = int(start_quality - (start_quality - end_quality) * (i / length))
        # Add some randomness
        quality = max(0, min(40, quality + random.randint(-2, 2)))
        # Convert to Phred+33 ASCII
        scores.append(chr(quality + 33))
    return ''.join(scores)

def generate_sequence(length, bases='ATCG'):
    """Generate a random DNA sequence."""
    return ''.join(random.choice(bases) for _ in range(length))

def reverse_complement(seq):
    """Generate reverse complement of a sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(base, 'N') for base in reversed(seq))

# Generate 10 paired-end reads
num_reads = 10
read_length = 75
adapter = "AGATCGGAAGAGC"

with open('sample_R1_realistic.fastq', 'w') as r1, open('sample_R2_realistic.fastq', 'w') as r2:
    for i in range(1, num_reads + 1):
        # Generate forward read
        fwd_seq = generate_sequence(read_length)
        # Add adapter to some reads (30% chance)
        if random.random() < 0.3:
            fwd_seq = fwd_seq + adapter  # Full adapter
        
        fwd_qual = generate_quality_scores(len(fwd_seq))
        
        # Generate reverse read (overlapping with forward)
        rev_seq = reverse_complement(fwd_seq[-30:]) + generate_sequence(read_length - 30)
        # Add adapter to some reverse reads
        if random.random() < 0.3:
            rev_seq = adapter + rev_seq
        
        rev_qual = generate_quality_scores(len(rev_seq))
        
        # Write R1
        r1.write(f"@SRR001666.{i} {i}/1\n")
        r1.write(f"{fwd_seq}\n")
        r1.write("+\n")
        r1.write(f"{fwd_qual}\n")
        
        # Write R2
        r2.write(f"@SRR001666.{i} {i}/2\n")
        r2.write(f"{rev_seq}\n")
        r2.write("+\n")
        r2.write(f"{rev_qual}\n")

print("Generated realistic FASTQ files:")
print("- sample_R1_realistic.fastq")
print("- sample_R2_realistic.fastq")
print("\nFeatures:")
print("- Quality scores degrade from Phred 40 to Phred 15")
print("- Some reads contain adapter contamination")
print("- Reads have overlapping regions for merging")


