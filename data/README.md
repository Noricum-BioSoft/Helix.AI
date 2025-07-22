# Helix.AI Data Directory

This directory contains all data files used by the Helix.AI platform.

## Directory Structure

### `/samples/`
Sample sequence files for testing and development:
- **msa_file.fa**: Multiple sequence alignment file

### `/phylogenetic/`
Phylogenetic tree datasets:
- **beautiful_phylogenetic_dataset.fasta**: Beautiful phylogenetic dataset
- **diverse_phylogenetic_dataset.fasta**: Diverse phylogenetic dataset  
- **realistic_phylogenetic_dataset.fasta**: Realistic phylogenetic dataset
- **diverse_sequences.fasta**: Diverse sequence dataset
- **phylogenetic_tree_dataset.fasta**: Standard phylogenetic tree dataset

## Data Formats

### FASTA Format
All sequence files use the standard FASTA format:
```
>sequence_name
ATGCGATCGATCGATCG...
```

### Usage
These datasets are used for:
- Testing phylogenetic tree generation
- Validating sequence alignment algorithms
- Demonstrating clustering functionality
- Development and debugging

## Adding New Data

When adding new data files:
1. Place in appropriate subdirectory
2. Use descriptive filenames
3. Include format documentation
4. Update this README if needed 