# 📝 Sequence Examples Guide

This document explains how sequences are handled in example commands throughout the Helix.AI frontend.

## Overview

All bioinformatics operations require DNA or protein sequences. To make examples functional and immediately usable, we've integrated sample sequences directly into example commands.

## How It Works

### 1. Sample Sequences Library

Located in `src/utils/sampleSequences.ts`, this file contains:
- **Short sequences**: For quick examples
- **Medium sequences**: For alignment examples  
- **FASTA format sequences**: Pre-formatted for alignment and phylogenetic tree operations
- **Single sequences**: For mutation examples
- **Plasmid inserts**: For plasmid visualization

### 2. Automatic Sequence Injection

The `getExampleWithSequences()` helper function automatically adds appropriate sequences to commands:

```typescript
// Example: Alignment command
getExampleWithSequences("align these sequences", "threeSequences")
// Returns: "align these sequences\n\n>seq1\nATGCGATCG...\n>seq2\n..."

// Example: Phylogenetic tree
getExampleWithSequences("visualize the phylogenetic tree", "phylogenetic")
// Returns: Command with 5 longer sequences for tree building
```

### 3. Example Categories

#### Commands WITH Sequences Included

These examples include complete, executable sequences:

- **Alignment Examples**: Include 2-3 sequences in FASTA format
- **Phylogenetic Tree**: Include 5 longer sequences (150+ bp each)
- **Mutation Examples**: Include a single sequence for mutation
- **Plasmid Visualization**: Include insert sequence

#### Commands WITHOUT Sequences (Context-Dependent)

These examples rely on sequences from previous steps or uploaded files:

- **Sequence Selection**: "select 10 representative sequences"
  - Requires: Previous alignment step
  - Uses: Sequences from workflow context

- **Plasmid Operations**: "insert each of the sequences into pUC19"
  - Requires: Sequences from previous step or uploaded file
  - Uses: Workflow context or file upload

- **Vendor Research**: "research DNA synthesis vendors for 1000bp sequences"
  - Requires: No sequences needed
  - Uses: Length parameter only

## Example Command Structure

### Complete Example (With Sequences)

```
visualize the phylogenetic tree

>Sequence_01
ACTCGATCACAAAGCTTAGGTCCGATCAATTTTGATAGTTACCCCCCACGGTCCAATCCGTTGGGTGAACACCGAGAAATTCGACAGATTTGCACTGCAAGTGCAGTCAGTAGGAGTTGCTGACTTACGGGCCGGGATGTCGTACGTCCACGG
>Sequence_02
GTGCCGAACTAAGGAGACGTTACAGTACGCACCAGCAGACTCTCACAAAGACTCTGGCTAGTCCGTCGAAACGGCCTGCTAGAACAATGAAAGAGCCACGTCAAAAGAAAACTTCGTTGTACCTAGCGTCAGGTTTCTGCTAGAAACAGCAAGATCGCAGTCGTATGATTGATGGGGTACTCAGCC
...
```

### Context-Dependent Example

```
select 10 representative sequences
```

*Note: This command will use sequences from a previous alignment step or uploaded file.*

## Sequence Types

### Short Sequences (10-50 bp)
- Used for: Quick demonstrations, simple mutations
- Example: `ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG`

### Medium Sequences (40-100 bp)
- Used for: Alignment examples, basic analysis
- Example: Multiple 40-50 bp sequences in FASTA format

### Long Sequences (150+ bp)
- Used for: Phylogenetic tree construction
- Example: 5 sequences of 150-200 bp each for realistic tree building

## User Workflows

### Workflow 1: Standalone Examples
1. User clicks example command
2. Command includes all necessary sequences
3. Command executes immediately
4. Results displayed

### Workflow 2: Multi-Step Workflow
1. User uploads file OR uses example with sequences
2. System stores sequences in workflow context
3. User runs subsequent commands
4. Commands use sequences from context

### Workflow 3: File Upload
1. User uploads FASTA file
2. File content stored
3. Commands reference uploaded sequences
4. No sequences needed in command text

## Implementation Details

### Component Updates

All example components now use the sequence helper:

- **WelcomeScreen**: Quick start examples with sequences
- **EmptyState**: Action cards with complete examples
- **ExampleCommandsPanel**: Categorized examples with sequences

### Helper Function Logic

```typescript
getExampleWithSequences(baseCommand, sequenceType)
```

**Parameters:**
- `baseCommand`: The natural language command
- `sequenceType`: Type of sequences to include ('threeSequences', 'phylogenetic', 'single', etc.)

**Returns:**
- Complete command with sequences appended in appropriate format

### Sequence Format

Sequences are appended in FASTA format:
```
>sequence_name
ATGCGATCGATCG...
```

Or inline for mutation commands:
```
mutate sequence ATGCGATCGATCG...
```

## Best Practices

### For Users

1. **Try Examples First**: All examples are complete and executable
2. **Upload Files**: For your own sequences, upload FASTA files
3. **Chain Commands**: Use workflow context to chain operations
4. **Check Descriptions**: Example descriptions indicate if sequences are included

### For Developers

1. **Use Helper Function**: Always use `getExampleWithSequences()` for examples
2. **Include Descriptions**: Clearly indicate if sequences are included
3. **Context Awareness**: Note when commands need previous steps
4. **Format Properly**: Use FASTA format for multiple sequences

## Example Descriptions

Each example includes a description that clarifies:

- ✅ **"with sample sequences"**: Sequences are included in the command
- ⚠️ **"requires previous step"**: Needs sequences from workflow context
- ℹ️ **"no sequences needed"**: Command doesn't require sequences (e.g., vendor research)

## Testing

All example commands have been tested to ensure:
- ✅ Sequences are properly formatted
- ✅ Commands are executable
- ✅ Results are meaningful
- ✅ Workflow context is preserved

---

**Last Updated:** $(date)  
**Location:** `frontend/src/utils/sampleSequences.ts`  
**Components:** WelcomeScreen, EmptyState, ExampleCommandsPanel



