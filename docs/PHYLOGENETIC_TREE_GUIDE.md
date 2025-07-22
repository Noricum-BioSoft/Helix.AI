# 🌳 Phylogenetic Tree Visualization Guide

## Overview

The Helix.AI platform includes an interactive phylogenetic tree visualization system that allows you to construct and visualize evolutionary relationships between DNA sequences using D3.js.

## Features

### 🎯 **Interactive Tree Rendering**
- **Zoom and Pan**: Navigate large trees with smooth zoom and pan controls
- **Node Highlighting**: Hover over nodes to highlight branches and relationships
- **Responsive Design**: Trees automatically adjust to container size
- **Real-time Updates**: Trees update immediately when new data is available

### 📊 **Tree Statistics**
- **Number of Sequences**: Total sequences in the analysis
- **Alignment Length**: Length of aligned sequences
- **Tree Height**: Maximum distance from root to leaves
- **Number of Clades**: Number of distinct groups in the tree
- **Distance Method**: Algorithm used for distance calculation (Identity)
- **Tree Construction Method**: Algorithm used (UPGMA)

### 🔧 **Technical Implementation**
- **Newick Format**: Standard phylogenetic tree format support
- **D3.js Integration**: Modern JavaScript visualization library
- **Error Handling**: Robust error handling for malformed data
- **TypeScript Support**: Full TypeScript integration for type safety

## Usage

### Basic Tree Construction

```bash
# Command: "visualize the variants in a phylogenetic tree"
```

This command will:
1. Use aligned sequences from the current session
2. Construct a phylogenetic tree using UPGMA method
3. Display the tree with interactive visualization
4. Show comprehensive tree statistics

### Tree Construction Process

1. **Sequence Input**: The system accepts aligned sequences in FASTA format
2. **Distance Calculation**: Uses identity-based distance calculation
3. **Tree Construction**: Applies UPGMA (Unweighted Pair Group Method with Arithmetic Mean) algorithm
4. **Newick Generation**: Converts tree to Newick format for visualization
5. **D3.js Rendering**: Renders interactive tree in the browser

### Example Workflow

```bash
# Step 1: Generate mutations
"generate 50 random mutations from this sequence"

# Step 2: Align sequences
"align the mutated sequences"

# Step 3: Visualize tree
"visualize the variants in a phylogenetic tree"
```

## Tree Interpretation

### Understanding the Visualization

- **Nodes**: Represent sequences or groups of sequences
- **Branches**: Show evolutionary relationships
- **Branch Lengths**: Indicate evolutionary distance
- **Root**: Common ancestor of all sequences
- **Leaves**: Individual sequences at the tips

### Reading Tree Statistics

- **Tree Height**: Lower values indicate more similar sequences
- **Number of Clades**: Higher values indicate more diverse sequences
- **Alignment Length**: Length of the aligned sequence region
- **Distance Method**: Currently uses identity-based distance

## Technical Details

### Newick Format

The system generates trees in Newick format:
```
((seq3:0.00000,seq2:0.00000)Inner1:0.00000,seq1:0.00000)Inner2:0.00000;
```

### D3.js Implementation

- **Tree Layout**: Uses D3's hierarchical layout
- **Coordinate System**: Cartesian coordinates with proper scaling
- **Interactivity**: Zoom, pan, and hover effects
- **Styling**: CSS-based styling for nodes and links

### Error Handling

- **Invalid Newick**: Graceful handling of malformed tree data
- **Empty Data**: Clear error messages for missing data
- **Rendering Errors**: Fallback display for rendering issues
- **Coordinate Validation**: Prevents NaN coordinate errors

## Customization

### Tree Appearance

The tree visualization can be customized by modifying the D3.js code in `frontend/src/components/PhylogeneticTree.tsx`:

- **Node Colors**: Change node and branch colors
- **Font Sizes**: Adjust text size for better readability
- **Tree Orientation**: Switch between horizontal and vertical layouts
- **Animation**: Add smooth transitions and animations

### Performance Optimization

For large trees (>100 sequences):
- **Node Clustering**: Group similar nodes for better performance
- **Lazy Loading**: Load tree data progressively
- **Simplified Rendering**: Use simplified graphics for large trees

## Troubleshooting

### Common Issues

1. **Tree Not Displaying**
   - Check browser console for JavaScript errors
   - Verify Newick string is valid
   - Ensure D3.js is properly loaded

2. **Poor Performance**
   - Reduce number of sequences for large trees
   - Check browser memory usage
   - Consider using simplified tree layout

3. **Incorrect Tree Structure**
   - Verify input sequences are properly aligned
   - Check distance calculation method
   - Validate Newick format

### Debug Information

The system provides debug information in the browser console:
- Newick string parsing
- Tree layout calculations
- Rendering performance metrics
- Error details for troubleshooting

## Future Enhancements

### Planned Features

- **Multiple Tree Methods**: NJ, ML, and Bayesian methods
- **Tree Comparison**: Side-by-side tree comparison
- **Export Options**: PNG, SVG, and PDF export
- **Advanced Statistics**: Bootstrap support and confidence intervals
- **Interactive Editing**: Manual tree editing capabilities

### Performance Improvements

- **WebGL Rendering**: GPU-accelerated rendering for large trees
- **Progressive Loading**: Load tree data in chunks
- **Caching**: Cache tree calculations for repeated analysis
- **Parallel Processing**: Multi-threaded tree construction

## API Reference

### Backend Endpoints

```
POST /execute
{
  "command": "visualize the variants in a phylogenetic tree",
  "session_id": "user_session"
}
```

### Response Format

```json
{
  "text": "Phylogenetic tree constructed successfully...",
  "tree_newick": "((seq3:0.00000,seq2:0.00000)Inner1:0.00000,seq1:0.00000)Inner2:0.00000;",
  "statistics": {
    "num_sequences": 96,
    "alignment_length": 10,
    "tree_height": 0,
    "num_clades": 96,
    "distance_method": "identity",
    "tree_method": "upgma"
  }
}
```

## Related Documentation

- [Natural Language Guide](NATURAL_LANGUAGE_GUIDE.md)
- [Session Management](HISTORY_TRACKING.md)
- [Command Router](COMMAND_ROUTER.md)
- [Frontend Components](../frontend/src/components/PhylogeneticTree.tsx) 