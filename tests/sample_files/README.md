# Sample Files for DataBloom.AI Testing

This directory contains sample files for testing the drag-and-drop upload functionality in the DataBloom.AI frontend.

## Files

### `sample_sequences.fasta`
- **Format**: FASTA format
- **Content**: 4 DNA sequences with headers
- **Use Case**: Test sequence alignment functionality
- **How to test**: Drag and drop into the input area, select "Align Sequences"

### `sample_sequences.csv`
- **Format**: CSV format
- **Content**: 4 DNA sequences with column headers
- **Use Case**: Test sequence analysis functionality
- **How to test**: Drag and drop into the input area, select "Analyze Data"

### `sample_sequences_with_data.csv`
- **Format**: CSV format with additional experimental data
- **Content**: 4 DNA sequences with OD600 measurements
- **Use Case**: Test data analysis with experimental measurements
- **How to test**: Drag and drop into the input area, select "Analyze Data"

### `single_sequence.txt`
- **Format**: Plain text
- **Content**: Single DNA sequence
- **Use Case**: Test mutation functionality
- **How to test**: Drag and drop into the input area, select "Mutate Sequence"

## Testing Instructions

1. Start the DataBloom.AI application using `./start-app.sh`
2. Open the frontend at `http://localhost:5173`
3. Drag and drop any of these files into the input area
4. Select the appropriate action from the modal that appears:
   - **Align Sequences**: For FASTA files or multiple sequences
   - **Mutate Sequence**: For single sequences
   - **Analyze Data**: For CSV files with data
   - **Cancel**: To abort the operation

## Expected Results

- **FASTA files**: Should trigger sequence alignment with multiple sequences
- **CSV files**: Should trigger data analysis with statistics and visualizations
- **Single sequences**: Should trigger mutation generation with variants 