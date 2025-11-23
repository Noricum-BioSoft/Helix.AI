# RNA-seq Demo Dataset

This directory contains sample FASTQ files for demonstrating RNA-seq preprocessing workflows in Helix.AI.

## Files

- **sample_R1.fastq**: Forward (R1) paired-end reads
- **sample_R2.fastq**: Reverse (R2) paired-end reads

## Dataset Details

- **Read Count**: 10 paired-end reads
- **Read Length**: ~70-75 bases per read
- **Format**: Standard FASTQ (4 lines per read)
- **Quality Encoding**: Phred+33 (Sanger/Illumina 1.8+)
- **Purpose**: Demonstration of read trimming, adapter removal, and read merging

## Usage

These files are used in the RNA-seq preprocessing workflow demo. See:
- [RNA-seq Preprocessing Demo](../../docs/demos/RNASEQ_PREPROCESSING_DEMO.md)
- [Command Script](../../docs/demos/RNASEQ_COMMANDS.txt)

## Workflow

1. Upload both FASTQ files to Helix.AI
2. Run quality trimming commands
3. Remove adapter sequences
4. Merge paired-end reads
5. Generate quality reports

## Notes

- These are synthetic sequences designed for demonstration
- Quality scores are uniform (all 'I' = Phred 40) for simplicity
- In real data, quality scores would vary and degrade toward read ends
- The sequences are designed to have overlapping regions for successful merging


