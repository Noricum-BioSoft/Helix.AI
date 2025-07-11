// Command Parser for MCP Integration

export interface ParsedCommand {
  type: 'sequence_alignment' | 'mutate_sequence' | 'analyze_sequence_data' | 'visualize_alignment' | 'dna_vendor_research' | 'general';
  parameters: Record<string, any>;
  originalCommand: string;
}

export class CommandParser {
  private static alignmentKeywords = [
    'align', 'alignment', 'sequence alignment', 'multiple sequence alignment',
    'clustal', 'muscle', 'mafft', 'align sequences'
  ];

  private static mutationKeywords = [
    'mutate', 'mutation', 'variant', 'variants', 'mutate sequence',
    'generate variants', 'create mutations'
  ];

  private static analysisKeywords = [
    'analyze', 'analysis', 'analyze sequence', 'sequence analysis',
    'phylogeny', 'composition', 'bioinformatics analysis'
  ];

  private static visualizationKeywords = [
    'visualize', 'visualization', 'plot', 'graph', 'chart',
    'visualize alignment', 'alignment visualization'
  ];

  private static vendorResearchKeywords = [
    'order', 'vendor', 'synthesis', 'test', 'assay', 'expression', 'function', 'binding',
    'dna synthesis', 'gene synthesis', 'find vendor', 'research vendor',
    'testing options', 'quality control', 'validation'
  ];

  static parseCommand(command: string): ParsedCommand {
    const lowerCommand = command.toLowerCase();
    
    // Check for sequence alignment
    if (this.alignmentKeywords.some(keyword => lowerCommand.includes(keyword))) {
      return this.parseAlignmentCommand(command);
    }
    
    // Check for mutation
    if (this.mutationKeywords.some(keyword => lowerCommand.includes(keyword))) {
      return this.parseMutationCommand(command);
    }
    
    // Check for analysis
    if (this.analysisKeywords.some(keyword => lowerCommand.includes(keyword))) {
      return this.parseAnalysisCommand(command);
    }
    
    // Check for visualization
    if (this.visualizationKeywords.some(keyword => lowerCommand.includes(keyword))) {
      return this.parseVisualizationCommand(command);
    }
    
    // Check for vendor research
    if (this.vendorResearchKeywords.some(keyword => lowerCommand.includes(keyword))) {
      return this.parseVendorResearchCommand(command);
    }
    
    // Default to general command
    return {
      type: 'general',
      parameters: { command },
      originalCommand: command
    };
  }

  private static parseAlignmentCommand(command: string): ParsedCommand {
    const lowerCommand = command.toLowerCase();
    
    // Extract algorithm if specified
    let algorithm = 'clustal';
    if (lowerCommand.includes('muscle')) algorithm = 'muscle';
    else if (lowerCommand.includes('mafft')) algorithm = 'mafft';
    
    // Try to extract sequences from the command
    // This is a simplified parser - in practice, you might want more sophisticated parsing
    const sequences = this.extractSequences(command);
    
    return {
      type: 'sequence_alignment',
      parameters: {
        sequences,
        algorithm
      },
      originalCommand: command
    };
  }

  private static parseMutationCommand(command: string): ParsedCommand {
    const lowerCommand = command.toLowerCase();
    
    // Extract number of variants if specified
    const variantMatch = command.match(/(\d+)\s*variants?/i);
    const numVariants = variantMatch ? parseInt(variantMatch[1]) : 96;
    
    // Extract mutation rate if specified
    const rateMatch = command.match(/(\d*\.?\d+)\s*mutation\s*rate/i);
    const mutationRate = rateMatch ? parseFloat(rateMatch[1]) : 0.1;
    
    // Try to extract sequence
    const sequence = this.extractSequence(command);
    
    return {
      type: 'mutate_sequence',
      parameters: {
        sequence,
        num_variants: numVariants,
        mutation_rate: mutationRate
      },
      originalCommand: command
    };
  }

  private static parseAnalysisCommand(command: string): ParsedCommand {
    const lowerCommand = command.toLowerCase();
    
    // Determine analysis type
    let analysisType = 'alignment';
    if (lowerCommand.includes('phylogeny')) analysisType = 'phylogeny';
    else if (lowerCommand.includes('composition')) analysisType = 'composition';
    
    // Try to extract data
    const data = this.extractData(command);
    
    return {
      type: 'analyze_sequence_data',
      parameters: {
        data,
        analysis_type: analysisType
      },
      originalCommand: command
    };
  }

  private static parseVisualizationCommand(command: string): ParsedCommand {
    const lowerCommand = command.toLowerCase();
    
    // Extract output format
    let outputFormat = 'png';
    if (lowerCommand.includes('svg')) outputFormat = 'svg';
    else if (lowerCommand.includes('pdf')) outputFormat = 'pdf';
    
    // Try to extract file path
    const alignmentFile = this.extractFilePath(command);
    
    return {
      type: 'visualize_alignment',
      parameters: {
        alignment_file: alignmentFile,
        output_format: outputFormat
      },
      originalCommand: command
    };
  }

  private static parseVendorResearchCommand(command: string): ParsedCommand {
    const lowerCommand = command.toLowerCase();
    
    // Extract vendor type
    let vendorType = 'dna_synthesis';
    if (lowerCommand.includes('gene_synthesis')) vendorType = 'gene_synthesis';
    else if (lowerCommand.includes('find_vendor')) vendorType = 'find_vendor';
    else if (lowerCommand.includes('research_vendor')) vendorType = 'research_vendor';
    
    // Extract specific details
    const details: Record<string, any> = {};
    const detailMatch = command.match(/(\w+)\s*:\s*(.+)/i);
    if (detailMatch) {
      details[detailMatch[1].toLowerCase()] = detailMatch[2];
    }

    return {
      type: 'dna_vendor_research',
      parameters: {
        vendor_type: vendorType,
        details
      },
      originalCommand: command
    };
  }

  private static extractSequences(command: string): string {
    // Look for FASTA format sequences in the command
    const fastaMatch = command.match(/(>[^\n]+\n[ATCGN\-\s]+)/gi);
    if (fastaMatch) {
      return fastaMatch.join('\n');
    }
    
    // Look for simple sequences
    const sequenceMatch = command.match(/([ATCGN]+)/gi);
    if (sequenceMatch) {
      return sequenceMatch.map((seq, i) => `>seq${i + 1}\n${seq}`).join('\n');
    }
    
    // Default sequences for testing
    return `>seq1
ACTGTTGAC
>seq2
ACTGCATCC
>seq3
ACTGCAATGAC`;
  }

  private static extractSequence(command: string): string {
    // Look for a single sequence
    const sequenceMatch = command.match(/([ATCGN]+)/i);
    if (sequenceMatch) {
      return sequenceMatch[1];
    }
    
    // Default sequence for testing
    return 'ACTGTTGAC';
  }

  private static extractData(command: string): string {
    // Look for file paths
    const fileMatch = command.match(/(\w+\.(csv|fasta|fa|txt))/i);
    if (fileMatch) {
      return fileMatch[1];
    }
    
    // Look for FASTA content
    const fastaMatch = command.match(/(>[^\n]+\n[ATCGN\-\s]+)/gi);
    if (fastaMatch) {
      return fastaMatch.join('\n');
    }
    
    // Default data for testing
    return `>seq1
ACTGTTGAC
>seq2
ACTGCATCC`;
  }

  private static extractFilePath(command: string): string {
    // Look for file paths
    const fileMatch = command.match(/(\w+\.(fasta|fa|aln))/i);
    if (fileMatch) {
      return fileMatch[1];
    }
    
    // Default file path for testing
    return 'tools/seqs.aln.fasta';
  }
} 