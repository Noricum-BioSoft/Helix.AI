// Sample DNA sequences for examples and testing
// These are realistic DNA sequences that can be used in example commands

export const sampleSequences = {
  // Short sequences for quick examples
  short: {
    seq1: "ATGCGATCG",
    seq2: "ATGCGATC",
    seq3: "ATGCGATCGATCG",
  },
  
  // Medium sequences for alignment examples
  medium: {
    seq1: "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
    seq2: "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
    seq3: "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
  },
  
  // FASTA format sequences for alignment
  fasta: {
    twoSequences: `>seq1
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq2
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG`,
    
    threeSequences: `>seq1
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq2
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq3
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG`,
    
    // Longer sequences for phylogenetic tree
    phylogenetic: `>Sequence_01
ACTCGATCACAAAGCTTAGGTCCGATCAATTTTGATAGTTACCCCCCACGGTCCAATCCGTTGGGTGAACACCGAGAAATTCGACAGATTTGCACTGCAAGTGCAGTCAGTAGGAGTTGCTGACTTACGGGCCGGGATGTCGTACGTCCACGG
>Sequence_02
GTGCCGAACTAAGGAGACGTTACAGTACGCACCAGCAGACTCTCACAAAGACTCTGGCTAGTCCGTCGAAACGGCCTGCTAGAACAATGAAAGAGCCACGTCAAAAGAAAACTTCGTTGTACCTAGCGTCAGGTTTCTGCTAGAAACAGCAAGATCGCAGTCGTATGATTGATGGGGTACTCAGCC
>Sequence_03
TTCTCACACTGTGTAAAAATTACACAAAAGATACGCCCAGTATTGGGGTGGGTATCCTCCGGGATGGGTAACTGGGGGTTCCCTTATGGTCAATGGAAAACCAGCCAAGATACATCTCATTGTTATAGGATGTTGAGCGCCATTAGCCTGCGATCACTGGGCGCCGTTTTTTCAACGTTTCTCCTCAC
>Sequence_04
TTTGTGTGTTCACCTGGTGTCCAACAATTCGATGGATCATTGGGCCGATCCGTTAGCGCCGAACGCGAGTGTTGGGAGTTTTCTGCCGTCGACGTCGTCGAGTGAAATATCAAGCCCTGCAGGTCGACTGCGGCGTGTTGACCGTTAGTGGTTTACAATGGCTGTTAACGTTAATTGCAGGTACCCTGCAG
>Sequence_05
TAAATGACGTCAGACTCTCTTAGTTATGCTCCGACTGGCTTTTACAGTTTCTTATAATAGGCTAGCAGCAAGAGGGTCCCGGTGTCCGTTTGGTTATCCTCGTTCTCAGTTGGTGAATAGGGACGCGGCATATGTACGGCAACGTATAATGT`,
  },
  
  // Single sequence for mutation examples
  single: "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
  
  // Plasmid insert sequence
  plasmidInsert: "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
};

// Helper function to format sequences for commands
export const formatSequencesForCommand = (sequences: string, command: string): string => {
  // If command already contains sequences, return as is
  if (command.includes('>seq') || command.includes('>Sequence')) {
    return command;
  }
  
  // Otherwise, append sequences to the command
  return `${command}\n\n${sequences}`;
};

// Helper function to get example command with sequences
export const getExampleWithSequences = (baseCommand: string, sequenceType: keyof typeof sampleSequences.fasta = 'threeSequences'): string => {
  if (baseCommand.includes('align') || baseCommand.includes('sequence alignment')) {
    return `${baseCommand}\n\n${sampleSequences.fasta[sequenceType]}`;
  }
  
  if (baseCommand.includes('mutate') || baseCommand.includes('variant')) {
    // For mutation commands, include sequence in the command text itself
    if (baseCommand.includes('sequence')) {
      // If command already mentions sequence, append it
      return `${baseCommand}\n\n${sampleSequences.single}`;
    } else {
      // Otherwise, add sequence parameter
      return `${baseCommand}\n\nSequence: ${sampleSequences.single}`;
    }
  }
  
  if (baseCommand.includes('visualize') && baseCommand.includes('tree')) {
    return `${baseCommand}\n\n${sampleSequences.fasta.phylogenetic}`;
  }
  
  if (baseCommand.includes('plasmid') || baseCommand.includes('insert')) {
    // For plasmid commands, include insert in the command if not already present
    if (baseCommand.includes('insert')) {
      return baseCommand; // Command already has insert
    } else {
      return `${baseCommand}\n\nInsert sequence: ${sampleSequences.plasmidInsert}`;
    }
  }
  
  // For commands that need sequences from context, return base command
  // (they'll use sequences from previous steps or uploaded files)
  return baseCommand;
};

