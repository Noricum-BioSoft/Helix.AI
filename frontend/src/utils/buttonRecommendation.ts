/**
 * Determines which button (Agent or Send Prompt) should be recommended based on the command.
 * 
 * @param command - The user's command text
 * @returns 'agent' | 'direct' | null - Which button to recommend, or null if no clear recommendation
 */
export function getRecommendedButton(command: string): 'agent' | 'direct' | null {
  if (!command || !command.trim()) {
    return null;
  }

  const trimmedCommand = command.trim().toLowerCase();

  // Questions and informational queries → Agent button
  const questionPatterns = [
    /^(what|how|why|when|where|who|which|can|could|should|would|is|are|does|do|will|tell me|explain|describe|help me|i need help|i want to know)/i,
    /\?$/, // Ends with question mark
  ];

  // Action commands → Direct/Send Prompt button
  const actionPatterns = [
    /^(visualize|align|create|generate|perform|mutate|select|insert|show|display|build|run|execute|analyze|compare|find|search|calculate|compute|trim|merge|assess|quality|phylogenetic|tree|plasmid|sequence|variant|mutation)/i,
    /^(visualize|align|create|generate|perform|mutate|select|insert|show|display|build|run|execute|analyze|compare|find|search|calculate|compute|trim|merge|assess|quality|phylogenetic|tree|plasmid|sequence|variant|mutation)\s+/i,
  ];

  // Check for questions first (higher priority)
  for (const pattern of questionPatterns) {
    if (pattern.test(trimmedCommand)) {
      return 'agent';
    }
  }

  // Check for action commands
  for (const pattern of actionPatterns) {
    if (pattern.test(trimmedCommand)) {
      return 'direct';
    }
  }

  // Default: if command is short and looks like a question, recommend agent
  // Otherwise, recommend direct for longer commands that seem like instructions
  if (trimmedCommand.length < 50 && /^(what|how|why|when|where|who|which|can|could|should|would|is|are|does|do|will)/i.test(trimmedCommand)) {
    return 'agent';
  }

  // Default to direct for action-like commands
  return 'direct';
}

