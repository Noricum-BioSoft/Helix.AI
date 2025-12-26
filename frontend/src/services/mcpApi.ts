import axios from 'axios';

// Use VITE_API_BASE_URL from environment for production builds
// Falls back to localhost for local development
const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:8001';

export const mcpApi = {
  // Health check
  healthCheck: async () => {
    const response = await axios.get(`${API_BASE_URL}/health`);
    return response.data;
  },

  // List available tools
  listTools: async () => {
    const response = await axios.get(`${API_BASE_URL}/mcp/tools`);
    return response.data;
  },

  // General command execution
  executeCommand: async (command: string, sessionId?: string) => {
    const response = await axios.post(`${API_BASE_URL}/execute`, { 
      command,
      session_id: sessionId 
    });
    return response.data;
  },

  // Agent orchestrator command
  agentCommand: async (payload: { prompt: string; session_id?: string | null; files?: Array<{ name: string; content: string }>; }) => {
    const response = await axios.post(`${API_BASE_URL}/agent`, payload);
    return response.data;
  },

  handleNaturalCommand: async (command: string, sessionId?: string) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/handle-natural-command`, {
      command,
      session_id: sessionId
    });
    return response.data;
  },

  // Parse command only
  parseCommand: async (command: string, sessionId?: string) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/parse-command`, {
      command,
      session_id: sessionId
    });
    return response.data;
  },

  // Execute parsed command
  executeParsedCommand: async (parsedCommand: any) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/execute-command`, {
      parsed_command: parsedCommand
    });
    return response.data;
  },

  // Session management
  createSession: async () => {
    const response = await axios.post(`${API_BASE_URL}/create_session`);
    return response.data;
  },

  getSessionInfo: async (sessionId: string) => {
    const response = await axios.get(`${API_BASE_URL}/session/${sessionId}`);
    return response.data;
  },

  // Jobs API (async jobs like EMR FastQC)
  getJob: async (jobId: string) => {
    const response = await axios.get(`${API_BASE_URL}/jobs/${jobId}`);
    return response.data;
  },

  getJobResults: async (jobId: string) => {
    const response = await axios.get(`${API_BASE_URL}/jobs/${jobId}/results`);
    return response.data;
  },

  // ⚠️ DEPRECATED: Direct tool calls - DO NOT USE
  // All commands should go through executeCommand() which routes through the agent
  // These methods are kept for backward compatibility but should not be used in new code
  
  // Sequence alignment (DEPRECATED - use executeCommand instead)
  sequenceAlignment: async (params: { sequences: string; algorithm?: string }, sessionId?: string) => {
    console.warn('⚠️ sequenceAlignment is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/mcp/sequence-alignment`, { ...params, session_id: sessionId });
    return response.data;
  },

  // Mutate sequence (DEPRECATED - use executeCommand instead)
  mutateSequence: async (params: { sequence: string; num_variants?: number; mutation_rate?: number }, sessionId?: string) => {
    console.warn('⚠️ mutateSequence is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/mcp/mutate-sequence`, { ...params, session_id: sessionId });
    return response.data;
  },

  // Analyze sequence data (DEPRECATED - use executeCommand instead)
  analyzeSequenceData: async (params: { data: string; analysis_type?: string }, sessionId?: string) => {
    console.warn('⚠️ analyzeSequenceData is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/mcp/analyze-sequence-data`, { ...params, session_id: sessionId });
    return response.data;
  },

  // Select variants (DEPRECATED - use executeCommand instead)
  selectVariants: async (params: { 
    session_id: string; 
    selection_criteria?: string; 
    num_variants?: number; 
    custom_filters?: any 
  }) => {
    console.warn('⚠️ selectVariants is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/mcp/select-variants`, params);
    return response.data;
  },

  // Visualize alignment (DEPRECATED - use executeCommand instead)
  visualizeAlignment: async (params: { alignment_file: string; output_format?: string }, sessionId?: string) => {
    console.warn('⚠️ visualizeAlignment is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/mcp/visualize-alignment`, { ...params, session_id: sessionId });
    return response.data;
  },

  // Plasmid visualization (DEPRECATED - use executeCommand instead)
  plasmidVisualization: async (params: { 
    vector_name: string; 
    cloning_sites: string; 
    insert_sequence: string 
  }, sessionId?: string) => {
    console.warn('⚠️ plasmidVisualization is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/mcp/plasmid-visualization`, {
      ...params,
      session_id: sessionId
    });
    return response.data;
  },

  // Plasmid for representatives (DEPRECATED - use executeCommand instead)
  plasmidForRepresentatives: async (params: { 
    representatives: string[]; 
    aligned_sequences: string; 
    vector_name?: string; 
    cloning_sites?: string 
  }, sessionId?: string) => {
    console.warn('⚠️ plasmidForRepresentatives is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/mcp/plasmid-for-representatives`, {
      ...params,
      session_id: sessionId
    });
    return response.data;
  },

  // Read trimming (DEPRECATED - use executeCommand instead)
  readTrimming: async (params: { reads: string; adapter?: string; quality_threshold?: number }, sessionId?: string) => {
    console.warn('⚠️ readTrimming is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/mcp/read-trimming`, { ...params, session_id: sessionId });
    return response.data;
  },

  // Read merging (DEPRECATED - use executeCommand instead)
  readMerging: async (params: { forward_reads: string; reverse_reads: string; min_overlap?: number }, sessionId?: string) => {
    console.warn('⚠️ readMerging is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/mcp/read-merging`, { ...params, session_id: sessionId });
    return response.data;
  },

  // Phylogenetic tree (DEPRECATED - use executeCommand instead)
  phylogeneticTree: async (params: { aligned_sequences: string }) => {
    console.warn('⚠️ phylogeneticTree is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/mcp/phylogenetic-tree`, params);
    return response.data;
  },

  // Sequence selection (DEPRECATED - use executeCommand instead)
  sequenceSelection: async (params: { 
    aligned_sequences: string; 
    selection_type?: string; 
    num_sequences?: number 
  }) => {
    console.warn('⚠️ sequenceSelection is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/mcp/sequence-selection`, params);
    return response.data;
  },

  // Synthesis submission (DEPRECATED - use executeCommand instead)
  synthesisSubmission: async (params: { 
    sequences: string; 
    vendor_preference?: string; 
    quantity?: string; 
    delivery_time?: string 
  }) => {
    console.warn('⚠️ synthesisSubmission is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/mcp/synthesis-submission`, params);
    return response.data;
  }
}; 