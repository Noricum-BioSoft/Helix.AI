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

  // Sequence alignment
  sequenceAlignment: async (params: { sequences: string; algorithm?: string }, sessionId?: string) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/sequence-alignment`, { ...params, session_id: sessionId });
    return response.data;
  },

  // Mutate sequence
  mutateSequence: async (params: { sequence: string; num_variants?: number; mutation_rate?: number }, sessionId?: string) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/mutate-sequence`, { ...params, session_id: sessionId });
    return response.data;
  },

  // Analyze sequence data
  analyzeSequenceData: async (params: { data: string; analysis_type?: string }, sessionId?: string) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/analyze-sequence-data`, { ...params, session_id: sessionId });
    return response.data;
  },

  // Select variants
  selectVariants: async (params: { 
    session_id: string; 
    selection_criteria?: string; 
    num_variants?: number; 
    custom_filters?: any 
  }) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/select-variants`, params);
    return response.data;
  },

  // Visualize alignment
  visualizeAlignment: async (params: { alignment_file: string; output_format?: string }, sessionId?: string) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/visualize-alignment`, { ...params, session_id: sessionId });
    return response.data;
  },

  // Plasmid visualization
  plasmidVisualization: async (params: { 
    vector_name: string; 
    cloning_sites: string; 
    insert_sequence: string 
  }, sessionId?: string) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/plasmid-visualization`, {
      ...params,
      session_id: sessionId
    });
    return response.data;
  },

  // Plasmid for representatives
  plasmidForRepresentatives: async (params: { 
    representatives: string[]; 
    aligned_sequences: string; 
    vector_name?: string; 
    cloning_sites?: string 
  }, sessionId?: string) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/plasmid-for-representatives`, {
      ...params,
      session_id: sessionId
    });
    return response.data;
  },

  // Phylogenetic tree

  readTrimming: async (params: { reads: string; adapter?: string; quality_threshold?: number }, sessionId?: string) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/read-trimming`, { ...params, session_id: sessionId });
    return response.data;
  },

  readMerging: async (params: { forward_reads: string; reverse_reads: string; min_overlap?: number }, sessionId?: string) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/read-merging`, { ...params, session_id: sessionId });
    return response.data;
  },

  phylogeneticTree: async (params: { aligned_sequences: string }) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/phylogenetic-tree`, params);
    return response.data;
  },

  // Sequence selection
  sequenceSelection: async (params: { 
    aligned_sequences: string; 
    selection_type?: string; 
    num_sequences?: number 
  }) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/sequence-selection`, params);
    return response.data;
  },

  // Synthesis submission
  synthesisSubmission: async (params: { 
    sequences: string; 
    vendor_preference?: string; 
    quantity?: string; 
    delivery_time?: string 
  }) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/synthesis-submission`, params);
    return response.data;
  }
}; 