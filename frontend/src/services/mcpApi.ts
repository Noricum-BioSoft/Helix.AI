import axios from 'axios';

const API_BASE_URL = 'http://localhost:8001';

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
  executeCommand: async (command: string) => {
    const response = await axios.post(`${API_BASE_URL}/execute`, { command });
    return response.data;
  },

  // Natural language command handling
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
  createSession: async (userId?: string) => {
    const response = await axios.post(`${API_BASE_URL}/session/create`, {
      user_id: userId
    });
    return response.data;
  },

  getSessionInfo: async (sessionId: string) => {
    const response = await axios.get(`${API_BASE_URL}/session/${sessionId}`);
    return response.data;
  },

  // Sequence alignment
  sequenceAlignment: async (params: { sequences: string; algorithm?: string }) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/sequence-alignment`, params);
    return response.data;
  },

  // Mutate sequence
  mutateSequence: async (params: { sequence: string; num_variants?: number; mutation_rate?: number }) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/mutate-sequence`, params);
    return response.data;
  },

  // Analyze sequence data
  analyzeSequenceData: async (params: { data: string; analysis_type?: string }) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/analyze-sequence-data`, params);
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
  visualizeAlignment: async (params: { alignment_file: string; output_format?: string }) => {
    const response = await axios.post(`${API_BASE_URL}/mcp/visualize-alignment`, params);
    return response.data;
  }
}; 