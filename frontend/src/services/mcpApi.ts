// MCP API Service for Bioinformatics Operations

const API_BASE_URL = 'http://localhost:8001';

export interface MCPToolRequest {
  tool_name: string;
  arguments: Record<string, any>;
}

export interface MCPResponse {
  success: boolean;
  result: Record<string, any>;
  error?: string;
}

export interface SequenceAlignmentRequest {
  sequences: string;
  algorithm?: string;
}

export interface MutationRequest {
  sequence: string;
  num_variants?: number;
  mutation_rate?: number;
}

export interface AnalysisRequest {
  data: string;
  analysis_type?: string;
}

export interface VisualizationRequest {
  alignment_file: string;
  output_format?: string;
}

class MCPApiService {
  private async makeRequest<T>(
    endpoint: string,
    method: 'GET' | 'POST' = 'POST',
    data?: any
  ): Promise<T> {
    const response = await fetch(`${API_BASE_URL}${endpoint}`, {
      method,
      headers: {
        'Content-Type': 'application/json',
      },
      body: data ? JSON.stringify(data) : undefined,
    });

    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }

    return response.json();
  }

  // Legacy endpoint for backward compatibility
  async executeCommand(command: string): Promise<any> {
    return this.makeRequest('/execute', 'POST', { command });
  }

  // MCP-specific endpoints
  async sequenceAlignment(request: SequenceAlignmentRequest): Promise<MCPResponse> {
    return this.makeRequest<MCPResponse>('/mcp/sequence-alignment', 'POST', request);
  }

  async mutateSequence(request: MutationRequest): Promise<MCPResponse> {
    return this.makeRequest<MCPResponse>('/mcp/mutate-sequence', 'POST', request);
  }

  async analyzeSequenceData(request: AnalysisRequest): Promise<MCPResponse> {
    return this.makeRequest<MCPResponse>('/mcp/analyze-sequence-data', 'POST', request);
  }

  async visualizeAlignment(request: VisualizationRequest): Promise<MCPResponse> {
    return this.makeRequest<MCPResponse>('/mcp/visualize-alignment', 'POST', request);
  }

  async listTools(): Promise<{ tools: Array<{ name: string; description: string; parameters: Record<string, string> }> }> {
    return this.makeRequest('/mcp/tools', 'GET');
  }

  async healthCheck(): Promise<{ status: string; service: string }> {
    return this.makeRequest('/health', 'GET');
  }
}

export const mcpApi = new MCPApiService(); 