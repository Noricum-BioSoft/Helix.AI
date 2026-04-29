import axios from 'axios';

export const normalizeBaseUrl = (value?: string): string | undefined => {
  const trimmed = value?.trim();
  if (!trimmed) return undefined;
  return trimmed.replace(/\/+$/, '');
};

// Resolution order:
// - Explicit build-time override via VITE_API_BASE_URL (for separate API domains)
// - Production default: same-origin (CloudFront can route API paths to the backend)
// - Dev default: local backend
export const API_BASE_URL =
  normalizeBaseUrl(import.meta.env.VITE_API_BASE_URL as string | undefined) ??
  (import.meta.env.PROD ? '' : 'http://localhost:8001');

export interface ColumnStat {
  name: string;
  dtype: string;
  n_missing?: number;
  pct_missing?: number;
  n_unique?: number;
  min?: number;
  max?: number;
  mean?: number;
  top_values?: Record<string, number>;
}

export interface SchemaPreview {
  format?: string;
  family?: string;
  n_records?: number | null;
  summary?: Record<string, unknown>;
  schema?: {
    columns?: ColumnStat[];
    obs_columns?: string[];
    var_columns?: string[];
    fields?: string[];
    cell_attributes?: string[];
    gene_attributes?: string[];
  };
  sample?: Record<string, unknown>[];
  available_sheets?: string[] | null;
  profiler_error?: string | null;
}

export interface UploadedFileResponse {
  file_id: string;
  filename: string;
  original_filename?: string;
  size: number;
  content_type?: string;
  local_path?: string;
  uploaded_at?: string;
  schema_preview?: SchemaPreview;
}

export const helixApi = {
  // Health check
  healthCheck: async () => {
    const response = await axios.get(`${API_BASE_URL}/health`);
    return response.data;
  },

  // List available tools
  listTools: async () => {
    const response = await axios.get(`${API_BASE_URL}/tools/list`);
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

  /**
   * Stream a command execution via SSE (POST /execute/stream).
   *
   * Immediately fires `onProgress` with the first "received" event so the UI
   * can show activity, then keeps calling `onProgress` with heartbeat phases
   * until the agent finishes.  `onResult` is called exactly once with the full
   * response payload (same shape as `executeCommand`).  `onError` is called on
   * network or server errors.
   *
   * Returns a cleanup function — call it on component unmount to abort.
   */
  executeCommandStream: (
    command: string,
    sessionId: string | undefined,
    onProgress: (phase: string, message: string) => void,
    onResult: (data: unknown) => void,
    onError: (err: Error) => void,
  ): (() => void) => {
    const controller = new AbortController();

    (async () => {
      try {
        const res = await fetch(`${API_BASE_URL}/execute/stream`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ command, session_id: sessionId }),
          signal: controller.signal,
        });

        if (!res.ok || !res.body) {
          onError(new Error(`HTTP ${res.status}: ${res.statusText}`));
          return;
        }

        const reader = res.body.getReader();
        const decoder = new TextDecoder();
        let buffer = '';

        while (true) {
          const { done, value } = await reader.read();
          if (done) break;

          buffer += decoder.decode(value, { stream: true });
          const lines = buffer.split('\n');
          buffer = lines.pop() ?? '';   // keep any partial line

          for (const line of lines) {
            if (!line.startsWith('data: ')) continue;
            try {
              const event = JSON.parse(line.slice(6));
              if (event.type === 'progress') {
                onProgress(event.phase ?? '', event.message ?? '');
              } else if (event.type === 'result') {
                onResult(event.data);
              } else if (event.type === 'error') {
                onError(new Error(event.detail ?? 'Server error'));
              }
            } catch { /* malformed line — skip */ }
          }
        }
      } catch (err) {
        if ((err as Error).name !== 'AbortError') {
          onError(err instanceof Error ? err : new Error(String(err)));
        }
      }
    })();

    return () => controller.abort();
  },

  // Dispatch a previously-planned pipeline (execute_plan=true flag)
  executePipelinePlan: async (command: string, sessionId?: string) => {
    const response = await axios.post(`${API_BASE_URL}/execute`, {
      command,
      session_id: sessionId,
      execute_plan: true,
    });
    return response.data;
  },

  // Agent orchestrator command
  agentCommand: async (payload: { prompt: string; session_id?: string | null; files?: Array<{ name: string; content: string }>; }) => {
    const response = await axios.post(`${API_BASE_URL}/agent`, payload);
    return response.data;
  },

  handleNaturalCommand: async (command: string, sessionId?: string) => {
    const response = await axios.post(`${API_BASE_URL}/tools/handle-natural-command`, {
      command,
      session_id: sessionId
    });
    return response.data;
  },

  // Parse command only
  parseCommand: async (command: string, sessionId?: string) => {
    const response = await axios.post(`${API_BASE_URL}/tools/parse-command`, {
      command,
      session_id: sessionId
    });
    return response.data;
  },

  // Execute parsed command
  executeParsedCommand: async (parsedCommand: any) => {
    const response = await axios.post(`${API_BASE_URL}/tools/execute-command`, {
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

  uploadSessionFiles: async (sessionId: string, files: File[]) => {
    const formData = new FormData();
    files.forEach((file) => formData.append('files', file));
    // Do NOT set Content-Type manually. When posting FormData the browser must
    // generate the multipart boundary itself; an explicit header without the
    // boundary causes a 400 "Missing boundary in multipart." error.
    const response = await axios.post(
      `${API_BASE_URL}/session/${sessionId}/uploads`,
      formData,
    );
    return response.data as {
      success: boolean;
      session_id: string;
      files: UploadedFileResponse[];
      uploaded_count: number;
      max_upload_mb?: number | null;
    };
  },

  // Jobs API (async jobs like EMR FastQC)
  getJob: async (jobId: string) => {
    const response = await axios.get(`${API_BASE_URL}/jobs/${jobId}`);
    // API returns {success: true, job: {...}}, extract the job object
    return response.data.job || response.data;
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
    const response = await axios.post(`${API_BASE_URL}/tools/sequence-alignment`, { ...params, session_id: sessionId });
    return response.data;
  },

  // Mutate sequence (DEPRECATED - use executeCommand instead)
  mutateSequence: async (params: { sequence: string; num_variants?: number; mutation_rate?: number }, sessionId?: string) => {
    console.warn('⚠️ mutateSequence is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/tools/mutate-sequence`, { ...params, session_id: sessionId });
    return response.data;
  },

  // Analyze sequence data (DEPRECATED - use executeCommand instead)
  analyzeSequenceData: async (params: { data: string; analysis_type?: string }, sessionId?: string) => {
    console.warn('⚠️ analyzeSequenceData is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/tools/analyze-sequence-data`, { ...params, session_id: sessionId });
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
    const response = await axios.post(`${API_BASE_URL}/tools/select-variants`, params);
    return response.data;
  },

  // Visualize alignment (DEPRECATED - use executeCommand instead)
  visualizeAlignment: async (params: { alignment_file: string; output_format?: string }, sessionId?: string) => {
    console.warn('⚠️ visualizeAlignment is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/tools/visualize-alignment`, { ...params, session_id: sessionId });
    return response.data;
  },

  // Plasmid visualization (DEPRECATED - use executeCommand instead)
  plasmidVisualization: async (params: { 
    vector_name: string; 
    cloning_sites: string; 
    insert_sequence: string 
  }, sessionId?: string) => {
    console.warn('⚠️ plasmidVisualization is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/tools/plasmid-visualization`, {
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
    const response = await axios.post(`${API_BASE_URL}/tools/plasmid-for-representatives`, {
      ...params,
      session_id: sessionId
    });
    return response.data;
  },

  // Read trimming (DEPRECATED - use executeCommand instead)
  readTrimming: async (params: { reads: string; adapter?: string; quality_threshold?: number }, sessionId?: string) => {
    console.warn('⚠️ readTrimming is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/tools/read-trimming`, { ...params, session_id: sessionId });
    return response.data;
  },

  // Read merging (DEPRECATED - use executeCommand instead)
  readMerging: async (params: { forward_reads: string; reverse_reads: string; min_overlap?: number }, sessionId?: string) => {
    console.warn('⚠️ readMerging is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/tools/read-merging`, { ...params, session_id: sessionId });
    return response.data;
  },

  // Phylogenetic tree (DEPRECATED - use executeCommand instead)
  phylogeneticTree: async (params: { aligned_sequences: string }) => {
    console.warn('⚠️ phylogeneticTree is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/tools/phylogenetic-tree`, params);
    return response.data;
  },

  // Sequence selection (DEPRECATED - use executeCommand instead)
  sequenceSelection: async (params: { 
    aligned_sequences: string; 
    selection_type?: string; 
    num_sequences?: number 
  }) => {
    console.warn('⚠️ sequenceSelection is deprecated. Use executeCommand() instead.');
    const response = await axios.post(`${API_BASE_URL}/tools/sequence-selection`, params);
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
    const response = await axios.post(`${API_BASE_URL}/tools/synthesis-submission`, params);
    return response.data;
  },

  /**
   * Subscribe to real-time Nextflow pipeline progress via Server-Sent Events.
   *
   * Opens an EventSource to GET /jobs/{jobId}/stream and fires callbacks as
   * Nextflow lifecycle events arrive.  Returns the EventSource so the caller
   * can close it manually if needed (e.g. component unmount).
   *
   * @param jobId        - The job_id returned by POST /execute for pipeline tools.
   * @param onEvent      - Called for each intermediate event (process_completed, etc.).
   * @param onDone       - Called once when the pipeline completes successfully.
   * @param onError      - Called when the pipeline fails or the connection errors.
   * @returns The underlying EventSource instance.
   */
  subscribeJobStream: (
    jobId: string,
    onEvent: (event: { type: string; process?: string; trace?: object }) => void,
    onDone: (result: { job_id: string; result_files: string[] }) => void,
    onError: (err: Error) => void,
  ): EventSource => {
    const source = new EventSource(`${API_BASE_URL}/jobs/${jobId}/stream`);

    // Intermediate progress events
    source.addEventListener('process_submitted',  (e: MessageEvent) => {
      try { onEvent({ type: 'process_submitted',  ...JSON.parse(e.data) }); } catch (_) { /* ignore */ }
    });
    source.addEventListener('process_completed',  (e: MessageEvent) => {
      try { onEvent({ type: 'process_completed',  ...JSON.parse(e.data) }); } catch (_) { /* ignore */ }
    });
    source.addEventListener('started', (e: MessageEvent) => {
      try { onEvent({ type: 'started', ...JSON.parse(e.data) }); } catch (_) { /* ignore */ }
    });
    source.addEventListener('heartbeat', () => {
      // keep-alive ping — no action needed
    });

    // Terminal events
    source.addEventListener('completed', (e: MessageEvent) => {
      try { onDone(JSON.parse(e.data)); } catch (_) { /* ignore */ }
      source.close();
    });
    source.addEventListener('failed', (e: MessageEvent) => {
      try {
        const data = JSON.parse(e.data);
        onError(new Error(data.error || 'Pipeline failed'));
      } catch (_) {
        onError(new Error('Pipeline failed'));
      }
      source.close();
    });
    source.addEventListener('error', (e: Event) => {
      onError(new Error('SSE connection error'));
      source.close();
    });

    return source;
  },
}; 