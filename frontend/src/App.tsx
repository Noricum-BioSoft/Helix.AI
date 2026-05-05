import React, { useState, useEffect, useMemo, useRef } from 'react';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import { MSAView } from 'react-msaview';
import Plot from 'react-plotly.js';
import { API_BASE_URL, helixApi } from './services/helixApi';
import type { UploadedFileResponse } from './services/helixApi';
import { SchemaPreviewPanel } from './components/SchemaPreviewPanel';
import { AnalysisPlanCard } from './components/AnalysisPlanCard';
import { WorkflowPlanCard } from './components/WorkflowPlanCard';
import { NeedsInputsCard } from './components/NeedsInputsCard';
import { CapabilityGrid } from './components/CapabilityGrid';
import { FollowUpChips } from './components/FollowUpChips';
import { getContextualPlaceholder } from './utils/followUpSuggestions';
import { CommandParser, ParsedCommand } from './utils/commandParser';
import { PlasmidDataVisualizer, PlasmidRepresentativesVisualizer } from './components/PlasmidVisualizer';
import { PhylogeneticTree } from './components/PhylogeneticTree';

import 'bootstrap/dist/css/bootstrap.min.css';
import './theme.css';
import Button from 'react-bootstrap/Button';
import Card from 'react-bootstrap/Card';
import Modal from 'react-bootstrap/Modal';
import { DesignOptionOne, DesignOptionTwo, DesignOptionThree } from './components/designs';
import type { PromptDesignProps, QuickExample, UploadedSessionFile } from './components/designs';
import { getExampleWithSequences, sampleSequences } from './utils/sampleSequences';
import { theme } from './theme';
import { JobsPanel } from './components/JobsPanel';
import { ExamplesPanel } from './components/ExamplesPanel';
import { DemoScenariosPanel } from './components/DemoScenariosPanel';
import { getDemoScenarioById, getDemoScenarioByTool, getDemoScenarioByCommandAndTool, getDemoScenarioByCommand, DataPreviewTable } from './data/demoScenarios';
import { ThinkingIndicator, ActivityIndicator } from './components/ThinkingIndicator';

const SESSION_STORAGE_KEY = 'helix_session_id';
import {
  MAX_UPLOAD_BYTES,
  MAX_UPLOAD_MB,
  ALLOWED_UPLOAD_EXTENSIONS,
  validateSelectedFiles,
} from './utils/uploadValidation';

interface HistoryItem {
  input: string;
  output: any;
  type: string;
  timestamp: Date;
  /** Set when the item was triggered by a demo scenario card — links back to followUpPrompt */
  scenarioId?: string;
  /** Bioinformatics run tracking IDs, populated from BioOrchestrator responses */
  run_id?: string;
  parent_run_id?: string;
  /** Temporary ID used to match a 'pending' placeholder with its resolved response */
  pendingId?: string;
}

const BIOINF_LOADING_MESSAGES = [
  'Analyzing your request…',
  'Consulting the bioinformatics knowledge base…',
  'Evaluating workflow requirements…',
  'Identifying required tools and inputs…',
  'Preparing a detailed plan…',
  'Checking tool dependencies…',
  'Almost ready…',
];

/** Map backend session history entries to frontend HistoryItem[] (newest first). */
function sessionHistoryToItems(session: { history?: Array<{ command?: string; tool?: string; result?: any; metadata?: any; run_id?: string; timestamp?: string }> }): HistoryItem[] {
  const raw = session?.history;
  if (!Array.isArray(raw) || raw.length === 0) return [];
  const items: HistoryItem[] = [];
  for (const entry of raw) {
    const cmd = entry.command ?? '';
    const out = entry.result ?? {};
    const meta = entry.metadata ?? {};
    const runId = entry.run_id ?? meta.run_id;
    const parentRunId = meta.parent_run_id;
    const ts = entry.timestamp || meta.timestamp;
    items.push({
      input: cmd,
      output: out,
      type: 'agent',
      timestamp: ts ? new Date(ts) : new Date(0),
      run_id: runId,
      parent_run_id: parentRunId,
    });
  }
  return items.reverse(); // backend order is oldest-first; UI shows newest first
}

interface WorkflowContext {
  alignedSequences?: string;
  selectedSequences?: string[];
  mutatedSequences?: string[];
  plasmidData?: any;
}

const DESIGN_OPTIONS = [
  { id: 'classic', label: 'Classic', component: null },
  { id: 'integrated', label: 'Integrated Panel', component: DesignOptionOne },
  { id: 'split', label: 'Split View', component: DesignOptionTwo },
  { id: 'workspace', label: 'Workspace Tabs', component: DesignOptionThree },
] as const;

type DesignOptionId = (typeof DESIGN_OPTIONS)[number]['id'];


function App() {
  const [command, setCommand] = useState('');
  const [history, setHistory] = useState<HistoryItem[]>([]);
  const [loading, setLoading] = useState(false);
  const [agentLoading, setAgentLoading] = useState(false);
  const [loadingMsg, setLoadingMsg] = useState(BIOINF_LOADING_MESSAGES[0]);
  const [serverStatus, setServerStatus] = useState<string>('unknown');
  const [dragActive, setDragActive] = useState(false);
  const [sessionId, setSessionId] = useState<string | null>(null);
  const [commandMode, setCommandMode] = useState<'structured' | 'natural'>('natural');
  const [activeTab, setActiveTab] = useState<string>('command');
  const [uploadedFiles, setUploadedFiles] = useState<UploadedSessionFile[]>([]);
  const [workflowContext, setWorkflowContext] = useState<WorkflowContext>({});
  const [isInitialized, setIsInitialized] = useState(false);
  // Track which pipeline commands have already been submitted so the Execute button
  // on the original plan card turns into "Submitted ✓" after the user clicks it.
  const [executedPipelineCommands, setExecutedPipelineCommands] = useState<Set<string>>(new Set());
  const historyTopRef = useRef<HTMLDivElement>(null);
  /** Scrolled-to ref for the pending loading card so it stays in view regardless of scroll position. */
  const pendingItemRef = useRef<HTMLDivElement>(null);
  const [selectedDesign, setSelectedDesign] = useState<DesignOptionId>('integrated');
  const [examplesOpen, setExamplesOpen] = useState(false);
  const [jobsOpen, setJobsOpen] = useState(false);
  const [demoOpen, setDemoOpen] = useState(false);
  /** Holds the scenario id of the most-recently loaded demo prompt, cleared after one submission. */
  const [pendingScenarioId, setPendingScenarioId] = useState<string | undefined>(undefined);
  /** Data-preview modal state */
  const [previewTables, setPreviewTables] = useState<DataPreviewTable[] | null>(null);
  const [sessionModalOpen, setSessionModalOpen] = useState(false);
  const [sessionInfo, setSessionInfo] = useState<any>(null);
  const [sessionInfoLoading, setSessionInfoLoading] = useState(false);
  const fileInputRef = useRef<HTMLInputElement>(null);
  const [activities, setActivities] = useState<Array<{
    id: string;
    message: string;
    status: 'active' | 'completed' | 'error';
    timestamp: Date;
  }>>([]);
  const SelectedDesignComponent = useMemo(() => {
    const option = DESIGN_OPTIONS.find(option => option.id === selectedDesign);
    return option?.component;
  }, [selectedDesign]);

  // Activity tracking helpers
  const addActivity = (message: string) => {
    const id = `activity-${Date.now()}-${Math.random()}`;
    setActivities(prev => [...prev, {
      id,
      message,
      status: 'active',
      timestamp: new Date()
    }]);
    return id;
  };

  const updateActivity = (id: string, status: 'completed' | 'error') => {
    setActivities(prev => prev.map(activity => 
      activity.id === id ? { ...activity, status } : activity
    ));
    
    // Auto-remove completed/error activities after 3 seconds
    setTimeout(() => {
      setActivities(prev => prev.filter(activity => activity.id !== id));
    }, 3000);
  };

  const handleExampleClick = (exampleCommand: string, scenarioId?: string) => {
    setCommand(exampleCommand);
    setPendingScenarioId(scenarioId);
  };

  const handleToggleExamples = () => {
    setExamplesOpen(prev => !prev);
  };

  const handleToggleJobs = () => {
    setJobsOpen(prev => !prev);
  };

  const handleSessionClick = async () => {
    if (!sessionId) return;
    setSessionModalOpen(true);
    setSessionInfoLoading(true);
    try {
      const info = await helixApi.getSessionInfo(sessionId);
      setSessionInfo(info);
    } catch (error) {
      console.error('Failed to fetch session info:', error);
      setSessionInfo({ error: 'Failed to fetch session information' });
    } finally {
      setSessionInfoLoading(false);
    }
  };

  const handleNewSession = async () => {
    try {
      const res = await helixApi.createSession() as { session_id: string };
      setSessionId(res.session_id);
      try {
        sessionStorage.setItem(SESSION_STORAGE_KEY, res.session_id);
      } catch {
        // ignore storage errors
      }
    } catch {
      setSessionId(null);
      try {
        sessionStorage.removeItem(SESSION_STORAGE_KEY);
      } catch {}
    }
    setHistory([]);
    setExecutedPipelineCommands(new Set());
    setCommand('');
    setUploadedFiles([]);
    setWorkflowContext({});
  };

  const extractJobIds = (obj: any, out: Set<string>) => {
    if (!obj) return;
    if (typeof obj === 'object') {
      if (typeof obj.job_id === 'string') out.add(obj.job_id);
      if (typeof obj.jobId === 'string') out.add(obj.jobId);
      for (const v of Array.isArray(obj) ? obj : Object.values(obj)) {
        extractJobIds(v, out);
      }
    }
  };

  const jobIds = useMemo(() => {
    const s = new Set<string>();
    for (const item of history) extractJobIds(item.output, s);
    return Array.from(s);
  }, [history]);
  // On mount: restore session from storage if valid (and restore history), otherwise create new; check server health
  // While loading: cycle through BIOINF_LOADING_MESSAGES client-side so the
  // user always sees activity even when the backend responds faster than the
  // 2.5 s SSE heartbeat interval.  SSE progress events can still override the
  // message with server-provided text via setLoadingMsg().
  // When loading ends: reset to the first message for the next request.
  useEffect(() => {
    if (!loading) {
      setLoadingMsg(BIOINF_LOADING_MESSAGES[0]);
      return;
    }
    let idx = 0;
    const timer = setInterval(() => {
      idx = (idx + 1) % BIOINF_LOADING_MESSAGES.length;
      setLoadingMsg(BIOINF_LOADING_MESSAGES[idx]);
    }, 1800);
    return () => clearInterval(timer);
  }, [loading]);

  useEffect(() => {
    const initializeApp = async () => {
      if (isInitialized) return; // Prevent duplicate initialization
      
      try {
        setIsInitialized(true);
        await checkServerHealth();
        let id: string | null = null;
        let sessionData: { session?: { history?: unknown[] } } | null = null;
        try {
          const stored = sessionStorage.getItem(SESSION_STORAGE_KEY);
          if (stored && stored.trim()) {
            const info = await helixApi.getSessionInfo(stored) as { session?: { history?: unknown[] } };
            id = stored;
            sessionData = info;
            console.log('Session restored from storage:', stored);
          }
        } catch {
          try {
            sessionStorage.removeItem(SESSION_STORAGE_KEY);
          } catch {}
        }
        if (!id) {
          const res = await helixApi.createSession() as { session_id: string };
          id = res.session_id;
          try {
            sessionStorage.setItem(SESSION_STORAGE_KEY, id);
          } catch {}
          console.log('Session created:', id);
        }
        setSessionId(id);
        // Restore conversation history when we re-open the same session (e.g. after refresh)
        if (sessionData?.session) {
          const restored = sessionHistoryToItems(sessionData.session as Parameters<typeof sessionHistoryToItems>[0]);
          if (restored.length > 0) {
            setHistory(restored);
            console.log('Restored', restored.length, 'history items from session');
          }
        }
      } catch (err) {
        console.error('Failed to initialize app:', err);
        setIsInitialized(false); // Reset on error
      }
    };
    initializeApp();
  }, [isInitialized]);

  const checkServerHealth = async () => {
    try {
      const health = await helixApi.healthCheck();
      const status = (health as any)?.status;
      if (typeof status !== 'string' || !status.trim()) {
        throw new Error('Invalid health response');
      }
      setServerStatus(status);
    } catch (error) {
      setServerStatus('error');
      console.error('Server health check failed:', error);
    }
  };

  const mapUploadedApiFile = (file: UploadedFileResponse): UploadedSessionFile => ({
    file_id: file.file_id,
    name: file.original_filename || file.filename,
    size: file.size,
    status: 'uploaded',
    local_path: file.local_path,
    uploaded_at: file.uploaded_at,
    schema_preview: file.schema_preview,
  });

  const ensureSessionId = async (): Promise<string> => {
    if (sessionId) return sessionId;
    const res = await helixApi.createSession() as { session_id: string };
    setSessionId(res.session_id);
    try {
      sessionStorage.setItem(SESSION_STORAGE_KEY, res.session_id);
    } catch {
      // ignore storage errors
    }
    return res.session_id;
  };

  const handleUploadFiles = async (files: File[]) => {
    if (!files.length) return;

    const { accepted, rejected } = validateSelectedFiles(files);
    if (rejected.length > 0) {
      setUploadedFiles((prev) => [...prev, ...rejected]);
    }
    if (accepted.length === 0) {
      return;
    }

    const uploading: UploadedSessionFile[] = accepted.map((file) => ({
      name: file.name,
      size: file.size,
      status: 'uploading',
    }));
    setUploadedFiles((prev) => [...prev, ...uploading]);

    try {
      const sid = await ensureSessionId();
      const response = await helixApi.uploadSessionFiles(sid, accepted);
      const uploaded = (response.files || []).map(mapUploadedApiFile);
      setUploadedFiles((prev) => {
        const withoutPending = prev.filter(
          (entry) => !(entry.status === 'uploading' && accepted.some((file) => file.name === entry.name && file.size === entry.size)),
        );
        return [...withoutPending, ...uploaded];
      });
    } catch (error: unknown) {
      const detail =
        (error as { response?: { data?: { detail?: string } } })?.response?.data?.detail
        ?? (error instanceof Error ? error.message : String(error));
      console.error('Upload failed:', error);
      setUploadedFiles((prev) => {
        const withoutPending = prev.filter(
          (entry) => !(entry.status === 'uploading' && accepted.some((file) => file.name === entry.name && file.size === entry.size)),
        );
        const failed = accepted.map((file) => ({
          name: file.name,
          size: file.size,
          status: 'failed' as const,
          error: `Upload failed: ${detail}`,
        }));
        return [...withoutPending, ...failed];
      });
    }
  };

  const updateWorkflowContext = (stepType: string, data: any) => {
    setWorkflowContext(prev => {
      switch (stepType) {
        case 'sequence_alignment':
          return { ...prev, alignedSequences: data };
        case 'sequence_selection':
          return { ...prev, selectedSequences: data };
        case 'mutate_sequence':
          return { ...prev, mutatedSequences: data };
        case 'plasmid_visualization':
          return { ...prev, plasmidData: data };
        default:
          return prev;
      }
    });
  };

  const enhanceCommandWithContext = (command: string): string => {
    const lowerCommand = command.toLowerCase();
    
    console.log('Enhancing command:', command);
    console.log('Current workflow context:', workflowContext);
    
    // If command mentions "aligned sequences" and we have them in context
    if ((lowerCommand.includes('aligned sequences') || lowerCommand.includes('select sequence')) && workflowContext.alignedSequences) {
      const enhancedCommand = `${command}\n\nAligned sequences from previous step:\n${workflowContext.alignedSequences}`;
      console.log('Enhanced command with aligned sequences:', enhancedCommand);
      return enhancedCommand;
    }
    
    // If command mentions "selected sequences" and we have them in context
    if (lowerCommand.includes('selected sequences') && workflowContext.selectedSequences) {
      const sequencesText = workflowContext.selectedSequences.map((seq, i) => `>selected_sequence_${i+1}\n${seq}`).join('\n');
      const enhancedCommand = `${command}\n\nSelected sequences from previous step:\n${sequencesText}`;
      console.log('Enhanced command with selected sequences:', enhancedCommand);
      return enhancedCommand;
    }
    
    // If command mentions "mutated sequences" and we have them in context
    if (lowerCommand.includes('mutated sequences') && workflowContext.mutatedSequences) {
      const sequencesText = workflowContext.mutatedSequences.map((seq, i) => `>mutant_${i+1}\n${seq}`).join('\n');
      const enhancedCommand = `${command}\n\nMutated sequences from previous step:\n${sequencesText}`;
      console.log('Enhanced command with mutated sequences:', enhancedCommand);
      return enhancedCommand;
    }
    
    console.log('No context enhancement applied');
    return command;
  };

  const executeCommand = async (
    commandText: string,
    scenarioIdOverride?: string,
    clearInputsAfter: boolean = true
  ) => {
    if (!commandText.trim()) return;

    // Immediately show user message + thinking card — don't wait for the backend
    const pendingId = `pending-${Date.now()}`;
    const pendingItem: HistoryItem = {
      input: commandText,
      output: null,
      type: 'pending',
      timestamp: new Date(),
      pendingId,
    };
    setHistory(prev => [pendingItem, ...prev]);
    // Scroll to the pending card so the loading dots are always visible,
    // regardless of where the user was in the conversation.
    setTimeout(() => pendingItemRef.current?.scrollIntoView({ behavior: 'smooth', block: 'start' }), 60);

    setLoading(true);
    const activityId = addActivity('Processing your request...');
    
    try {
      let response;
      let parsedCommand: ParsedCommand | undefined;
      
      console.log('Command mode:', commandMode);
      console.log('Session ID:', sessionId);
      console.log('Command:', commandText);
      console.log('Workflow context:', workflowContext);
      
      // Enhance command with workflow context
      let finalCommand = enhanceCommandWithContext(commandText);
      
      console.log('Original command:', commandText);
      console.log('Enhanced command:', finalCommand);
      console.log('Workflow context before sending:', workflowContext);
      
      // Uploaded files are now persisted server-side under the session directory.
      // Commands should reference session context rather than embedding file contents.
      
      // Use SSE streaming endpoint for improved perceived latency.
      // Progress events update the loading message; the final "result" event
      // resolves this promise with the same payload as the REST endpoint.
      console.log('Calling agent via /execute/stream endpoint…');
      response = await new Promise<unknown>((resolve, reject) => {
        const cleanup = helixApi.executeCommandStream(
          finalCommand,
          sessionId || undefined,
          (phase, message) => {
            // Update the cycling loading message from SSE progress events
            setLoadingMsg(message || phase);
          },
          (data) => {
            cleanup();
            resolve(data);
          },
          (err) => {
            cleanup();
            reject(err);
          },
        );
      });
      console.log('Agent response:', response);

      // ── Async pipeline job — switch to SSE mode ────────────────────────
      // When a pipeline tool (chip_seq_analysis, etc.) is dispatched, the
      // backend returns {job_id, status: 'submitted'} immediately instead of
      // a synchronous result.  Open an EventSource to receive live progress.
      if ((response as any).job_id && (response as any).status === 'submitted') {
        const jobId    = (response as any).job_id as string;
        const pipeline = (response as any).pipeline || (response as any).tool || 'pipeline';

        // Augment the raw response so the renderer can show a "submitted" card
        const submittedResponse = {
          ...(response as any),
          _pipelineJob: {
            jobId,
            pipeline,
            status: 'submitted' as string,
            steps:  [] as Array<{ type: string; process?: string }>,
            streamUrl: (response as any).stream_url || `/jobs/${jobId}/stream`,
          },
        };

        const historyItem: HistoryItem = {
          input:     commandText,
          output:    submittedResponse,
          type:      'agent',
          timestamp: new Date(),
          scenarioId: scenarioIdOverride ?? pendingScenarioId,
        };

        setPendingScenarioId(undefined);
        setHistory(prev => prev.map(item => item.pendingId === pendingId ? historyItem : item));
        setTimeout(() => historyTopRef.current?.scrollIntoView({ behavior: 'smooth' }), 100);
        updateActivity(activityId, 'completed');

        // Open the SSE stream and update the history item in-place as events arrive
        helixApi.subscribeJobStream(
          jobId,
          (evt) => {
            setHistory(prev => prev.map(item => {
              if (!item.output?._pipelineJob || item.output._pipelineJob.jobId !== jobId) return item;
              return {
                ...item,
                output: {
                  ...item.output,
                  _pipelineJob: {
                    ...item.output._pipelineJob,
                    status: 'running',
                    steps:  [...(item.output._pipelineJob.steps || []), evt],
                  },
                },
              };
            }));
          },
          (result) => {
            setHistory(prev => prev.map(item => {
              if (!item.output?._pipelineJob || item.output._pipelineJob.jobId !== jobId) return item;
              return {
                ...item,
                output: {
                  ...item.output,
                  status: 'completed',
                  text: `Pipeline **${pipeline}** completed. ${result.result_files?.length ?? 0} output file(s) ready.`,
                  _pipelineJob: {
                    ...item.output._pipelineJob,
                    status: 'completed',
                    resultFiles: result.result_files || [],
                  },
                },
              };
            }));
          },
          (err) => {
            setHistory(prev => prev.map(item => {
              if (!item.output?._pipelineJob || item.output._pipelineJob.jobId !== jobId) return item;
              return {
                ...item,
                output: {
                  ...item.output,
                  status: 'error',
                  text: `Pipeline **${pipeline}** failed: ${err.message}`,
                  _pipelineJob: {
                    ...item.output._pipelineJob,
                    status: 'failed',
                    error: err.message,
                  },
                },
              };
            }));
          },
        );

        setLoading(false);
        if (clearInputsAfter) setCommand('');
        return;   // skip the normal synchronous history push below
      }
      // ── end async pipeline job ─────────────────────────────────────────

      // Update session ID if the backend created one automatically
      if ((response as any).session_id && !sessionId) {
        const sid = (response as any).session_id;
        setSessionId(sid);
        try {
          sessionStorage.setItem(SESSION_STORAGE_KEY, sid);
        } catch {}
        console.log('Session created automatically by backend:', sid);
      }
      
      // Update workflow context based on the response
      if ((response as any).success && (response as any).result) {
        const result = (response as any).result.result || (response as any).result;
        
        console.log('Processing response result:', result);
        
        // Extract aligned sequences from alignment response
        // Check for the tool response structure
        if (result.messages && Array.isArray(result.messages)) {
          // Look for tool messages with sequence alignment output
          for (const message of result.messages) {
            if (message.type === 'tool' && message.name === 'sequence_alignment') {
              try {
                const toolResult = JSON.parse(message.content);
                if (toolResult.output && Array.isArray(toolResult.output)) {
                  const alignedSequences = toolResult.output.map((seq: any) => 
                    `>${seq.name}\n${seq.sequence}`
                  ).join('\n');
                  console.log('Updating workflow context with aligned sequences:', alignedSequences);
                  updateWorkflowContext('sequence_alignment', alignedSequences);
                  break;
                }
              } catch (e) {
                console.log('Could not parse tool result:', e);
              }
            }
            
            // Look for tool messages with sequence selection output
            if (message.type === 'tool' && message.name === 'sequence_selection') {
              try {
                const toolResult = JSON.parse(message.content);
                if (toolResult.output && Array.isArray(toolResult.output)) {
                  const selectedSequences = toolResult.output.map((seq: any) => seq.sequence);
                  console.log('Updating workflow context with selected sequences:', selectedSequences);
                  updateWorkflowContext('sequence_selection', selectedSequences);
                  break;
                }
              } catch (e) {
                console.log('Could not parse tool result:', e);
              }
            }
          }
        }
        
        // Fallback: try direct result.output structure
        if (result.output && Array.isArray(result.output) && result.output.length > 0 && result.output[0].name && result.output[0].sequence) {
          const alignedSequences = result.output.map((seq: any) => 
            `>${seq.name}\n${seq.sequence}`
          ).join('\n');
          console.log('Updating workflow context with aligned sequences (fallback):', alignedSequences);
          updateWorkflowContext('sequence_alignment', alignedSequences);
        }
        
        // Extract selected sequences from selection response (fallback)
        if (result.output && Array.isArray(result.output) && result.output.length > 0 && result.output[0].sequence) {
          const selectedSequences = result.output.map((seq: any) => seq.sequence);
          console.log('Updating workflow context with selected sequences (fallback):', selectedSequences);
          updateWorkflowContext('sequence_selection', selectedSequences);
        }
        
        // Extract mutated sequences from mutation response
        if (result.output && result.output.variants && Array.isArray(result.output.variants)) {
          console.log('Updating workflow context with mutated sequences:', result.output.variants);
          updateWorkflowContext('mutate_sequence', result.output.variants);
        }
        
        // Extract plasmid data from plasmid visualization response
        if (result.output && result.output.features) {
          console.log('Updating workflow context with plasmid data:', result.output);
          updateWorkflowContext('plasmid_visualization', result.output);
        }
      }
      
      // Add to history
      // Since /execute always routes through the agent, all responses are agent responses
      const _r = response as any;
      const historyItem: HistoryItem = {
        input: commandText,
        output: response,
        type: 'agent',
        timestamp: new Date(),
        scenarioId: scenarioIdOverride ?? pendingScenarioId,
        run_id: _r?.run_id || _r?.result?.run_id,
        parent_run_id: _r?.parent_run_id || _r?.result?.parent_run_id,
      };
      
      console.log('🔍 Adding to history:', historyItem);
      console.log('🔍 Response structure:', JSON.stringify(response, null, 2));
      console.log('🔍 Response success:', (response as any).success);
      console.log('🔍 Response result keys:', (response as any).result ? Object.keys((response as any).result) : 'No result');
      
      setPendingScenarioId(undefined);
      // Replace the pending placeholder with the resolved item
      setHistory(prev => prev.map(item => item.pendingId === pendingId ? historyItem : item));
      setTimeout(() => historyTopRef.current?.scrollIntoView({ behavior: 'smooth' }), 100);
      updateActivity(activityId, 'completed');
      
    } catch (error) {
      console.error('Error executing command:', error);
      updateActivity(activityId, 'error');
      
      // Add error to history
      // Since /execute always routes through the agent, errors are agent errors
      const historyItem: HistoryItem = {
        input: commandText,
        output: { error: error instanceof Error ? error.message : 'Unknown error' },
        type: 'agent_error',
        timestamp: new Date(),
        scenarioId: scenarioIdOverride ?? pendingScenarioId,
      };
      
      // Replace the pending placeholder with the error item
      setHistory(prev => prev.map(item => item.pendingId === pendingId ? historyItem : item));
      setTimeout(() => historyTopRef.current?.scrollIntoView({ behavior: 'smooth' }), 100);
    } finally {
      if (clearInputsAfter) {
        setCommand('');
      }
      setLoading(false);
    }
  };

  const handleSubmit = async () => {
    if (!command.trim()) return;
    return executeCommand(command, pendingScenarioId, true);
  };

  const handleAgentSubmit = async () => {
    if (!command.trim()) return;

    setAgentLoading(true);

    try {
      let finalCommand = enhanceCommandWithContext(command);

      // Uploaded files are persisted by session and should not be inlined into prompts.

      // ALWAYS use /execute endpoint - the agent handles everything
      // If files are needed, they should be included in the command text or handled by the agent
      const response = await helixApi.executeCommand(finalCommand, sessionId || undefined);

      if ((response as any).session_id && !sessionId) {
        const sid = (response as any).session_id;
        setSessionId(sid);
        try {
          sessionStorage.setItem(SESSION_STORAGE_KEY, sid);
        } catch {}
      }

      const historyItem: HistoryItem = {
        input: command,
        output: (response && response.result) ? response.result : response,
        type: 'agent',
        timestamp: new Date(),
      };

      setHistory(prev => [historyItem, ...prev]);
    } catch (error) {
      console.error('Error executing agent command:', error);
      const historyItem: HistoryItem = {
        input: command,
        output: { error: error instanceof Error ? error.message : 'Unknown error' },
        type: 'agent_error',
        timestamp: new Date(),
      };
      setHistory(prev => [historyItem, ...prev]);
    } finally {
      setCommand('');
      setAgentLoading(false);
    }
  };



  // Execute a previously-planned pipeline when the user clicks "Execute Pipeline"
  const handleExecutePipeline = async (originalCommand: string) => {
    setLoading(true);
    // Immediately mark the plan as submitted so its Execute button turns into "Submitted ✓"
    setExecutedPipelineCommands(prev => new Set([...prev, originalCommand]));
    const activityId = addActivity('Submitting pipeline jobs...');
    try {
      const response = await helixApi.executePipelinePlan(originalCommand, sessionId || undefined);
      if ((response as any).session_id && !sessionId) {
        const sid = (response as any).session_id;
        setSessionId(sid);
        try {
          sessionStorage.setItem(SESSION_STORAGE_KEY, sid);
        } catch {}
      }
      const historyItem: HistoryItem = {
        input: `▶ Execute Pipeline`,
        output: response,
        type: 'agent',
        timestamp: new Date(),
      };
      setHistory(prev => [historyItem, ...prev]);
      updateActivity(activityId, 'completed');
      // Scroll to the top of the history so the user sees the new execution result
      setTimeout(() => historyTopRef.current?.scrollIntoView({ behavior: 'smooth' }), 100);
    } catch (error) {
      console.error('Error executing pipeline:', error);
      // Unmark the command so the user can retry
      setExecutedPipelineCommands(prev => {
        const next = new Set(prev);
        next.delete(originalCommand);
        return next;
      });
      updateActivity(activityId, 'error');
      const historyItem: HistoryItem = {
        input: `▶ Execute Pipeline`,
        output: { error: error instanceof Error ? error.message : 'Failed to submit pipeline' },
        type: 'agent_error',
        timestamp: new Date(),
      };
      setHistory(prev => [historyItem, ...prev]);
    } finally {
      setLoading(false);
    }
  };

  // Drag-and-drop handlers
  const handleDragOver = (e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setDragActive(true);
  };
  const handleDragLeave = (e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setDragActive(false);
  };
  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setDragActive(false);
    if (e.dataTransfer.files && e.dataTransfer.files.length > 0) {
      const files = Array.from(e.dataTransfer.files);
      void handleUploadFiles(files);
    }
  };

  const handleBrowseClick = () => {
    fileInputRef.current?.click();
  };

  const handleFileRemove = (index?: number) => {
    if (index !== undefined) {
      setUploadedFiles(prev => prev.filter((_, i) => i !== index));
    } else {
      setUploadedFiles([]);
    }
  };

  const handleFileInput = (e: React.ChangeEvent<HTMLInputElement>) => {
    if (e.target.files && e.target.files.length > 0) {
      const files = Array.from(e.target.files);
      void handleUploadFiles(files);
      e.target.value = '';
    }
  };

  const extractAgentMessageText = (content: any): string => {
    if (!content) return '';
    if (typeof content === 'string') return content;
    if (Array.isArray(content)) {
      return content
        .map((part) => {
          if (!part) return '';
          if (typeof part === 'string') return part;
          if (typeof part === 'object') {
            if (typeof part.text === 'string') return part.text;
            if (typeof part.value === 'string') return part.value;
            if (part.content) return extractAgentMessageText(part.content);
          }
          return typeof part === 'object' ? JSON.stringify(part) : String(part);
        })
        .filter(Boolean)
        .join('\n');
    }
    if (typeof content === 'object') {
      if (typeof content.text === 'string') return content.text;
      if (typeof content.value === 'string') return content.value;
      if (content.content) return extractAgentMessageText(content.content);
    }
    return typeof content === 'object' ? JSON.stringify(content, null, 2) : String(content);
  };

  const formatLabel = (key: string) => {
    return key
      .replace(/_/g, ' ')
      .replace(/([a-z0-9])([A-Z])/g, '$1 $2')
      .replace(/^./, (char) => char.toUpperCase());
  };

  const formatNewickForDisplay = (newick: string): string => {
    const source = newick.trim();
    if (!source) return source;

    let depth = 0;
    let out = '';

    for (let i = 0; i < source.length; i++) {
      const char = source[i];

      if (char === '(') {
        out += '(\n' + '  '.repeat(depth + 1);
        depth += 1;
      } else if (char === ',') {
        out += ',\n' + '  '.repeat(depth);
      } else if (char === ')') {
        depth = Math.max(0, depth - 1);
        out += '\n' + '  '.repeat(depth) + ')';
      } else if (char === ';') {
        out += ';\n';
      } else {
        out += char;
      }
    }

    return out.trim();
  };

  const renderStructuredData = (data: any, depth = 0, visited = new WeakSet()): React.ReactNode => {
    if (data === null || data === undefined) return null;

    if (typeof data === 'string' || typeof data === 'number' || typeof data === 'boolean') {
      return <span>{String(data)}</span>;
    }

    if (Array.isArray(data)) {
      if (data.length === 0) return null;
      return (
        <ul className="structured-list">
          {data.map((item, index) => (
            <li key={index}>{renderStructuredData(item, depth + 1, visited)}</li>
          ))}
        </ul>
      );
    }

    if (typeof data === 'object') {
      if (visited.has(data)) {
        return <span>[circular]</span>;
      }
      visited.add(data);

      if (typeof data.text === 'string' && Object.keys(data).length === 1) {
        return <span>{data.text}</span>;
      }

      const entries = Object.entries(data).filter(([_, value]) => value !== undefined && value !== null);
      if (entries.length === 0) return null;

      return (
        <div className={`structured-object structured-object-depth-${depth}`}>
          {entries.map(([key, value]) => (
            <div key={key} className="structured-row">
              <div className="structured-label">{formatLabel(key)}</div>
              <div className="structured-value">
                {/newick/i.test(key) && typeof value === 'string' ? (
                  <pre className="newick-display">{formatNewickForDisplay(value)}</pre>
                ) : (
                  renderStructuredData(value, depth + 1, visited)
                )}
              </div>
            </div>
          ))}
        </div>
      );
    }

    return null;
  };


  const renderAgentResponse = (agentOutput: any) => {
    const agentResult = agentOutput?.result ?? agentOutput;
    
    console.log('🔍 renderAgentResponse called');
    console.log('🔍 agentOutput keys:', Object.keys(agentOutput || {}));
    console.log('🔍 agentResult keys:', agentResult && typeof agentResult === 'object' ? Object.keys(agentResult) : 'not an object');
    
    // Extract data from multiple possible locations (matching renderOutput logic)
    const rawResult = agentOutput?.raw_result || agentResult?.raw_result || (agentOutput?.result && agentOutput.result.raw_result);
    const actualResult = (agentOutput?.result && agentOutput.result.result) 
      ? agentOutput.result.result 
      : (agentOutput?.result || agentOutput?.raw_result || agentOutput || agentResult);
    
    // Debug: Log the structure to find tree_newick
    console.log('🔍 ========== RESPONSE STRUCTURE DEBUG ==========');
    console.log('🔍 agentOutput.raw_result keys:', rawResult && typeof rawResult === 'object' ? Object.keys(rawResult) : 'not an object');
    console.log('🔍 agentOutput.raw_result.result keys:', rawResult?.result && typeof rawResult.result === 'object' ? Object.keys(rawResult.result) : 'not an object');
    console.log('🔍 agentOutput.raw_result.tree_newick:', rawResult?.tree_newick ? `FOUND (${rawResult.tree_newick.length} chars)` : 'NOT FOUND');
    console.log('🔍 agentOutput.raw_result.result.tree_newick:', rawResult?.result?.tree_newick ? `FOUND (${rawResult.result.tree_newick.length} chars)` : 'NOT FOUND');
    console.log('🔍 agentOutput.tree_newick:', agentOutput?.tree_newick ? `FOUND (${agentOutput.tree_newick.length} chars)` : 'NOT FOUND');
    console.log('🔍 agentOutput.data:', agentOutput?.data ? 'EXISTS' : 'NOT FOUND');
    console.log('🔍 agentOutput.data keys:', agentOutput?.data && typeof agentOutput.data === 'object' ? Object.keys(agentOutput.data) : 'N/A');
    console.log('🔍 agentOutput.visualization_type:', agentOutput?.visualization_type);
    console.log('🔍 agentOutput.tool:', agentOutput?.tool);
    console.log('🔍 ==============================================');
    
    // Defined here (before any specialized renderers) so downloadLinksSection can use it.
    const normalizeAssetUrl = (url?: string) => {
      const u = (url || '').trim();
      if (!u) return u;
      if (/^https?:\/\//i.test(u)) return u;
      if (u.startsWith('/')) return `${API_BASE_URL}${u}`;
      return u;
    };

    // ── Download links (computed early so every renderer can include them) ──────
    // Prefer normalized backend field data.downloadable_artifacts.
    // Fallback to legacy links only when this is an executed/completed response.
    const statusCandidates = [
      agentOutput?.status,
      agentOutput?.result?.status,
      agentOutput?.raw_result?.status,
      agentOutput?.raw_result?.result?.status,
      actualResult?.status,
    ]
      .map((s) => (s == null ? '' : String(s).toLowerCase()))
      .filter(Boolean);
    const planningStatuses = new Set(['workflow_planned', 'needs_inputs', 'tool_mapped', 'workflow_needs_clarification']);
    // Prefer planning statuses when mixed envelopes disagree (e.g., top-level success + inner plan).
    const responseStatus =
      statusCandidates.find((s) => planningStatuses.has(s)) ||
      statusCandidates[0] ||
      '';
    const nonExecutedStatuses = new Set([
      'workflow_planned',
      'needs_inputs',
      'tool_mapped',
      'pipeline_submitted',
      'workflow_needs_clarification',
      'submitted',
      'running',
      'queued',
      'pending',
    ]);
    const isExecutedForDownloads =
      !!responseStatus &&
      !nonExecutedStatuses.has(responseStatus) &&
      !['error', 'failed', 'workflow_failed'].includes(responseStatus);
    const executeReadyFlag = (
      agentOutput?.execute_ready ??
      agentOutput?.result?.execute_ready ??
      agentOutput?.raw_result?.execute_ready ??
      agentOutput?.raw_result?.result?.execute_ready
    );
    const topText = String(agentOutput?.text || agentOutput?.result?.text || '');
    const looksLikePlanText = /##\s*Pipeline Plan/i.test(topText);
    const isNonExecutablePlanResponse =
      responseStatus === 'workflow_planned' ||
      responseStatus === 'needs_inputs' ||
      (executeReadyFlag === false && looksLikePlanText);

    const _downloadArtifactsRaw = agentOutput?.data?.downloadable_artifacts;
    const _hasDownloadArtifacts =
      Array.isArray(_downloadArtifactsRaw) && _downloadArtifactsRaw.length > 0;
    const _legacyLinksRaw: any[] =
      agentOutput?.data?.links || actualResult?.links || rawResult?.links || [];
    const _earlyLinksRaw: any[] = isNonExecutablePlanResponse
      ? []
      : (_hasDownloadArtifacts
        ? (_downloadArtifactsRaw as any[])
        : (isExecutedForDownloads ? _legacyLinksRaw : []));
    const _earlyLinks: any[] = Array.isArray(_earlyLinksRaw)
      ? _earlyLinksRaw
          .map((l: any) => (l && typeof l === 'object' ? { ...l, url: normalizeAssetUrl(l.url) } : l))
          .filter((l: any) => l && typeof l === 'object' && typeof l.url === 'string' && l.url.trim().length > 0)
      : [];

    const downloadLinksSection =
      _earlyLinks.length > 0 ? (
        <div className="mt-3 pt-3" style={{ borderTop: '1px solid #E2E8F0' }}>
          <div style={{ fontSize: '0.8rem', fontWeight: 600, color: '#475569', marginBottom: '6px' }}>
            Download results
          </div>
          <div style={{ display: 'flex', flexWrap: 'wrap', gap: '8px' }}>
            {_earlyLinks.map((l: any, idx: number) => {
              const label: string = l?.label || l?.s3_uri || 'artifact';
              const isScript = label.endsWith('.py');
              const isBundle = label.endsWith('.zip');
              const isDownloadable = (
                isScript ||
                isBundle ||
                /\.(newick|nwk|csv|tsv|txt|json|fasta|fa)$/i.test(label)
              );
              const icon = isBundle ? '📦' : isScript ? '📄' : '⬇';
              const styleMap: React.CSSProperties = isBundle
                ? { background: '#0f4c81', color: '#e0f0ff', border: '1px solid #1e6bb8', fontWeight: 600 }
                : isScript
                ? { background: '#1e293b', color: '#f1f5f9', border: '1px solid #334155', fontFamily: 'monospace' }
                : { background: 'transparent', color: '#2563eb', border: '1px solid #2563eb' };
              return (
                <a
                  key={`dl-${label}-${idx}`}
                  href={l.url}
                  download={isDownloadable ? label : undefined}
                  target={isDownloadable ? undefined : '_blank'}
                  rel="noreferrer"
                  className="btn btn-sm"
                  style={{ fontSize: '0.78rem', ...styleMap }}
                >
                  {icon} {label}
                </a>
              );
            })}
          </div>
        </div>
      ) : null;

    // ========== PHYLOGENETIC TREE VISUALIZATIONS ==========
    const extractNewickFromText = (text?: string): string | undefined => {
      if (!text || typeof text !== 'string') return undefined;

      // Prefer fenced code blocks if present
      const fenced = text.match(/```(?:newick)?\s*([\s\S]*?)```/i);
      if (fenced?.[1]) {
        const candidate = fenced[1].trim();
        if (candidate.includes('(') && candidate.includes(')') && candidate.includes(';')) {
          return candidate;
        }
      }

      // Fallback: grab the first Newick-like expression ending in semicolon
      const inline = text.match(/(\([^;]*\)[^;]*;)/s);
      if (inline?.[1]) {
        const candidate = inline[1].trim();
        if (candidate.includes('(') && candidate.includes(')') && candidate.includes(';')) {
          return candidate;
        }
      }

      return undefined;
    };

    // Check ALL possible locations including top-level, raw_result, data, and nested structures
    // Priority: raw_result.result first (most likely location based on agent execution structure), then raw_result, then top-level
    const treeNewick = rawResult?.result?.tree_newick ||  // Agent execution result structure
                       rawResult?.result?.newick_tree ||
                       rawResult?.tree_newick ||
                       rawResult?.newick_tree ||
                       agentOutput?.tree_newick ||
                       agentOutput?.newick_tree ||
                       agentResult?.tree_newick ||
                       agentResult?.newick_tree ||
                       actualResult?.tree_newick ||
                       actualResult?.newick_tree ||
                       agentOutput?.data?.tree_newick ||
                       agentOutput?.data?.newick_tree ||
                       (agentOutput?.result && agentOutput.result.tree_newick) ||
                       (agentOutput?.result && agentOutput.result.newick_tree) ||
                       (agentOutput?.raw_result && agentOutput.raw_result.tree_newick) ||
                       (agentOutput?.raw_result && agentOutput.raw_result.newick_tree) ||
                       (agentOutput?.raw_result?.result && agentOutput.raw_result.result.tree_newick) ||
                       (agentOutput?.raw_result?.result && agentOutput.raw_result.result.newick_tree) ||
                       extractNewickFromText(actualResult?.text) ||
                       extractNewickFromText(rawResult?.text) ||
                       extractNewickFromText(agentOutput?.text);
    
    const eteVisualization = rawResult?.result?.ete_visualization ||  // Agent execution result structure
                             agentOutput?.ete_visualization ||
                             agentResult?.ete_visualization ||
                             actualResult?.ete_visualization ||
                             rawResult?.ete_visualization ||
                             agentOutput?.data?.ete_visualization ||
                             (agentOutput?.result && agentOutput.result.ete_visualization) ||
                             (agentOutput?.raw_result?.result && agentOutput.raw_result.result.ete_visualization);
    
    const clusteringResult = rawResult?.result?.clustering_result ||  // Agent execution result structure
                             actualResult?.clustering_result || 
                             rawResult?.clustering_result ||
                             agentOutput?.data?.clustering_result ||
                             (agentOutput?.result && agentOutput.result.clustering_result) ||
                             (agentOutput?.raw_result?.result && agentOutput.raw_result.result.clustering_result);
    
    const clusteredVisualization = rawResult?.result?.clustered_visualization ||  // Agent execution result structure
                                   actualResult?.clustered_visualization || 
                                   rawResult?.clustered_visualization ||
                                   agentOutput?.data?.clustered_visualization ||
                                   (agentOutput?.result && agentOutput.result.clustered_visualization) ||
                                   (agentOutput?.raw_result?.result && agentOutput.raw_result.result.clustered_visualization);
    
    console.log('🔍 Final treeNewick check:', treeNewick ? `FOUND (${treeNewick.length} chars)` : 'NOT FOUND');
    
    if (treeNewick) {
      console.log('🔍 Rendering phylogenetic tree visualization');
      const parseStatsFromText = (text?: string): Record<string, string> => {
        if (!text || typeof text !== 'string') return {};
        const stats: Record<string, string> = {};
        const lines = text.split('\n').map((l) => l.trim()).filter(Boolean);

        for (const line of lines) {
          // Matches markdown bold form, e.g. "**Sequences analysed:** 3"
          const boldMatch = line.match(/^\*\*([^*]+)\*\*:\s*(.+)$/);
          if (boldMatch) {
            stats[boldMatch[1].trim()] = boldMatch[2].trim();
            continue;
          }
          // Matches plain key:value form
          const plainMatch = line.match(/^([A-Za-z][A-Za-z0-9 _%()-]+):\s*(.+)$/);
          if (plainMatch) {
            const key = plainMatch[1].trim();
            if (!/^(step|tool|status)$/i.test(key)) {
              stats[key] = plainMatch[2].trim();
            }
          }
        }

        return stats;
      };

      const phyloStatsRaw =
        actualResult?.statistics ||
        rawResult?.result?.statistics ||
        rawResult?.statistics ||
        agentOutput?.data?.results?.result?.statistics ||
        agentOutput?.data?.results?.statistics ||
        agentOutput?.result?.statistics ||
        agentResult?.statistics;

      const parsedTextStats = {
        ...parseStatsFromText(actualResult?.text),
        ...parseStatsFromText(rawResult?.text),
        ...parseStatsFromText(agentOutput?.text),
      };

      const normalizedPhyloStats: Record<string, string> = {};
      if (phyloStatsRaw && typeof phyloStatsRaw === 'object') {
        Object.entries(phyloStatsRaw).forEach(([key, value]) => {
          normalizedPhyloStats[key] = String(value);
        });
      }
      Object.entries(parsedTextStats).forEach(([key, value]) => {
        if (!(key in normalizedPhyloStats)) {
          normalizedPhyloStats[key] = value;
        }
      });

      const rawVisuals: any[] = [
        ...(Array.isArray(rawResult?.result?.visuals) ? rawResult.result.visuals : []),
        ...(Array.isArray(rawResult?.visuals) ? rawResult.visuals : []),
        ...(Array.isArray(actualResult?.visuals) ? actualResult.visuals : []),
        ...(Array.isArray(agentOutput?.data?.visuals) ? agentOutput.data.visuals : []),
      ];

      const staticImageVisuals: Array<{ title: string; src: string }> = [];
      const seenSources = new Set<string>();
      const addImageVisual = (title: string, src?: string) => {
        if (!src || typeof src !== 'string') return;
        const trimmed = src.trim();
        if (!trimmed) return;
        if (seenSources.has(trimmed)) return;
        seenSources.add(trimmed);
        staticImageVisuals.push({ title, src: trimmed });
      };

      // Visuals payload from backend (preferred)
      rawVisuals.forEach((v: any, idx: number) => {
        if (!v || typeof v !== 'object') return;
        const title = v.title || `Plot ${idx + 1}`;
        if (v.type === 'image_b64' && typeof v.data === 'string') {
          addImageVisual(title, `data:image/png;base64,${v.data}`);
          return;
        }
        if (v.type === 'image' && typeof v.url === 'string') {
          addImageVisual(title, v.url);
          return;
        }
        if (typeof v.svg === 'string' && v.svg.trim().startsWith('<svg')) {
          addImageVisual(title, `data:image/svg+xml;utf8,${encodeURIComponent(v.svg)}`);
        }
      });

      // Direct phylo image fields for payloads that omit `visuals`
      const maybeB64 = (value: any) =>
        typeof value === 'string' && value.trim() && !value.trim().startsWith('<')
          ? `data:image/png;base64,${value.trim()}`
          : undefined;
      addImageVisual(
        'Phylogenetic Tree',
        maybeB64(rawResult?.result?.tree_plot_b64 || actualResult?.tree_plot_b64 || rawResult?.tree_plot_b64),
      );
      addImageVisual(
        'Pairwise Distance Matrix',
        maybeB64(rawResult?.result?.dist_matrix_b64 || actualResult?.dist_matrix_b64 || rawResult?.dist_matrix_b64),
      );
      if (typeof clusteredVisualization === 'string') {
        if (clusteredVisualization.trim().startsWith('<svg')) {
          addImageVisual('Clustered Tree View', `data:image/svg+xml;utf8,${encodeURIComponent(clusteredVisualization)}`);
        } else {
          addImageVisual('Clustered Tree View', `data:image/png;base64,${clusteredVisualization}`);
        }
      } else if (clusteredVisualization?.svg) {
        addImageVisual('Clustered Tree View', `data:image/svg+xml;utf8,${encodeURIComponent(clusteredVisualization.svg)}`);
      } else if (clusteredVisualization?.data) {
        addImageVisual('Clustered Tree View', `data:image/png;base64,${clusteredVisualization.data}`);
      }

      if (typeof eteVisualization === 'string') {
        if (eteVisualization.trim().startsWith('<svg')) {
          addImageVisual('Annotated Tree (ETE3)', `data:image/svg+xml;utf8,${encodeURIComponent(eteVisualization)}`);
        } else {
          addImageVisual('Annotated Tree (ETE3)', `data:image/png;base64,${eteVisualization}`);
        }
      } else if (eteVisualization?.svg) {
        addImageVisual('Annotated Tree (ETE3)', `data:image/svg+xml;utf8,${encodeURIComponent(eteVisualization.svg)}`);
      } else if (eteVisualization?.data) {
        addImageVisual('Annotated Tree (ETE3)', `data:image/png;base64,${eteVisualization.data}`);
      }

      return (
        <div>
          {actualResult?.text && (
            <div className="agent-response-markdown mb-3">
              <ReactMarkdown remarkPlugins={[remarkGfm]}>{actualResult.text}</ReactMarkdown>
            </div>
          )}

          <div className="row g-3 mb-3">
            <div className="col-12 col-xl-8">
              <div className="bg-light p-3 border rounded h-100">
                <h5 className="mb-3">🧬 Interactive Phylogenetic Tree</h5>
                <PhylogeneticTree newick={treeNewick} />
              </div>
            </div>
            <div className="col-12 col-xl-4">
              <div className="d-flex flex-column gap-3 h-100">
                <div className="bg-light p-3 border rounded">
                  <h5 className="mb-3">📈 Tree Statistics</h5>
                  {Object.keys(normalizedPhyloStats).length > 0 ? (
                    <div className="row">
                      {Object.entries(normalizedPhyloStats).map(([key, value]) => (
                        <div key={key} className="col-12 mb-2">
                          <strong>{key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}:</strong> {String(value)}
                        </div>
                      ))}
                    </div>
                  ) : (
                    <div className="text-muted">Summary metrics unavailable in this response.</div>
                  )}
                </div>
                <div className="bg-light p-3 border rounded flex-grow-1">
                  <h5 className="mb-2">🌿 Newick Tree</h5>
                  <pre className="newick-display mb-0">{formatNewickForDisplay(String(treeNewick))}</pre>
                </div>
              </div>
            </div>
          </div>

          {staticImageVisuals.length > 0 && (
            <div className="bg-light p-3 border rounded mb-3">
              <h5 className="mb-3">🖼️ Analysis Visuals</h5>
              <div style={{ display: 'flex', flexWrap: 'wrap', gap: '16px' }}>
                {staticImageVisuals.map((v, idx) => (
                  <div key={`${v.title}-${idx}`} style={{ flex: '1 1 420px', minWidth: '280px', maxWidth: '700px' }}>
                    <div style={{ fontSize: '0.9rem', fontWeight: 600, color: '#334155', marginBottom: '6px' }}>
                      {v.title}
                    </div>
                    <a href={v.src} target="_blank" rel="noreferrer" title="Open full-size">
                      <img
                        src={v.src}
                        alt={v.title}
                        style={{ width: '100%', borderRadius: '8px', border: '1px solid #E2E8F0', background: '#fff', cursor: 'pointer' }}
                      />
                    </a>
                  </div>
                ))}
              </div>
            </div>
          )}

          {/* Clustering Results */}
          {clusteringResult && (
            <div className="bg-light p-3 border rounded mb-3">
              <h5>🧬 Clustering Analysis Results</h5>
              <div className="row">
                <div className="col-md-6">
                  <h6>Summary</h6>
                  <ul className="list-unstyled">
                    <li><strong>Total Sequences:</strong> {clusteringResult.total_sequences}</li>
                    <li><strong>Number of Clusters:</strong> {clusteringResult.num_clusters}</li>
                    <li><strong>Representatives:</strong> {clusteringResult.representatives?.join(', ')}</li>
                  </ul>
                </div>
                <div className="col-md-6">
                  <h6>Cluster Details</h6>
                  {clusteringResult.clusters && (
                    <div>
                      {clusteringResult.clusters.map((cluster: any, index: number) => (
                        <div key={index} className="mb-2">
                          <strong>Cluster {cluster.cluster_id}:</strong>
                          <ul className="list-unstyled ml-3">
                            <li>Size: {cluster.size} sequences</li>
                            <li>Representative: {cluster.representative}</li>
                            <li>Average Distance: {cluster.average_distance?.toFixed(4)}</li>
                          </ul>
                        </div>
                      ))}
                    </div>
                  )}
                </div>
              </div>
            </div>
          )}

          {downloadLinksSection}

          {/* Debug Section - Consistent across all commands */}
          {renderDebugInfo(agentOutput, actualResult)}
        </div>
      );
    }
    
    // ========== PLASMID VISUALIZATIONS ==========
    const plasmidData = rawResult?.result?.plasmid_data ||  // Agent execution result structure
                        actualResult?.plasmid_data || 
                        rawResult?.plasmid_data ||
                        (agentOutput?.result && agentOutput.result.plasmid_data) ||
                        (agentOutput?.raw_result?.result && agentOutput.raw_result.result.plasmid_data);
    
    const plasmidResults = rawResult?.result?.plasmid_results ||  // Agent execution result structure
                           actualResult?.plasmid_results ||
                           rawResult?.plasmid_results ||
                           (agentOutput?.result && agentOutput.result.plasmid_results) ||
                           (agentOutput?.raw_result?.result && agentOutput.raw_result.result.plasmid_results);
    
    if (plasmidData) {
      console.log('🔍 Rendering plasmid visualization');
      console.log('🔍 plasmidData:', plasmidData);
      console.log('🔍 plasmidData.sequence:', plasmidData?.sequence);
      console.log('🔍 plasmidData type:', typeof plasmidData);
      
      // Validate plasmid data before rendering
      if (!plasmidData.sequence || typeof plasmidData.sequence !== 'string') {
        console.error('❌ Invalid plasmid data: missing or invalid sequence');
        return (
          <div>
            {actualResult?.text && (
              <pre className="bg-light p-3 border rounded mb-3">{actualResult.text}</pre>
            )}
            <div className="alert alert-warning">
              <h5>🧬 Plasmid Visualization</h5>
              <p>Plasmid data received but sequence is missing or invalid.</p>
              <details>
                <summary>Debug Info</summary>
                <pre>{JSON.stringify(plasmidData, null, 2)}</pre>
              </details>
            </div>
            
            {/* Debug Section - Consistent across all commands */}
            {renderDebugInfo(agentOutput, actualResult)}
          </div>
        );
      }
      
      return (
        <div>
          {actualResult?.text && (
            <pre className="bg-light p-3 border rounded mb-3">{actualResult.text}</pre>
          )}
          <div className="bg-light p-3 border rounded mb-3">
            <h5>🧬 Plasmid Visualization</h5>
            <PlasmidDataVisualizer data={plasmidData} />
          </div>
          
          {/* Debug Section - Consistent across all commands */}
          {renderDebugInfo(agentOutput, actualResult)}
        </div>
      );
    }
    
    if (plasmidResults) {
      console.log('🔍 Rendering plasmid representatives visualization');
      console.log('🔍 plasmidResults:', plasmidResults);
      
      // Validate plasmid results
      if (!Array.isArray(plasmidResults) || plasmidResults.length === 0) {
        console.error('❌ Invalid plasmid results: not an array or empty');
        return (
          <div>
            {actualResult?.text && (
              <pre className="bg-light p-3 border rounded mb-3">{actualResult.text}</pre>
            )}
            <div className="alert alert-warning">
              <h5>🧬 Plasmid Representatives Visualization</h5>
              <p>Plasmid results received but data is invalid.</p>
            </div>
            
            {/* Debug Section - Consistent across all commands */}
            {renderDebugInfo(agentOutput, actualResult)}
          </div>
        );
      }
      
      return (
        <div>
          {actualResult?.text && (
            <pre className="bg-light p-3 border rounded mb-3">{actualResult.text}</pre>
          )}
          <PlasmidRepresentativesVisualizer 
            plasmidResults={plasmidResults}
            vectorName={rawResult?.result?.vector_name || actualResult?.vector_name || "pUC19"}
            cloningSites={rawResult?.result?.cloning_sites || actualResult?.cloning_sites || "EcoRI, BamHI, HindIII"}
          />
          
          {/* Debug Section - Consistent across all commands */}
          {renderDebugInfo(agentOutput, actualResult)}
        </div>
      );
    }
    
    // ========== QUALITY ASSESSMENT VISUALIZATIONS ==========
    const qualityMetrics = rawResult?.result?.quality_metrics ||  // Agent execution result structure
                          actualResult?.quality_metrics ||
                          rawResult?.quality_metrics ||
                          (agentOutput?.result && agentOutput.result.quality_metrics) ||
                          (agentOutput?.raw_result?.result && agentOutput.raw_result.result.quality_metrics);
    
    const qualityPlotData = rawResult?.result?.plot_data ||  // Agent execution result structure
                           actualResult?.plot_data || 
                           rawResult?.plot_data ||
                           (agentOutput?.result && agentOutput.result.plot_data) ||
                           (agentOutput?.result && agentOutput.result.result && agentOutput.result.result.plot_data) ||
                           (agentOutput?.raw_result?.result && agentOutput.raw_result.result.plot_data);
    
    const qualitySummary = rawResult?.result?.summary ||  // Agent execution result structure
                          actualResult?.summary ||
                          rawResult?.summary ||
                          (agentOutput?.result && agentOutput.result.summary) ||
                          (agentOutput?.raw_result?.result && agentOutput.raw_result.result.summary);
    
    // Skip legacy quality-assessment renderer when the response already uses the
    // modern `results_viewer` path (which handles FastQC via BioOrchestrator).
    const isResultsViewer = agentOutput?.visualization_type === 'results_viewer';
    if (!isResultsViewer && (qualityMetrics || qualityPlotData || qualitySummary)) {
      console.log('🔍 Rendering quality assessment visualization');
      return (
        <div>
          {actualResult?.text && (
            <pre className="bg-light p-3 border rounded mb-3">{actualResult.text}</pre>
          )}
          
          {qualitySummary && (
            <div className="bg-light p-3 border rounded mb-3">
              <h5>Quality Assessment Summary</h5>
              <div className="row">
                {Object.entries(qualitySummary).map(([key, value]) => (
                  <div key={key} className="col-md-6 mb-2">
                    <strong>{key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}:</strong>{' '}
                    {typeof value === 'object' && value !== null
                      ? <code style={{fontSize: '0.78rem'}}>{JSON.stringify(value)}</code>
                      : String(value)}
                  </div>
                ))}
              </div>
            </div>
          )}
          
          {qualityPlotData && (
            <div className="bg-light p-3 border rounded mb-3">
              <h5>Quality Assessment Plots</h5>
              
              {qualityPlotData.length_distribution && (
                <div className="mb-4">
                  <h6>Length Distribution</h6>
                  <Plot
                    data={[{
                      x: qualityPlotData.length_distribution.x,
                      y: qualityPlotData.length_distribution.y,
                      type: 'bar',
                      marker: { color: '#007bff' }
                    }]}
                    layout={{
                      title: 'Length Distribution',
                      xaxis: { title: 'Length (bp)' },
                      yaxis: { title: 'Count' },
                      height: 400,
                      margin: { l: 60, r: 20, t: 60, b: 60 }
                    }}
                    style={{ width: '100%' }}
                  />
                </div>
              )}
              
              {qualityPlotData.gc_content_distribution && (
                <div className="mb-4">
                  <h6>GC Content Distribution</h6>
                  <Plot
                    data={[{
                      x: qualityPlotData.gc_content_distribution.x,
                      y: qualityPlotData.gc_content_distribution.y,
                      type: 'bar',
                      marker: { color: '#28a745' }
                    }]}
                    layout={{
                      title: 'GC Content Distribution',
                      xaxis: { title: 'GC Content (%)' },
                      yaxis: { title: 'Count' },
                      height: 400,
                      margin: { l: 60, r: 20, t: 60, b: 60 }
                    }}
                    style={{ width: '100%' }}
                  />
                </div>
              )}
              
              {qualityPlotData.base_composition && (
                <div className="mb-4">
                  <h6>Base Composition</h6>
                  <Plot
                    data={[{
                      x: qualityPlotData.base_composition.x,
                      y: qualityPlotData.base_composition.y,
                      type: 'bar',
                      marker: { color: '#ffc107' }
                    }]}
                    layout={{
                      title: 'Base Composition',
                      xaxis: { title: 'Base' },
                      yaxis: { title: 'Count' },
                      height: 400,
                      margin: { l: 60, r: 20, t: 60, b: 60 }
                    }}
                    style={{ width: '100%' }}
                  />
                </div>
              )}
            </div>
          )}
          
          {downloadLinksSection}

          {/* Debug Section - Consistent across all commands */}
          {renderDebugInfo(agentOutput, actualResult)}
        </div>
      );
    }
    
    // ========== SEQUENCE ALIGNMENT VISUALIZATIONS ==========
    // Check for aligned sequences in various formats
    const alignedSequences = rawResult?.result?.aligned_sequences ||  // Agent execution result structure
                            rawResult?.result?.output ||
                            actualResult?.output ||
                            actualResult?.aligned_sequences ||
                            rawResult?.aligned_sequences ||
                            (agentOutput?.result && agentOutput.result.aligned_sequences) ||
                            (agentOutput?.raw_result?.result && agentOutput.raw_result.result.aligned_sequences);
    
    if (alignedSequences && Array.isArray(alignedSequences)) {
      console.log('🔍 Rendering sequence alignment visualization');
      return (
        <div>
          {actualResult?.text && (
            <pre className="bg-light p-3 border rounded mb-3">{actualResult.text}</pre>
          )}
          <div className="bg-light p-3 border rounded mb-3">
            <h5>Sequence Alignment Result</h5>
            <div className="table-responsive">
              <table className="table table-sm">
                <thead>
                  <tr>
                    <th>Name</th>
                    <th>Sequence</th>
                  </tr>
                </thead>
                <tbody>
                  {alignedSequences.map((seq: any, index: number) => (
                    <tr key={index}>
                      <td>{seq.name || `Sequence_${index + 1}`}</td>
                      <td className="font-monospace">{seq.sequence || seq}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
          {actualResult?.statistics && (
            <div className="bg-light p-3 border rounded mb-3">
              <h6>Alignment Statistics</h6>
              <div className="row">
                {Object.entries(actualResult.statistics).map(([key, value]) => (
                  <div key={key} className="col-md-6 mb-2">
                    <strong>{key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}:</strong> {String(value)}
                  </div>
                ))}
              </div>
            </div>
          )}
          
          {/* Debug Section - Consistent across all commands */}
          {renderDebugInfo(agentOutput, actualResult)}
        </div>
      );
    }
    
    // ========== PLOTLY PLOTS ==========
    const plotData = rawResult?.result?.plot ||  // Agent execution result structure
                    actualResult?.plot ||
                    rawResult?.plot ||
                    (agentOutput?.result && agentOutput.result.plot) ||
                    (agentOutput?.raw_result?.result && agentOutput.raw_result.result.plot);
    
    if (plotData && plotData.data) {
      console.log('🔍 Rendering Plotly visualization');
      return (
        <div>
          {actualResult?.text && (
            <pre className="bg-light p-3 border rounded mb-3">{actualResult.text}</pre>
          )}
          <div className="bg-light p-3 border rounded mb-3">
            <Plot
              data={plotData.data}
              layout={plotData.layout || { title: 'Visualization' }}
              style={{ width: '100%' }}
            />
          </div>
          
          {/* Debug Section - Consistent across all commands */}
          {renderDebugInfo(agentOutput, actualResult)}
        </div>
      );
    }
    
    // ========== VENDOR RESEARCH ==========
    const vendors = rawResult?.result?.vendors ||  // Agent execution result structure
                   actualResult?.vendors ||
                   rawResult?.vendors ||
                   (agentOutput?.result && agentOutput.result.vendors) ||
                   (agentOutput?.raw_result?.result && agentOutput.raw_result.result.vendors);
    
    const recommendations = rawResult?.result?.recommendations ||  // Agent execution result structure
                          actualResult?.recommendations ||
                          rawResult?.recommendations ||
                          (agentOutput?.result && agentOutput.result.recommendations) ||
                          (agentOutput?.raw_result?.result && agentOutput.raw_result.result.recommendations) ||
                          [];
    
    if (vendors) {
      console.log('🔍 Rendering vendor research results');
      return (
        <div>
          {actualResult?.text && (
            <pre className="bg-light p-3 border rounded mb-3">{actualResult.text}</pre>
          )}
          <div className="bg-light p-3 border rounded mb-3">
            <h5>DNA Synthesis Vendors</h5>
            <div className="row">
              {Object.entries(vendors).map(([vendorId, vendor]: [string, any]) => (
                <div key={vendorId} className="col-md-6 mb-3">
                  <div className="card h-100">
                    <div className="card-body">
                      <h6 className="card-title">{vendor.name}</h6>
                      <p className="card-text">
                        <strong>Services:</strong> {vendor.services.join(', ')}
                      </p>
                      <p className="card-text">
                        <strong>Pricing:</strong> {vendor.pricing_range}
                      </p>
                      <p className="card-text">
                        <strong>Turnaround:</strong> {vendor.turnaround_time}
                      </p>
                      <p className="card-text">
                        <strong>Max Length:</strong> {vendor.max_length}
                      </p>
                      <p className="card-text">
                        <strong>Specialties:</strong> {vendor.specialties.join(', ')}
                      </p>
                      <a href={vendor.website} target="_blank" rel="noopener noreferrer" className="btn btn-sm btn-outline-primary">
                        Visit Website
                      </a>
                    </div>
                  </div>
                </div>
              ))}
            </div>
          </div>
          {recommendations.length > 0 && (
            <div className="bg-light p-3 border rounded mb-3">
              <h5>Recommendations</h5>
              <ul className="list-unstyled">
                {recommendations.map((rec: string, index: number) => (
                  <li key={index} className="mb-2">
                    <i className="bi bi-lightbulb text-warning me-2"></i>
                    {rec}
                  </li>
                ))}
              </ul>
            </div>
          )}
          
          {/* Debug Section - Consistent across all commands */}
          {renderDebugInfo(agentOutput, actualResult)}
        </div>
      );
    }
    
    // Extract messages from multiple possible locations
    const messages: any[] = Array.isArray(agentResult?.messages)
      ? agentResult.messages
      : Array.isArray(agentResult?.output?.messages)
        ? agentResult.output.messages
        : Array.isArray(agentOutput?.raw_result?.messages)
          ? agentOutput.raw_result.messages
          : Array.isArray(actualResult?.messages)
            ? actualResult.messages
            : Array.isArray(rawResult?.messages)
              ? rawResult.messages
              : [];

    // Debug logging for questions
    if (agentOutput?.tool === 'agent' || agentOutput?.raw_result?.task_type === 'qa') {
      console.log('🔍 [Question Detection] Detected question response');
      console.log('🔍 [Question Detection] Messages count:', messages.length);
      console.log('🔍 [Question Detection] Raw result keys:', agentOutput?.raw_result ? Object.keys(agentOutput.raw_result) : 'none');
    }

    // Filter for assistant messages only (skip system messages)
    const assistantMessages = messages.filter((msg) => {
      const role = msg?.role?.toLowerCase();
      const type = msg?.type?.toLowerCase();
      return (role === 'assistant' || type === 'ai') && role !== 'system' && type !== 'system';
    });
    
    // Get the last assistant message (most recent answer)
    const finalAssistant = assistantMessages.length > 0 ? assistantMessages[assistantMessages.length - 1] : null;
    let finalText = finalAssistant
      ? extractAgentMessageText((finalAssistant as any).content ?? (finalAssistant as any).text ?? (finalAssistant as any).value)
      : '';
    
    // CRITICAL FIX: If the text is JSON wrapped in code blocks, parse it and extract the answer
    // The agent returns JSON with details_markdown or user_friendly_summary containing the actual markdown answer
    if (finalText && (finalText.includes('```json') || finalText.includes('```')) && (finalText.includes('user_friendly') || finalText.includes('details_markdown') || finalText.includes('task_type'))) {
      try {
        // Extract JSON from code block (handle both ```json and ```)
        const jsonMatch = finalText.match(/```(?:json)?\s*([\s\S]*?)\s*```/);
        if (jsonMatch && jsonMatch[1]) {
          const jsonString = jsonMatch[1].trim();
          console.log('🔍 [Question Detection] Extracted JSON string length:', jsonString.length);
          console.log('🔍 [Question Detection] JSON string preview:', jsonString.substring(0, 200));
          
          const jsonContent = JSON.parse(jsonString);
          console.log('🔍 [Question Detection] Parsed JSON keys:', Object.keys(jsonContent));
          
          // Extract the answer field - check for details_markdown first (full answer), then user_friendly_summary
          // The backend returns details_markdown (full markdown answer) and user_friendly_summary (brief summary)
          if (jsonContent.details_markdown && typeof jsonContent.details_markdown === 'string') {
            finalText = jsonContent.details_markdown;
            console.log('✅ [Question Detection] Extracted details_markdown from JSON, length:', finalText.length);
          } else if (jsonContent.user_friendly_summary && typeof jsonContent.user_friendly_summary === 'string') {
            finalText = jsonContent.user_friendly_summary;
            console.log('✅ [Question Detection] Extracted user_friendly_summary from JSON, length:', finalText.length);
          } else if (jsonContent.user_friendly && typeof jsonContent.user_friendly === 'string') {
            finalText = jsonContent.user_friendly;
            console.log('✅ [Question Detection] Extracted user_friendly from JSON, length:', finalText.length);
          } else if (jsonContent.answer && typeof jsonContent.answer === 'string') {
            finalText = jsonContent.answer;
            console.log('✅ [Question Detection] Extracted answer from JSON, length:', finalText.length);
          } else if (jsonContent.text && typeof jsonContent.text === 'string') {
            finalText = jsonContent.text;
            console.log('✅ [Question Detection] Extracted text from JSON, length:', finalText.length);
          } else if (jsonContent.response && typeof jsonContent.response === 'string') {
            finalText = jsonContent.response;
            console.log('✅ [Question Detection] Extracted response from JSON, length:', finalText.length);
          } else if (jsonContent.content && typeof jsonContent.content === 'string') {
            finalText = jsonContent.content;
            console.log('✅ [Question Detection] Extracted content from JSON, length:', finalText.length);
          } else {
            console.warn('🔍 [Question Detection] JSON found but no details_markdown/user_friendly_summary/user_friendly/answer/text/response/content field');
            console.warn('🔍 [Question Detection] Available fields:', Object.keys(jsonContent));
            console.warn('🔍 [Question Detection] Full JSON content:', JSON.stringify(jsonContent, null, 2).substring(0, 500));
            // If JSON parsing worked but no text field found, try regex fallback on original text
            // Try details_markdown first, then user_friendly_summary, then user_friendly
            let regexMatch = null;
            const detailsMarkdownRegex = /"details_markdown"\s*:\s*"((?:[^"\\]|\\.|\\n|\\r|\\t)*)"/;
            regexMatch = finalText.match(detailsMarkdownRegex);
            if (!regexMatch) {
              const userFriendlySummaryRegex = /"user_friendly_summary"\s*:\s*"((?:[^"\\]|\\.|\\n|\\r|\\t)*)"/;
              regexMatch = finalText.match(userFriendlySummaryRegex);
            }
            if (!regexMatch) {
              const userFriendlyRegex = /"user_friendly"\s*:\s*"((?:[^"\\]|\\.|\\n|\\r|\\t)*)"/;
              regexMatch = finalText.match(userFriendlyRegex);
            }
            if (regexMatch && regexMatch[1]) {
              finalText = regexMatch[1]
                .replace(/\\"/g, '"')
                .replace(/\\n/g, '\n')
                .replace(/\\t/g, '\t')
                .replace(/\\r/g, '\r')
                .replace(/\\\\/g, '\\');
              console.log('✅ [Question Detection] Extracted text via regex fallback, length:', finalText.length);
            } else {
              console.warn('🔍 [Question Detection] Regex fallback also found no match');
            }
          }
        }
      } catch (e) {
        console.warn('🔍 [Question Detection] Failed to parse JSON from assistant message:', e);
        // Try to extract the answer value directly from the string as fallback
        try {
          // Try details_markdown first, then user_friendly_summary, then user_friendly
          let regexMatch = null;
          const detailsMarkdownRegex = /"details_markdown"\s*:\s*"((?:[^"\\]|\\.|\\n|\\r|\\t)*)"/s;
          regexMatch = finalText.match(detailsMarkdownRegex);
          if (!regexMatch) {
            const userFriendlySummaryRegex = /"user_friendly_summary"\s*:\s*"((?:[^"\\]|\\.|\\n|\\r|\\t)*)"/s;
            regexMatch = finalText.match(userFriendlySummaryRegex);
          }
          if (!regexMatch) {
            const userFriendlyRegex = /"user_friendly"\s*:\s*"((?:[^"\\]|\\.|\\n|\\r|\\t)*)"/s;
            regexMatch = finalText.match(userFriendlyRegex);
          }
          if (regexMatch && regexMatch[1]) {
            finalText = regexMatch[1]
              .replace(/\\"/g, '"')
              .replace(/\\n/g, '\n')
              .replace(/\\t/g, '\t')
              .replace(/\\r/g, '\r')
              .replace(/\\\\/g, '\\');
            console.log('✅ [Question Detection] Extracted text via regex fallback, length:', finalText.length);
          } else {
            console.warn('🔍 [Question Detection] Regex fallback found no match');
          }
        } catch (regexError) {
          console.warn('🔍 [Question Detection] Regex fallback also failed:', regexError);
        }
      }
    }
    
    // If we got text but it looks like a system prompt, try the previous assistant message
    if (finalText && (finalText.includes('BioAgent') || finalText.includes('autonomous bioinformatics assistant') || finalText.length > 5000)) {
      if (assistantMessages.length > 1) {
        const prevAssistant = assistantMessages[assistantMessages.length - 2];
        const prevText = extractAgentMessageText((prevAssistant as any).content ?? (prevAssistant as any).text ?? (prevAssistant as any).value);
        if (prevText && !prevText.includes('BioAgent') && prevText.length < 5000) {
          finalText = prevText;
        }
      }
    }
    
    // Debug logging
    if (agentOutput?.tool === 'agent' || agentOutput?.raw_result?.task_type === 'qa') {
      console.log('🔍 [Question Detection] Assistant messages count:', assistantMessages.length);
      console.log('🔍 [Question Detection] Final text length:', finalText?.length || 0);
      console.log('🔍 [Question Detection] Final text preview:', finalText?.substring(0, 100) || 'empty');
    }

    if (!finalText && typeof agentResult?.final_output === 'string') {
      finalText = agentResult.final_output;
    }
    if (!finalText && typeof agentResult?.response === 'string') {
      finalText = agentResult.response;
    }
    if (!finalText && typeof agentResult?.text === 'string') {
      finalText = agentResult.text;
    }
    // Also check top-level text from standard response
    if (!finalText && typeof agentOutput?.text === 'string') {
      finalText = agentOutput.text;
    }
    // Check raw_result for text fields
    if (!finalText && typeof agentOutput?.raw_result?.text === 'string') {
      finalText = agentOutput.raw_result.text;
    }
    if (!finalText && typeof actualResult?.text === 'string') {
      finalText = actualResult.text;
    }

    // Results viewer (S3 / FastQC HTML)
    // resultsLinks already computed above as _earlyLinks (before specialized renderers)
    const resultsLinks = _earlyLinks;
    const resultsVisualsRaw = agentOutput?.data?.visuals || actualResult?.visuals || rawResult?.visuals || [];

    const resultsVisuals = Array.isArray(resultsVisualsRaw)
      ? resultsVisualsRaw.map((v: any) => (v && typeof v === 'object' ? { ...v, url: normalizeAssetUrl(v.url) } : v))
      : [];
    const mainResultsIframe = Array.isArray(resultsVisuals)
      ? resultsVisuals.find((v: any) => v && v.type === 'iframe' && typeof v.url === 'string' && v.url.length > 0)
      : null;

    // Check if there are special visualizations (which would indicate it's not just a question)
    // This needs to be defined before isQuestion check
    const hasSpecialVisualizations = 
      treeNewick || 
      plasmidData || 
      plasmidResults || 
      qualityMetrics || 
      qualityPlotData || 
      alignedSequences || 
      plotData || 
      vendors ||
      agentOutput?.visualization_type === 'results_viewer' ||
      (Array.isArray(resultsVisuals) && resultsVisuals.length > 0);
    
    // Check if this is a question response (QA intent) - do this before checking finalText
    // Questions should only show the markdown answer, not debug info or full JSON
    const prompt = agentOutput?.prompt || '';
    const promptLower = prompt.toLowerCase().trim();
    const looksLikeQuestion = promptLower.endsWith('?') || 
      /^(what|how|why|when|where|who|which|can|could|should|would|is|are|does|do|will|tell me|explain|describe|help me|i need help|i want to know)/.test(promptLower);
    
    // Check for task_type: "qa" in various locations (this is a key indicator from the backend)
    let taskType = agentOutput?.raw_result?.task_type || 
                   actualResult?.task_type || 
                   agentResult?.task_type ||
                   rawResult?.task_type ||
                   rawResult?.result?.task_type;
    
    // If taskType is still undefined, try to extract it from the assistant message JSON
    // This handles the case where the agent returns JSON in a code block
    if (!taskType && finalText && (finalText.includes('```json') || finalText.includes('```')) && finalText.includes('task_type')) {
      try {
        const jsonMatch = finalText.match(/```(?:json)?\s*([\s\S]*?)\s*```/);
        if (jsonMatch && jsonMatch[1]) {
          const jsonContent = JSON.parse(jsonMatch[1]);
          taskType = jsonContent.task_type;
          console.log('✅ [Question Detection] Extracted task_type from JSON:', taskType);
        }
      } catch (e) {
        // Try regex fallback
        const taskTypeMatch = finalText.match(/"task_type"\s*:\s*"([^"]+)"/);
        if (taskTypeMatch && taskTypeMatch[1]) {
          taskType = taskTypeMatch[1];
          console.log('✅ [Question Detection] Extracted task_type via regex:', taskType);
        }
      }
    }
    
    // If taskType is still undefined, try to extract it from the assistant message JSON
    // This handles the case where the agent returns JSON in a code block
    if (!taskType && finalText && (finalText.includes('```json') || finalText.includes('```')) && finalText.includes('task_type')) {
      try {
        const jsonMatch = finalText.match(/```(?:json)?\s*([\s\S]*?)\s*```/);
        if (jsonMatch && jsonMatch[1]) {
          const jsonContent = JSON.parse(jsonMatch[1]);
          taskType = jsonContent.task_type;
          console.log('✅ [Question Detection] Extracted task_type from JSON:', taskType);
        }
      } catch (e) {
        // Try regex fallback
        const taskTypeMatch = finalText.match(/"task_type"\s*:\s*"([^"]+)"/);
        if (taskTypeMatch && taskTypeMatch[1]) {
          taskType = taskTypeMatch[1];
          console.log('✅ [Question Detection] Extracted task_type via regex:', taskType);
        }
      }
    }
    
    // Enhanced question detection with better logging
    const isQuestion = 
      taskType === 'qa' ||
      agentOutput?.tool === 'agent' || 
      agentOutput?.visualization_type === 'markdown' ||
      (agentOutput?.raw_result?.intent === 'qa') ||
      (agentOutput?.raw_result?.intent_reason && agentOutput.raw_result.intent_reason.includes('question')) ||
      (actualResult?.intent === 'qa') ||
      (actualResult?.intent_reason && actualResult.intent_reason.includes('question')) ||
      (rawResult?.intent === 'qa') ||
      (rawResult?.result?.intent === 'qa') ||
      (looksLikeQuestion && agentOutput?.tool === 'agent' && !hasSpecialVisualizations);
    
    // Debug logging for question detection
    if (agentOutput?.tool === 'agent' || looksLikeQuestion) {
      console.log('🔍 [Question Detection] taskType:', taskType);
      console.log('🔍 [Question Detection] isQuestion:', isQuestion);
      console.log('🔍 [Question Detection] tool:', agentOutput?.tool);
      console.log('🔍 [Question Detection] visualization_type:', agentOutput?.visualization_type);
      console.log('🔍 [Question Detection] hasSpecialVisualizations:', hasSpecialVisualizations);
    }

    // ── Tabular analysis execution result ────────────────────────────────────
    // Must be checked BEFORE the generic finalText block so the rich renderer
    // (title + stat tiles + plot + interpretation) isn't short-circuited.
    const taResult =
      agentOutput?.raw_result ?? agentOutput?.result?.raw_result ?? agentOutput;
    const isTaResult =
      (agentOutput?.tool === 'tabular_analysis' || agentOutput?.result?.tool === 'tabular_analysis') &&
      taResult?.status === 'success';

    if (isTaResult) {
      const interpretation: string = taResult?.text ?? '';
      const plotB64: string | undefined = taResult?.plot_base64;
      const planTitle: string | undefined = taResult?.plan_title;
      const resultData = taResult?.result_data;

      // Flatten result_data dict into displayable key-value stat tiles,
      // skipping internal keys like plot_path
      const SKIP_KEYS = new Set(['plot_path', 'plot_tmp']);

      /** Convert any result value to a concise human-readable string for a stat tile. */
      const formatStatValue = (v: unknown): string => {
        if (v === null || v === undefined) return '—';
        if (typeof v === 'boolean') return v ? 'true' : 'false';
        if (typeof v === 'string') return v.length > 60 ? v.slice(0, 57) + '…' : v;
        if (typeof v === 'number') {
          return Math.abs(v) < 0.001 || Math.abs(v) >= 1e6
            ? v.toExponential(3)
            : v.toPrecision(4);
        }
        if (Array.isArray(v)) {
          if (v.length === 0) return '(empty)';
          // If it's a flat array of numbers/strings, show first few values
          const allScalar = v.every((x) => typeof x !== 'object' || x === null);
          if (allScalar && v.length <= 5) return v.join(', ');
          if (allScalar) return `${v.slice(0, 3).join(', ')} … (${v.length} items)`;
          return `${v.length} items`;
        }
        if (typeof v === 'object') {
          const keys = Object.keys(v as object);
          if (keys.length === 0) return '{}';
          // Show key: value pairs for small flat dicts (like missing-value counts)
          const entries = Object.entries(v as Record<string, unknown>);
          const allNumeric = entries.every(([, val]) => typeof val === 'number');
          if (allNumeric && keys.length <= 4) {
            return entries.map(([k, val]) => `${k}: ${val}`).join(', ');
          }
          if (allNumeric) {
            return entries.slice(0, 3).map(([k, val]) => `${k}: ${val}`).join(', ') + ` … (${keys.length} cols)`;
          }
          return `${keys.length} entries`;
        }
        return String(v);
      };

      const resultEntries: Array<[string, string]> = (() => {
        if (!resultData) return [];
        const raw: Record<string, unknown> =
          resultData.type === 'dict' ? (resultData.value as Record<string, unknown>) ?? {} : {};
        return Object.entries(raw)
          .filter(([k]) => !SKIP_KEYS.has(k))
          .map(([k, v]) => {
            const label = k.replace(/_/g, ' ').replace(/\b\w/g, (c) => c.toUpperCase());
            return [label, formatStatValue(v)] as [string, string];
          });
      })();

      return (
        <div>
          {planTitle && (
            <div style={{ fontWeight: 700, fontSize: '0.95rem', color: '#1e293b', marginBottom: 10 }}>
              🧬 {planTitle}
            </div>
          )}

          {/* Computed numeric values (e.g. Pearson r, p-value) */}
          {resultEntries.length > 0 && (
            <div style={{ display: 'flex', flexWrap: 'wrap', gap: 8, marginBottom: 14 }}>
              {resultEntries.map(([label, value]) => (
                <div
                  key={label}
                  style={{
                    background: '#f0fdf4',
                    border: '1px solid #bbf7d0',
                    borderRadius: 7,
                    padding: '6px 12px',
                    minWidth: 110,
                  }}
                >
                  <div style={{ fontSize: '0.68rem', color: '#166534', fontWeight: 600, textTransform: 'uppercase', letterSpacing: '0.04em' }}>
                    {label}
                  </div>
                  <div style={{ fontSize: '0.95rem', fontWeight: 700, color: '#14532d', fontFamily: 'monospace' }}>
                    {value}
                  </div>
                </div>
              ))}
            </div>
          )}

          {plotB64 && (
            <div style={{ marginBottom: 14, textAlign: 'center' }}>
              <img
                src={plotB64}
                alt="Analysis plot"
                style={{ maxWidth: '100%', borderRadius: 6, border: '1px solid #e2e8f0' }}
              />
            </div>
          )}
          {interpretation && (
            <div style={{ lineHeight: 1.7, color: '#334155', fontSize: '0.88rem' }}>
              <ReactMarkdown remarkPlugins={[remarkGfm]}>{interpretation}</ReactMarkdown>
            </div>
          )}
          {downloadLinksSection}
        </div>
      );
    }

    if (finalText) {
      // ── Advisory JSON detection ──────────────────────────────────────────
      // The backend normalises advisory responses to { helix_type: "advisory", ... }.
      // Also handle legacy / cached variants with or without code-fence wrapping.
      {
        let parseCandidate = finalText.trimStart();
        // Strip fenced code blocks (```json ... ``` or ``` ... ```)
        if (parseCandidate.startsWith('```')) {
          const fenceMatch = parseCandidate.match(/```(?:json)?\s*([\s\S]*?)\s*```/);
          if (fenceMatch) parseCandidate = fenceMatch[1].trim();
        }
        if (parseCandidate.startsWith('{')) {
          try {
            const maybeAdvisory = JSON.parse(parseCandidate);
            if (maybeAdvisory && typeof maybeAdvisory === 'object') {
              const isAdvisory = (
                // Canonical (backend-normalised)
                maybeAdvisory.helix_type === 'advisory' ||
                // Legacy shape 1: planning advisor
                maybeAdvisory.classification ||
                maybeAdvisory.recommended_workflow ||
                (maybeAdvisory.requirements && !maybeAdvisory.domain) ||
                maybeAdvisory.minimum_information_needed_from_you ||
                // Legacy shape 2: explanation advisor
                (maybeAdvisory.type === 'answer' && Array.isArray(maybeAdvisory.sections)) ||
                (Array.isArray(maybeAdvisory.sections) && maybeAdvisory.summary)
              );
              if (isAdvisory) {
                return <div>{renderAdvisoryJSON(maybeAdvisory)}</div>;
              }
            }
          } catch {
            // Not valid JSON — fall through to markdown rendering
          }
        }
      }

      // If it's a question and has no special visualizations, render only the markdown
      if (isQuestion && !hasSpecialVisualizations) {
        return (
          <div>
            <div className="agent-response-markdown">
              <ReactMarkdown remarkPlugins={[remarkGfm]}>{finalText}</ReactMarkdown>
            </div>
          </div>
        );
      }
      
      // Otherwise, render with debug info (for commands/actions)
      return (
        <div>
          <div className="agent-response-markdown">
            <ReactMarkdown remarkPlugins={[remarkGfm]}>{finalText}</ReactMarkdown>
          </div>

          {/* Results Viewer (embed main HTML report + links) */}
          {agentOutput?.visualization_type === 'results_viewer' && (
            <div className="p-3 border rounded mb-3" style={{ background: '#F8FAFC' }}>
              {/* iFrame reports (e.g. FastQC HTML) */}
              {mainResultsIframe?.url && (
                <div style={{ width: '100%', marginBottom: '16px' }}>
                  <iframe
                    src={mainResultsIframe.url}
                    title={mainResultsIframe.title || 'Results report'}
                    style={{ width: '100%', height: '720px', border: 0, background: '#fff' }}
                  />
                </div>
              )}

              {/* Image visuals — supports both url-based ('image') and inline base64 ('image_b64') */}
              {Array.isArray(resultsVisuals) && resultsVisuals.some((v: any) => v?.type === 'image' || v?.type === 'image_b64') && (
                <div style={{ display: 'flex', flexWrap: 'wrap', gap: '16px', marginBottom: '16px' }}>
                  {resultsVisuals
                    .filter((v: any) => (v?.type === 'image' && v?.url) || (v?.type === 'image_b64' && v?.data))
                    .map((v: any, idx: number) => {
                      const imgSrc = v.type === 'image_b64'
                        ? `data:image/png;base64,${v.data}`
                        : v.url;
                      return (
                        <div key={idx} style={{ flex: '1 1 420px', minWidth: '280px', maxWidth: '600px' }}>
                          <div style={{ fontSize: '0.8rem', fontWeight: 600, color: '#475569', marginBottom: '4px' }}>
                            {v.title || `Plot ${idx + 1}`}
                          </div>
                          <a href={imgSrc} target="_blank" rel="noreferrer" title="Open full-size">
                            <img
                              src={imgSrc}
                              alt={v.title || `Plot ${idx + 1}`}
                              style={{ width: '100%', borderRadius: '6px', border: '1px solid #E2E8F0', background: '#fff', cursor: 'pointer' }}
                            />
                          </a>
                        </div>
                      );
                    })}
                </div>
              )}

              {/* Downloadable artifacts */}
              {downloadLinksSection}
            </div>
          )}
          
          {/* Debug Section - Consistent across all commands */}
          {renderDebugInfo(agentOutput, actualResult)}
        </div>
      );
    }

    // If no finalText but it's a question, try to extract from messages more aggressively
    if (isQuestion && !hasSpecialVisualizations && assistantMessages.length > 0) {
      // Try to find any assistant message with content (skip system prompts)
      for (let i = assistantMessages.length - 1; i >= 0; i--) {
        const msg = assistantMessages[i];
        const extracted = extractAgentMessageText(msg.content ?? msg.text ?? msg.value);
        if (extracted && extracted.trim().length > 0) {
          // Skip system prompts - they're usually very long and contain specific keywords
          const isSystemPrompt = extracted.includes('BioAgent') || 
                                 extracted.includes('autonomous bioinformatics assistant') ||
                                 extracted.includes('Core Capabilities:') ||
                                 extracted.includes('Non-Negotiable Principles') ||
                                 extracted.length > 5000;
          
          if (!isSystemPrompt) {
            return (
              <div>
                <div className="agent-response-markdown">
                  <ReactMarkdown remarkPlugins={[remarkGfm]}>{extracted}</ReactMarkdown>
                </div>
              </div>
            );
          }
        }
      }
    }

    // CRITICAL: Don't render structured data for questions - this is what's causing the JSON to show
    // If it's a question, we must return early and never call renderStructuredData
    if (isQuestion && !hasSpecialVisualizations) {
      // First, try one more time to extract from all messages
      if (messages.length > 0) {
        // Look through all messages (not just assistant) to find the answer
        for (let i = messages.length - 1; i >= 0; i--) {
          const msg = messages[i];
          // Skip system messages
          if (msg?.role === 'system' || msg?.type === 'system') continue;
          
          const extracted = extractAgentMessageText(msg.content ?? msg.text ?? msg.value);
          if (extracted && extracted.trim().length > 0) {
            const isSystemPrompt = extracted.includes('BioAgent') || 
                                   extracted.includes('autonomous bioinformatics assistant') ||
                                   extracted.includes('Core Capabilities:') ||
                                   extracted.includes('Non-Negotiable Principles') ||
                                   extracted.length > 5000;
            
            if (!isSystemPrompt) {
              return (
                <div>
                  <div className="agent-response-markdown">
                    <ReactMarkdown remarkPlugins={[remarkGfm]}>{extracted}</ReactMarkdown>
                  </div>
                </div>
              );
            }
          }
        }
      }
      // Last resort: try to find any text content in the response
      const allTextFields = [
        agentOutput?.text,
        agentOutput?.raw_result?.text,
        actualResult?.text,
        agentResult?.text,
        agentResult?.response,
        agentResult?.final_output,
      ].filter(Boolean);
      
      if (allTextFields.length > 0) {
        const foundText = allTextFields[0];
        if (typeof foundText === 'string' && foundText.trim().length > 0) {
          return (
            <div>
              <div className="agent-response-markdown">
                <ReactMarkdown remarkPlugins={[remarkGfm]}>{foundText}</ReactMarkdown>
              </div>
            </div>
          );
        }
      }
      
      // If we still can't find text, show a message instead of structured data
      // This prevents renderStructuredData from being called
      return (
        <div>
          <div className="alert alert-info">
            <p>Question detected but answer could not be extracted from the response structure.</p>
          </div>
        </div>
      );
    }

    // CRITICAL: Never render structured data for questions - this is the root cause
    // Even if isQuestion check somehow fails, we should check again here
    const finalIsQuestion = isQuestion || 
                            taskType === 'qa' || 
                            (agentOutput?.tool === 'agent' && !hasSpecialVisualizations && looksLikeQuestion);
    
    if (finalIsQuestion && !hasSpecialVisualizations) {
      console.log('🔍 [Question Detection] Final check: This is a question, preventing structured data render');
      // Return a simple message instead of structured data
      return (
        <div>
          <div className="alert alert-warning">
            <p>Question detected but answer text could not be extracted. Check console for response structure.</p>
            <details className="mt-2">
              <summary className="text-muted small" style={{ cursor: 'pointer' }}>Debug: Response Structure</summary>
              <pre className="small mt-2 bg-white p-2 border rounded" style={{maxHeight: '200px', overflow: 'auto'}}>
                {JSON.stringify({taskType, isQuestion, finalIsQuestion, messagesCount: messages.length, assistantMessagesCount: assistantMessages.length}, null, 2)}
              </pre>
            </details>
          </div>
        </div>
      );
    }
    
    // Only render structured data if it's NOT a question
    // This check prevents questions from showing JSON structure
    if (!finalIsQuestion && renderStructuredData(agentResult)) {
      return (
        <div>
          <div className="bg-light p-3 border rounded">
            {renderStructuredData(agentResult)}
          </div>
          
          {/* Debug Section - Consistent across all commands */}
          {renderDebugInfo(agentOutput, actualResult)}
        </div>
      );
    }
    
    // Final fallback - but only if NOT a question
    if (!finalIsQuestion) {
      return (
        <div>
          <details className="bg-light p-2 border rounded">
            <summary className="text-muted" style={{ cursor: 'pointer' }}>View Raw Agent Output</summary>
            <pre className="bg-white border rounded p-2 mt-2">{JSON.stringify(agentResult, null, 2)}</pre>
          </details>
          
          {/* Debug Section - Consistent across all commands */}
          {renderDebugInfo(agentOutput, actualResult)}
        </div>
      );
    }
    
    // If we get here and it's a question, something went wrong
    return (
      <div>
        <div className="alert alert-warning">
          <p>Unable to render question response. Please check the browser console for details.</p>
        </div>
      </div>
    );

    return (
      <div>
        <details className="bg-light p-2 border rounded">
          <summary className="text-muted" style={{ cursor: 'pointer' }}>View Raw Agent Output</summary>
          <pre className="bg-white border rounded p-2 mt-2">{JSON.stringify(agentResult, null, 2)}</pre>
        </details>
        
        {/* Debug Section - Consistent across all commands */}
        {renderDebugInfo(agentOutput, actualResult)}
      </div>
    );
  };


  // Helper function to render consistent Debug Info section
  const renderDebugInfo = (output: any, actualResult?: any) => {
    const outputKeys = output && typeof output === 'object' ? Object.keys(output) : [];
    const resultKeys = output?.result && typeof output.result === 'object' ? Object.keys(output.result) : [];
    const actualResultKeys = actualResult && typeof actualResult === 'object' ? Object.keys(actualResult) : [];
    
    return (
      <details className="bg-light p-2 border rounded mb-3">
        <summary className="text-muted" style={{cursor: 'pointer'}}>
          🔍 Debug Info (click to expand)
        </summary>
        <div className="mt-2">
          <p><strong>Output keys:</strong> {outputKeys.length > 0 ? outputKeys.join(', ') : 'No output'}</p>
          <p><strong>Result keys:</strong> {resultKeys.length > 0 ? resultKeys.join(', ') : 'No result'}</p>
          {actualResultKeys.length > 0 && (
            <p><strong>ActualResult keys:</strong> {actualResultKeys.join(', ')}</p>
          )}
          <details className="mt-2">
            <summary className="text-muted small" style={{cursor: 'pointer'}}>Full Output (JSON)</summary>
            <pre className="small mt-2 bg-white p-2 border rounded" style={{maxHeight: '400px', overflow: 'auto'}}>
              {JSON.stringify(output, null, 2)}
            </pre>
          </details>
          {actualResult && (
            <details className="mt-2">
              <summary className="text-muted small" style={{cursor: 'pointer'}}>Full ActualResult (JSON)</summary>
              <pre className="small mt-2 bg-white p-2 border rounded" style={{maxHeight: '400px', overflow: 'auto'}}>
                {JSON.stringify(actualResult, null, 2)}
              </pre>
            </details>
          )}
        </div>
      </details>
    );
  };

  // ── Advisory JSON renderer ───────────────────────────────────────────────
  // The LLM sometimes returns a rich structured JSON advisory plan (classification,
  // answer summary, requirements, workflow steps, questions for the user).
  // This helper detects that shape and renders it as a clean, sectioned card.
  // ── Canonical HelixAdvisory renderer ────────────────────────────────────
  // Renders a normalised `{ helix_type: "advisory", ... }` object produced by
  // the backend advisory_normalizer.  Also handles legacy ad-hoc shapes as a
  // fallback so old cached responses continue to display correctly.
  const renderAdvisoryJSON = (json: Record<string, any>): React.ReactNode => {
    // ── Normalise to canonical shape on the client (legacy / cached responses)
    // The backend normalises new responses, but old history entries may still
    // carry the raw LLM JSON.  Re-normalise here for safety.
    const isCanonical = json.helix_type === 'advisory';

    // Pull fields from either the canonical schema or any legacy variant.
    const cls = json.classification as { domain?: string; task_type?: string; feasible?: boolean } | undefined;

    const title: string = json.title || cls?.task_type || '';

    // Summary: canonical → json.summary; legacy shape 1 → json.answer?.summary; legacy shape 2 → json.summary
    const summary: string =
      json.summary ||
      (typeof json.answer === 'string' ? json.answer : json.answer?.summary) ||
      '';

    // Sections (canonical + shape 2)
    const sections: Array<{
      heading: string;
      content?: string;
      items?: Array<{ label?: string; description?: string; examples?: string[]; tools?: string[]; term?: string; meaning?: string; metric?: string; interpretation?: string } | string>;
    }> = json.sections || [];

    // Workflow steps (canonical: workflow_steps; legacy shape 1: recommended_workflow)
    const workflowSteps: Array<{ step: number; name: string; description?: string }> =
      json.workflow_steps || json.recommended_workflow || [];

    // Requirements (canonical: requirements[]; legacy shape 1: requirements.{input_data,...})
    // Flatten all requirement groups into a single list for rendering
    const requirementsList: Array<{ label: string; description?: string; examples?: string[]; tools?: string[] }> = (() => {
      if (Array.isArray(json.requirements)) return json.requirements;
      if (json.requirements && typeof json.requirements === 'object') {
        const { input_data = [], reference_files = [], software_or_pipeline_components = [] } = json.requirements as any;
        return [
          ...input_data.map((d: any) => ({ label: d.item || d.label || '', description: d.details || d.description, examples: d.examples })),
          ...reference_files.map((f: any) => ({ label: f.item || f.label || '', description: f.details || f.description, examples: f.examples })),
          ...software_or_pipeline_components.map((s: any) => ({ label: s.step || s.label || '', description: s.purpose || s.description, tools: s.tools })),
        ];
      }
      return [];
    })();

    // Questions for user (canonical; legacy: minimum_information_needed_from_you)
    const questions: Array<{ label?: string; question?: string; description?: string; examples?: string[] }> =
      json.questions_for_user || json.minimum_information_needed_from_you || [];

    // Next steps (canonical: next_steps[]; legacy: next_step?.message)
    const nextSteps: string[] = Array.isArray(json.next_steps) ? json.next_steps : (
      json.next_step?.message ? [json.next_step.message] : []
    );

    // ── Helpers ──────────────────────────────────────────────────────────
    const sectionLabel = (text: string, color: string) => (
      <div style={{ fontWeight: 600, color, marginBottom: 8, textTransform: 'uppercase' as const, fontSize: '0.72rem', letterSpacing: '0.05em' }}>
        {text}
      </div>
    );

    const badge = (label: string, bg: string, fg: string) => (
      <span style={{ background: bg, color: fg, borderRadius: 6, padding: '2px 10px', fontSize: '0.75rem', fontWeight: 600, marginRight: 6 }}>
        {label}
      </span>
    );

    const renderItem = (item: any, i: number, bullet: string, bulletColor: string) => {
      const label = item.label || item.item || item.step || item.metric || item.term || item.question || (typeof item === 'string' ? item : '');
      const desc = item.description || item.details || item.purpose || item.content || item.meaning || item.interpretation || '';
      const examples: string[] = item.examples || [];
      const tools: string[] = item.tools || [];
      return (
        <div key={i} style={{ display: 'flex', gap: 10, alignItems: 'flex-start' }}>
          <span style={{ color: bulletColor, fontWeight: 700, minWidth: 14, marginTop: 1 }}>{bullet}</span>
          <div>
            <span style={{ fontWeight: 600, color: '#0F172A' }}>{label}</span>
            {tools.length > 0 && (
              <span style={{ marginLeft: 6 }}>
                {tools.map((t: string) => (
                  <span key={t} style={{ background: '#E0F2FE', color: '#0369A1', borderRadius: 4, padding: '1px 6px', fontSize: '0.75rem', fontFamily: 'monospace', marginRight: 4 }}>{t}</span>
                ))}
              </span>
            )}
            {examples.length > 0 && <span style={{ color: '#9CA3AF', fontSize: '0.8rem' }}> ({examples.join(', ')})</span>}
            {desc && <div style={{ color: '#64748B', marginTop: 1, lineHeight: 1.5 }}>{desc}</div>}
          </div>
        </div>
      );
    };

    // ── Render ────────────────────────────────────────────────────────────
    return (
      <div style={{ display: 'flex', flexDirection: 'column', gap: 18, fontSize: '0.88rem' }}>

        {/* Title + classification badges */}
        <div style={{ display: 'flex', alignItems: 'center', gap: 10, flexWrap: 'wrap' }}>
          {title && <span style={{ fontWeight: 700, fontSize: '1rem', color: '#0F172A' }}>{title}</span>}
          {cls?.domain && badge(cls.domain, '#EFF6FF', '#1D4ED8')}
          {cls?.task_type && !title && badge(cls.task_type, '#F0FDF4', '#166534')}
          {cls?.feasible === true && badge('✓ Feasible', '#ECFDF5', '#059669')}
          {cls?.feasible === false && badge('⚠ Needs Review', '#FFF7ED', '#C2410C')}
        </div>

        {/* Summary */}
        {summary && (
          <div style={{ background: '#F8FAFC', border: '1px solid #E2E8F0', borderRadius: 8, padding: '12px 16px', lineHeight: 1.7, color: '#334155' }}>
            {summary}
          </div>
        )}

        {/* Sections (canonical + shape 2 explanation sections) */}
        {sections.map((sec, si) => (
          <div key={si}>
            {sectionLabel(sec.heading, '#0369A1')}
            {sec.content && (
              <div style={{ color: '#334155', lineHeight: 1.7, marginBottom: sec.items ? 8 : 0 }}>
                {sec.content}
              </div>
            )}
            {Array.isArray(sec.items) && sec.items.length > 0 && (
              <div style={{ display: 'flex', flexDirection: 'column', gap: 6 }}>
                {sec.items.map((item, ii) => renderItem(item, ii, '•', '#3B82F6'))}
              </div>
            )}
          </div>
        ))}

        {/* Flat requirements list (legacy shape 1 or canonical) */}
        {!isCanonical && requirementsList.length > 0 && (
          <div>
            {sectionLabel('Requirements', '#1E40AF')}
            <div style={{ display: 'flex', flexDirection: 'column', gap: 6 }}>
              {requirementsList.map((r, i) => renderItem(r, i, '•', '#3B82F6'))}
            </div>
          </div>
        )}
        {isCanonical && requirementsList.length > 0 && (
          <div>
            {sectionLabel('Requirements', '#1E40AF')}
            <div style={{ display: 'flex', flexDirection: 'column', gap: 6 }}>
              {requirementsList.map((r, i) => renderItem(r, i, '•', '#3B82F6'))}
            </div>
          </div>
        )}

        {/* Workflow steps */}
        {workflowSteps.length > 0 && (
          <div>
            {sectionLabel('Recommended Workflow', '#047857')}
            <div style={{ display: 'flex', flexDirection: 'column', gap: 8 }}>
              {workflowSteps.map((w: any) => (
                <div key={w.step} style={{ display: 'flex', gap: 12, alignItems: 'flex-start' }}>
                  <div style={{ minWidth: 26, height: 26, borderRadius: '50%', background: '#D1FAE5', color: '#065F46', fontWeight: 700, fontSize: '0.78rem', display: 'flex', alignItems: 'center', justifyContent: 'center', flexShrink: 0 }}>
                    {w.step}
                  </div>
                  <div>
                    <span style={{ fontWeight: 600, color: '#0F172A' }}>{w.name}</span>
                    {w.description && <div style={{ color: '#64748B', marginTop: 2, lineHeight: 1.5 }}>{w.description}</div>}
                  </div>
                </div>
              ))}
            </div>
          </div>
        )}

        {/* Questions for the user */}
        {questions.length > 0 && (
          <div style={{ background: '#FFFBEB', border: '1px solid #FDE68A', borderRadius: 8, padding: '12px 16px' }}>
            {sectionLabel('To proceed, I need the following information', '#92400E')}
            <div style={{ display: 'flex', flexDirection: 'column', gap: 6 }}>
              {questions.map((q: any, i) => renderItem({ ...q, label: q.label || q.question || '' }, i, '?', '#D97706'))}
            </div>
          </div>
        )}

        {/* Next steps */}
        {nextSteps.length > 0 && (
          <div style={{ borderTop: '1px solid #E2E8F0', paddingTop: 12 }}>
            {sectionLabel('Next Steps', '#475569')}
            <div style={{ display: 'flex', flexDirection: 'column', gap: 4 }}>
              {nextSteps.map((ns, i) => (
                <div key={i} style={{ display: 'flex', gap: 8, color: '#475569', lineHeight: 1.6 }}>
                  <span style={{ color: '#94A3B8' }}>→</span>
                  <span>{ns}</span>
                </div>
              ))}
            </div>
          </div>
        )}
      </div>
    );
  };

  const renderOutput = (item: HistoryItem) => {
    const { output, type } = item;

    // ── Async Nextflow pipeline job card ──────────────────────────────────
    const pj = output?._pipelineJob;
    if (pj) {
      const statusColors: Record<string, string> = {
        submitted: '#6366f1',
        running:   '#f59e0b',
        completed: '#10b981',
        failed:    '#ef4444',
      };
      const color = statusColors[pj.status] ?? '#6b7280';

      const statusIcon: Record<string, string> = {
        submitted: '⏳',
        running:   '⚙️',
        completed: '✅',
        failed:    '❌',
      };
      const icon = statusIcon[pj.status] ?? '🔄';

      const steps: Array<{ type: string; process?: string }> = pj.steps ?? [];

      return (
        <div style={{
          border: `1.5px solid ${color}`,
          borderRadius: '8px',
          padding: '14px 18px',
          background: '#fafafa',
          maxWidth: '640px',
        }}>
          {/* Header */}
          <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '10px' }}>
            <span style={{ fontSize: '20px' }}>{icon}</span>
            <div>
              <div style={{ fontWeight: 600, color: '#111' }}>
                Pipeline: {pj.pipeline || output?.tool || 'Nextflow'}
              </div>
              <div style={{ fontSize: '12px', color: '#6b7280' }}>
                Job <code style={{ background: '#f3f4f6', padding: '1px 4px', borderRadius: '3px' }}>
                  {(pj.jobId ?? '').slice(0, 8)}…
                </code>
                {' · '}
                <span style={{ color, fontWeight: 500, textTransform: 'capitalize' }}>{pj.status}</span>
              </div>
            </div>
            {pj.status === 'running' && (
              <div style={{
                marginLeft: 'auto',
                width: '18px', height: '18px',
                border: '2px solid #f59e0b',
                borderTopColor: 'transparent',
                borderRadius: '50%',
                animation: 'spin 0.8s linear infinite',
              }} />
            )}
          </div>

          {/* Step timeline */}
          {steps.length > 0 && (
            <div style={{ marginTop: '8px', borderTop: '1px solid #e5e7eb', paddingTop: '8px' }}>
              <div style={{ fontSize: '12px', color: '#6b7280', marginBottom: '4px' }}>Steps</div>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '3px' }}>
                {steps.map((s, i) => (
                  <div key={i} style={{ display: 'flex', alignItems: 'center', gap: '6px', fontSize: '12px' }}>
                    <span style={{ color: s.type === 'process_completed' ? '#10b981' : '#6b7280' }}>
                      {s.type === 'process_completed' ? '✓' : s.type === 'failed' ? '✗' : '·'}
                    </span>
                    <span style={{ color: '#374151' }}>{s.process || s.type}</span>
                  </div>
                ))}
              </div>
            </div>
          )}

          {/* Result files */}
          {pj.resultFiles && pj.resultFiles.length > 0 && (
            <div style={{ marginTop: '8px', borderTop: '1px solid #e5e7eb', paddingTop: '8px' }}>
              <div style={{ fontSize: '12px', color: '#6b7280', marginBottom: '4px' }}>Output files</div>
              {pj.resultFiles.map((f: string, i: number) => (
                <div key={i} style={{ fontSize: '12px', color: '#374151', fontFamily: 'monospace' }}>
                  📄 {f}
                </div>
              ))}
            </div>
          )}

          {/* Error detail */}
          {pj.error && (
            <div style={{
              marginTop: '8px', padding: '8px', borderRadius: '4px',
              background: '#fef2f2', color: '#dc2626', fontSize: '12px',
              fontFamily: 'monospace', whiteSpace: 'pre-wrap', maxHeight: '120px', overflow: 'auto',
            }}>
              {pj.error}
            </div>
          )}

          {/* Stream URL hint while running */}
          {pj.status === 'submitted' || pj.status === 'running' ? (
            <div style={{ marginTop: '8px', fontSize: '11px', color: '#9ca3af' }}>
              Live stream: {pj.streamUrl}
            </div>
          ) : null}
        </div>
      );
    }
    // ── end pipeline job card ─────────────────────────────────────────────

    if (type === 'agent_error') {
      const message = output?.error || 'Agent encountered an error while processing the request.';
      return <div className="alert alert-danger mb-0">{message}</div>;
    }

    // Resolve scenario by explicit scenarioId first, then by command/tool match.
    // When backend returns tool "__plan__", no scenario has that tool — match by command so "Load & run" still shows.
    const responseTool =
      output?.tool || output?.tool_name || output?.result?.tool;
    const planFirstStepTool =
      output?.raw_result?.result?.steps?.[0]?.tool_name ||
      output?.result?.result?.steps?.[0]?.tool_name ||
      output?.result?.steps?.[0]?.tool_name;
    const commandForMatch = item.input ?? output?.prompt ?? '';
    const scenario =
      (item.scenarioId ? getDemoScenarioById(item.scenarioId) : undefined) ??
      (responseTool === '__plan__'
        ? getDemoScenarioByCommand(commandForMatch) ||
          (planFirstStepTool ? getDemoScenarioByTool(planFirstStepTool) : undefined)
        : responseTool
          ? getDemoScenarioByCommandAndTool(commandForMatch, responseTool) ||
            getDemoScenarioByTool(responseTool)
          : undefined);

    // ── Single source of truth: backend workflow_state ───────────────────────
    // Canonical field emitted by build_standard_response. Fall back to heuristics
    // only if talking to an older backend that doesn't emit it yet.
    const workflowState: string = (
      output?.workflow_state ||
      output?.result?.workflow_state ||
      output?.raw_result?.workflow_state ||
      ''
    ).toUpperCase();

    const resolvedStatus =
      output?.status ||
      output?.result?.status ||
      output?.raw_result?.status ||
      output?.raw_result?.result?.status ||
      output?.data?.results?.status ||
      output?.data?.results?.result?.status;

    // If backend sends workflow_state, trust it completely.
    // Otherwise fall back to the old flag/status heuristics.
    const planAction: 'approve' | 'execute' | 'none' = (() => {
      if (workflowState === 'WAITING_FOR_APPROVAL') return 'approve';
      if (workflowState === 'PLANNING') {
        // PLANNING means execute_ready=true + no approval gate
        const er = output?.execute_ready ?? output?.result?.execute_ready ?? false;
        return er ? 'execute' : 'none';
      }
      if (workflowState === 'WAITING_FOR_INPUTS' || workflowState === 'IDLE' ||
          workflowState === 'COMPLETED' || workflowState === 'FAILED' ||
          workflowState === 'EXECUTING') return 'none';
      // Legacy fallback for older backend responses without workflow_state
      const status = (resolvedStatus || '').toLowerCase();
      if (status !== 'workflow_planned') return 'none';
      const approvalRequired =
        output?.approval_required === true ||
        output?.result?.approval_required === true ||
        output?.raw_result?.approval_required === true;
      const executeReady =
        output?.execute_ready === true ||
        output?.result?.execute_ready === true;
      if (approvalRequired) return 'approve';
      if (executeReady) return 'execute';
      return 'none';
    })();

    const isApprovalPlan = planAction === 'approve';
    const showExecutePipeline = planAction === 'execute';
    const isWorkflowPlan = planAction !== 'none' || workflowState === 'PLANNING';

    // needs_inputs: waiting for data before planning/execution can proceed
    const explicitNeedsInputs =
      workflowState === 'WAITING_FOR_INPUTS' ||
      (resolvedStatus || '').toLowerCase() === 'needs_inputs';

    const getAssistantTextForDetection = (o: any): string => {
      const msgs =
        o?.raw_result?.messages ||
        o?.result?.raw_result?.messages ||
        o?.data?.results?.messages ||
        o?.data?.results?.result?.messages ||
        [];
      if (!Array.isArray(msgs) || msgs.length === 0) return '';
      const lastAssistant = [...msgs].reverse().find((m: any) => {
        const t = (m?.type || '').toString().toLowerCase();
        const r = (m?.role || '').toString().toLowerCase();
        return t === 'ai' || t === 'assistant' || r === 'assistant';
      });
      return extractAgentMessageText(
        (lastAssistant as any)?.content ?? (lastAssistant as any)?.text ?? (lastAssistant as any)?.value
      );
    };

    const agentText = getAssistantTextForDetection(output);
    const looksLikeNeedsInputsFromAgent =
      (output?.tool === 'agent' || output?.result?.tool === 'agent' || output?.visualization_type === 'markdown') &&
      /please provide|required information|required inputs|before .* can .* proceed/i.test(agentText) &&
      /count matrix|sample metadata|design formula/i.test(agentText);

    // Some demo scenarios are explicitly authored as "needs_inputs". Historically
    // we forced the needs-inputs UX whenever the backend emitted workflow_planned
    // for these demos, because the router used to produce only thin/empty plans.
    // Now the router emits rich `router_reasoning.suggested_steps`; when those
    // exist we trust the backend and render the WorkflowPlanCard so the user can
    // see the proposed pipeline (still alongside the "Load & run" shortcut).
    const planSuggestedSteps =
      ((output?.data as any)?.results?.steps?.[0]?.arguments?.router_reasoning?.suggested_steps as
        | string[]
        | undefined) ?? [];
    const hasRichPlan = Array.isArray(planSuggestedSteps) && planSuggestedSteps.length > 0;

    const scenarioForcesNeedsInputs =
      isWorkflowPlan &&
      !!scenario?.followUpPrompt &&
      scenario?.expectedBehavior === 'needs_inputs' &&
      !hasRichPlan;

    const isNeedsInputs = explicitNeedsInputs || looksLikeNeedsInputsFromAgent || scenarioForcesNeedsInputs;

    const renderedResponse = renderAgentResponse(output);

    // ── Tabular analysis plan — rendered as AnalysisPlanCard ─────────────────
    // Detect the structured plan object produced by analysis_planner.py
    const analysisPlan =
      output?.result?.analysis_plan ??
      output?.analysis_plan ??
      output?.raw_result?.analysis_plan ??
      output?.data?.analysis_plan;

    if (analysisPlan && isWorkflowPlan && !scenarioForcesNeedsInputs) {
      const alreadyApproved = executedPipelineCommands.has(item.input);
      const introText: string =
        output?.result?.text ?? output?.text ?? output?.raw_result?.text ?? '';
      return (
        <div>
          {introText && (
            <p style={{ color: '#475569', fontSize: '0.88rem', marginBottom: 12 }}>
              {introText}
            </p>
          )}
          <AnalysisPlanCard
            plan={analysisPlan}
            onApprove={() => executeCommand('Approve.', item.scenarioId, false)}
            onCancel={() => {/* no-op — user can scroll away */}}
            loading={loading}
            approved={alreadyApproved}
          />
        </div>
      );
    }

    // Workflow plan: __plan__ tool → rich WorkflowPlanCard; other plans → generic renderer.
    if (isWorkflowPlan && !scenarioForcesNeedsInputs) {
      const alreadyExecuted = executedPipelineCommands.has(item.input);
      const isRouterPlan = output?.tool === '__plan__' || output?.result?.tool === '__plan__';
      return (
        <div>
          {isRouterPlan ? (
            <WorkflowPlanCard
              output={output as Record<string, unknown>}
              planAction={planAction}
              alreadyApproved={alreadyExecuted}
              loading={loading}
              onApprove={() => executeCommand('Approve.', item.scenarioId, false)}
              onExecute={() => handleExecutePipeline(item.input)}
            />
          ) : (
            <>
              {renderedResponse}
              <div className="mt-3 d-flex align-items-center gap-3">
                {alreadyExecuted ? (
                  <Button variant="outline-success" disabled className="px-4">
                    ✓ Submitted
                  </Button>
                ) : (
                  <>
                    {planAction === 'approve' && (
                      <Button
                        variant="primary"
                        onClick={() => executeCommand('Approve.', item.scenarioId, false)}
                        disabled={loading}
                        className="px-4"
                      >
                        {loading ? 'Approving…' : '✅ I approve'}
                      </Button>
                    )}
                    {showExecutePipeline && (
                      <Button
                        variant="success"
                        onClick={() => handleExecutePipeline(item.input)}
                        disabled={loading}
                        className="px-4"
                      >
                        {loading ? 'Submitting…' : '▶ Execute Pipeline'}
                      </Button>
                    )}
                  </>
                )}
                <span className="text-muted small">
                  {alreadyExecuted
                    ? 'Pipeline submitted — see results above.'
                    : planAction === 'approve'
                      ? 'Approval confirms execution of the pending reviewed plan.'
                      : planAction === 'execute'
                        ? 'Confirms the plan above and queues all steps for execution.'
                        : 'This is a design/preview plan. Provide inputs or approval context to continue.'}
                </span>
              </div>
            </>
          )}
          {scenario?.followUpPrompt && (
            <div
              className="mt-3 p-3 rounded-3 d-flex align-items-start gap-3"
              style={{
                background: 'linear-gradient(135deg, #EFF6FF, #F0FDF4)',
                border: '1px solid #BFDBFE',
              }}
            >
              <span style={{ fontSize: '1.4rem', lineHeight: 1, flexShrink: 0 }}>📋</span>
              <div className="flex-grow-1">
                <div style={{ fontWeight: 600, fontSize: '0.88rem', color: '#1E40AF', marginBottom: '4px' }}>
                  Example data available for this demo
                </div>
                <div style={{ fontSize: '0.8rem', color: '#475569', marginBottom: '10px' }}>
                  {(() => {
                    const domainPhrases: Record<string, string> = {
                      'Bulk RNA-seq':   'bulk RNA-seq count matrix + sample metadata',
                      'Single-Cell':    'single-cell gene-expression matrix (10x HDF5)',
                      'Sequencing QC':  'paired-end FASTQ files',
                      'Phylogenetics':  'aligned FASTA sequences',
                    };
                    const dataDesc = domainPhrases[scenario.domain] ?? `${scenario.domain} dataset`;
                    return `Example ${dataDesc} hosted on S3 — ready to load and run.`;
                  })()}
                </div>
                <div className="d-flex align-items-center gap-2 flex-wrap">
                  {scenario.dataPreview && (
                    <Button
                      size="sm"
                      variant="outline-secondary"
                      onClick={() => setPreviewTables(scenario.dataPreview!)}
                      style={{ fontSize: '0.8rem', fontWeight: 600 }}
                    >
                      🔍 Preview data
                    </Button>
                  )}
                  <Button
                    size="sm"
                    variant="success"
                    disabled={loading}
                    onClick={() => executeCommand(scenario.followUpPrompt!, item.scenarioId ?? scenario.id, true)}
                    style={{ fontSize: '0.8rem', fontWeight: 600 }}
                  >
                    {loading ? 'Running…' : '▶ Load & run'}
                  </Button>
                </div>
              </div>
            </div>
          )}
        </div>
      );
    }

    // needs_inputs only (no workflow plan, or scenario forces needs_inputs): show Load & run only.
    if (isNeedsInputs) {
      const hasStructuredInputs =
        Array.isArray((output?.data as any)?.results?.required_inputs) &&
        ((output?.data as any)?.results?.required_inputs as unknown[]).length > 0;

      if (scenario?.followUpPrompt) {
        // Demo scenario: show NeedsInputsCard + "Load & run" shortcut
        return (
          <div>
            <NeedsInputsCard
              output={output as Record<string, unknown>}
            />
            <div
              className="mt-3 p-3 rounded-3 d-flex align-items-start gap-3"
              style={{
                background: 'linear-gradient(135deg, #EFF6FF, #F0FDF4)',
                border: '1px solid #BFDBFE',
              }}
            >
              <span style={{ fontSize: '1.4rem', lineHeight: 1, flexShrink: 0 }}>📋</span>
              <div className="flex-grow-1">
                <div style={{ fontWeight: 600, fontSize: '0.88rem', color: '#1E40AF', marginBottom: '4px' }}>
                  Example data available for this demo
                </div>
                <div style={{ fontSize: '0.8rem', color: '#475569', marginBottom: '10px' }}>
                  {(() => {
                    const domainPhrases: Record<string, string> = {
                      'Bulk RNA-seq':   'bulk RNA-seq count matrix + sample metadata',
                      'Single-Cell':    'single-cell gene-expression matrix (10x HDF5)',
                      'Sequencing QC':  'paired-end FASTQ files',
                      'Phylogenetics':  'aligned FASTA sequences',
                    };
                    const dataDesc = domainPhrases[scenario.domain] ?? `${scenario.domain} dataset`;
                    return `Example ${dataDesc} hosted on S3 — ready to load and run.`;
                  })()}
                </div>
                <div className="d-flex align-items-center gap-2 flex-wrap">
                  {scenario.dataPreview && (
                    <Button
                      size="sm"
                      variant="outline-secondary"
                      onClick={() => setPreviewTables(scenario.dataPreview!)}
                      style={{ fontSize: '0.8rem', fontWeight: 600 }}
                    >
                      🔍 Preview data
                    </Button>
                  )}
                  <Button
                    size="sm"
                    variant="success"
                    disabled={loading}
                    onClick={() => executeCommand(scenario.followUpPrompt!, item.scenarioId ?? scenario.id, true)}
                    style={{ fontSize: '0.8rem', fontWeight: 600 }}
                  >
                    {loading ? 'Running…' : '▶ Load & run'}
                  </Button>
                </div>
              </div>
            </div>
          </div>
        );
      }

      // Non-demo path: show NeedsInputsCard with upload CTA when structured data present,
      // otherwise fall through to plain markdown.
      if (hasStructuredInputs) {
        return (
          <NeedsInputsCard
            output={output as Record<string, unknown>}
            onUpload={() => fileInputRef.current?.click()}
          />
        );
      }
    }

    return renderedResponse;
  };


  const renderWorkflowContextSection = () => {
    const hasContext = Boolean(
      workflowContext.alignedSequences ||
      workflowContext.selectedSequences ||
      workflowContext.mutatedSequences ||
      workflowContext.plasmidData
    );

    if (!hasContext) {
      return null;
    }

    return (
      <section className="alert alert-success mb-4">
        <div className="d-flex justify-content-between align-items-center mb-2">
          <h2 className="h5 mb-0">🔄 Workflow Context Available</h2>
          <Button
            variant="outline-secondary"
            size="sm"
            onClick={() => setWorkflowContext({})}
          >
            Clear Context
          </Button>
        </div>
        <div className="row gy-2">
          {workflowContext.alignedSequences && (
            <div className="col-md-6">
              <strong>📊 Aligned Sequences</strong>
              <div className="text-muted small">
                {workflowContext.alignedSequences.split('\n').length} lines available
              </div>
            </div>
          )}
          {workflowContext.selectedSequences && (
            <div className="col-md-6">
              <strong>🎯 Selected Sequences</strong>
              <div className="text-muted small">
                {workflowContext.selectedSequences.length} sequences available
              </div>
            </div>
          )}
          {workflowContext.mutatedSequences && (
            <div className="col-md-6">
              <strong>🧬 Mutated Sequences</strong>
              <div className="text-muted small">
                {workflowContext.mutatedSequences.length} variants available
              </div>
            </div>
          )}
          {workflowContext.plasmidData && (
            <div className="col-md-6">
              <strong>🧪 Plasmid Data</strong>
              <div className="text-muted small">
                Vector visualization available
              </div>
            </div>
          )}
        </div>
      </section>
    );
  };

  const renderHistorySection = () => {
    if (history.length === 0) {
      return (
        <section className="mt-4">
          <h2 className="h5 mb-3">Conversation</h2>
          <CapabilityGrid onSelectPrompt={(prompt) => { setCommand(prompt); }} />
        </section>
      );
    }

    return (
      <section className="mt-4">
        <h2 className="h5 mb-3">Conversation</h2>
        <div ref={historyTopRef} className="d-flex flex-column gap-4">
          {history.map((item, index) => {
            const timestamp = item.timestamp instanceof Date ? item.timestamp : new Date(item.timestamp);

            // ── Pending placeholder while backend is processing ──────────
            if (item.type === 'pending') {
              return (
                <div key={item.pendingId || index} ref={pendingItemRef} className="d-flex flex-column gap-3">
                  <div
                    className="align-self-end bg-blue-subtle border border-brand-blue rounded-4 shadow-sm px-3 py-2"
                    style={{ maxWidth: '50%', wordBreak: 'break-word', overflowWrap: 'anywhere' }}
                  >
                    <div className="text-muted small text-uppercase fw-semibold mb-1">Prompt</div>
                    <div className="prompt-bubble-text">{item.input}</div>
                    <div className="text-muted small mt-1">{timestamp.toLocaleTimeString()}</div>
                  </div>
                  <Card className="border-0 shadow-sm">
                    <Card.Body>
                      <div style={{ display: 'flex', alignItems: 'center', gap: '12px', padding: '6px 0' }}>
                        <div style={{ display: 'flex', gap: '5px', alignItems: 'center' }}>
                          {[0, 1, 2].map(i => (
                            <div
                              key={i}
                              style={{
                                width: 9, height: 9, borderRadius: '50%',
                                background: '#3B82F6',
                                animation: `helix-bounce 1.2s ease-in-out ${i * 0.2}s infinite`,
                              }}
                            />
                          ))}
                        </div>
                        <span style={{ fontSize: '0.88rem', color: '#64748B', fontStyle: 'italic' }}>
                          {loadingMsg}
                        </span>
                      </div>
                    </Card.Body>
                  </Card>
                </div>
              );
            }

            return (
              <div key={index} className="d-flex flex-column gap-3">
                <div
                  className="align-self-end bg-blue-subtle border border-brand-blue rounded-4 shadow-sm px-3 py-2"
                  style={{ maxWidth: '50%', wordBreak: 'break-word', overflowWrap: 'anywhere' }}
                >
                  <div className="text-muted small text-uppercase fw-semibold mb-1">Prompt</div>
                  <div className="prompt-bubble-text">{item.input}</div>
                  <div className="text-muted small mt-1">{timestamp.toLocaleTimeString()} • {item.type}</div>
                </div>
                <Card className="border-0 shadow-sm">
                  <Card.Body>
                    <div className="text-muted small text-uppercase fw-semibold mb-2">Response</div>
                    {renderOutput(item)}
                    {item.run_id && (
                      <div className="mt-3 d-flex align-items-center gap-2 flex-wrap" style={{ borderTop: '1px solid #e9ecef', paddingTop: '0.6rem' }}>
                        <span
                          className="badge rounded-pill"
                          style={{ background: '#EFF6FF', color: '#1D4ED8', fontFamily: 'monospace', fontSize: '0.72rem', fontWeight: 500 }}
                          title={`Run ID: ${item.run_id}`}
                        >
                          🔬 run:{item.run_id.slice(0, 8)}
                        </span>
                        {item.parent_run_id && (
                          <>
                            <span className="text-muted" style={{ fontSize: '0.72rem' }}>←</span>
                            <span
                              className="badge rounded-pill"
                              style={{ background: '#F0FDF4', color: '#166534', fontFamily: 'monospace', fontSize: '0.72rem', fontWeight: 500 }}
                              title={`Parent Run ID: ${item.parent_run_id}`}
                            >
                              parent:{item.parent_run_id.slice(0, 8)}
                            </span>
                            <button
                              className="btn btn-link btn-sm p-0"
                              style={{ fontSize: '0.72rem', color: '#6B7280' }}
                              onClick={() => executeCommand(`what changed between the runs?`)}
                            >
                              diff ↗
                            </button>
                          </>
                        )}
                      </div>
                    )}
                    {/* Contextual follow-up chips — only on the most recent item */}
                    {index === 0 && (
                      <FollowUpChips
                        tool={item.type}
                        onSelect={(prompt) => executeCommand(prompt)}
                      />
                    )}
                  </Card.Body>
                </Card>
              </div>
            );
          })}
        </div>
      </section>
    );
  };

  const quickExamples = useMemo<QuickExample[]>(() => [
    {
      title: 'Visualize Phylogenetic Tree',
      description: 'Build a tree with sample sequences',
      command: getExampleWithSequences('visualize the phylogenetic tree', 'phylogenetic'),
    },
    {
      title: 'Align Sample Sequences',
      description: 'Perform alignment with sample FASTA data',
      command: getExampleWithSequences('align these sequences', 'threeSequences'),
    },
    {
      title: 'Create 96 Variants',
      description: 'Generate mutations from a sample sequence',
      command: getExampleWithSequences('create 96 variants', 'single'),
    },
    {
      title: 'Research DNA Vendors',
      description: 'Find synthesis vendors for 1000bp sequences',
      command: 'research DNA synthesis vendors for 1000bp sequences',
    },
  ], []);

  const workflowContextContent = renderWorkflowContextSection();
  const historyContent = renderHistorySection();

  const designProps: PromptDesignProps = {
    command,
    onCommandChange: setCommand,
    onSubmit: handleSubmit,
    onSuggestSubmit: (q: string) => void executeCommand(q),
    loading,
    agentLoading,
    onAgentSubmit: handleAgentSubmit,
    placeholder: getContextualPlaceholder(history[0]?.type),
    dragActive,
    uploadedFiles,
    onFileRemove: handleFileRemove,
    onDropZoneDragOver: handleDragOver,
    onDropZoneDragLeave: handleDragLeave,
    onDropZoneDrop: handleDrop,
    onBrowseClick: handleBrowseClick,
    onToggleExamples: handleToggleExamples,
    jobsOpen,
    onToggleJobs: handleToggleJobs,
    jobsCount: jobIds.length,
    workflowContextContent,
    historyContent,
  };

  if (SelectedDesignComponent) {
    return (
      <div className="container-fluid py-4 px-5">
        <input
          type="file"
          ref={fileInputRef}
          style={{ display: 'none' }}
          accept=".fasta,.fa,.fas,.fastq,.fq,.gz,.csv,.tsv,.xlsx,.xls,.txt"
          multiple
          onChange={handleFileInput}
        />

        <div className="d-flex justify-content-between align-items-start flex-wrap mb-4">
          <div className="mb-3">
            <h1 className="mb-1">Helix.AI</h1>
            <div className="d-flex flex-column gap-2">
              <span className={`badge ${serverStatus === 'healthy' ? 'bg-success' : 'bg-danger'}`} style={{ width: '150px' }}>
                Server: {serverStatus}
              </span>
              {sessionId && (
                <span 
                  className="badge" 
                  style={{ width: '150px', backgroundColor: theme.colors.blue, color: theme.colors.white, cursor: 'pointer' }}
                  onClick={handleSessionClick}
                >
                  Session Info
                </span>
              )}
            </div>
          </div>

          <div className="helix-header-actions">
            <Button
              variant="outline-secondary"
              className="prompt-toolbar-button"
              onClick={handleToggleExamples}
              aria-label="Toggle examples"
              style={{
                background: 'linear-gradient(135deg, #3A60A8, #7B3FA8)',
                color: '#FFFFFF',
                border: 'none',
                fontWeight: 700,
                letterSpacing: '0.02em',
              }}
            >
              📚 Examples
            </Button>
            <Button
              variant="outline-secondary"
              className="prompt-toolbar-button"
              onClick={() => setDemoOpen(true)}
              aria-label="Open demo scenarios"
              style={{
                background: 'linear-gradient(135deg, #3A60A8, #7B3FA8)',
                color: '#FFFFFF',
                border: 'none',
                fontWeight: 700,
                letterSpacing: '0.02em',
              }}
            >
              🧬 Demo
            </Button>
            <Button
              variant="outline-secondary"
              className="prompt-toolbar-button"
              onClick={handleToggleJobs}
              aria-label="Toggle jobs"
              style={{
                background: 'linear-gradient(135deg, #3A60A8, #7B3FA8)',
                color: '#FFFFFF',
                border: 'none',
                fontWeight: 700,
                letterSpacing: '0.02em',
              }}
            >
              🧾 Jobs{jobIds.length > 0 ? ` (${jobIds.length})` : ''}
            </Button>
          </div>
        </div>

        <SelectedDesignComponent {...designProps} />
        <JobsPanel show={jobsOpen} onHide={() => setJobsOpen(false)} jobIds={jobIds} />
        <ExamplesPanel show={examplesOpen} onHide={() => setExamplesOpen(false)} onSelect={handleExampleClick} />
        <DemoScenariosPanel show={demoOpen} onHide={() => setDemoOpen(false)} onSelect={handleExampleClick} />

        {/* Data Preview Modal */}
        <Modal show={!!previewTables} onHide={() => setPreviewTables(null)} size="xl" centered>
          <Modal.Header closeButton style={{ background: '#F8FAFC', borderBottom: '1px solid #E2E8F0' }}>
            <Modal.Title style={{ fontSize: '1rem', fontWeight: 700, color: '#1E293B' }}>
              🔍 Example Data Preview
            </Modal.Title>
          </Modal.Header>
          <Modal.Body style={{ background: '#F8FAFC', maxHeight: '70vh', overflowY: 'auto' }}>
            {(previewTables ?? []).map((table, ti) => (
              <div key={ti} className={ti > 0 ? 'mt-4' : ''}>
                <div style={{ fontWeight: 600, fontSize: '0.9rem', color: '#1E293B', marginBottom: '6px' }}>
                  {table.title}
                </div>
                <div style={{ overflowX: 'auto' }}>
                  <table className="table table-sm table-bordered mb-1" style={{ fontSize: '0.8rem', background: '#fff' }}>
                    <thead style={{ background: '#EFF6FF' }}>
                      <tr>
                        {table.headers.map((h, hi) => (
                          <th key={hi} style={{ whiteSpace: 'nowrap', color: '#1E40AF', fontWeight: 600 }}>{h}</th>
                        ))}
                      </tr>
                    </thead>
                    <tbody>
                      {table.rows.map((row, ri) => (
                        <tr key={ri}>
                          {row.map((cell, ci) => (
                            <td key={ci} style={{ whiteSpace: 'nowrap', fontFamily: 'monospace', color: ci === 0 ? '#7C3AED' : '#334155' }}>
                              {cell}
                            </td>
                          ))}
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
                {table.note && (
                  <div style={{ fontSize: '0.75rem', color: '#64748B', marginTop: '4px' }}>{table.note}</div>
                )}
              </div>
            ))}
          </Modal.Body>
          <Modal.Footer style={{ background: '#F8FAFC', borderTop: '1px solid #E2E8F0' }}>
            <Button variant="outline-secondary" size="sm" onClick={() => setPreviewTables(null)}>
              Close
            </Button>
          </Modal.Footer>
        </Modal>

        {/* Session Info Modal */}
        <Modal show={sessionModalOpen} onHide={() => setSessionModalOpen(false)} size="lg">
          <Modal.Header closeButton>
            <Modal.Title>Session Information</Modal.Title>
          </Modal.Header>
          <Modal.Body>
            {sessionInfoLoading ? (
              <div className="text-center py-4">
                <div className="spinner-border" role="status">
                  <span className="visually-hidden">Loading...</span>
                </div>
              </div>
            ) : (
              <pre style={{ 
                backgroundColor: '#f8f9fa', 
                padding: '1rem', 
                borderRadius: '0.25rem',
                overflow: 'auto',
                maxHeight: '70vh',
                fontSize: '0.875rem'
              }}>
                {JSON.stringify(sessionInfo, null, 2)}
              </pre>
            )}
          </Modal.Body>
          <Modal.Footer>
            <Button
              variant="outline-primary"
              onClick={() => {
                handleNewSession();
                setSessionModalOpen(false);
              }}
              aria-label="Start a new session"
            >
              ＋ New Session
            </Button>
            <Button variant="secondary" onClick={() => setSessionModalOpen(false)}>
              Close
            </Button>
          </Modal.Footer>
        </Modal>
      </div>
    );
  }


  return (
    <div className="container-fluid py-4 px-5">
      <div className="row">
        <div className="col-md-9 pe-5">
          {/* Header Row with Tips */}
          <div className="row mb-4">
            <div className="col-md-8">
              <div className="row">
                <div className="col-12">
                  <h1 className="mb-0">Helix.AI</h1>
                </div>
                <div className="col-12 mt-2">
                  <div className="d-flex flex-column gap-2" style={{ width: 'fit-content' }}>
                    <span className={`badge ${serverStatus === 'healthy' ? 'bg-success' : 'bg-danger'}`} style={{ width: '150px' }}>
                      Server: {serverStatus}
                    </span>
                    {sessionId && (
                      <span 
                        className="badge bg-info" 
                        style={{ width: '150px', cursor: 'pointer' }}
                        onClick={handleSessionClick}
                      >
                        Session Info
                      </span>
                    )}
                  </div>
                </div>
              </div>
            </div>
            <div className="col-md-4">
              {/* Tips Section - Spans both rows */}
              <div className="small text-muted p-3 bg-light rounded border h-100">
                <strong className="text-primary">💡 Tips:</strong>
                <ul className="mb-0 mt-1">
                  <li>Upload FASTA files by dragging and dropping</li>
                  <li>Use natural language for complex workflows</li>
                  <li>Combine multiple steps in one command</li>
                  <li>Ask for vendor research and testing options</li>
                </ul>
              </div>
            </div>
          </div>
          
          {/* Workflow Context Display */}
          {(workflowContext.alignedSequences || workflowContext.selectedSequences || workflowContext.mutatedSequences || workflowContext.plasmidData) && (
            <div className="alert alert-success mb-3">
              <h6 className="mb-2">🔄 Workflow Context Available:</h6>
              <div className="row">
                {workflowContext.alignedSequences && (
                  <div className="col-md-6 mb-2">
                    <strong>📊 Aligned Sequences:</strong>
                    <br />
                    <small className="text-muted">
                      {workflowContext.alignedSequences.split('\n').length} lines available
                    </small>
                  </div>
                )}
                {workflowContext.selectedSequences && (
                  <div className="col-md-6 mb-2">
                    <strong>🎯 Selected Sequences:</strong>
                    <br />
                    <small className="text-muted">
                      {workflowContext.selectedSequences.length} sequences available
                    </small>
                  </div>
                )}
                {workflowContext.mutatedSequences && (
                  <div className="col-md-6 mb-2">
                    <strong>🧬 Mutated Sequences:</strong>
                    <br />
                    <small className="text-muted">
                      {workflowContext.mutatedSequences.length} variants available
                    </small>
                  </div>
                )}
                {workflowContext.plasmidData && (
                  <div className="col-md-6 mb-2">
                    <strong>🧪 Plasmid Data:</strong>
                    <br />
                    <small className="text-muted">
                      Vector visualization available
                    </small>
                  </div>
                )}
              </div>
              <button 
                className="btn btn-sm btn-outline-secondary mt-2"
                onClick={() => setWorkflowContext({})}
              >
                🗑️ Clear Context
              </button>
            </div>
          )}

          {/* Uploaded Files Display */}
          {uploadedFiles.length > 0 && (
            <div className="alert alert-info mb-3">
              <div className="d-flex justify-content-between align-items-center mb-2">
                <strong>📁 Uploaded Files ({uploadedFiles.length})</strong>
                <button
                  className="btn btn-sm btn-outline-secondary"
                  onClick={() => setUploadedFiles([])}
                >
                  Clear All
                </button>
              </div>
              {uploadedFiles.map((file, index) => (
                <div key={index} className="mb-2 pb-2 border-bottom">
                  <div className="d-flex justify-content-between align-items-center">
                    <div>
                      <strong>{file.name}</strong>
                      <br />
                      <small className="text-muted">
                        {(file.size / (1024 * 1024)).toFixed(2)} MB
                        {file.status === 'uploading' && ' • uploading…'}
                        {file.status === 'uploaded' && ' • ready'}
                        {file.status === 'failed' && ` • ${file.error ?? 'failed'}`}
                      </small>
                    </div>
                    <button
                      className="btn btn-sm btn-outline-secondary"
                      onClick={() => handleFileRemove(index)}
                    >
                      ✕
                    </button>
                  </div>
                  {file.schema_preview && file.status === 'uploaded' && (
                    <SchemaPreviewPanel
                      filename={file.name}
                      preview={file.schema_preview}
                      onSuggest={(q) => setCommand(q)}
                    />
                  )}
                </div>
              ))}
            </div>
          )}
          
          {/* Command Input with Drag-and-Drop */}
      <div
        className={`mb-4${dragActive ? ' border border-primary border-3' : ''}`}
        onDragOver={handleDragOver}
        onDragLeave={handleDragLeave}
        onDrop={handleDrop}
        style={{ position: 'relative', background: dragActive ? '#e6f0ff' : undefined }}
      >
        <textarea
          className="form-control mb-3"
          value={command}
          onChange={(e) => setCommand(e.target.value)}
          rows={command.split('\n').length < 6 ? 6 : command.split('\n').length}
          style={{ 
            resize: 'vertical', 
            minHeight: 120,
            fontSize: '1rem',
            fontFamily: 'monospace',
            lineHeight: '1.4'
          }}
          onKeyDown={(e) => {
            if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') {
              e.preventDefault();
              handleSubmit();
            }
          }}
          placeholder={getContextualPlaceholder(history[0]?.type)}
        />
        <div className="text-center">
          <button 
            className="btn btn-primary btn-lg" 
            onClick={handleSubmit}
            disabled={loading}
            style={{ 
              padding: '15px 40px', 
              fontSize: '1.2rem',
              fontWeight: '600',
              borderRadius: '8px',
              boxShadow: '0 4px 8px rgba(0,123,255,0.3)'
            }}
          >
            {loading ? (
              <>
                <span className="spinner-border spinner-border-sm me-2" role="status" aria-hidden="true"></span>
                Processing...
              </>
            ) : (
              'Submit'
            )}
          </button>
        </div>
        {dragActive && (
          <div
            style={{
              position: 'absolute',
              top: 0,
              left: 0,
              width: '100%',
              height: '100%',
              background: 'rgba(0,123,255,0.1)',
              zIndex: 10,
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              color: '#007bff',
              fontWeight: 'bold',
              fontSize: '1.2rem',
              pointerEvents: 'none',
            }}
          >
            Drop file to upload
          </div>
        )}
      </div>

      {/* History */}
          {history.map((item, index) => (
        <div className="mt-4" key={index}>
          <div className="mb-2">
                <strong>Command:</strong> <span style={{ whiteSpace: 'pre-wrap', wordBreak: 'break-word' }}>{item.input}</span>
                <small className="text-muted ms-2">
                  ({item.type}) - {item.timestamp.toLocaleTimeString()}
                </small>
              </div>
              {renderOutput(item)}
              {item.run_id && (
                <div className="mt-2 d-flex align-items-center gap-2 flex-wrap">
                  <span
                    className="badge rounded-pill"
                    style={{ background: '#EFF6FF', color: '#1D4ED8', fontFamily: 'monospace', fontSize: '0.72rem', fontWeight: 500 }}
                    title={`Run ID: ${item.run_id}`}
                  >
                    🔬 run:{item.run_id.slice(0, 8)}
                  </span>
                  {item.parent_run_id && (
                    <>
                      <span className="text-muted" style={{ fontSize: '0.72rem' }}>←</span>
                      <span
                        className="badge rounded-pill"
                        style={{ background: '#F0FDF4', color: '#166534', fontFamily: 'monospace', fontSize: '0.72rem', fontWeight: 500 }}
                        title={`Parent Run ID: ${item.parent_run_id}`}
                      >
                        parent:{item.parent_run_id.slice(0, 8)}
                      </span>
                    </>
                  )}
                </div>
              )}
            </div>
          ))}

          </div>

        {/* Sidebar with Available Tools */}
        <div className="col-md-3 ps-4">
          {/* Example Commands */}
          <div className="card mb-3">
            <div className="card-header bg-primary text-white">
              <h6 className="mb-0">💡 Example Commands</h6>
            </div>
            <div className="card-body p-3" style={{ maxHeight: '400px', overflowY: 'auto' }}>
              {commandMode === 'natural' ? (
                <>
                  <div className="mb-2">
                    <strong className="text-primary">🧬 Sequence Analysis:</strong>
                    <div className="small text-muted mb-1 cursor-pointer p-1 rounded" onClick={() => handleExampleClick("align these sequences: >seq1 ATGCGATCG >seq2 ATGCGATC")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "align these sequences: {'>'}seq1 ATGCGATCG {'>'}seq2 ATGCGATC"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer p-1 rounded" onClick={() => handleExampleClick("perform multiple sequence alignment on the uploaded sequences")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "perform multiple sequence alignment on the uploaded sequences"
                </div>
                    <div className="small text-muted cursor-pointer p-1 rounded" onClick={() => handleExampleClick("show me the alignment of these DNA sequences")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "show me the alignment of these DNA sequences"
            </div>
          </div>
          
                  <div className="mb-2">
                    <strong className="text-primary">🎯 Sequence Selection:</strong>
                    <div className="small text-muted mb-1 cursor-pointer p-1 rounded" onClick={() => handleExampleClick("from the sequence variants, pick 10 sequences randomly")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "from the sequence variants, pick 10 sequences randomly"
            </div>
                    <div className="small text-muted mb-1 cursor-pointer p-1 rounded" onClick={() => handleExampleClick("select 5 sequences with the highest mutation rate")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "select 5 sequences with the highest mutation rate"
                  </div>
                    <div className="small text-muted cursor-pointer p-1 rounded" onClick={() => handleExampleClick("choose the most diverse sequences from the alignment")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "choose the most diverse sequences from the alignment"
                  </div>
                  </div>
                  
                  <div className="mb-2">
                    <strong className="text-primary">🧪 Mutation Generation:</strong>
                    <div className="small text-muted mb-1 cursor-pointer p-1 rounded" onClick={() => handleExampleClick("mutate the sequence ATGCGATCG to create 96 variants")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "mutate the sequence ATGCGATCG to create 96 variants"
                  </div>
                    <div className="small text-muted mb-1 cursor-pointer p-1 rounded" onClick={() => handleExampleClick("generate variants of this DNA sequence")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "generate variants of this DNA sequence"
                    </div>
                    <div className="small text-muted cursor-pointer p-1 rounded" onClick={() => handleExampleClick("create mutations with 0.2 mutation rate")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "create mutations with 0.2 mutation rate"
                    </div>
                  </div>
                  
                  <div className="mb-2">
                    <strong className="text-primary">🔬 Data Analysis:</strong>
                    <div className="small text-muted mb-1 cursor-pointer p-1 rounded" onClick={() => handleExampleClick("analyze the alignment and show me the most conserved regions")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "analyze the alignment and show me the most conserved regions"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer p-1 rounded" onClick={() => handleExampleClick("create visualizations of the sequence data")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "create visualizations of the sequence data"
                    </div>
                    <div className="small text-muted cursor-pointer p-1 rounded" onClick={() => handleExampleClick("show me statistics for these sequences")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "show me statistics for these sequences"
                    </div>
                  </div>
                  
                  <div className="mb-2">
                    <strong className="text-primary">🧬 DNA Synthesis & Testing:</strong>
                    <div className="small text-muted mb-1 cursor-pointer p-1 rounded" onClick={() => handleExampleClick("I want to order these sequences from a DNA synthesis vendor")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "I want to order these sequences from a DNA synthesis vendor"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer p-1 rounded" onClick={() => handleExampleClick("find DNA synthesis companies for my sequences")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "find DNA synthesis companies for my sequences"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer p-1 rounded" onClick={() => handleExampleClick("what testing options are available for my sequences?")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "what testing options are available for my sequences?"
                    </div>
                    <div className="small text-muted cursor-pointer p-1 rounded" onClick={() => handleExampleClick("research vendors for gene synthesis and validation")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "research vendors for gene synthesis and validation"
                    </div>
                  </div>
                  
                  <div className="mb-2">
                    <strong className="text-primary">🔄 Multi-step Workflows:</strong>
                    <div className="small text-muted mb-1 cursor-pointer p-1 rounded" onClick={() => handleExampleClick("mutate this sequence, then align the variants and pick the best ones")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "mutate this sequence, then align the variants and pick the best ones"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer p-1 rounded" onClick={() => handleExampleClick("analyze these sequences and then find vendors to synthesize them")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "analyze these sequences and then find vendors to synthesize them"
                    </div>
                    <div className="small text-muted cursor-pointer p-1 rounded" onClick={() => handleExampleClick("align sequences and find synthesis vendors with testing options")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "align sequences and find synthesis vendors with testing options"
                    </div>
                  </div>
                  
                  <div className="mb-2">
                    <strong className="text-primary">🔬 Plasmid Visualization:</strong>
                    <div className="small text-muted mb-1 cursor-pointer p-1 rounded" onClick={() => handleExampleClick("Visualize plasmid pUC19 with EcoRI site and insert ATGCGATCG")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "Visualize plasmid pUC19 with EcoRI site and insert ATGCGATCG"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer p-1 rounded" onClick={() => handleExampleClick("Create plasmid visualization for pBR322 with BamHI site")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "Create plasmid visualization for pBR322 with BamHI site"
                    </div>
                    <div className="small text-muted cursor-pointer p-1 rounded" onClick={() => handleExampleClick("Show plasmid map with cloning sites and insert sequence")} style={{cursor: 'pointer', fontSize: '0.8rem'}}>
                      "Show plasmid map with cloning sites and insert sequence"
                    </div>
                  </div>
                </>
              ) : (
                <>
                  <div className="mb-3">
                    <strong>🧬 Sequence Alignment:</strong>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("align sequences ACTGTTGAC ACTGCATCC")} style={{cursor: 'pointer'}}>
                      "align sequences ACTGTTGAC ACTGCATCC"
                  </div>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("align with clustal algorithm")} style={{cursor: 'pointer'}}>
                      "align with clustal algorithm"
                  </div>
                    <div className="small text-muted cursor-pointer" onClick={() => handleExampleClick("multiple sequence alignment")} style={{cursor: 'pointer'}}>
                      "multiple sequence alignment"
                  </div>
                  </div>
                  
                  <div className="mb-3">
                    <strong>🧪 Mutation Analysis:</strong>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("mutate sequence ACTGTTGAC with 96 variants")} style={{cursor: 'pointer'}}>
                      "mutate sequence ACTGTTGAC with 96 variants"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("generate 10 variants of ACTGTTGAC")} style={{cursor: 'pointer'}}>
                      "generate 10 variants of ACTGTTGAC"
                    </div>
                    <div className="small text-muted cursor-pointer" onClick={() => handleExampleClick("create mutations with 0.1 mutation rate")} style={{cursor: 'pointer'}}>
                      "create mutations with 0.1 mutation rate"
                    </div>
                  </div>
                  
                  <div className="mb-3">
                    <strong>📊 Data Analysis:</strong>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("analyze sequence data for phylogeny")} style={{cursor: 'pointer'}}>
                      "analyze sequence data for phylogeny"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("sequence composition analysis")} style={{cursor: 'pointer'}}>
                      "sequence composition analysis"
                    </div>
                    <div className="small text-muted cursor-pointer" onClick={() => handleExampleClick("statistical analysis of alignment")} style={{cursor: 'pointer'}}>
                      "statistical analysis of alignment"
                    </div>
                  </div>
                  
                  <div className="mb-3">
                    <strong>📈 Visualization:</strong>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("visualize alignment in PNG format")} style={{cursor: 'pointer'}}>
                      "visualize alignment in PNG format"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("create phylogenetic tree plot")} style={{cursor: 'pointer'}}>
                      "create phylogenetic tree plot"
                    </div>
                    <div className="small text-muted cursor-pointer" onClick={() => handleExampleClick("generate sequence composition chart")} style={{cursor: 'pointer'}}>
                      "generate sequence composition chart"
                    </div>
                  </div>
                  
                  <div className="mb-3">
                    <strong>🔬 Plasmid Visualization:</strong>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("insert sequence ATGGTGCACCTGACTGATGCTGAGAAGTCTGCGGTACTGCCTGCTGGGGGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA into pUC19 vector")} style={{cursor: 'pointer'}}>
                      "Insert GFP gene into pUC19 vector"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("express sequence ATGGTGCACCTGACTGATGCTGAGAAGTCTGCGGTACTGCCTGCTGGGGGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA in pUC19 vector")} style={{cursor: 'pointer'}}>
                      "Express beta-globin gene in pUC19"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("insert representatives into pUC19 vector")} style={{cursor: 'pointer'}}>
                      "Insert representative sequences into pUC19 vector"
                    </div>
                  </div>
                </>
              )}
          </div>
        </div>
          
        </div>
      </div>
      
      {/* Session Info Modal */}
      <Modal show={sessionModalOpen} onHide={() => setSessionModalOpen(false)} size="lg">
        <Modal.Header closeButton>
          <Modal.Title>Session Information</Modal.Title>
        </Modal.Header>
        <Modal.Body>
          {sessionInfoLoading ? (
            <div className="text-center py-4">
              <div className="spinner-border" role="status">
                <span className="visually-hidden">Loading...</span>
              </div>
            </div>
          ) : (
            <pre style={{ 
              backgroundColor: '#f8f9fa', 
              padding: '1rem', 
              borderRadius: '0.25rem',
              overflow: 'auto',
              maxHeight: '70vh',
              fontSize: '0.875rem'
            }}>
              {JSON.stringify(sessionInfo, null, 2)}
            </pre>
          )}
        </Modal.Body>
        <Modal.Footer>
          <Button
            variant="outline-primary"
            onClick={() => {
              handleNewSession();
              setSessionModalOpen(false);
            }}
            aria-label="Start a new session"
          >
            ＋ New Session
          </Button>
          <Button variant="secondary" onClick={() => setSessionModalOpen(false)}>
            Close
          </Button>
        </Modal.Footer>
      </Modal>

    </div>
  );
}

export default App;
