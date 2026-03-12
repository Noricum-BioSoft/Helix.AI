import React, { useState, useEffect, useMemo, useRef } from 'react';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import { MSAView } from 'react-msaview';
import Plot from 'react-plotly.js';
import { API_BASE_URL, mcpApi } from './services/mcpApi';
import { CommandParser, ParsedCommand } from './utils/commandParser';
import { PlasmidDataVisualizer, PlasmidRepresentativesVisualizer } from './components/PlasmidVisualizer';
import { PhylogeneticTree } from './components/PhylogeneticTree';

import 'bootstrap/dist/css/bootstrap.min.css';
import './theme.css';
import Button from 'react-bootstrap/Button';
import Card from 'react-bootstrap/Card';
import Modal from 'react-bootstrap/Modal';
import { DesignOptionOne, DesignOptionTwo, DesignOptionThree } from './components/designs';
import type { PromptDesignProps, QuickExample } from './components/designs';
import { getExampleWithSequences, sampleSequences } from './utils/sampleSequences';
import { theme } from './theme';
import { JobsPanel } from './components/JobsPanel';
import { ExamplesPanel } from './components/ExamplesPanel';
import { DemoScenariosPanel } from './components/DemoScenariosPanel';
import { getDemoScenarioById, getDemoScenarioByTool, getDemoScenarioByCommandAndTool, DataPreviewTable } from './data/demoScenarios';
import { ThinkingIndicator, ActivityIndicator } from './components/ThinkingIndicator';

const SESSION_STORAGE_KEY = 'helix_session_id';

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
}

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
  const [serverStatus, setServerStatus] = useState<string>('unknown');
  const [dragActive, setDragActive] = useState(false);
  const [sessionId, setSessionId] = useState<string | null>(null);
  const [commandMode, setCommandMode] = useState<'structured' | 'natural'>('natural');
  const [activeTab, setActiveTab] = useState<string>('command');
  const [uploadedFiles, setUploadedFiles] = useState<Array<{ name: string; content: string }>>([]);
  const [workflowContext, setWorkflowContext] = useState<WorkflowContext>({});
  const [isInitialized, setIsInitialized] = useState(false);
  // Track which pipeline commands have already been submitted so the Execute button
  // on the original plan card turns into "Submitted ✓" after the user clicks it.
  const [executedPipelineCommands, setExecutedPipelineCommands] = useState<Set<string>>(new Set());
  const historyTopRef = useRef<HTMLDivElement>(null);
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
      const info = await mcpApi.getSessionInfo(sessionId);
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
      const res = await mcpApi.createSession() as { session_id: string };
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
            const info = await mcpApi.getSessionInfo(stored) as { session?: { history?: unknown[] } };
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
          const res = await mcpApi.createSession() as { session_id: string };
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
      const health = await mcpApi.healthCheck();
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
      
      // If there are uploaded files, append their content to the command
      if (uploadedFiles.length > 0) {
        uploadedFiles.forEach((file) => {
          // Detect R1/R2 pattern for paired-end FASTQ files
          if (file.name.includes('R1') || file.name.includes('_1') || file.name.endsWith('_1.fastq') || file.name.includes('_R1')) {
            finalCommand = `${finalCommand}\n\nForward reads (${file.name}):\n${file.content}`;
          } else if (file.name.includes('R2') || file.name.includes('_2') || file.name.endsWith('_2.fastq') || file.name.includes('_R2')) {
            finalCommand = `${finalCommand}\n\nReverse reads (${file.name}):\n${file.content}`;
          } else {
            finalCommand = `${finalCommand}\n\nFile content (${file.name}):\n${file.content}`;
          }
        });
      }
      
      // ALWAYS use the agent endpoint - never bypass it
      // The agent will handle routing, tool selection, and execution
      console.log('Calling agent via /execute endpoint...');
      response = await mcpApi.executeCommand(finalCommand, sessionId || undefined);
      console.log('Agent response:', response);
      
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
      setHistory(prev => [historyItem, ...prev]);
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
      
      setHistory(prev => [historyItem, ...prev]);
    } finally {
      if (clearInputsAfter) {
        setCommand('');
        setUploadedFiles([]); // Clear uploaded files after processing
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

      if (uploadedFiles.length > 0) {
        uploadedFiles.forEach((file) => {
          if (file.name.includes('R1') || file.name.includes('_1') || file.name.endsWith('_1.fastq') || file.name.includes('_R1')) {
            finalCommand = `${finalCommand}\n\nForward reads (${file.name}):\n${file.content}`;
          } else if (file.name.includes('R2') || file.name.includes('_2') || file.name.endsWith('_2.fastq') || file.name.includes('_R2')) {
            finalCommand = `${finalCommand}\n\nReverse reads (${file.name}):\n${file.content}`;
          } else {
            finalCommand = `${finalCommand}\n\nFile content (${file.name}):\n${file.content}`;
          }
        });
      }

      // ALWAYS use /execute endpoint - the agent handles everything
      // If files are needed, they should be included in the command text or handled by the agent
      const response = await mcpApi.executeCommand(finalCommand, sessionId || undefined);

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
      setUploadedFiles([]);
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
      const response = await mcpApi.executePipelinePlan(originalCommand, sessionId || undefined);
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
      const newFiles: Array<{ name: string; content: string }> = [];
      let filesRead = 0;
      
      files.forEach((file) => {
        const reader = new FileReader();
        reader.onload = (event) => {
          const content = event.target?.result as string;
          newFiles.push({ name: file.name, content });
          filesRead++;
          
          // When all files are read, update state
          if (filesRead === files.length) {
            setUploadedFiles(prev => [...prev, ...newFiles]);
          }
        };
        reader.readAsText(file);
      });
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
      const newFiles: Array<{ name: string; content: string }> = [];
      let filesRead = 0;
      
      files.forEach((file) => {
        const reader = new FileReader();
        reader.onload = (event) => {
          const content = event.target?.result as string;
          newFiles.push({ name: file.name, content });
          filesRead++;
          
          // When all files are read, update state
          if (filesRead === files.length) {
            setUploadedFiles(prev => [...prev, ...newFiles]);
          }
        };
        reader.readAsText(file);
      });
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
    // Links come from data.links (build_standard_response) or raw_result.links
    // (patch_and_rerun / legacy paths).
    const _earlyLinksRaw: any[] =
      agentOutput?.data?.links || actualResult?.links || rawResult?.links || [];
    const _earlyLinks: any[] = Array.isArray(_earlyLinksRaw)
      ? _earlyLinksRaw.map((l: any) =>
          l && typeof l === 'object' ? { ...l, url: normalizeAssetUrl(l.url) } : l
        )
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
              const isDownloadable = isScript || isBundle;
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

    if (finalText) {
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

  const renderOutput = (item: HistoryItem) => {
    const { output, type } = item;

    if (type === 'agent_error') {
      const message = output?.error || 'Agent encountered an error while processing the request.';
      return <div className="alert alert-danger mb-0">{message}</div>;
    }

    // Resolve scenario by explicit scenarioId first, then by command/tool match.
    // We do this early so rendering can stay deterministic for demo cards.
    const responseTool =
      output?.tool || output?.tool_name || output?.result?.tool;
    const scenario =
      (item.scenarioId ? getDemoScenarioById(item.scenarioId) : undefined) ??
      (responseTool
        ? getDemoScenarioByCommandAndTool(item.input ?? '', responseTool)
        : undefined);

    // Detect a workflow plan response so we can append the "Execute Pipeline" button
    const isWorkflowPlan =
      output?.execute_ready === true ||
      output?.result?.execute_ready === true ||
      output?.status === 'workflow_planned' ||
      output?.result?.status === 'workflow_planned';

    // Detect a needs_inputs response so we can show the "Use example data" button.
    // Some demos come back as plain agent markdown (status=success) but still clearly ask
    // for count matrix / sample metadata. Handle that case too.
    const explicitNeedsInputs =
      output?.status === 'needs_inputs' ||
      output?.result?.status === 'needs_inputs';

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

    // Some demo scenarios are explicitly authored as "needs_inputs". If backend
    // emits workflow_planned for these, keep the demo UX consistent and prefer
    // the scenario's follow-up data path.
    const scenarioForcesNeedsInputs =
      isWorkflowPlan &&
      !!scenario?.followUpPrompt &&
      scenario?.expectedBehavior === 'needs_inputs';

    const isNeedsInputs = explicitNeedsInputs || looksLikeNeedsInputsFromAgent || scenarioForcesNeedsInputs;

    const renderedResponse = renderAgentResponse(output);

    if (isWorkflowPlan && !scenarioForcesNeedsInputs) {
      const alreadyExecuted = executedPipelineCommands.has(item.input);
      return (
        <div>
          {renderedResponse}
          <div className="mt-3 d-flex align-items-center gap-3">
            {alreadyExecuted ? (
              <Button variant="outline-success" disabled className="px-4">
                ✓ Submitted
              </Button>
            ) : (
              <Button
                variant="success"
                onClick={() => handleExecutePipeline(item.input)}
                disabled={loading}
                className="px-4"
              >
                {loading ? 'Submitting…' : '▶ Execute Pipeline'}
              </Button>
            )}
            <span className="text-muted small">
              {alreadyExecuted
                ? 'Pipeline submitted — see results above.'
                : 'Confirms the plan above and queues all steps for execution.'}
            </span>
          </div>
        </div>
      );
    }

    if (isNeedsInputs) {
      if (scenario?.followUpPrompt) {
        return (
          <div>
            {renderedResponse}
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
          <div className="text-muted">Run a command to see your prompts and responses here.</div>
        </section>
      );
    }

    return (
      <section className="mt-4">
        <h2 className="h5 mb-3">Conversation</h2>
        <div ref={historyTopRef} className="d-flex flex-column gap-4">
          {history.map((item, index) => {
            const timestamp = item.timestamp instanceof Date ? item.timestamp : new Date(item.timestamp);
            return (
              <div key={index} className="d-flex flex-column gap-3">
                <div
                  className="align-self-end bg-blue-subtle border border-brand-blue rounded-4 shadow-sm px-3 py-2"
                  style={{ maxWidth: '50%', wordBreak: 'break-word', overflowWrap: 'anywhere' }}
                >
                  <div className="text-muted small text-uppercase fw-semibold mb-1">Prompt</div>
                  <div style={{ whiteSpace: 'pre-wrap', wordBreak: 'break-word', overflowWrap: 'anywhere' }}>{item.input}</div>
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
    loading,
    agentLoading,
    onAgentSubmit: handleAgentSubmit,
    placeholder: 'Ask anything or upload a FASTA file...',
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
          accept=".fasta,.fa,.fas,.fastq,.csv,.txt"
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
                <div key={index} className="d-flex justify-content-between align-items-center mb-2 pb-2 border-bottom">
                  <div>
                    <strong>{file.name}</strong>
                    <br />
                    <small className="text-muted">
                      Content length: {file.content.length.toLocaleString()} characters
                    </small>
                  </div>
                  <button 
                    className="btn btn-sm btn-outline-secondary"
                    onClick={() => handleFileRemove(index)}
                  >
                    ✕ Remove
                  </button>
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
          placeholder="Ask anything or upload a FASTA file..."
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
