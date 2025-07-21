import React, { useState, useEffect } from 'react';
import { MSAView } from 'react-msaview';
import Plot from 'react-plotly.js';
import { mcpApi } from './services/mcpApi';
import { CommandParser, ParsedCommand } from './utils/commandParser';
import { PlasmidDataVisualizer } from './components/PlasmidVisualizer';
import { PhylogeneticTree } from './components/PhylogeneticTree';

import 'bootstrap/dist/css/bootstrap.min.css';
import Modal from 'react-bootstrap/Modal';
import Button from 'react-bootstrap/Button';
import Nav from 'react-bootstrap/Nav';
import Tab from 'react-bootstrap/Tab';

interface HistoryItem {
  input: string;
  output: any;
  type: string;
  timestamp: Date;
}

interface WorkflowContext {
  alignedSequences?: string;
  selectedSequences?: string[];
  mutatedSequences?: string[];
  plasmidData?: any;
}

function App() {
  const [command, setCommand] = useState('');
  const [history, setHistory] = useState<HistoryItem[]>([]);
  const [loading, setLoading] = useState(false);
  const [availableTools, setAvailableTools] = useState<any[]>([]);
  const [serverStatus, setServerStatus] = useState<string>('unknown');
  const [dragActive, setDragActive] = useState(false);
  const [sessionId, setSessionId] = useState<string | null>(null);
  const [commandMode, setCommandMode] = useState<'structured' | 'natural'>('natural');
  const [activeTab, setActiveTab] = useState<string>('command');
  const [uploadedFile, setUploadedFile] = useState<{ name: string; content: string } | null>(null);
  const [workflowContext, setWorkflowContext] = useState<WorkflowContext>({});
  const [isInitialized, setIsInitialized] = useState(false);

  const handleExampleClick = (exampleCommand: string) => {
    setCommand(exampleCommand);
  };

  // Create a session on mount and check server health
  useEffect(() => {
    const initializeApp = async () => {
      if (isInitialized) return; // Prevent duplicate initialization
      
      try {
        setIsInitialized(true);
        
        // Check server health first
        await checkServerHealth();
        
        // Load available tools
        await loadAvailableTools();
        
        // Create session
        const res = await mcpApi.createSession() as { session_id: string };
        setSessionId(res.session_id);
        console.log('Session created:', res.session_id);
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
      setServerStatus((health as any).status);
    } catch (error) {
      setServerStatus('error');
      console.error('Server health check failed:', error);
    }
  };

  const loadAvailableTools = async () => {
    try {
      const toolsResponse = await mcpApi.listTools();
      setAvailableTools((toolsResponse as any).tools);
    } catch (error) {
      console.error('Failed to load available tools:', error);
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

  const handleSubmit = async () => {
    if (!command.trim()) return;
    
    setLoading(true);
    
    try {
      let response;
      let parsedCommand: ParsedCommand | undefined;
      
      console.log('Command mode:', commandMode);
      console.log('Session ID:', sessionId);
      console.log('Command:', command);
      console.log('Workflow context:', workflowContext);
      
      // Enhance command with workflow context
      let finalCommand = enhanceCommandWithContext(command);
      
      console.log('Original command:', command);
      console.log('Enhanced command:', finalCommand);
      console.log('Workflow context before sending:', workflowContext);
      
      // If there's an uploaded file, append its content to the command
      if (uploadedFile) {
        finalCommand = `${finalCommand}\n\nFile content:\n${uploadedFile.content}`;
      }
      
      if (commandMode === 'natural') {
        // Use natural language command handling
        console.log('Calling handleNaturalCommand...');
        response = await mcpApi.executeCommand(finalCommand, sessionId || undefined);
        console.log('Natural command response:', response);
      } else {
        // Use structured command parsing
        parsedCommand = CommandParser.parseCommand(command);
        
        switch (parsedCommand.type) {
          case 'sequence_alignment':
            response = await mcpApi.sequenceAlignment(parsedCommand.parameters as any, sessionId || undefined);
            break;
            
          case 'mutate_sequence':
            response = await mcpApi.mutateSequence(parsedCommand.parameters as any, sessionId || undefined);
            break;
            
          case 'analyze_sequence_data':
            response = await mcpApi.analyzeSequenceData(parsedCommand.parameters as any, sessionId || undefined);
            break;
            
          case 'visualize_alignment':
            response = await mcpApi.visualizeAlignment(parsedCommand.parameters as any, sessionId || undefined);
            break;
            
          case 'dna_vendor_research':
            response = await mcpApi.executeCommand(command, sessionId || undefined);
            break;
            
          default:
            // Fallback to general command execution
            response = await mcpApi.executeCommand(command, sessionId || undefined);
            break;
        }
      }
      
      // Update session ID if the backend created one automatically
      if ((response as any).session_id && !sessionId) {
        setSessionId((response as any).session_id);
        console.log('Session created automatically by backend:', (response as any).session_id);
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
      const historyItem: HistoryItem = {
        input: command,
        output: response,
        type: commandMode === 'natural' ? 'natural_command' : parsedCommand?.type || 'structured_command',
        timestamp: new Date()
      };
      
      setHistory(prev => [historyItem, ...prev]);
      
    } catch (error) {
      console.error('Error executing command:', error);
      
      // Add error to history
      const historyItem: HistoryItem = {
        input: command,
        output: { error: error instanceof Error ? error.message : 'Unknown error' },
        type: 'error',
        timestamp: new Date()
      };
      
      setHistory(prev => [historyItem, ...prev]);
    } finally {
      setCommand('');
      setUploadedFile(null); // Clear uploaded file after processing
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
      const file = e.dataTransfer.files[0];
      const reader = new FileReader();
      reader.onload = (event) => {
        const content = event.target?.result as string;
        setUploadedFile({ name: file.name, content });
        // Update command placeholder to suggest what to do with the file
        setCommand(`analyze the uploaded file ${file.name}`);
      };
      reader.readAsText(file);
    }
  };



  const renderOutput = (item: HistoryItem) => {
    const { output, type } = item;
    
    // Add debugging
    console.log('renderOutput called with:', { output, type });
    console.log('output structure:', JSON.stringify(output, null, 2));
    console.log('output keys:', Object.keys(output));
    
    // --- NEW: Render phylogenetic tree if present ---
    if (output && output.tree_newick) {
      return (
        <div>
          <PhylogeneticTree newick={output.tree_newick} />
          {output.text && (
            <pre className="bg-light p-3 border rounded mb-3">{output.text}</pre>
          )}
          {output.statistics && (
            <div className="bg-light p-3 border rounded mb-3">
              <h6>Tree Statistics</h6>
              <div className="row">
                {Object.entries(output.statistics).map(([key, value]) => (
                  <div key={key} className="col-md-6 mb-2">
                    <strong>{key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}:</strong> {String(value)}
                  </div>
                ))}
              </div>
            </div>
          )}
        </div>
      );
    }

    // --- NEW: Render vendor research results if present ---
    if (output && output.output && output.output.result && output.output.result.vendors) {
      const vendors = output.output.result.vendors;
      const recommendations = output.output.result.recommendations || [];
      
      return (
        <div>
          {output.text && (
            <pre className="bg-light p-3 border rounded mb-3">{output.text}</pre>
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
        </div>
      );
    }
    
    if (type === 'error') {
      return (
        <div className="alert alert-danger">
          <strong>Error:</strong> {output.error}
        </div>
      );
    }
    
    // Handle responses with messages array directly in output
    if (output.messages && Array.isArray(output.messages)) {
      console.log('Found messages array directly in output with length:', output.messages.length);
      console.log('Messages:', output.messages);
      
      // Find the last tool message (most recent tool result)
      const toolMessages = output.messages.filter(
        (msg: any) => msg.type === 'tool' && msg.content
      );
      console.log('Tool messages found:', toolMessages.length);
      
      if (toolMessages.length > 0) {
        const lastToolMsg = toolMessages[toolMessages.length - 1];
        console.log('Last tool message:', lastToolMsg);
        
        try {
          const toolResult = JSON.parse(lastToolMsg.content);
          console.log('Parsed last tool result:', toolResult);
          
          // Handle sequence selection results
          if (lastToolMsg.name === 'sequence_selection' && toolResult.output && Array.isArray(toolResult.output)) {
            console.log('Found selected sequences array:', toolResult.output);
            return (
              <div>
                {toolResult.text && (
                  <pre className="bg-light p-3 border rounded mb-3">{toolResult.text}</pre>
                )}
                <div className="bg-light p-3 border rounded mb-3">
                  <h5>Selected Sequences</h5>
                  <div className="table-responsive">
                    <table className="table table-sm">
                      <thead>
                        <tr>
                          <th>Name</th>
                          <th>Sequence</th>
                          <th>Length</th>
                          <th>GC Content</th>
                          <th>Gaps</th>
                          <th>Conservation Score</th>
                        </tr>
                      </thead>
                      <tbody>
                        {toolResult.output.map((seq: any, index: number) => (
                          <tr key={index}>
                            <td>{seq.name}</td>
                            <td className="font-monospace">{seq.sequence}</td>
                            <td>{seq.statistics?.length || 'N/A'} bp</td>
                            <td>{seq.statistics?.gc_content || 'N/A'}%</td>
                            <td>{seq.statistics?.gap_count || 0} ({seq.statistics?.gap_percentage || 0}%)</td>
                            <td>{seq.statistics?.conservation_score || 'N/A'}%</td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                </div>
                {toolResult.selection_criteria && (
                  <div className="bg-light p-3 border rounded mb-3">
                    <h6>Selection Criteria</h6>
                    <div className="row">
                      {Object.entries(toolResult.selection_criteria).map(([key, value]) => (
                        <div key={key} className="col-md-6 mb-2">
                          <strong>{key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}:</strong> {String(value)}
                        </div>
                      ))}
                    </div>
                  </div>
                )}
              </div>
            );
          }
          
          // Handle mutation results
          if (lastToolMsg.name === 'mutate_sequence' && toolResult.output && toolResult.output.variants) {
            console.log('Found mutation variants:', toolResult.output.variants);
            return (
              <div>
                {toolResult.text && (
                  <pre className="bg-light p-3 border rounded mb-3">{toolResult.text}</pre>
                )}
                <div className="bg-light p-3 border rounded mb-3">
                  <h5>Mutation Results</h5>
                  <div className="row">
                    <div className="col-md-6">
                      <h6>Statistics</h6>
                      <ul>
                        <li><strong>Total Variants:</strong> {toolResult.output.total_variants}</li>
                        <li><strong>Substitutions:</strong> {toolResult.output.substitution_count}</li>
                        <li><strong>Insertions:</strong> {toolResult.output.insertion_count}</li>
                        <li><strong>Deletions:</strong> {toolResult.output.deletion_count}</li>
                        <li><strong>Mutation Rate:</strong> {(toolResult.output.mutation_rate * 100).toFixed(2)}%</li>
                      </ul>
                    </div>
                    <div className="col-md-6">
                      <h6>First 10 Variants</h6>
                      <div className="table-responsive">
                        <table className="table table-sm">
                          <thead>
                            <tr>
                              <th>#</th>
                              <th>Sequence</th>
                            </tr>
                          </thead>
                          <tbody>
                            {toolResult.output.variants.slice(0, 10).map((variant: string, index: number) => (
                              <tr key={index}>
                                <td>{index + 1}</td>
                                <td className="font-monospace">{variant}</td>
                              </tr>
                            ))}
                          </tbody>
                        </table>
                      </div>
                      {toolResult.output.variants.length > 10 && (
                        <p className="text-muted">... and {toolResult.output.variants.length - 10} more variants</p>
                      )}
                    </div>
                  </div>
                </div>
              </div>
            );
          }
          
          // Handle alignment results
          if (lastToolMsg.name === 'sequence_alignment' && toolResult.output && Array.isArray(toolResult.output)) {
            console.log('Found alignment array:', toolResult.output);
            return (
              <div>
                {toolResult.text && (
                  <pre className="bg-light p-3 border rounded mb-3">{toolResult.text}</pre>
                )}
                <div className="bg-light p-3 border rounded mb-3">
                  <h5>Alignment Result</h5>
                  <div className="table-responsive">
                    <table className="table table-sm">
                      <thead>
                        <tr>
                          <th>Name</th>
                          <th>Sequence</th>
                        </tr>
                      </thead>
                      <tbody>
                        {toolResult.output.map((seq: any, index: number) => (
                          <tr key={index}>
                            <td>{seq.name}</td>
                            <td className="font-monospace">{seq.sequence}</td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                </div>
                {toolResult.statistics && (
                  <div className="bg-light p-3 border rounded mb-3">
                    <h6>Alignment Statistics</h6>
                    <div className="row">
                      {Object.entries(toolResult.statistics).map(([key, value]) => (
                        <div key={key} className="col-md-6 mb-2">
                          <strong>{key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}:</strong> {String(value)}
                        </div>
                      ))}
                    </div>
                  </div>
                )}
              </div>
            );
          }
          
          // Handle plasmid visualization results
          if (lastToolMsg.name === 'plasmid_visualization' && toolResult.plasmid_data) {
            return (
              <div>
                {toolResult.text && (
                  <pre className="bg-light p-3 border rounded mb-3">{toolResult.text}</pre>
                )}
                <div className="bg-light p-3 border rounded mb-3">
                  <h5>Plasmid Visualization</h5>
                  <div className="row">
                    <div className="col-md-6">
                      <h6>Plasmid Details</h6>
                      <ul className="list-unstyled">
                        <li><strong>Name:</strong> {toolResult.plasmid_data.name}</li>
                        <li><strong>Size:</strong> {toolResult.plasmid_data.size} bp</li>
                        <li><strong>Description:</strong> {toolResult.plasmid_data.description}</li>
                        <li><strong>Visualization Type:</strong> {toolResult.visualization_type}</li>
                      </ul>
                    </div>
                    <div className="col-md-6">
                      <h6>Features</h6>
                      {toolResult.plasmid_data.features && toolResult.plasmid_data.features.length > 0 ? (
                        <ul className="list-unstyled">
                          {toolResult.plasmid_data.features.map((feature: any, index: number) => (
                            <li key={index}>
                              <strong>{feature.name}:</strong> {feature.description} 
                              (position {feature.start}-{feature.end})
                            </li>
                          ))}
                        </ul>
                      ) : (
                        <p className="text-muted">No features defined</p>
                      )}
                    </div>
                  </div>
                  {toolResult.metadata && (
                    <div className="mt-3">
                      <h6>Input Parameters</h6>
                      <div className="row">
                        {Object.entries(toolResult.metadata).map(([key, value]) => (
                          <div key={key} className="col-md-6 mb-2">
                            <strong>{key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}:</strong> {String(value)}
                          </div>
                        ))}
                      </div>
                    </div>
                  )}
                </div>
                <div className="my-3">
                  <PlasmidDataVisualizer data={toolResult.plasmid_data} />
                </div>
              </div>
            );
          }
          
          // Fallback: show any tool result as text
          if (toolResult.text) {
            return (
              <div>
                <pre className="bg-light p-3 border rounded">{toolResult.text}</pre>
              </div>
            );
          }
          
        } catch (e) {
          console.log('Error parsing tool result:', e);
        }
      }
    }
    
    // Handle MCP responses
    if (output.success !== undefined) {
      if (!output.success) {
        return (
          <div className="alert alert-danger">
            <strong>MCP Error:</strong> {output.error}
          </div>
        );
      }
      
      const result = output.result;
      console.log('result structure:', JSON.stringify(result, null, 2));
      
      // Handle nested result structure from natural language commands
      const actualResult = result.result || result;
      console.log('actualResult structure:', JSON.stringify(actualResult, null, 2));

      // --- NEW: Render phylogenetic tree if present in actualResult ---
      if (actualResult && actualResult.tree_newick) {
        return (
          <div>
            <PhylogeneticTree newick={actualResult.tree_newick} />
            {actualResult.text && (
              <pre className="bg-light p-3 border rounded mb-3">{actualResult.text}</pre>
            )}
            {actualResult.statistics && (
              <div className="bg-light p-3 border rounded mb-3">
                <h6>Tree Statistics</h6>
                <div className="row">
                  {Object.entries(actualResult.statistics).map(([key, value]) => (
                    <div key={key} className="col-md-6 mb-2">
                      <strong>{key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}:</strong> {String(value)}
                    </div>
                  ))}
                </div>
              </div>
            )}
          </div>
        );
      }
      // --- END NEW ---
      
      // --- NEW: Handle plasmid visualization results directly in actualResult ---
      if (actualResult && actualResult.plasmid_data) {
        return (
          <div>
            {actualResult.text && (
              <pre className="bg-light p-3 border rounded mb-3">{actualResult.text}</pre>
            )}
            <div className="bg-light p-3 border rounded mb-3">
              <h5>Plasmid Visualization</h5>
              <div className="row">
                <div className="col-md-6">
                  <h6>Plasmid Details</h6>
                  <ul className="list-unstyled">
                    <li><strong>Name:</strong> {actualResult.plasmid_data.name}</li>
                    <li><strong>Size:</strong> {actualResult.plasmid_data.size} bp</li>
                    <li><strong>Description:</strong> {actualResult.plasmid_data.description}</li>
                    <li><strong>Visualization Type:</strong> {actualResult.visualization_type}</li>
                  </ul>
                </div>
                <div className="col-md-6">
                  <h6>Features</h6>
                  {actualResult.plasmid_data.features && actualResult.plasmid_data.features.length > 0 ? (
                    <ul className="list-unstyled">
                      {actualResult.plasmid_data.features.map((feature: any, index: number) => (
                        <li key={index}>
                          <strong>{feature.name}:</strong> {feature.description} 
                          (position {feature.start}-{feature.end})
                        </li>
                      ))}
                    </ul>
                  ) : (
                    <p className="text-muted">No features defined</p>
                  )}
                </div>
              </div>
              {actualResult.metadata && (
                <div className="mt-3">
                  <h6>Input Parameters</h6>
                  <div className="row">
                    {Object.entries(actualResult.metadata).map(([key, value]) => (
                      <div key={key} className="col-md-6 mb-2">
                        <strong>{key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}:</strong> {String(value)}
                      </div>
                    ))}
                  </div>
                </div>
              )}
            </div>
            <div className="my-3">
              <PlasmidDataVisualizer data={actualResult.plasmid_data} />
            </div>
          </div>
        );
      }
      // --- END NEW ---

      // --- NEW: Handle alignment in ToolMessage for natural language commands ---
      // Check if messages is in actualResult (nested) or directly in output
      const messages = actualResult.messages || output.messages;
      console.log('Looking for messages in:', { actualResultMessages: actualResult.messages, outputMessages: output.messages, finalMessages: messages });
      
      if (messages && Array.isArray(messages)) {
        console.log('Found messages array with length:', messages.length);
        console.log('Messages:', messages);
        
        // Find the last tool message (most recent tool result)
        const toolMessages = messages.filter(
          (msg: any) => msg.type === 'tool' && msg.content
        );
        console.log('Tool messages found:', toolMessages.length);
        
        if (toolMessages.length > 0) {
          const lastToolMsg = toolMessages[toolMessages.length - 1];
          console.log('Last tool message:', lastToolMsg);
          
          try {
            const toolResult = JSON.parse(lastToolMsg.content);
            console.log('Parsed last tool result:', toolResult);
            
            // Handle sequence selection results
            if (lastToolMsg.name === 'sequence_selection' && toolResult.output && Array.isArray(toolResult.output)) {
              console.log('Found selected sequences array:', toolResult.output);
              return (
                <div>
                  {toolResult.text && (
                    <pre className="bg-light p-3 border rounded mb-3">{toolResult.text}</pre>
                  )}
                  <div className="bg-light p-3 border rounded mb-3">
                    <h5>Selected Sequences</h5>
                    <div className="table-responsive">
                      <table className="table table-sm">
                        <thead>
                          <tr>
                            <th>Name</th>
                            <th>Sequence</th>
                            <th>Length</th>
                            <th>GC Content</th>
                            <th>Gaps</th>
                            <th>Conservation Score</th>
                          </tr>
                        </thead>
                        <tbody>
                          {toolResult.output.map((seq: any, index: number) => (
                            <tr key={index}>
                              <td>{seq.name}</td>
                              <td className="font-monospace">{seq.sequence}</td>
                              <td>{seq.statistics?.length || 'N/A'} bp</td>
                              <td>{seq.statistics?.gc_content || 'N/A'}%</td>
                              <td>{seq.statistics?.gap_count || 0} ({seq.statistics?.gap_percentage || 0}%)</td>
                              <td>{seq.statistics?.conservation_score || 'N/A'}%</td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  </div>
                  {toolResult.selection_criteria && (
                    <div className="bg-light p-3 border rounded mb-3">
                      <h6>Selection Criteria</h6>
                      <div className="row">
                        {Object.entries(toolResult.selection_criteria).map(([key, value]) => (
                          <div key={key} className="col-md-6 mb-2">
                            <strong>{key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}:</strong> {String(value)}
                          </div>
                        ))}
                      </div>
                    </div>
                  )}
                </div>
              );
            }
            
            // Handle mutation results
            if (lastToolMsg.name === 'mutate_sequence' && toolResult.output && toolResult.output.variants) {
              console.log('Found mutation variants:', toolResult.output.variants);
              return (
                <div>
                  {toolResult.text && (
                    <pre className="bg-light p-3 border rounded mb-3">{toolResult.text}</pre>
                  )}
                  <div className="bg-light p-3 border rounded mb-3">
                    <h5>Mutation Results</h5>
                    <div className="row">
                      <div className="col-md-6">
                        <h6>Statistics</h6>
                        <ul>
                          <li><strong>Total Variants:</strong> {toolResult.output.total_variants}</li>
                          <li><strong>Substitutions:</strong> {toolResult.output.substitution_count}</li>
                          <li><strong>Insertions:</strong> {toolResult.output.insertion_count}</li>
                          <li><strong>Deletions:</strong> {toolResult.output.deletion_count}</li>
                          <li><strong>Mutation Rate:</strong> {(toolResult.output.mutation_rate * 100).toFixed(2)}%</li>
                        </ul>
                      </div>
                      <div className="col-md-6">
                        <h6>First 10 Variants</h6>
                        <div className="table-responsive">
                          <table className="table table-sm">
                            <thead>
                              <tr>
                                <th>#</th>
                                <th>Sequence</th>
                              </tr>
                            </thead>
                            <tbody>
                              {toolResult.output.variants.slice(0, 10).map((variant: string, index: number) => (
                                <tr key={index}>
                                  <td>{index + 1}</td>
                                  <td className="font-monospace">{variant}</td>
                                </tr>
                              ))}
                            </tbody>
                          </table>
                        </div>
                        {toolResult.output.variants.length > 10 && (
                          <p className="text-muted">... and {toolResult.output.variants.length - 10} more variants</p>
                        )}
                      </div>
                    </div>
                  </div>
                </div>
              );
            }
            
            // Handle alignment results
            if (lastToolMsg.name === 'sequence_alignment' && toolResult.output && Array.isArray(toolResult.output)) {
              console.log('Found alignment array:', toolResult.output);
              return (
                <div>
                  {toolResult.text && (
                    <pre className="bg-light p-3 border rounded mb-3">{toolResult.text}</pre>
                  )}
                  <div className="bg-light p-3 border rounded mb-3">
                    <h5>Alignment Result</h5>
                    <div className="table-responsive">
                      <table className="table table-sm">
                        <thead>
                          <tr>
                            <th>Name</th>
                            <th>Sequence</th>
                          </tr>
                        </thead>
                        <tbody>
                          {toolResult.output.map((seq: any, index: number) => (
                            <tr key={index}>
                              <td>{seq.name}</td>
                              <td className="font-monospace">{seq.sequence}</td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  </div>
                  {toolResult.statistics && (
                    <div className="bg-light p-3 border rounded mb-3">
                      <h6>Alignment Statistics</h6>
                      <div className="row">
                        {Object.entries(toolResult.statistics).map(([key, value]) => (
                          <div key={key} className="col-md-6 mb-2">
                            <strong>{key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}:</strong> {String(value)}
                          </div>
                        ))}
                      </div>
                    </div>
                  )}
                </div>
              );
            }
            
            // Handle plasmid visualization results
            if (lastToolMsg.name === 'plasmid_visualization' && toolResult.plasmid_data) {
              return (
                <div>
                  {toolResult.text && (
                    <pre className="bg-light p-3 border rounded mb-3">{toolResult.text}</pre>
                  )}
                  <div className="bg-light p-3 border rounded mb-3">
                    <h5>Plasmid Visualization</h5>
                    <div className="row">
                      <div className="col-md-6">
                        <h6>Plasmid Details</h6>
                        <ul className="list-unstyled">
                          <li><strong>Name:</strong> {toolResult.plasmid_data.name}</li>
                          <li><strong>Size:</strong> {toolResult.plasmid_data.size} bp</li>
                          <li><strong>Description:</strong> {toolResult.plasmid_data.description}</li>
                          <li><strong>Visualization Type:</strong> {toolResult.visualization_type}</li>
                        </ul>
                      </div>
                      <div className="col-md-6">
                        <h6>Features</h6>
                        {toolResult.plasmid_data.features && toolResult.plasmid_data.features.length > 0 ? (
                          <ul className="list-unstyled">
                            {toolResult.plasmid_data.features.map((feature: any, index: number) => (
                              <li key={index}>
                                <strong>{feature.name}:</strong> {feature.description} 
                                (position {feature.start}-{feature.end})
                              </li>
                            ))}
                          </ul>
                        ) : (
                          <p className="text-muted">No features defined</p>
                        )}
                      </div>
                    </div>
                    {toolResult.metadata && (
                      <div className="mt-3">
                        <h6>Input Parameters</h6>
                        <div className="row">
                          {Object.entries(toolResult.metadata).map(([key, value]) => (
                            <div key={key} className="col-md-6 mb-2">
                              <strong>{key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}:</strong> {String(value)}
                            </div>
                          ))}
                        </div>
                      </div>
                    )}
                  </div>
                  <div className="my-3">
                    <PlasmidDataVisualizer data={toolResult.plasmid_data} />
                  </div>
                </div>
              );
            }
            
            // Fallback: show any tool result as text
            if (toolResult.text) {
              return (
                <div>
                  <pre className="bg-light p-3 border rounded">{toolResult.text}</pre>
                </div>
              );
            }
            
          } catch (e) {
            console.log('Error parsing tool result:', e);
            // Fallback to default rendering
          }
        } else {
          console.log('No alignment tool message found');
        }
      } else {
        console.log('No messages array found');
      }
      // --- END NEW ---
      
      return (
        <div>
          {actualResult.text && (
            <pre className="bg-light p-3 border rounded mb-3">
              {actualResult.text}
            </pre>
          )}
          
          {actualResult.selected_variants && Array.isArray(actualResult.selected_variants) && (
            <div className="bg-light p-3 border rounded mb-3">
              <h5>Selected Variants ({actualResult.num_variants_selected}/{actualResult.num_variants_requested})</h5>
              <div className="table-responsive">
                <table className="table table-sm">
                  <thead>
                    <tr>
                      <th>Name</th>
                      <th>Sequence</th>
                    </tr>
                  </thead>
                  <tbody>
                    {actualResult.selected_variants.map((variant: any, index: number) => (
                      <tr key={index}>
                        <td>{variant.name}</td>
                        <td className="font-monospace">{variant.sequence}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          )}
          
          {actualResult.alignment && Array.isArray(actualResult.alignment) && (
            <div className="bg-light p-3 border rounded mb-3">
              <h5>Alignment Result</h5>
              <div className="table-responsive">
                <table className="table table-sm">
                  <thead>
                    <tr>
                      <th>Name</th>
                      <th>Sequence</th>
                    </tr>
                  </thead>
                  <tbody>
                    {actualResult.alignment.map((seq: any, index: number) => (
                      <tr key={index}>
                        <td>{seq.name}</td>
                        <td className="font-monospace">{seq.sequence}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          )}
          
          {actualResult.analysis && (
            <div className="bg-light p-3 border rounded mb-3">
              <h5>Analysis Results</h5>
              <div className="row">
                <div className="col-md-6">
                  <h6>Selection Statistics</h6>
                  <ul className="list-unstyled">
                    <li><strong>Criteria:</strong> {actualResult.analysis.criteria_used}</li>
                    <li><strong>Selection Ratio:</strong> {actualResult.analysis.selection_ratio}</li>
                    <li><strong>Total Variants:</strong> {actualResult.analysis.total_variants}</li>
                    <li><strong>Selected Variants:</strong> {actualResult.analysis.selected_variants}</li>
                  </ul>
                </div>
                {actualResult.analysis.length_stats && (
                  <div className="col-md-6">
                    <h6>Length Statistics</h6>
                    <ul className="list-unstyled">
                      <li><strong>Mean Length:</strong> {actualResult.analysis.length_stats.selected_variants.mean.toFixed(2)}</li>
                      <li><strong>Length Range:</strong> {actualResult.analysis.length_stats.selected_variants.min} - {actualResult.analysis.length_stats.selected_variants.max}</li>
                    </ul>
                  </div>
                )}
              </div>
            </div>
          )}
          
          {actualResult.output && (
            <div className="bg-light p-3 border rounded mb-3">
              <h5>Output Data</h5>
              <pre className="mb-0">{JSON.stringify(actualResult.output, null, 2)}</pre>
            </div>
          )}
          
          {actualResult.plot && (
            <Plot
              data={actualResult.plot.data}
              layout={actualResult.plot.layout}
            />
          )}
          
          {/* DNA Vendor Research Results */}
          {(actualResult.result && actualResult.result.vendors) || (actualResult.output && actualResult.output.result && actualResult.output.result.vendors) ? (
            <div className="bg-light p-3 border rounded mb-3">
              <h5>DNA Synthesis Vendors</h5>
              <div className="row">
                {Object.entries((actualResult.result && actualResult.result.vendors) || (actualResult.output && actualResult.output.result && actualResult.output.result.vendors) || {}).map(([vendorId, vendor]: [string, any]) => (
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
          ) : null}
          
          {/* Testing Options Results */}
          {actualResult.result && actualResult.result.testing_options && (
            <div className="bg-light p-3 border rounded mb-3">
              <h5>Testing Options</h5>
              <div className="row">
                {Object.entries(actualResult.result.testing_options).map(([testId, test]: [string, any]) => (
                  <div key={testId} className="col-md-6 mb-3">
                    <div className="card h-100">
                      <div className="card-body">
                        <h6 className="card-title">{test.name}</h6>
                        <p className="card-text">
                          <strong>Services:</strong> {test.services.join(', ')}
                        </p>
                        <p className="card-text">
                          <strong>Vendors:</strong> {test.vendors.join(', ')}
                        </p>
                        <p className="card-text">
                          <strong>Pricing:</strong> {test.pricing_range}
                        </p>
                        <p className="card-text">
                          <strong>Turnaround:</strong> {test.turnaround_time}
                        </p>
                        <p className="card-text">
                          <strong>Description:</strong> {test.description}
                        </p>
                      </div>
                    </div>
                  </div>
                ))}
              </div>
            </div>
          )}
          
          {/* Recommendations */}
          {(actualResult.result && actualResult.result.recommendations) || (actualResult.output && actualResult.output.result && actualResult.output.result.recommendations) ? (
            <div className="bg-light p-3 border rounded mb-3">
              <h5>Recommendations</h5>
              <ul className="list-unstyled">
                {((actualResult.result && actualResult.result.recommendations) || (actualResult.output && actualResult.output.result && actualResult.output.result.recommendations) || []).map((rec: string, index: number) => (
                  <li key={index} className="mb-2">
                    <i className="bi bi-lightbulb text-warning me-2"></i>
                    {rec}
                  </li>
                ))}
              </ul>
            </div>
          ) : null}
          
          {/* --- NEW: Handle plasmid visualization results directly in actualResult --- */}
          {(() => {
            console.log('Checking for plasmid data in actualResult:', actualResult);
            if (actualResult && actualResult.plasmid_data) {
              console.log('Found plasmid visualization data directly in actualResult:', actualResult.plasmid_data);
              console.log('About to render PlasmidDataVisualizer component');
              console.log('DEBUG: ActualResult structure:', JSON.stringify(actualResult, null, 2));
              return true;
            }
            return false;
          })() && (
            <div>
              {actualResult.text && (
                <pre className="bg-light p-3 border rounded mb-3">{actualResult.text}</pre>
              )}
              <div className="bg-light p-3 border rounded mb-3">
                <h5>Plasmid Visualization</h5>
                <div className="row">
                  <div className="col-md-6">
                    <h6>Plasmid Details</h6>
                    <ul className="list-unstyled">
                      <li><strong>Name:</strong> {actualResult.plasmid_data.name}</li>
                      <li><strong>Size:</strong> {actualResult.plasmid_data.size} bp</li>
                      <li><strong>Description:</strong> {actualResult.plasmid_data.description}</li>
                      <li><strong>Visualization Type:</strong> {actualResult.visualization_type}</li>
                    </ul>
                  </div>
                  <div className="col-md-6">
                    <h6>Features</h6>
                    {actualResult.plasmid_data.features && actualResult.plasmid_data.features.length > 0 ? (
                      <ul className="list-unstyled">
                        {actualResult.plasmid_data.features.map((feature: any, index: number) => (
                          <li key={index}>
                            <strong>{feature.name}:</strong> {feature.description} 
                            (position {feature.start}-{feature.end})
                          </li>
                        ))}
                      </ul>
                    ) : (
                      <p className="text-muted">No features defined</p>
                    )}
                  </div>
                </div>
                {actualResult.metadata && (
                  <div className="mt-3">
                    <h6>Input Parameters</h6>
                    <div className="row">
                      {Object.entries(actualResult.metadata).map(([key, value]) => (
                        <div key={key} className="col-md-6 mb-2">
                          <strong>{key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}:</strong> {String(value)}
                        </div>
                      ))}
                    </div>
                  </div>
                )}
              </div>
              <div className="my-3">
                <PlasmidDataVisualizer data={actualResult.plasmid_data} />
              </div>
            </div>
          )}
          {/* --- END NEW --- */}
        </div>
      );
    }
    
    // Handle legacy responses
    if (output.text) {
      return (
        <pre className="bg-light p-3 border rounded mb-3">
          {output.text}
        </pre>
      );
    }
    
    if (output.plot) {
      return (
        <Plot
          data={output.plot.data}
          layout={output.plot.layout}
        />
      );
    }

    if (output.tree_newick) {
      console.log('🔍 Detected tree_newick in output:', output.tree_newick);
      return (
        <div>
          {output.text && (
            <pre className="bg-light p-3 border rounded mb-3">{output.text}</pre>
          )}
          <PhylogeneticTree newick={output.tree_newick} />
          {output.statistics && (
            <div className="bg-light p-3 border rounded mb-3">
              <h6>Tree Statistics</h6>
              <div className="row">
                {Object.entries(output.statistics).map(([key, value]) => (
                  <div key={key} className="col-md-6 mb-2">
                    <strong>{key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}:</strong> {String(value)}
                  </div>
                ))}
              </div>
            </div>
          )}
        </div>
      );
    }
    
    // Intelligent fallback for unknown response types
    const formatResponse = (data: any): React.ReactNode => {
      // If it's a simple string or has a message
      if (typeof data === 'string') {
    return (
          <div className="bg-light p-3 border rounded mb-3">
            <p className="mb-0">{data}</p>
          </div>
        );
      }
      
      // If it has a status and result
      if (data.status && data.result) {
        return (
          <div className="bg-light p-3 border rounded mb-3">
            <div className={`alert ${data.status === 'success' ? 'alert-success' : 'alert-warning'} mb-3`}>
              <strong>Status:</strong> {data.status}
            </div>
            {typeof data.result === 'string' ? (
              <p className="mb-0">{data.result}</p>
            ) : (
              <div>
                <h6>Results:</h6>
                <pre className="small">{JSON.stringify(data.result, null, 2)}</pre>
              </div>
            )}
          </div>
        );
      }
      
      // If it has a text field
      if (data.text) {
        return (
          <div className="bg-light p-3 border rounded mb-3">
            <p className="mb-0">{data.text}</p>
          </div>
        );
      }
      
      // If it has statistics
      if (data.statistics) {
        return (
          <div className="bg-light p-3 border rounded mb-3">
            <h6>Analysis Results</h6>
            <div className="row">
              {Object.entries(data.statistics).map(([key, value]) => (
                <div key={key} className="col-md-6 mb-2">
                  <strong>{key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}:</strong> {String(value)}
                </div>
              ))}
            </div>
            {data.variants && (
              <div className="mt-3">
                <h6>Generated Variants:</h6>
                <div className="small">
                  {data.variants.slice(0, 5).map((variant: string, index: number) => (
                    <div key={index} className="mb-1">
                      <code>{variant}</code>
                    </div>
                  ))}
                  {data.variants.length > 5 && (
                    <p className="text-muted">... and {data.variants.length - 5} more variants</p>
                  )}
                </div>
              </div>
            )}
          </div>
        );
      }
      
      // If it has sequences
      if (data.sequences) {
        return (
          <div className="bg-light p-3 border rounded mb-3">
            <h6>Sequences</h6>
            <div className="small">
              {data.sequences.map((seq: any, index: number) => (
                <div key={index} className="mb-2">
                  <strong>{seq.name || `Sequence ${index + 1}`}:</strong>
                  <br />
                  <code>{seq.sequence}</code>
                </div>
              ))}
            </div>
          </div>
        );
      }
      
      // If it has vendors
      if (data.vendors) {
        return (
          <div className="bg-light p-3 border rounded mb-3">
            <h6>DNA Synthesis Vendors</h6>
            <div className="row">
              {Object.entries(data.vendors).map(([vendorId, vendor]: [string, any]) => (
                <div key={vendorId} className="col-md-6 mb-3">
                  <div className="card h-100">
                    <div className="card-body">
                      <h6 className="card-title">{vendor.name}</h6>
                      <p className="card-text">
                        <strong>Services:</strong> {vendor.services?.join(', ') || 'N/A'}
                      </p>
                      <p className="card-text">
                        <strong>Pricing:</strong> {vendor.pricing_range || 'N/A'}
                      </p>
                      <p className="card-text">
                        <strong>Turnaround:</strong> {vendor.turnaround_time || 'N/A'}
                      </p>
                    </div>
                  </div>
                </div>
              ))}
            </div>
          </div>
        );
      }
      
      // If it has testing options
      if (data.testing_options) {
        return (
          <div className="bg-light p-3 border rounded mb-3">
            <h6>Testing Options</h6>
            <div className="row">
              {Object.entries(data.testing_options).map(([testId, test]: [string, any]) => (
                <div key={testId} className="col-md-6 mb-3">
                  <div className="card h-100">
                    <div className="card-body">
                      <h6 className="card-title">{test.name}</h6>
                      <p className="card-text">
                        <strong>Services:</strong> {test.services?.join(', ') || 'N/A'}
                      </p>
                      <p className="card-text">
                        <strong>Vendors:</strong> {test.vendors?.join(', ') || 'N/A'}
                      </p>
                      <p className="card-text">
                        <strong>Pricing:</strong> {test.pricing_range || 'N/A'}
                      </p>
                    </div>
                  </div>
                </div>
              ))}
            </div>
          </div>
        );
      }
      
      // If it has recommendations
      if (data.recommendations) {
        return (
          <div className="bg-light p-3 border rounded mb-3">
            <h6>Recommendations</h6>
            <ul className="list-unstyled">
              {data.recommendations.map((rec: string, index: number) => (
                <li key={index} className="mb-2">
                  <i className="bi bi-lightbulb text-warning me-2"></i>
                  {rec}
                </li>
              ))}
            </ul>
          </div>
        );
      }
      
      // If it has error information
      if (data.error) {
        return (
          <div className="alert alert-danger">
            <h6>Error</h6>
            <p className="mb-0">{data.error}</p>
          </div>
        );
      }
      
      // If it has success information
      if (data.success !== undefined) {
        return (
          <div className={`alert ${data.success ? 'alert-success' : 'alert-warning'}`}>
            <h6>{data.success ? 'Success' : 'Warning'}</h6>
            {data.message && <p className="mb-0">{data.message}</p>}
            {data.result && (
              <div className="mt-2">
                <strong>Result:</strong>
                <pre className="small mt-1">{JSON.stringify(data.result, null, 2)}</pre>
              </div>
            )}
          </div>
        );
      }
      
      // Last resort: show formatted JSON with a note
      return (
        <div className="bg-light p-3 border rounded mb-3">
          <div className="alert alert-info mb-3">
            <i className="bi bi-info-circle me-2"></i>
            Raw response data (for debugging)
          </div>
          <details>
            <summary>View Raw Data</summary>
            <pre className="small mt-2">{JSON.stringify(data, null, 2)}</pre>
          </details>
        </div>
      );
    };
    
    return formatResponse(output);
  };

  return (
    <div className="container py-5">
      <div className="row">
        <div className="col-md-8">
          <h1 className="mb-4">Helix.AI Bioinformatics Interface</h1>
          
          {/* Navigation Tabs */}
          <Nav variant="tabs" activeKey={activeTab} onSelect={(k) => setActiveTab(k || 'command')} className="mb-4">
            <Nav.Item>
              <Nav.Link eventKey="command">Command Interface</Nav.Link>
            </Nav.Item>
          </Nav>
          
          {/* Tab Content */}
          <Tab.Content>
            <Tab.Pane eventKey="command" active={activeTab === 'command'}>
          
          {/* Server Status */}
          <div className="mb-3">
            <span className={`badge ${serverStatus === 'healthy' ? 'bg-success' : 'bg-danger'}`}>
              Server: {serverStatus}
            </span>
            {sessionId && (
              <span className="badge bg-info ms-2">
                Session: {sessionId.substring(0, 8)}...
              </span>
            )}
          </div>
          
          {/* Command Mode Toggle */}
          <div className="mb-3">
            <div className="btn-group" role="group">
              <input
                type="radio"
                className="btn-check"
                name="commandMode"
                id="naturalMode"
                checked={commandMode === 'natural'}
                onChange={() => setCommandMode('natural')}
              />
              <label className="btn btn-outline-primary" htmlFor="naturalMode">
                Natural Language
              </label>
              
              <input
                type="radio"
                className="btn-check"
                name="commandMode"
                id="structuredMode"
                checked={commandMode === 'structured'}
                onChange={() => setCommandMode('structured')}
              />
              <label className="btn btn-outline-primary" htmlFor="structuredMode">
                Structured Commands
              </label>
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
              <button 
                className="btn btn-sm btn-outline-primary mt-2 ms-2"
                onClick={() => {
                  const testSequences = ">sequence1\nATGCGATCGATGCGATCG\n>sequence2\nATGCGATCGATGCGATC-\n>sequence3\nATGCGATCGATGCGATC-";
                  setWorkflowContext(prev => ({ ...prev, alignedSequences: testSequences }));
                  console.log('Test: Set aligned sequences in workflow context');
                }}
              >
                🧪 Test Context
              </button>
            </div>
          )}

          {/* Uploaded File Display */}
          {uploadedFile && (
            <div className="alert alert-info mb-3">
              <div className="d-flex justify-content-between align-items-center">
                <div>
                  <strong>📁 File uploaded:</strong> {uploadedFile.name}
                  <br />
                  <small className="text-muted">
                    Content length: {uploadedFile.content.length} characters
                  </small>
                </div>
                <button 
                  className="btn btn-sm btn-outline-secondary"
                  onClick={() => setUploadedFile(null)}
                >
                  ✕ Remove
                </button>
              </div>
            </div>
          )}
          
          {/* Command Input with Drag-and-Drop */}
      <div
        className={`input-group mb-4${dragActive ? ' border border-primary border-3' : ''}`}
        onDragOver={handleDragOver}
        onDragLeave={handleDragLeave}
        onDrop={handleDrop}
        style={{ position: 'relative', background: dragActive ? '#e6f0ff' : undefined }}
      >
        <textarea
          className="form-control"
          value={command}
          onChange={(e) => setCommand(e.target.value)}
          rows={command.split('\n').length < 4 ? 4 : command.split('\n').length}
          style={{ resize: 'vertical', minHeight: 80 }}
          onKeyDown={(e) => {
            if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') {
              e.preventDefault();
              handleSubmit();
            }
          }}
          placeholder={commandMode === 'natural' 
            ? `Enter natural language command (multi-line supported, Ctrl+Enter to submit)\nExample for FASTA:\n>seq1\nATCGATCGATCG\n>seq2\nATCGATCGATCG`
            : `Enter your bioinformatics command (multi-line supported, Ctrl+Enter to submit)\nExample for FASTA:\n>seq1\nATCGATCGATCG\n>seq2\nATCGATCGATCG`
          }
        />
        <button 
          className="btn btn-primary" 
          onClick={handleSubmit}
          disabled={loading}
          style={{ marginLeft: 8 }}
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
                <strong>Command:</strong> {item.input}
                <small className="text-muted ms-2">
                  ({item.type}) - {item.timestamp.toLocaleTimeString()}
                </small>
              </div>
              {renderOutput(item)}
            </div>
          ))}
            </Tab.Pane>
            

          </Tab.Content>
          </div>

        {/* Sidebar with Available Tools */}
        <div className="col-md-4">
          {/* Example Commands */}
          <div className="card">
            <div className="card-header">
              <h5 className="mb-0">Example Commands</h5>
            </div>
            <div className="card-body">
              {commandMode === 'natural' ? (
                <>
                  <div className="mb-3">
                    <strong>🧬 Sequence Analysis:</strong>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("align these sequences: >seq1 ATGCGATCG >seq2 ATGCGATC")} style={{cursor: 'pointer'}}>
                      "align these sequences: {'>'}seq1 ATGCGATCG {'>'}seq2 ATGCGATC"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("perform multiple sequence alignment on the uploaded sequences")} style={{cursor: 'pointer'}}>
                      "perform multiple sequence alignment on the uploaded sequences"
                </div>
                    <div className="small text-muted cursor-pointer" onClick={() => handleExampleClick("show me the alignment of these DNA sequences")} style={{cursor: 'pointer'}}>
                      "show me the alignment of these DNA sequences"
            </div>
          </div>
          
                  <div className="mb-3">
                    <strong>🎯 Sequence Selection:</strong>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("from the sequence variants, pick 10 sequences randomly")} style={{cursor: 'pointer'}}>
                      "from the sequence variants, pick 10 sequences randomly"
            </div>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("select 5 sequences with the highest mutation rate")} style={{cursor: 'pointer'}}>
                      "select 5 sequences with the highest mutation rate"
                  </div>
                    <div className="small text-muted cursor-pointer" onClick={() => handleExampleClick("choose the most diverse sequences from the alignment")} style={{cursor: 'pointer'}}>
                      "choose the most diverse sequences from the alignment"
                  </div>
                  </div>
                  
                  <div className="mb-3">
                    <strong>🧪 Mutation Generation:</strong>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("mutate the sequence ATGCGATCG to create 96 variants")} style={{cursor: 'pointer'}}>
                      "mutate the sequence ATGCGATCG to create 96 variants"
                  </div>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("generate variants of this DNA sequence")} style={{cursor: 'pointer'}}>
                      "generate variants of this DNA sequence"
                    </div>
                    <div className="small text-muted cursor-pointer" onClick={() => handleExampleClick("create mutations with 0.2 mutation rate")} style={{cursor: 'pointer'}}>
                      "create mutations with 0.2 mutation rate"
                    </div>
                  </div>
                  
                  <div className="mb-3">
                    <strong>🔬 Data Analysis:</strong>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("analyze the alignment and show me the most conserved regions")} style={{cursor: 'pointer'}}>
                      "analyze the alignment and show me the most conserved regions"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("create visualizations of the sequence data")} style={{cursor: 'pointer'}}>
                      "create visualizations of the sequence data"
                    </div>
                    <div className="small text-muted cursor-pointer" onClick={() => handleExampleClick("show me statistics for these sequences")} style={{cursor: 'pointer'}}>
                      "show me statistics for these sequences"
                    </div>
                  </div>
                  
                  <div className="mb-3">
                    <strong>🧬 DNA Synthesis & Testing:</strong>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("I want to order these sequences from a DNA synthesis vendor")} style={{cursor: 'pointer'}}>
                      "I want to order these sequences from a DNA synthesis vendor"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("find DNA synthesis companies for my sequences")} style={{cursor: 'pointer'}}>
                      "find DNA synthesis companies for my sequences"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("what testing options are available for my sequences?")} style={{cursor: 'pointer'}}>
                      "what testing options are available for my sequences?"
                    </div>
                    <div className="small text-muted cursor-pointer" onClick={() => handleExampleClick("research vendors for gene synthesis and validation")} style={{cursor: 'pointer'}}>
                      "research vendors for gene synthesis and validation"
                    </div>
                  </div>
                  
                  <div className="mb-3">
                    <strong>🔄 Multi-step Workflows:</strong>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("mutate this sequence, then align the variants and pick the best ones")} style={{cursor: 'pointer'}}>
                      "mutate this sequence, then align the variants and pick the best ones"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("analyze these sequences and then find vendors to synthesize them")} style={{cursor: 'pointer'}}>
                      "analyze these sequences and then find vendors to synthesize them"
                    </div>
                    <div className="small text-muted cursor-pointer" onClick={() => handleExampleClick("align sequences and find synthesis vendors with testing options")} style={{cursor: 'pointer'}}>
                      "align sequences and find synthesis vendors with testing options"
                    </div>
                  </div>
                  
                  <div className="mb-3">
                    <strong>🔬 Plasmid Visualization:</strong>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("Visualize plasmid pUC19 with EcoRI site and insert ATGCGATCG")} style={{cursor: 'pointer'}}>
                      "Visualize plasmid pUC19 with EcoRI site and insert ATGCGATCG"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("Create plasmid visualization for pBR322 with BamHI site")} style={{cursor: 'pointer'}}>
                      "Create plasmid visualization for pBR322 with BamHI site"
                    </div>
                    <div className="small text-muted cursor-pointer" onClick={() => handleExampleClick("Show plasmid map with cloning sites and insert sequence")} style={{cursor: 'pointer'}}>
                      "Show plasmid map with cloning sites and insert sequence"
                    </div>
                  </div>
                  
                  <div className="mb-3">
                    <strong>🔬 Plasmid Visualization:</strong>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("Visualize plasmid pUC19 with EcoRI site and insert ATGCGATCG")} style={{cursor: 'pointer'}}>
                      "Visualize plasmid pUC19 with EcoRI site and insert ATGCGATCG"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("Create plasmid visualization for pBR322 with BamHI site")} style={{cursor: 'pointer'}}>
                      "Create plasmid visualization for pBR322 with BamHI site"
                    </div>
                    <div className="small text-muted cursor-pointer" onClick={() => handleExampleClick("Show plasmid map with cloning sites and insert sequence")} style={{cursor: 'pointer'}}>
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
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("express ATGCGATCG in pTet vector")} style={{cursor: 'pointer'}}>
                      "express ATGCGATCG in pTet vector"
                    </div>
                    <div className="small text-muted mb-1 cursor-pointer" onClick={() => handleExampleClick("create plasmid visualization for the variants")} style={{cursor: 'pointer'}}>
                      "create plasmid visualization for the variants"
                    </div>
                    <div className="small text-muted cursor-pointer" onClick={() => handleExampleClick("show cloning sites in pUC19 vector")} style={{cursor: 'pointer'}}>
                      "show cloning sites in pUC19 vector"
                    </div>
                    <div className="small text-muted cursor-pointer" onClick={() => handleExampleClick("Visualize plasmid pUC19 with EcoRI site and insert ATGCGATCG")} style={{cursor: 'pointer'}}>
                      "Visualize plasmid pUC19 with EcoRI site and insert ATGCGATCG"
                    </div>
                  </div>
                </>
              )}
              
              <hr className="my-3" />
              <div className="small text-muted">
                <strong>💡 Tips:</strong>
                <ul className="mb-0 mt-1">
                  <li>Upload FASTA files by dragging and dropping</li>
                  <li>Use natural language for complex workflows</li>
                  <li>Combine multiple steps in one command</li>
                  <li>Ask for vendor research and testing options</li>
                </ul>
            </div>
          </div>
        </div>
          
          <div className="card mt-3">
            <div className="card-header">
              <h5 className="mb-0">Available MCP Tools</h5>
            </div>
            <div className="card-body">
              {availableTools.map((tool, index) => (
                <div key={index} className="mb-3">
                  <h6 className="text-primary">{tool.name}</h6>
                  <p className="small text-muted mb-1">{tool.description}</p>
                  {tool.parameters && (
                    <div className="small">
                      <strong>Parameters:</strong>
                      <ul className="list-unstyled mt-1">
                        {Object.entries(tool.parameters).map(([key, value]) => (
                          <li key={key}><code>{key}</code>: {String(value)}</li>
                        ))}
                      </ul>
                    </div>
                  )}
                </div>
              ))}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}

export default App;
