import React, { useState, useEffect } from 'react';
import { MSAView } from 'react-msaview';
import Plot from 'react-plotly.js';
import { mcpApi } from './services/mcpApi';
import { CommandParser, ParsedCommand } from './utils/commandParser';

import 'bootstrap/dist/css/bootstrap.min.css';
import Modal from 'react-bootstrap/Modal';
import Button from 'react-bootstrap/Button';

interface HistoryItem {
  input: string;
  output: any;
  type: string;
  timestamp: Date;
}

function App() {
  const [command, setCommand] = useState('');
  const [history, setHistory] = useState<HistoryItem[]>([]);
  const [loading, setLoading] = useState(false);
  const [availableTools, setAvailableTools] = useState<any[]>([]);
  const [serverStatus, setServerStatus] = useState<string>('unknown');
  const [dragActive, setDragActive] = useState(false);
  const [droppedFile, setDroppedFile] = useState<File | null>(null);
  const [fileContent, setFileContent] = useState<string>('');
  const [showFileModal, setShowFileModal] = useState(false);
  const [fileAction, setFileAction] = useState<string>('');
  const [sessionId, setSessionId] = useState<string>('');
  const [commandMode, setCommandMode] = useState<'structured' | 'natural'>('natural');

  useEffect(() => {
    // Check server health and load available tools on component mount
    checkServerHealth();
    loadAvailableTools();
    createSession();
  }, []);

  const createSession = async () => {
    try {
      const session = await mcpApi.createSession();
      setSessionId((session as any).session_id);
    } catch (error) {
      console.error('Failed to create session:', error);
    }
  };

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

  const handleSubmit = async () => {
    if (!command.trim()) return;
    
    setLoading(true);
    
    try {
      let response;
      let parsedCommand: ParsedCommand | undefined;
      
      console.log('Command mode:', commandMode);
      console.log('Session ID:', sessionId);
      console.log('Command:', command);
      
      if (commandMode === 'natural') {
        // Use natural language command handling
        console.log('Calling handleNaturalCommand...');
        response = await mcpApi.handleNaturalCommand(command, sessionId);
        console.log('Natural command response:', response);
      } else {
        // Use structured command parsing
        parsedCommand = CommandParser.parseCommand(command);
        
        switch (parsedCommand.type) {
          case 'sequence_alignment':
            response = await mcpApi.sequenceAlignment(parsedCommand.parameters as any);
            break;
            
          case 'mutate_sequence':
            response = await mcpApi.mutateSequence(parsedCommand.parameters as any);
            break;
            
          case 'analyze_sequence_data':
            response = await mcpApi.analyzeSequenceData(parsedCommand.parameters as any);
            break;
            
          case 'visualize_alignment':
            response = await mcpApi.visualizeAlignment(parsedCommand.parameters as any);
            break;
            
          default:
            // Fallback to general command execution
            response = await mcpApi.executeCommand(command);
            break;
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
      setDroppedFile(file);
      const reader = new FileReader();
      reader.onload = (event) => {
        setFileContent(event.target?.result as string);
        setShowFileModal(true);
      };
      reader.readAsText(file);
    }
  };

  // Handle file action selection
  const handleFileAction = async (action: string) => {
    setFileAction(action);
    setShowFileModal(false);
    if (!fileContent) return;
    setLoading(true);
    let response;
    let inputText = '';
    try {
      if (action === 'align') {
        inputText = `align sequences ${fileContent.trim()}`;
        response = await mcpApi.sequenceAlignment({ sequences: fileContent });
      } else if (action === 'mutate') {
        inputText = `mutate sequence ${fileContent.trim()}`;
        response = await mcpApi.mutateSequence({ sequence: fileContent });
      } else if (action === 'analyze') {
        inputText = `analyze sequence data`;
        response = await mcpApi.analyzeSequenceData({ data: fileContent });
      } else {
        // Cancel
        setDroppedFile(null);
        setFileContent('');
        setFileAction('');
        setLoading(false);
        return;
      }
      // Add to history
      const historyItem: HistoryItem = {
        input: inputText,
        output: response,
        type: action,
        timestamp: new Date()
      };
      setHistory(prev => [historyItem, ...prev]);
    } catch (error) {
      const historyItem: HistoryItem = {
        input: inputText,
        output: { error: error instanceof Error ? error.message : 'Unknown error' },
        type: 'error',
        timestamp: new Date()
      };
      setHistory(prev => [historyItem, ...prev]);
    } finally {
      setDroppedFile(null);
      setFileContent('');
      setFileAction('');
      setLoading(false);
    }
  };

  const renderOutput = (item: HistoryItem) => {
    const { output, type } = item;
    
    if (type === 'error') {
      return (
        <div className="alert alert-danger">
          <strong>Error:</strong> {output.error}
        </div>
      );
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
      
      // Handle nested result structure from natural language commands
      const actualResult = result.result || result;
      
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
    
    return (
      <pre className="bg-light p-3 border rounded mb-3">
        {JSON.stringify(output, null, 2)}
      </pre>
    );
  };

  return (
    <div className="container py-5">
      <div className="row">
        <div className="col-md-8">
          <h1 className="mb-4">DataBloom.AI Bioinformatics Interface</h1>
          
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
          
          {/* Command Input with Drag-and-Drop */}
      <div
        className={`input-group mb-4${dragActive ? ' border border-primary border-3' : ''}`}
        onDragOver={handleDragOver}
        onDragLeave={handleDragLeave}
        onDrop={handleDrop}
        style={{ position: 'relative', background: dragActive ? '#e6f0ff' : undefined }}
      >
        <input
          type="text"
          className="form-control"
          value={command}
          onChange={(e) => setCommand(e.target.value)}
          onKeyDown={(e) => {
            if (e.key === 'Enter') {
              e.preventDefault();
              handleSubmit();
            }
          }}
          placeholder={commandMode === 'natural' 
            ? "Enter natural language command (e.g., 'from the sequence variants, pick 10 sequences randomly') or drop a file here"
            : "Enter your bioinformatics command (e.g., 'align sequences', 'mutate ACTGTTGAC', 'analyze sequence data') or drop a file here"
          }
        />
        <button 
          className="btn btn-primary" 
          onClick={handleSubmit}
          disabled={loading}
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

      {/* File Action Modal */}
      <Modal show={showFileModal} onHide={() => setShowFileModal(false)}>
        <Modal.Header closeButton>
          <Modal.Title>File Uploaded</Modal.Title>
        </Modal.Header>
        <Modal.Body>
          <p>What would you like to do with the file <strong>{droppedFile?.name}</strong>?</p>
          <div className="d-grid gap-2">
            <Button variant="primary" onClick={() => handleFileAction('align')}>
              Align Sequences
            </Button>
            <Button variant="secondary" onClick={() => handleFileAction('mutate')}>
              Mutate Sequence
            </Button>
            <Button variant="success" onClick={() => handleFileAction('analyze')}>
              Analyze Data
            </Button>
            <Button variant="outline-secondary" onClick={() => handleFileAction('cancel')}>
              Cancel
            </Button>
          </div>
        </Modal.Body>
      </Modal>

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
          </div>

        {/* Sidebar with Available Tools */}
        <div className="col-md-4">
          <div className="card">
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
                          <li key={key}><code>{key}</code>: {value}</li>
                        ))}
                      </ul>
                    </div>
                  )}
                </div>
              ))}
            </div>
          </div>
          
          {/* Example Commands */}
          <div className="card mt-3">
            <div className="card-header">
              <h5 className="mb-0">Example Commands</h5>
            </div>
            <div className="card-body">
              {commandMode === 'natural' ? (
                <>
                  <div className="mb-2">
                    <strong>Natural Language:</strong>
                    <div className="small text-muted">"from the sequence variants, pick 10 sequences randomly and output them"</div>
                  </div>
                  <div className="mb-2">
                    <strong>Variant Selection:</strong>
                    <div className="small text-muted">"select 5 sequences with the highest mutation rate"</div>
                  </div>
                  <div className="mb-2">
                    <strong>Analysis:</strong>
                    <div className="small text-muted">"analyze the alignment and show me the most conserved regions"</div>
                  </div>
                  <div className="mb-2">
                    <strong>Multi-step:</strong>
                    <div className="small text-muted">"mutate this sequence, then align the variants and pick the best ones"</div>
                  </div>
                </>
              ) : (
                <>
                  <div className="mb-2">
                    <strong>Sequence Alignment:</strong>
                    <div className="small text-muted">"align sequences ACTGTTGAC ACTGCATCC"</div>
                  </div>
                  <div className="mb-2">
                    <strong>Mutation:</strong>
                    <div className="small text-muted">"mutate sequence ACTGTTGAC with 10 variants"</div>
                  </div>
                  <div className="mb-2">
                    <strong>Analysis:</strong>
                    <div className="small text-muted">"analyze sequence data for phylogeny"</div>
                  </div>
                  <div className="mb-2">
                    <strong>Visualization:</strong>
                    <div className="small text-muted">"visualize alignment in PNG format"</div>
                  </div>
                </>
              )}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}

export default App;
