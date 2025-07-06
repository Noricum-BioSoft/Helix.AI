import React, { useState, useEffect } from 'react';
import { MSAView } from 'react-msaview';
import Plot from 'react-plotly.js';
import { mcpApi } from './services/mcpApi';
import { CommandParser, ParsedCommand } from './utils/commandParser';

import 'bootstrap/dist/css/bootstrap.min.css';

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

  useEffect(() => {
    // Check server health and load available tools on component mount
    checkServerHealth();
    loadAvailableTools();
  }, []);

  const checkServerHealth = async () => {
    try {
      const health = await mcpApi.healthCheck();
      setServerStatus(health.status);
    } catch (error) {
      setServerStatus('error');
      console.error('Server health check failed:', error);
    }
  };

  const loadAvailableTools = async () => {
    try {
      const toolsResponse = await mcpApi.listTools();
      setAvailableTools(toolsResponse.tools);
    } catch (error) {
      console.error('Failed to load available tools:', error);
    }
  };

  const handleSubmit = async () => {
  if (!command.trim()) return;
    
  setLoading(true);
    
    try {
      // Parse the command to determine the appropriate MCP endpoint
      const parsedCommand: ParsedCommand = CommandParser.parseCommand(command);
      
      let response;
      
      switch (parsedCommand.type) {
        case 'sequence_alignment':
          response = await mcpApi.sequenceAlignment(parsedCommand.parameters);
          break;
          
        case 'mutate_sequence':
          response = await mcpApi.mutateSequence(parsedCommand.parameters);
          break;
          
        case 'analyze_sequence_data':
          response = await mcpApi.analyzeSequenceData(parsedCommand.parameters);
          break;
          
        case 'visualize_alignment':
          response = await mcpApi.visualizeAlignment(parsedCommand.parameters);
          break;
          
        default:
          // Fallback to general command execution
          response = await mcpApi.executeCommand(command);
          break;
      }
      
      // Add to history
      const historyItem: HistoryItem = {
        input: command,
        output: response,
        type: parsedCommand.type,
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
      
      return (
        <div>
          {result.text && (
            <pre className="bg-light p-3 border rounded mb-3">
              {result.text}
            </pre>
          )}
          
          {result.alignment && Array.isArray(result.alignment) && (
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
                    {result.alignment.map((seq: any, index: number) => (
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
          
          {result.output && (
            <div className="bg-light p-3 border rounded mb-3">
              <h5>Output Data</h5>
              <pre className="mb-0">{JSON.stringify(result.output, null, 2)}</pre>
            </div>
          )}
          
          {result.plot && (
            <Plot
              data={result.plot.data}
              layout={result.plot.layout}
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
          </div>
          
          {/* Command Input */}
      <div className="input-group mb-4">
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
              placeholder="Enter your bioinformatics command (e.g., 'align sequences', 'mutate ACTGTTGAC', 'analyze sequence data')"
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
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}

export default App;
