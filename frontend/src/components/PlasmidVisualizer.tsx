import React from 'react';
import { Card } from 'react-bootstrap';
import { SeqViz } from 'seqviz';

// Types for plasmid data
export interface PlasmidFeature {
  name: string;
  start: number;
  end: number;
  type: string;
  color?: string;
  description?: string;
}

export interface PlasmidData {
  name: string;
  sequence: string;
  features: PlasmidFeature[];
  size: number;
  description?: string;
}

export interface PlasmidDataVisualizerProps {
  data: PlasmidData;
}

// Convert our plasmid data format to SeqViz annotations format
function convertFeaturesToAnnotations(features: PlasmidFeature[]) {
  return features.map(feature => ({
    start: feature.start,
    end: feature.end,
    name: feature.name,
    direction: 1, // Default direction
    color: feature.color || "#ff0000"
  }));
}

// Component for displaying plasmid data directly (used in command results)
export const PlasmidDataVisualizer: React.FC<PlasmidDataVisualizerProps> = ({ data }) => {
  // Validate data before rendering
  if (!data) {
    return (
      <div className="alert alert-warning">
        <p>No plasmid data provided.</p>
      </div>
    );
  }
  
  // Validate sequence - SeqViz requires a non-empty string
  const sequence = data.sequence;
  if (!sequence || typeof sequence !== 'string' || sequence.length === 0) {
    return (
      <Card>
        <Card.Header>
          <h6>{data.name || 'Plasmid'} - {data.description || 'No description'}</h6>
        </Card.Header>
        <Card.Body>
          <div className="alert alert-warning">
            <p><strong>Invalid sequence data:</strong> Sequence is missing or invalid.</p>
            <details>
              <summary>Debug Info</summary>
              <pre>{JSON.stringify(data, null, 2)}</pre>
            </details>
          </div>
        </Card.Body>
      </Card>
    );
  }
  
  const features = data.features || [];
  const annotations = convertFeaturesToAnnotations(features);
  const [viewerType, setViewerType] = React.useState<'circular' | 'linear' | 'both'>('circular');
  
  return (
    <Card>
      <Card.Header>
        <h6>{data.name || 'Plasmid'} - {data.description || 'No description'}</h6>
      </Card.Header>
      <Card.Body>
        <div style={{ marginBottom: 20 }}>
          <strong>Name:</strong> {data.name || 'N/A'}<br />
          <strong>Size:</strong> {data.size || sequence.length} bp<br />
          <strong>Sequence Length:</strong> {sequence.length} bp<br />
          <strong>Features:</strong> {features.length}
          {features.length > 0 && (
            <ul>
              {features.map((feature, idx) => (
                <li key={idx}>
                  <strong>{feature.name}</strong> ({feature.type}): {feature.description || ''} [<span>start: {feature.start}, end: {feature.end}</span>]
                </li>
              ))}
            </ul>
          )}
        </div>
        
        {/* Viewer Type Selector */}
        <div className="mb-3">
          <div className="d-flex align-items-center">
            <label htmlFor="viewer-select" className="me-2">View:</label>
            <select 
              id="viewer-select"
              className="form-select form-select-sm" 
              style={{width: 'auto'}}
              value={viewerType}
              onChange={(e) => setViewerType(e.target.value as 'linear' | 'circular' | 'both')}
            >
              <option value="circular">Circular</option>
              <option value="linear">Linear</option>
              <option value="both">Both</option>
            </select>
          </div>
        </div>
        
        <div style={{ height: '500px', border: '1px solid #ccc', borderRadius: '8px', overflow: 'hidden' }}>
          <SeqViz
            name={data.name || 'Plasmid'}
            seq={sequence}
            seqType="dna"
            annotations={annotations}
            viewer={viewerType}
            style={{ width: '100%', height: '100%' }}
          />
        </div>
      </Card.Body>
    </Card>
  );
};

// Component for displaying multiple plasmid visualizations for representatives
export const PlasmidRepresentativesVisualizer: React.FC<{ 
  plasmidResults: any[]; 
  vectorName: string; 
  cloningSites: string;
}> = ({ plasmidResults, vectorName, cloningSites }) => {
  const [selectedPlasmid, setSelectedPlasmid] = React.useState<number>(0);
  const [viewerType, setViewerType] = React.useState<'circular' | 'linear' | 'both'>('circular');

  if (!plasmidResults || plasmidResults.length === 0) {
    return <div className="alert alert-warning">No plasmid results available.</div>;
  }

  return (
    <div className="bg-light p-3 border rounded mb-3">
      <h5>🧬 Plasmid Visualizations for Representative Sequences</h5>
      <div className="mb-3">
        <p><strong>Vector:</strong> {vectorName}</p>
        <p><strong>Cloning Sites:</strong> {cloningSites}</p>
        <p><strong>Total Representatives:</strong> {plasmidResults.length}</p>
      </div>
      
      {/* Plasmid Selector */}
      <div className="mb-3">
        <label htmlFor="plasmid-select" className="form-label">Select Representative:</label>
        <select 
          id="plasmid-select"
          className="form-select"
          value={selectedPlasmid}
          onChange={(e) => setSelectedPlasmid(parseInt(e.target.value))}
        >
          {plasmidResults.map((result, index) => (
            <option key={index} value={index}>
              {result.representative_name} ({result.sequence_length} bp)
            </option>
          ))}
        </select>
      </div>
      
      {/* Viewer Type Selector */}
      <div className="mb-3">
        <div className="d-flex align-items-center">
          <label htmlFor="viewer-type-select" className="me-2">View:</label>
          <select 
            id="viewer-type-select"
            className="form-select form-select-sm" 
            style={{width: 'auto'}}
            value={viewerType}
            onChange={(e) => setViewerType(e.target.value as 'linear' | 'circular' | 'both')}
          >
            <option value="circular">Circular</option>
            <option value="linear">Linear</option>
            <option value="both">Both</option>
          </select>
        </div>
      </div>
      
      {/* Selected Plasmid Visualization */}
      {plasmidResults[selectedPlasmid] && (() => {
        const selectedResult = plasmidResults[selectedPlasmid];
        const sequence = selectedResult?.plasmid_data?.sequence;
        const isValidSequence = sequence && typeof sequence === 'string' && sequence.length > 0;
        
        if (!isValidSequence) {
          return (
            <div className="alert alert-warning">
              <p><strong>Invalid sequence data:</strong> Sequence is missing or invalid for {selectedResult.representative_name}.</p>
            </div>
          );
        }
        
        return (
          <div>
            <h6>Plasmid: {selectedResult.representative_name}</h6>
            <div style={{ height: '500px', border: '1px solid #ccc', borderRadius: '8px', overflow: 'hidden' }}>
              <SeqViz
                name={selectedResult.plasmid_data.name}
                seq={sequence}
                seqType="dna"
                annotations={convertFeaturesToAnnotations(selectedResult.plasmid_data.features || [])}
                viewer={viewerType}
                style={{ width: '100%', height: '100%' }}
              />
            </div>
          </div>
        );
      })()}
      
      {/* Summary Table */}
      <div className="mt-3">
        <h6>All Representatives Summary</h6>
        <div className="table-responsive">
          <table className="table table-sm">
            <thead>
              <tr>
                <th>Name</th>
                <th>Sequence Length</th>
                <th>Insert Sequence (first 50 bp)</th>
              </tr>
            </thead>
            <tbody>
              {plasmidResults.map((result, index) => (
                <tr key={index} className={index === selectedPlasmid ? 'table-primary' : ''}>
                  <td>{result.representative_name}</td>
                  <td>{result.sequence_length} bp</td>
                  <td className="font-monospace">
                    {result.insert_sequence.substring(0, 50)}
                    {result.insert_sequence.length > 50 ? '...' : ''}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  );
}; 