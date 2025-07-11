import React, { useState } from 'react';
import { SeqViz } from 'seqviz';
import { mcpApi } from '../services/mcpApi';
import Button from 'react-bootstrap/Button';
import Form from 'react-bootstrap/Form';
import Card from 'react-bootstrap/Card';
import Alert from 'react-bootstrap/Alert';

interface PlasmidVisualizerProps {
  sessionId?: string;
}

interface PlasmidData {
  name: string;
  sequence: string;
  features: Array<{
    name: string;
    start: number;
    end: number;
    type: string;
    color: string;
    description: string;
  }>;
  size: number;
  description: string;
}

export const PlasmidVisualizer: React.FC<PlasmidVisualizerProps> = ({ sessionId }) => {
  const [vectorName, setVectorName] = useState('pTet');
  const [cloningSites, setCloningSites] = useState('BsaI:1-100');
  const [insertSequence, setInsertSequence] = useState('ACGT');
  const [plasmidData, setPlasmidData] = useState<PlasmidData | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleSubmit = async () => {
    setLoading(true);
    setError(null);
    
    try {
      const response = await mcpApi.plasmidVisualization({
        vector_name: vectorName,
        cloning_sites: cloningSites,
        insert_sequence: insertSequence
      }, sessionId);
      
      if (response.success && response.result.plasmid_data) {
        setPlasmidData(response.result.plasmid_data);
      } else {
        setError(response.result.text || 'Failed to generate plasmid visualization');
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : 'An error occurred');
    } finally {
      setLoading(false);
    }
  };

  const convertToSeqVizFormat = (data: PlasmidData) => {
    // Convert our plasmid data format to SeqViz format
    const features = data.features.map(feature => ({
      name: feature.name,
      start: feature.start,
      end: feature.end,
      color: feature.color,
      type: feature.type === 'restriction_site' ? 'enzyme' : 'feature'
    }));

    return {
      name: data.name,
      sequence: data.sequence,
      features: features,
      size: data.size
    };
  };

  return (
    <div className="plasmid-visualizer">
      <Card className="mb-3">
        <Card.Header>
          <h5>Plasmid Visualization</h5>
        </Card.Header>
        <Card.Body>
          <Form>
            <Form.Group className="mb-3">
              <Form.Label>Vector Name</Form.Label>
              <Form.Control
                type="text"
                value={vectorName}
                onChange={(e) => setVectorName(e.target.value)}
                placeholder="e.g., pTet"
              />
              <Form.Text className="text-muted">
                Enter the name of your plasmid vector
              </Form.Text>
            </Form.Group>

            <Form.Group className="mb-3">
              <Form.Label>Cloning Sites</Form.Label>
              <Form.Control
                type="text"
                value={cloningSites}
                onChange={(e) => setCloningSites(e.target.value)}
                placeholder="e.g., BsaI:123-456, EcoRI:789-1012"
              />
              <Form.Text className="text-muted">
                Enter restriction sites in format: Enzyme:start-end
              </Form.Text>
            </Form.Group>

            <Form.Group className="mb-3">
              <Form.Label>Insert Sequence</Form.Label>
              <Form.Control
                as="textarea"
                rows={3}
                value={insertSequence}
                onChange={(e) => setInsertSequence(e.target.value)}
                placeholder="e.g., ACGT"
              />
              <Form.Text className="text-muted">
                Enter the DNA sequence to insert (A, T, C, G only)
              </Form.Text>
            </Form.Group>

            <Button 
              variant="primary" 
              onClick={handleSubmit}
              disabled={loading}
            >
              {loading ? 'Generating...' : 'Generate Plasmid Visualization'}
            </Button>
          </Form>
        </Card.Body>
      </Card>

      {error && (
        <Alert variant="danger" className="mb-3">
          {error}
        </Alert>
      )}

      {plasmidData && (
        <Card>
          <Card.Header>
            <h6>{plasmidData.name} - {plasmidData.description}</h6>
          </Card.Header>
          <Card.Body>
            <div style={{ height: '400px', width: '100%' }}>
              <SeqViz
                {...convertToSeqVizFormat(plasmidData)}
                viewer="circular"
                showFeatures={true}
                showPrimers={false}
                showAnnotations={true}
                showGC={false}
                showORFs={false}
                showRestrictionSites={true}
                showAxis={true}
                showSequence={false}
                showComplement={false}
                showReverse={false}
                zoom={{ linear: 50, circular: 50 }}
                colors={{
                  features: plasmidData.features.map(f => f.color),
                  primers: "#ff0000",
                  orfs: "#00ff00",
                  restrictionSites: "#0000ff",
                  gc: "#ff00ff",
                  axis: "#000000",
                  sequence: "#000000",
                  complement: "#666666",
                  reverse: "#999999"
                }}
              />
            </div>
          </Card.Body>
        </Card>
      )}
    </div>
  );
}; 