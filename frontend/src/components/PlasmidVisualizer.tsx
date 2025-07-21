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

    // Convert our plasmid data format to SeqViz format
function convertToSeqVizFormat(data: PlasmidData) {
    const features = data.features.map(feature => ({
      name: feature.name,
      start: feature.start,
      end: feature.end,
      color: feature.color,
    type: feature.type === 'restriction_site' ? 'enzyme' : 'feature',
    description: feature.description,
    }));

    return {
      name: data.name,
      sequence: data.sequence,
      features: features,
    size: data.size,
  };
}

// Component for displaying plasmid data directly (used in command results)
export const PlasmidDataVisualizer: React.FC<PlasmidDataVisualizerProps> = ({ data }) => {
  const features = data.features || [];
  const seqvizData = convertToSeqVizFormat(data);
  return (
        <Card>
          <Card.Header>
        <h6>{data.name} - {data.description}</h6>
          </Card.Header>
          <Card.Body>
        <div style={{ marginBottom: 20 }}>
          <strong>Name:</strong> {data.name}<br />
          <strong>Size:</strong> {data.size} bp<br />
          <strong>Sequence:</strong> {data.sequence}<br />
          <strong>Features:</strong> {features.length}
          <ul>
            {features.map((feature, idx) => (
              <li key={idx}>
                <strong>{feature.name}</strong> ({feature.type}): {feature.description || ''} [<span>start: {feature.start}, end: {feature.end}</span>]
              </li>
            ))}
          </ul>
        </div>
        <div>
              <SeqViz
            name={seqvizData.name}
            seq={seqvizData.sequence}
            features={seqvizData.features}
            circular={true}
            zoom={{ linear: 1, circular: 1 }}
            showAnnotations={true}
                showPrimers={false}
            showCutsites={false}
            showComplement={true}
            showIndex={true}
            style={{ width: '100%', minHeight: 400 }}
              />
            </div>
          </Card.Body>
        </Card>
  );
}; 