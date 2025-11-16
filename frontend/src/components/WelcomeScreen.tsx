import React, { useState } from 'react';
import Button from 'react-bootstrap/Button';
import Card from 'react-bootstrap/Card';
import { theme } from '../theme';
import { getExampleWithSequences, sampleSequences } from '../utils/sampleSequences';

interface WelcomeScreenProps {
  onDismiss: () => void;
  onExampleClick: (command: string) => void;
}

export const WelcomeScreen: React.FC<WelcomeScreenProps> = ({ 
  onDismiss, 
  onExampleClick 
}) => {
  const [dontShowAgain, setDontShowAgain] = useState(false);

  const handleDismiss = () => {
    if (dontShowAgain) {
      localStorage.setItem('helix-ai-welcome-dismissed', 'true');
    }
    onDismiss();
  };

  const quickStartExamples = [
    {
      command: getExampleWithSequences("visualize the phylogenetic tree", "phylogenetic"),
      label: "visualize the phylogenetic tree",
      description: "Create an interactive phylogenetic tree"
    },
    {
      command: "research DNA synthesis vendors for 1000bp sequences",
      label: "research DNA synthesis vendors",
      description: "Find DNA synthesis and testing options"
    },
    {
      command: getExampleWithSequences("align these sequences", "threeSequences"),
      label: "align sequences",
      description: "Perform multiple sequence alignment"
    },
    {
      command: getExampleWithSequences("create 96 variants", "single"),
      label: "create 96 variants",
      description: "Generate sequence mutations"
    }
  ];

  return (
    <div className="welcome-screen" style={{
      position: 'fixed',
      top: 0,
      left: 0,
      right: 0,
      bottom: 0,
      backgroundColor: 'rgba(0, 0, 0, 0.5)',
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      zIndex: 1000,
      padding: '20px'
    }}>
      <Card style={{ maxWidth: '700px', width: '100%' }}>
        <Card.Body style={{ padding: '2rem' }}>
          <div className="text-center mb-4">
            <h1 className="mb-3">🧬 Welcome to Helix.AI</h1>
            <p className="lead text-muted">
              Your AI-powered bioinformatics assistant for managing biotechnology workflows
            </p>
            <p className="text-muted">
              Use natural language commands to align sequences, build phylogenetic trees, 
              research DNA synthesis vendors, and more.
            </p>
          </div>

          <div className="mb-4">
            <h5 className="mb-3">🚀 Quick Start</h5>
            <div className="row g-2">
              {quickStartExamples.map((example, idx) => (
                <div key={idx} className="col-md-6">
                  <Button
                    variant="outline-primary"
                    className="w-100 text-start"
                    onClick={() => {
                      onExampleClick(example.command);
                      handleDismiss();
                    }}
                    style={{ 
                      whiteSpace: 'normal',
                      wordBreak: 'break-word',
                      height: 'auto',
                      padding: '10px',
                      borderColor: theme.colors.gold,
                      color: theme.colors.gold
                    }}
                    title={example.description}
                  >
                    <div>
                      <strong>{example.label}</strong>
                      <br />
                      <small className="text-muted" style={{ fontSize: '0.85rem' }}>
                        {example.description}
                      </small>
                    </div>
                  </Button>
                </div>
              ))}
            </div>
          </div>

          <div className="mb-4">
            <h6>What you can do:</h6>
            <ul className="list-unstyled">
              <li className="mb-2">📊 <strong>Analyze sequences</strong> - Align, compare, and visualize DNA/RNA sequences</li>
              <li className="mb-2">🌳 <strong>Build phylogenetic trees</strong> - Create interactive tree visualizations</li>
              <li className="mb-2">🧬 <strong>Generate variants</strong> - Create and analyze sequence mutations</li>
              <li className="mb-2">🔬 <strong>Research vendors</strong> - Find DNA synthesis and testing options</li>
              <li className="mb-2">🧪 <strong>Visualize plasmids</strong> - Create plasmid and vector visualizations</li>
            </ul>
          </div>

          <div className="d-flex justify-content-between align-items-center">
            <div className="form-check">
              <input
                className="form-check-input"
                type="checkbox"
                id="dontShowAgain"
                checked={dontShowAgain}
                onChange={(e) => setDontShowAgain(e.target.checked)}
              />
              <label className="form-check-label" htmlFor="dontShowAgain">
                Don't show this again
              </label>
            </div>
            <Button 
              variant="primary" 
              onClick={handleDismiss}
              style={{ 
                backgroundColor: theme.colors.gold,
                borderColor: theme.colors.gold,
                color: theme.colors.textOnGold
              }}
            >
              Get Started
            </Button>
          </div>
        </Card.Body>
      </Card>
    </div>
  );
};

