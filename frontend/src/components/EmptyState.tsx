import React from 'react';
import Button from 'react-bootstrap/Button';
import Card from 'react-bootstrap/Card';
import { theme } from '../theme';
import { getExampleWithSequences, sampleSequences } from '../utils/sampleSequences';

interface EmptyStateProps {
  onExampleClick: (command: string) => void;
  onFileUpload: () => void;
}

export const EmptyState: React.FC<EmptyStateProps> = ({ 
  onExampleClick, 
  onFileUpload 
}) => {
  const quickActions = [
    {
      icon: "🌳",
      title: "Visualize Tree",
      command: getExampleWithSequences("visualize the phylogenetic tree", "phylogenetic"),
      description: "Create an interactive phylogenetic tree with sample sequences"
    },
    {
      icon: "📊",
      title: "Align Sequences",
      command: getExampleWithSequences("align these sequences", "threeSequences"),
      description: "Perform multiple sequence alignment with sample sequences"
    },
    {
      icon: "🔬",
      title: "Research Vendors",
      command: "research DNA synthesis vendors for 1000bp sequences",
      description: "Find DNA synthesis and testing options (no sequences needed)"
    },
    {
      icon: "🧬",
      title: "Create Variants",
      command: getExampleWithSequences("create 96 variants", "single"),
      description: "Generate 96 sequence mutations from a sample sequence"
    }
  ];

  return (
    <div className="empty-state text-center py-5" style={{ minHeight: '400px' }}>
      <div className="mb-4">
        <div style={{ fontSize: '4rem', marginBottom: '1rem' }}>🧬</div>
        <h2 className="mb-3">Get Started with Helix.AI</h2>
        <p className="text-muted lead">
          Start by trying an example command or uploading a file
        </p>
      </div>

      <div className="row g-3 mb-4">
        {quickActions.map((action, idx) => (
          <div key={idx} className="col-md-6">
            <Card 
              className="h-100 hover-shadow"
              style={{ 
                cursor: 'pointer',
                transition: 'all 0.2s',
                border: '1px solid #dee2e6'
              }}
              onClick={() => onExampleClick(action.command)}
              onMouseEnter={(e) => {
                e.currentTarget.style.transform = 'translateY(-2px)';
                e.currentTarget.style.boxShadow = theme.shadows.gold;
                e.currentTarget.style.borderColor = theme.colors.gold;
              }}
              onMouseLeave={(e) => {
                e.currentTarget.style.transform = 'translateY(0)';
                e.currentTarget.style.boxShadow = 'none';
                e.currentTarget.style.borderColor = '#dee2e6';
              }}
            >
              <Card.Body>
                <div style={{ fontSize: '2rem', marginBottom: '0.5rem' }}>
                  {action.icon}
                </div>
                <h6 className="mb-2">{action.title}</h6>
                <p className="text-muted small mb-0">{action.description}</p>
                <small className="text-primary mt-2 d-block">Click to try →</small>
              </Card.Body>
            </Card>
          </div>
        ))}
      </div>

      <div className="mb-4">
        <Button 
          variant="outline-primary" 
          size="lg"
          onClick={onFileUpload}
          style={{ 
            marginRight: '10px',
            borderColor: theme.colors.gold,
            color: theme.colors.gold
          }}
          onMouseEnter={(e) => {
            e.currentTarget.style.backgroundColor = theme.colors.gold;
            e.currentTarget.style.color = theme.colors.textOnGold;
          }}
          onMouseLeave={(e) => {
            e.currentTarget.style.backgroundColor = 'transparent';
            e.currentTarget.style.color = theme.colors.gold;
          }}
        >
          📁 Upload a FASTA File
        </Button>
        <Button 
          variant="outline-secondary" 
          size="lg"
          onClick={() => onExampleClick(getExampleWithSequences("visualize the phylogenetic tree", "phylogenetic"))}
          style={{
            borderColor: theme.colors.blue,
            color: theme.colors.blue
          }}
          onMouseEnter={(e) => {
            e.currentTarget.style.backgroundColor = theme.colors.blue;
            e.currentTarget.style.color = theme.colors.white;
          }}
          onMouseLeave={(e) => {
            e.currentTarget.style.backgroundColor = 'transparent';
            e.currentTarget.style.color = theme.colors.blue;
          }}
        >
          🚀 Try Example Command
        </Button>
      </div>

      <div className="text-muted small">
        <p className="mb-1">
          💡 Tip: You can use natural language commands like "align sequences" or "create a phylogenetic tree"
        </p>
        <p className="mb-0">
          📖 Need help? Check out the example commands panel for more ideas
        </p>
      </div>
    </div>
  );
};

