import React, { useState } from 'react';
import Button from 'react-bootstrap/Button';
import Card from 'react-bootstrap/Card';
import { theme } from '../theme';
import { getExampleWithSequences, sampleSequences } from '../utils/sampleSequences';

interface ExampleCommand {
  command: string;
  description?: string;
}

interface ExampleCategory {
  category: string;
  icon: string;
  commands: ExampleCommand[];
}

interface ExampleCommandsPanelProps {
  onCommandSelect: (command: string) => void;
}

export const ExampleCommandsPanel: React.FC<ExampleCommandsPanelProps> = ({ 
  onCommandSelect 
}) => {
  const [expandedCategories, setExpandedCategories] = useState<Set<string>>(new Set());

  const exampleCategories: ExampleCategory[] = [
    {
      category: "Getting Started",
      icon: "🚀",
      commands: [
        {
          command: "What is bioinformatics?",
          description: "Ask the agent for a plain-language overview of bioinformatics"
        },
        {
          command: "What is a sequence alignment?",
          description: "Get a short explanation of sequence alignments and when to use them"
        },
        {
          command: "How do you choose a DNA sequencing workflow?",
          description: "Request a suggested workflow for a sequencing or analysis question"
        }
      ]
    },
    {
      category: "Sequence Analysis",
      icon: "📊",
      commands: [
        { 
          command: getExampleWithSequences("perform multiple sequence alignment", "threeSequences"), 
          description: "Align multiple sequences (requires sequences in command or uploaded file)" 
        },
        { 
          command: "select 10 representative sequences", 
          description: "Choose representative sequences from alignment (requires previous alignment step)" 
        },
        { 
          command: "select 5 sequences with the highest mutation rate", 
          description: "Select sequences based on criteria (requires previous analysis)" 
        }
      ]
    },
    {
      category: "Mutation & Variants",
      icon: "🧬",
      commands: [
        { 
          command: getExampleWithSequences("create 96 variants", "single"), 
          description: "Generate 96 sequence variants from a sample sequence" 
        },
        { 
          command: getExampleWithSequences("mutate sequence with 0.1 mutation rate", "single"), 
          description: "Create mutations with specific rate (includes full sample sequence)" 
        }
      ]
    },
    {
      category: "Visualization",
      icon: "🌳",
      commands: [
        { 
          command: getExampleWithSequences("visualize the phylogenetic tree", "phylogenetic"), 
          description: "Display phylogenetic tree with sample sequences" 
        },
        { 
          command: `show me the plasmid visualization with insert ${sampleSequences.plasmidInsert}`, 
          description: "Visualize plasmid with sample insert sequence" 
        },
        { 
          command: "insert each of the sequences into pUC19", 
          description: "Create plasmid visualizations (requires sequences from previous step or uploaded file)" 
        }
      ]
    },
    {
      category: "Vendor Research",
      icon: "🔬",
      commands: [
        { 
          command: "research DNA synthesis vendors for 1000bp sequences", 
          description: "Find DNA synthesis vendors (no sequences needed - uses length parameter)" 
        },
        { 
          command: "find testing options for protein expression", 
          description: "Research testing options (no sequences needed)" 
        },
        { 
          command: "compare DNA synthesis vendors for large quantities", 
          description: "Compare vendor options (no sequences needed)" 
        }
      ]
    },
    {
      category: "Plasmid Operations",
      icon: "🧪",
      commands: [
        { 
          command: "insert each of the sequences into pUC19", 
          description: "Insert sequences into plasmid (requires sequences from previous step)" 
        },
        { 
          command: `visualize plasmid pUC19 with insert ${sampleSequences.plasmidInsert}`, 
          description: "Visualize specific plasmid with sample insert sequence" 
        }
      ]
    }
  ];

  const toggleCategory = (category: string) => {
    const newExpanded = new Set(expandedCategories);
    if (newExpanded.has(category)) {
      newExpanded.delete(category);
    } else {
      newExpanded.add(category);
    }
    setExpandedCategories(newExpanded);
  };

  const copyToClipboard = (text: string) => {
    navigator.clipboard.writeText(text);
    // Could add a toast notification here
  };

  return (
    <Card className="mb-4">
      <Card.Body style={{ maxHeight: '500px', overflowY: 'auto', overflowX: 'hidden' }}>
        {exampleCategories.map((category, idx) => (
          <div key={idx} className="mb-3">
            <div
              style={{
                cursor: 'pointer',
                padding: '8px',
                backgroundColor: expandedCategories.has(category.category) ? theme.colors.bgBlueSubtle : theme.colors.bgGoldSubtle,
                borderRadius: '4px',
                display: 'flex',
                justifyContent: 'space-between',
                alignItems: 'center',
                marginBottom: '8px',
                border: expandedCategories.has(category.category) ? `1px solid ${theme.colors.blue}` : `1px solid ${theme.colors.gold}`
              }}
              onClick={() => toggleCategory(category.category)}
            >
              <strong>
                {category.icon} {category.category}
              </strong>
              <span className="text-muted">
                {expandedCategories.has(category.category) ? '▼' : '▶'}
              </span>
            </div>
            {expandedCategories.has(category.category) && (
              <div>
                {category.commands.map((cmd, cmdIdx) => (
                  <div
                    key={cmdIdx}
                    className="mb-2 p-2 border rounded"
                    style={{
                      backgroundColor: '#fff',
                      cursor: 'pointer',
                      transition: 'all 0.2s'
                    }}
                    onMouseEnter={(e) => {
                      e.currentTarget.style.backgroundColor = theme.colors.bgBlueSubtle;
                      e.currentTarget.style.borderColor = theme.colors.blue;
                    }}
                    onMouseLeave={(e) => {
                      e.currentTarget.style.backgroundColor = '#fff';
                      e.currentTarget.style.borderColor = '#dee2e6';
                    }}
                    onClick={() => onCommandSelect(cmd.command)}
                  >
                    <div className="d-flex justify-content-between align-items-start">
                      <div style={{ flex: 1 }}>
                        <div className="fw-semibold" style={{ color: theme.colors.blue, lineHeight: 1.4 }}>
                          {cmd.description || cmd.command}
                        </div>
                      </div>
                      <Button
                        variant="link"
                        size="sm"
                        className="ms-2"
                        onClick={(e) => {
                          e.stopPropagation();
                          copyToClipboard(cmd.command);
                        }}
                        title="Copy to clipboard"
                      >
                        📋
                      </Button>
                    </div>
                  </div>
                ))}
              </div>
            )}
          </div>
        ))}
      </Card.Body>
    </Card>
  );
};

