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
          command: "What are common bioinformatics challenges and how can AI help overcome them?",
          description: "Learn about typical Bioinformatics problems and how AI can help."
        },
        {
          command: "List all tools Helix.AI has access to (registered MCP tools, discovered @tool functions, and local/EC2 CLI tools)",
          description: "Discover available tools and capabilities"
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
          command: "select 2 representative sequences", 
          description: "Choose representative sequences from alignment (requires previous alignment step)" 
        }
      ]
    },
    {
      category: "Mutation & Variants",
      icon: "🧬",
      commands: [
        { 
          command: `create 96 variants Sequence: ${sampleSequences.single}`, 
          description: "Generate 96 sequence variants from a sample sequence" 
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
          command: `visualize the plasmid sequence ${sampleSequences.plasmidInsert}`, 
          description: "Visualize complete plasmid from full sequence" 
        }
      ]
    },
    {
      category: "Plasmid Operations",
      icon: "🧪",
      commands: [
        { 
          command: `insert sequence ${sampleSequences.plasmidInsert} into pUC19 at position 1500`, 
          description: "Insert sequence into plasmid at specific position" 
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

