import React, { useState } from 'react';
import Modal from 'react-bootstrap/Modal';
import Button from 'react-bootstrap/Button';
import Badge from 'react-bootstrap/Badge';
import { demoScenarios, DemoScenario, ExpectedBehavior } from '../data/demoScenarios';
import { theme } from '../theme';

interface DemoScenariosPanelProps {
  show: boolean;
  onHide: () => void;
  onSelect: (prompt: string, scenarioId?: string) => void;
}

// ── Behaviour badge config ─────────────────────────────────────────────────────
const BEHAVIOR_CONFIG: Record<
  ExpectedBehavior,
  { label: string; bg: string; color: string; icon: string; description: string }
> = {
  needs_inputs: {
    label: 'Requests Data',
    bg: '#FFF3CD',
    color: '#856404',
    icon: '📋',
    description: 'Helix identifies the analysis and lists exactly what data files are needed.',
  },
  executes_pipeline: {
    label: 'Executes Pipeline',
    bg: '#D1FAE5',
    color: '#065F46',
    icon: '▶',
    description: 'S3 paths are included — Helix runs the full pipeline immediately.',
  },
  multi_step_plan: {
    label: 'Multi-Step Plan',
    bg: '#EDE9FE',
    color: '#5B21B6',
    icon: '🔗',
    description: 'Helix orchestrates multiple tools in sequence — fetch → align → tree.',
  },
};

// ── Output type icons ─────────────────────────────────────────────────────────
const OUTPUT_ICONS: Record<string, string> = {
  plot:   '📊',
  csv:    '📄',
  fasta:  '🧬',
  report: '📋',
  newick: '🌳',
};

// ── Single scenario card ───────────────────────────────────────────────────────
const ScenarioCard: React.FC<{
  scenario: DemoScenario;
  isSelected: boolean;
  onPreview: () => void;
  onLoad: () => void;
}> = ({ scenario, isSelected, onPreview, onLoad }) => {
  const behavior = BEHAVIOR_CONFIG[scenario.expectedBehavior];

  return (
    <div
      onClick={onPreview}
      style={{
        border: isSelected
          ? `2px solid ${scenario.domainColor}`
          : '1px solid #E2E8F0',
        borderRadius: '12px',
        padding: '20px',
        cursor: 'pointer',
        backgroundColor: isSelected ? `${scenario.domainColor}08` : '#FFFFFF',
        transition: 'all 0.18s ease',
        display: 'flex',
        flexDirection: 'column',
        gap: '12px',
        boxShadow: isSelected
          ? `0 4px 16px ${scenario.domainColor}30`
          : '0 1px 4px rgba(0,0,0,0.06)',
        height: '100%',
      }}
      onMouseEnter={(e) => {
        if (!isSelected) {
          e.currentTarget.style.borderColor = scenario.domainColor;
          e.currentTarget.style.boxShadow = `0 4px 12px ${scenario.domainColor}25`;
        }
      }}
      onMouseLeave={(e) => {
        if (!isSelected) {
          e.currentTarget.style.borderColor = '#E2E8F0';
          e.currentTarget.style.boxShadow = '0 1px 4px rgba(0,0,0,0.06)';
        }
      }}
    >
      {/* Header row */}
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
        <span style={{ fontSize: '2rem', lineHeight: 1 }}>{scenario.icon}</span>
        <span
          style={{
            fontSize: '0.7rem',
            fontWeight: 700,
            letterSpacing: '0.04em',
            padding: '3px 8px',
            borderRadius: '9999px',
            backgroundColor: `${scenario.domainColor}18`,
            color: scenario.domainColor,
            whiteSpace: 'nowrap',
          }}
        >
          {scenario.domain}
        </span>
      </div>

      {/* Title + subtitle */}
      <div>
        <div style={{ fontWeight: 700, fontSize: '0.95rem', color: '#1A202C', lineHeight: 1.3 }}>
          {scenario.title}
        </div>
        <div style={{ fontSize: '0.78rem', color: '#718096', marginTop: '4px', lineHeight: 1.4 }}>
          {scenario.subtitle}
        </div>
      </div>

      {/* Behaviour badge */}
      <div
        style={{
          display: 'inline-flex',
          alignItems: 'center',
          gap: '5px',
          padding: '4px 10px',
          borderRadius: '6px',
          backgroundColor: behavior.bg,
          color: behavior.color,
          fontSize: '0.72rem',
          fontWeight: 700,
          width: 'fit-content',
        }}
      >
        <span>{behavior.icon}</span>
        <span>{behavior.label}</span>
      </div>

      {/* Tags */}
      <div style={{ display: 'flex', flexWrap: 'wrap', gap: '4px' }}>
        {scenario.tags.slice(0, 4).map((tag) => (
          <span
            key={tag}
            style={{
              fontSize: '0.65rem',
              padding: '2px 7px',
              borderRadius: '9999px',
              backgroundColor: '#F1F5F9',
              color: '#475569',
              border: '1px solid #E2E8F0',
            }}
          >
            {tag}
          </span>
        ))}
      </div>

      {/* Inputs */}
      {scenario.inputs && scenario.inputs.length > 0 && (
        <div>
          <div style={{ fontSize: '0.68rem', color: '#94A3B8', marginBottom: '4px', fontWeight: 600, letterSpacing: '0.04em' }}>
            INPUTS
          </div>
          <div style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
            {scenario.inputs.slice(0, 3).map((inp) => (
              <span key={inp.label} style={{ fontSize: '0.68rem', color: '#64748B', lineHeight: 1.3 }}>
                • {inp.label}
              </span>
            ))}
            {scenario.inputs.length > 3 && (
              <span style={{ fontSize: '0.65rem', color: '#94A3B8' }}>+{scenario.inputs.length - 3} more</span>
            )}
          </div>
        </div>
      )}

      {/* Outputs row */}
      <div style={{ marginTop: 'auto' }}>
        <div style={{ fontSize: '0.68rem', color: '#94A3B8', marginBottom: '4px', fontWeight: 600, letterSpacing: '0.04em' }}>
          OUTPUTS
        </div>
        <div style={{ display: 'flex', flexWrap: 'wrap', gap: '4px' }}>
          {scenario.outputs.map((o) => (
            <span key={o.label} style={{ fontSize: '0.7rem', color: '#64748B' }}>
              {OUTPUT_ICONS[o.type]} {o.label}
              {scenario.outputs.indexOf(o) < scenario.outputs.length - 1 ? '' : ''}
            </span>
          ))}
        </div>
      </div>

      {/* Load button */}
      <button
        onClick={(e) => { e.stopPropagation(); onLoad(); }}
        style={{
          marginTop: '4px',
          padding: '8px 0',
          borderRadius: '8px',
          border: 'none',
          backgroundColor: scenario.domainColor,
          color: '#FFFFFF',
          fontSize: '0.8rem',
          fontWeight: 700,
          cursor: 'pointer',
          width: '100%',
          transition: 'opacity 0.15s',
          letterSpacing: '0.02em',
        }}
        onMouseEnter={(e) => (e.currentTarget.style.opacity = '0.88')}
        onMouseLeave={(e) => (e.currentTarget.style.opacity = '1')}
      >
        Try this scenario →
      </button>
    </div>
  );
};

// ── Main modal ─────────────────────────────────────────────────────────────────
export const DemoScenariosPanel: React.FC<DemoScenariosPanelProps> = ({
  show,
  onHide,
  onSelect,
}) => {
  const [selectedId, setSelectedId] = useState<string>(demoScenarios[0].id);

  const handleLoad = (scenario: DemoScenario) => {
    onSelect(scenario.prompt, scenario.id);
    onHide();
  };

  return (
    <Modal
      show={show}
      onHide={onHide}
      size="xl"
      centered
      dialogClassName="demo-scenarios-modal"
    >
      <Modal.Header
        closeButton
        style={{
          background: 'linear-gradient(135deg, #1A202C, #2D3748)',
          borderBottom: '1px solid #2D3748',
          padding: '20px 28px',
        }}
      >
        <div style={{ display: 'flex', alignItems: 'center', gap: '14px' }}>
          <div
            style={{
              width: '40px',
              height: '40px',
              borderRadius: '10px',
              background: 'linear-gradient(135deg, #3A60A8, #7B3FA8)',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              fontSize: '1.2rem',
            }}
          >
            🧬
          </div>
          <div>
            <div style={{ color: '#FFFFFF', fontWeight: 800, fontSize: '1.15rem', lineHeight: 1.2 }}>
              Helix.AI Demo Scenarios
            </div>
            <div style={{ color: '#94A3B8', fontSize: '0.8rem', marginTop: '2px' }}>
              5 real bioinformatics workflows — click any scenario to load it
            </div>
          </div>
          <Badge
            style={{
              marginLeft: '8px',
              backgroundColor: '#3A60A8',
              color: '#FFFFFF',
              fontSize: '0.7rem',
              padding: '4px 10px',
              borderRadius: '9999px',
            }}
          >
            5 scenarios
          </Badge>
        </div>
      </Modal.Header>

      <Modal.Body style={{ padding: 0, backgroundColor: '#F8FAFC' }}>
        <div style={{ padding: '20px', overflowX: 'auto', overflowY: 'auto' }}>
            {/* Legend — full descriptions, no truncation */}
            <div
              style={{
                display: 'flex',
                gap: '16px',
                flexWrap: 'wrap',
                marginBottom: '16px',
                padding: '12px 14px',
                borderRadius: '8px',
                backgroundColor: '#FFFFFF',
                border: '1px solid #E2E8F0',
              }}
            >
              {(Object.entries(BEHAVIOR_CONFIG) as [ExpectedBehavior, typeof BEHAVIOR_CONFIG[ExpectedBehavior]][]).map(
                ([key, cfg]) => (
                  <div
                    key={key}
                    style={{
                      display: 'flex',
                      alignItems: 'flex-start',
                      gap: '8px',
                      fontSize: '0.72rem',
                      flex: '1 1 200px',
                      minWidth: 0,
                    }}
                  >
                    <span
                      style={{
                        flexShrink: 0,
                        display: 'inline-block',
                        width: '10px',
                        height: '10px',
                        borderRadius: '3px',
                        marginTop: '3px',
                        backgroundColor: cfg.bg,
                        border: `1px solid ${cfg.color}50`,
                      }}
                    />
                    <div style={{ minWidth: 0 }}>
                      <span style={{ color: cfg.color, fontWeight: 700 }}>{cfg.label}</span>
                      <span style={{ color: '#64748B', display: 'block', marginTop: '2px', lineHeight: 1.4 }}>
                        {cfg.description}
                      </span>
                    </div>
                  </div>
                )
              )}
            </div>

            {/* Cards — all 5 in one row; minWidth so they don't collapse on narrow modals */}
            <div
              style={{
                display: 'grid',
                gridTemplateColumns: 'repeat(5, minmax(168px, 1fr))',
                gap: '12px',
                minWidth: '880px',
              }}
            >
              {demoScenarios.map((scenario) => (
                <ScenarioCard
                  key={scenario.id}
                  scenario={scenario}
                  isSelected={scenario.id === selectedId}
                  onPreview={() => setSelectedId(scenario.id)}
                  onLoad={() => handleLoad(scenario)}
                />
              ))}
            </div>
        </div>
      </Modal.Body>

      <Modal.Footer
        style={{
          backgroundColor: '#F8FAFC',
          borderTop: '1px solid #E2E8F0',
          padding: '12px 24px',
          display: 'flex',
          justifyContent: 'space-between',
          alignItems: 'center',
        }}
      >
        <span style={{ fontSize: '0.78rem', color: '#94A3B8' }}>
          "Try this scenario →" loads the prompt into Helix
        </span>
        <Button variant="outline-secondary" size="sm" onClick={onHide}>
          Close
        </Button>
      </Modal.Footer>
    </Modal>
  );
};

export default DemoScenariosPanel;
