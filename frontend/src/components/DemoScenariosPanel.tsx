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

// ── Preview pane ───────────────────────────────────────────────────────────────
const ScenarioPreview: React.FC<{
  scenario: DemoScenario;
  onLoad: () => void;
}> = ({ scenario, onLoad }) => {
  const behavior = BEHAVIOR_CONFIG[scenario.expectedBehavior];

  return (
    <div style={{ height: '100%', display: 'flex', flexDirection: 'column', gap: '16px' }}>
      {/* Scenario header */}
      <div
        style={{
          padding: '20px',
          borderRadius: '12px',
          background: `linear-gradient(135deg, ${scenario.domainColor}15, ${scenario.domainColor}05)`,
          border: `1px solid ${scenario.domainColor}30`,
        }}
      >
        <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '8px' }}>
          <span style={{ fontSize: '2rem' }}>{scenario.icon}</span>
          <div>
            <div style={{ fontWeight: 800, fontSize: '1.05rem', color: '#1A202C' }}>
              {scenario.title}
            </div>
            <div style={{ fontSize: '0.82rem', color: '#718096' }}>{scenario.subtitle}</div>
          </div>
        </div>

        <div style={{ display: 'flex', flexWrap: 'wrap', gap: '6px', marginTop: '10px' }}>
          <span
            style={{
              padding: '4px 10px',
              borderRadius: '6px',
              backgroundColor: behavior.bg,
              color: behavior.color,
              fontSize: '0.75rem',
              fontWeight: 700,
            }}
          >
            {behavior.icon} {behavior.label}
          </span>
          <span
            style={{
              padding: '4px 10px',
              borderRadius: '6px',
              backgroundColor: '#F1F5F9',
              color: '#475569',
              fontSize: '0.75rem',
            }}
          >
            ⏱ {scenario.estimatedRuntime}
          </span>
          <span
            style={{
              padding: '4px 10px',
              borderRadius: '6px',
              backgroundColor: '#F1F5F9',
              color: '#475569',
              fontSize: '0.75rem',
            }}
          >
            🔧 {scenario.tool}
          </span>
        </div>
      </div>

      {/* What Helix will do */}
      <div
        style={{
          padding: '14px 16px',
          borderRadius: '8px',
          backgroundColor: behavior.bg,
          border: `1px solid ${
            scenario.expectedBehavior === 'needs_inputs' ? '#FDE68A' :
            scenario.expectedBehavior === 'executes_pipeline' ? '#A7F3D0' : '#DDD6FE'
          }`,
        }}
      >
        <div style={{ fontSize: '0.7rem', fontWeight: 700, letterSpacing: '0.06em', color: behavior.color, marginBottom: '4px' }}>
          WHAT HELIX DOES WITH THIS PROMPT
        </div>
        <div style={{ fontSize: '0.82rem', color: behavior.color }}>{behavior.description}</div>
      </div>

      {/* Expected outputs */}
      <div>
        <div style={{ fontSize: '0.7rem', fontWeight: 700, letterSpacing: '0.06em', color: '#94A3B8', marginBottom: '8px' }}>
          EXPECTED OUTPUTS
        </div>
        <div style={{ display: 'flex', flexDirection: 'column', gap: '6px' }}>
          {scenario.outputs.map((o) => (
            <div
              key={o.label}
              style={{
                display: 'flex',
                alignItems: 'center',
                gap: '8px',
                padding: '7px 10px',
                borderRadius: '6px',
                backgroundColor: '#F8FAFC',
                border: '1px solid #E2E8F0',
                fontSize: '0.82rem',
                color: '#374151',
              }}
            >
              <span>{OUTPUT_ICONS[o.type]}</span>
              <span>{o.label}</span>
            </div>
          ))}
        </div>
      </div>

      {/* Prompt preview */}
      <div style={{ flex: 1, minHeight: 0 }}>
        <div style={{ fontSize: '0.7rem', fontWeight: 700, letterSpacing: '0.06em', color: '#94A3B8', marginBottom: '8px' }}>
          PROMPT PREVIEW
        </div>
        <div
          style={{
            padding: '12px',
            borderRadius: '8px',
            backgroundColor: '#0F172A',
            color: '#E2E8F0',
            fontSize: '0.72rem',
            fontFamily: 'Consolas, "Liberation Mono", Menlo, Courier, monospace',
            whiteSpace: 'pre-wrap',
            overflowY: 'auto',
            maxHeight: '220px',
            lineHeight: 1.6,
            border: '1px solid #1E293B',
          }}
        >
          {scenario.prompt}
        </div>
      </div>

      {/* CTA */}
      <button
        onClick={onLoad}
        style={{
          padding: '12px 0',
          borderRadius: '10px',
          border: 'none',
          backgroundColor: scenario.domainColor,
          color: '#FFFFFF',
          fontSize: '0.9rem',
          fontWeight: 700,
          cursor: 'pointer',
          width: '100%',
          letterSpacing: '0.02em',
          transition: 'opacity 0.15s',
        }}
        onMouseEnter={(e) => (e.currentTarget.style.opacity = '0.88')}
        onMouseLeave={(e) => (e.currentTarget.style.opacity = '1')}
      >
        Load into Helix ↗
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
  const selected = demoScenarios.find((s) => s.id === selectedId)!;

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
        <div style={{ display: 'flex', height: '580px' }}>

          {/* Left: scenario grid */}
          <div
            style={{
              flex: '0 0 62%',
              padding: '20px',
              overflowY: 'auto',
              borderRight: '1px solid #E2E8F0',
            }}
          >
            {/* Legend */}
            <div
              style={{
                display: 'flex',
                gap: '12px',
                flexWrap: 'wrap',
                marginBottom: '16px',
                padding: '10px 14px',
                borderRadius: '8px',
                backgroundColor: '#FFFFFF',
                border: '1px solid #E2E8F0',
              }}
            >
              {(Object.entries(BEHAVIOR_CONFIG) as [ExpectedBehavior, typeof BEHAVIOR_CONFIG[ExpectedBehavior]][]).map(
                ([key, cfg]) => (
                  <div key={key} style={{ display: 'flex', alignItems: 'center', gap: '5px', fontSize: '0.72rem' }}>
                    <span
                      style={{
                        display: 'inline-block',
                        width: '10px',
                        height: '10px',
                        borderRadius: '3px',
                        backgroundColor: cfg.bg,
                        border: `1px solid ${cfg.color}50`,
                      }}
                    />
                    <span style={{ color: cfg.color, fontWeight: 700 }}>{cfg.label}</span>
                    <span style={{ color: '#94A3B8' }}>— {cfg.description.split(' ').slice(0, 8).join(' ')}…</span>
                  </div>
                )
              )}
            </div>

            {/* Cards grid */}
            <div
              style={{
                display: 'grid',
                gridTemplateColumns: 'repeat(2, 1fr)',
                gap: '14px',
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

          {/* Right: preview pane */}
          <div
            style={{
              flex: '0 0 38%',
              padding: '20px',
              overflowY: 'auto',
              backgroundColor: '#FFFFFF',
            }}
          >
            <ScenarioPreview
              scenario={selected}
              onLoad={() => handleLoad(selected)}
            />
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
          Click any card to preview · "Try this scenario →" loads the prompt directly
        </span>
        <Button variant="outline-secondary" size="sm" onClick={onHide}>
          Close
        </Button>
      </Modal.Footer>
    </Modal>
  );
};

export default DemoScenariosPanel;
