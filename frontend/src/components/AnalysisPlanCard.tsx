/**
 * AnalysisPlanCard
 *
 * Renders a structured tabular analysis plan returned by the backend
 * analysis_planner when the user submits an analytical request.
 *
 * Shows:
 *  - Goal statement
 *  - Numbered step list with type badges and descriptions
 *  - Expected outputs
 *  - "Approve & Run" button  (calls onApprove)
 *  - "Cancel" link
 */

import React, { useState } from 'react';

// ---------------------------------------------------------------------------
// Types (mirrors backend analysis_planner output)
// ---------------------------------------------------------------------------

export interface AnalysisPlanStep {
  id: number;
  name: string;
  description: string;
  type: 'load' | 'filter' | 'compute' | 'visualize' | 'interpret' | string;
  operation?: string;
}

export interface AnalysisPlan {
  title: string;
  goal: string;
  steps: AnalysisPlanStep[];
  expected_outputs?: string[];
  type?: string;
}

interface AnalysisPlanCardProps {
  plan: AnalysisPlan;
  /** Called when the user clicks "Approve & Run" */
  onApprove: () => void;
  /** Called when the user clicks "Cancel" */
  onCancel?: () => void;
  /** True while the execution is in flight */
  loading?: boolean;
  /** True after the user already approved (prevents double-submit) */
  approved?: boolean;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

const STEP_TYPE_META: Record<string, { label: string; color: string; icon: string }> = {
  load:        { label: 'Load',       color: '#6366f1', icon: '📂' },
  filter:      { label: 'Filter',     color: '#0ea5e9', icon: '🔍' },
  compute:     { label: 'Compute',    color: '#10b981', icon: '⚙️' },
  visualize:   { label: 'Visualize',  color: '#f59e0b', icon: '📊' },
  interpret:   { label: 'Interpret',  color: '#8b5cf6', icon: '🔬' },
};

function stepMeta(type: string) {
  return STEP_TYPE_META[type] ?? { label: type, color: '#64748b', icon: '▶' };
}

// ---------------------------------------------------------------------------
// Component
// ---------------------------------------------------------------------------

export const AnalysisPlanCard: React.FC<AnalysisPlanCardProps> = ({
  plan,
  onApprove,
  onCancel,
  loading = false,
  approved = false,
}) => {
  const [expanded, setExpanded] = useState<Set<number>>(new Set());

  const toggleStep = (id: number) =>
    setExpanded((prev) => {
      const next = new Set(prev);
      next.has(id) ? next.delete(id) : next.add(id);
      return next;
    });

  return (
    <div
      style={{
        background: '#f8fafc',
        border: '1px solid #e2e8f0',
        borderRadius: 10,
        padding: '20px 22px',
        marginTop: 8,
        fontFamily: 'inherit',
      }}
    >
      {/* ── Header ─────────────────────────────────────────────────────── */}
      <div style={{ display: 'flex', alignItems: 'flex-start', gap: 10, marginBottom: 12 }}>
        <span style={{ fontSize: '1.4rem', lineHeight: 1 }}>🧬</span>
        <div>
          <div style={{ fontWeight: 700, fontSize: '1rem', color: '#1e293b' }}>
            {plan.title}
          </div>
          <div style={{ fontSize: '0.85rem', color: '#475569', marginTop: 3 }}>
            {plan.goal}
          </div>
        </div>
      </div>

      {/* ── Steps ──────────────────────────────────────────────────────── */}
      <div style={{ marginBottom: 14 }}>
        {(plan.steps ?? []).map((step) => {
          const meta = stepMeta(step.type);
          const isOpen = expanded.has(step.id);
          const detail = step.operation || step.description;
          return (
            <div
              key={step.id}
              style={{
                display: 'flex',
                gap: 10,
                padding: '8px 0',
                borderBottom: '1px solid #e2e8f0',
              }}
            >
              {/* Step number */}
              <div
                style={{
                  width: 24,
                  height: 24,
                  borderRadius: '50%',
                  background: meta.color,
                  color: '#fff',
                  fontSize: '0.72rem',
                  fontWeight: 700,
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'center',
                  flexShrink: 0,
                  marginTop: 2,
                }}
              >
                {step.id}
              </div>

              {/* Content */}
              <div style={{ flex: 1 }}>
                <div
                  style={{
                    display: 'flex',
                    alignItems: 'center',
                    gap: 6,
                    cursor: detail ? 'pointer' : 'default',
                  }}
                  onClick={() => detail && toggleStep(step.id)}
                >
                  <span style={{ fontSize: '0.88rem', fontWeight: 600, color: '#1e293b' }}>
                    {step.name}
                  </span>
                  <span
                    style={{
                      fontSize: '0.68rem',
                      fontWeight: 600,
                      color: meta.color,
                      background: `${meta.color}18`,
                      border: `1px solid ${meta.color}44`,
                      borderRadius: 4,
                      padding: '1px 5px',
                      textTransform: 'uppercase',
                      letterSpacing: '0.04em',
                    }}
                  >
                    {meta.icon} {meta.label}
                  </span>
                  {detail && (
                    <span style={{ fontSize: '0.72rem', color: '#94a3b8', marginLeft: 'auto' }}>
                      {isOpen ? '▲' : '▼'}
                    </span>
                  )}
                </div>

                {isOpen && detail && (
                  <div
                    style={{
                      marginTop: 4,
                      fontSize: '0.8rem',
                      color: '#64748b',
                      lineHeight: 1.5,
                    }}
                  >
                    {detail}
                  </div>
                )}
              </div>
            </div>
          );
        })}
      </div>

      {/* ── Expected outputs ───────────────────────────────────────────── */}
      {plan.expected_outputs && plan.expected_outputs.length > 0 && (
        <div style={{ marginBottom: 16 }}>
          <div
            style={{
              fontSize: '0.75rem',
              fontWeight: 600,
              color: '#64748b',
              textTransform: 'uppercase',
              letterSpacing: '0.05em',
              marginBottom: 5,
            }}
          >
            Expected outputs
          </div>
          <div style={{ display: 'flex', flexWrap: 'wrap', gap: 6 }}>
            {plan.expected_outputs.map((o, i) => (
              <span
                key={i}
                style={{
                  fontSize: '0.78rem',
                  background: '#e0f2fe',
                  color: '#0369a1',
                  border: '1px solid #bae6fd',
                  borderRadius: 5,
                  padding: '2px 8px',
                }}
              >
                {o}
              </span>
            ))}
          </div>
        </div>
      )}

      {/* ── Actions ────────────────────────────────────────────────────── */}
      <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
        {approved ? (
          <span
            style={{
              fontSize: '0.85rem',
              fontWeight: 600,
              color: '#059669',
            }}
          >
            ✓ Approved — running analysis…
          </span>
        ) : (
          <button
            onClick={onApprove}
            disabled={loading}
            style={{
              background: loading ? '#94a3b8' : '#1d4ed8',
              color: '#fff',
              border: 'none',
              borderRadius: 6,
              padding: '8px 20px',
              fontWeight: 700,
              fontSize: '0.88rem',
              cursor: loading ? 'not-allowed' : 'pointer',
              display: 'flex',
              alignItems: 'center',
              gap: 6,
            }}
          >
            {loading ? (
              <>
                <span
                  style={{
                    display: 'inline-block',
                    width: 12,
                    height: 12,
                    border: '2px solid #fff',
                    borderTopColor: 'transparent',
                    borderRadius: '50%',
                    animation: 'spin 0.8s linear infinite',
                  }}
                />
                Running…
              </>
            ) : (
              '▶ Approve & Run'
            )}
          </button>
        )}

        {!approved && onCancel && (
          <button
            onClick={onCancel}
            style={{
              background: 'transparent',
              color: '#64748b',
              border: 'none',
              padding: '8px 4px',
              fontSize: '0.82rem',
              cursor: 'pointer',
              textDecoration: 'underline',
            }}
          >
            Cancel
          </button>
        )}

        {!approved && (
          <span style={{ fontSize: '0.75rem', color: '#94a3b8', marginLeft: 'auto' }}>
            {plan.steps?.length ?? 0} steps · no code runs until you approve
          </span>
        )}
      </div>

      {/* CSS for spinner */}
      <style>{`
        @keyframes spin { to { transform: rotate(360deg); } }
      `}</style>
    </div>
  );
};

export default AnalysisPlanCard;
