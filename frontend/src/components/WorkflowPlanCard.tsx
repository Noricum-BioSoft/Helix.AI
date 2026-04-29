/**
 * WorkflowPlanCard
 *
 * Rich card renderer for __plan__ / workflow_planned responses.
 *
 * Extracts the structured step list from the router_reasoning payload embedded
 * in `data.results.steps[0].arguments.router_reasoning.suggested_steps` and
 * renders a clean numbered step list with an approve/execute CTA.
 */

import React from 'react';

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

interface RouterReasoning {
  problem_description?: string;
  intent?: string;
  suggested_steps?: string[];
}

interface PlanStep {
  id?: string;
  tool_name?: string;
  description?: string;
  arguments?: {
    router_reasoning?: RouterReasoning;
    command?: string;
  };
}

interface PlanResults {
  steps?: PlanStep[];
  execute_ready?: boolean;
  approval_required?: boolean;
}

export interface WorkflowPlanCardProps {
  /** Full output object from the backend /execute response */
  output: Record<string, unknown>;
  /** "approve" = show "✅ I approve"; "execute" = show "▶ Execute Pipeline" */
  planAction: 'approve' | 'execute' | 'none';
  alreadyApproved: boolean;
  loading: boolean;
  onApprove: () => void;
  onExecute: () => void;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

const STEP_COLORS = [
  '#6366f1', '#0ea5e9', '#10b981', '#f59e0b',
  '#ef4444', '#8b5cf6', '#ec4899', '#14b8a6',
];

function extractPlanData(output: Record<string, unknown>): {
  title: string;
  steps: string[];
  intent: string;
} {
  const results =
    (output?.data as any)?.results as PlanResults | undefined ??
    (output?.result as any)?.result as PlanResults | undefined ??
    (output?.result as any) as PlanResults | undefined;

  const firstStep = results?.steps?.[0];
  const reasoning = firstStep?.arguments?.router_reasoning;

  const problem =
    reasoning?.problem_description ||
    firstStep?.description ||
    (output?.prompt as string) ||
    'Workflow plan';

  const suggested = reasoning?.suggested_steps ?? [];
  const intent = reasoning?.intent ?? '';

  // If no suggested_steps, surface the tool + description as a single step
  const steps =
    suggested.length > 0
      ? suggested
      : firstStep
        ? [`${firstStep.tool_name ?? 'execute'}: ${firstStep.description ?? problem}`]
        : [problem];

  return { title: problem, steps, intent };
}

// ---------------------------------------------------------------------------
// Component
// ---------------------------------------------------------------------------

export const WorkflowPlanCard: React.FC<WorkflowPlanCardProps> = ({
  output,
  planAction,
  alreadyApproved,
  loading,
  onApprove,
  onExecute,
}) => {
  const { title, steps, intent } = extractPlanData(output);
  const isMultiStep = intent === 'multi-step' || steps.length > 1;

  return (
    <div
      style={{
        background: '#f8fafc',
        border: '1.5px solid #e2e8f0',
        borderRadius: 12,
        padding: '18px 22px 16px',
        marginTop: 6,
        fontFamily: 'inherit',
      }}
    >
      {/* ── Header ──────────────────────────────────────────────────────────── */}
      <div style={{ display: 'flex', alignItems: 'flex-start', gap: 10, marginBottom: 14 }}>
        <span style={{ fontSize: '1.3rem', lineHeight: 1, flexShrink: 0 }}>🧬</span>
        <div style={{ flex: 1 }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: 8, flexWrap: 'wrap' }}>
            <span style={{ fontWeight: 700, fontSize: '0.97rem', color: '#1e293b' }}>
              Proposed Workflow
            </span>
            {isMultiStep && (
              <span
                style={{
                  fontSize: '0.68rem',
                  fontWeight: 600,
                  color: '#6366f1',
                  background: '#eef2ff',
                  border: '1px solid #c7d2fe',
                  borderRadius: 4,
                  padding: '1px 6px',
                  textTransform: 'uppercase',
                  letterSpacing: '0.04em',
                }}
              >
                Multi-step
              </span>
            )}
          </div>
          <div
            style={{
              fontSize: '0.83rem',
              color: '#475569',
              marginTop: 3,
              lineHeight: 1.45,
            }}
          >
            {title}
          </div>
        </div>
      </div>

      {/* ── Step list ───────────────────────────────────────────────────────── */}
      <div
        style={{
          marginBottom: 16,
          borderLeft: '2px solid #e2e8f0',
          paddingLeft: 14,
        }}
      >
        {steps.map((step, i) => (
          <div
            key={i}
            style={{
              display: 'flex',
              alignItems: 'flex-start',
              gap: 10,
              padding: '6px 0',
              borderBottom: i < steps.length - 1 ? '1px dashed #f1f5f9' : 'none',
            }}
          >
            {/* Numbered circle */}
            <div
              style={{
                width: 22,
                height: 22,
                borderRadius: '50%',
                background: STEP_COLORS[i % STEP_COLORS.length],
                color: '#fff',
                fontSize: '0.7rem',
                fontWeight: 700,
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                flexShrink: 0,
                marginTop: 1,
              }}
            >
              {i + 1}
            </div>
            <span
              style={{
                fontSize: '0.86rem',
                color: '#334155',
                lineHeight: 1.4,
                paddingTop: 2,
              }}
            >
              {step}
            </span>
          </div>
        ))}
      </div>

      {/* ── Actions ─────────────────────────────────────────────────────────── */}
      <div style={{ display: 'flex', alignItems: 'center', gap: 12, flexWrap: 'wrap' }}>
        {alreadyApproved ? (
          <span style={{ fontSize: '0.85rem', fontWeight: 600, color: '#059669' }}>
            ✓ Approved — workflow queued
          </span>
        ) : (
          <>
            {planAction === 'approve' && (
              <button
                onClick={onApprove}
                disabled={loading}
                style={{
                  background: loading ? '#94a3b8' : '#1d4ed8',
                  color: '#fff',
                  border: 'none',
                  borderRadius: 7,
                  padding: '8px 20px',
                  fontWeight: 700,
                  fontSize: '0.88rem',
                  cursor: loading ? 'not-allowed' : 'pointer',
                  display: 'flex',
                  alignItems: 'center',
                  gap: 6,
                  transition: 'background 0.15s',
                }}
              >
                {loading ? (
                  <>
                    <span
                      style={{
                        display: 'inline-block',
                        width: 11,
                        height: 11,
                        border: '2px solid #fff',
                        borderTopColor: 'transparent',
                        borderRadius: '50%',
                        animation: 'wpc-spin 0.8s linear infinite',
                      }}
                    />
                    Approving…
                  </>
                ) : (
                  '✅ I approve'
                )}
              </button>
            )}
            {planAction === 'execute' && (
              <button
                onClick={onExecute}
                disabled={loading}
                style={{
                  background: loading ? '#94a3b8' : '#059669',
                  color: '#fff',
                  border: 'none',
                  borderRadius: 7,
                  padding: '8px 20px',
                  fontWeight: 700,
                  fontSize: '0.88rem',
                  cursor: loading ? 'not-allowed' : 'pointer',
                }}
              >
                {loading ? 'Submitting…' : '▶ Execute Pipeline'}
              </button>
            )}
          </>
        )}
        <span
          style={{
            fontSize: '0.73rem',
            color: '#94a3b8',
            marginLeft: alreadyApproved ? 0 : 'auto',
          }}
        >
          {alreadyApproved
            ? ''
            : `${steps.length} step${steps.length !== 1 ? 's' : ''} · no code runs until you approve`}
        </span>
      </div>

      <style>{`@keyframes wpc-spin { to { transform: rotate(360deg); } }`}</style>
    </div>
  );
};

export default WorkflowPlanCard;
