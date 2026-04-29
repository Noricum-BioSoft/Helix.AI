/**
 * NeedsInputsCard
 *
 * Rich card renderer for needs_inputs / WAITING_FOR_INPUTS responses.
 *
 * Replaces the plain markdown dump with a structured card that shows:
 *  - Tool display name and description
 *  - Required inputs as a clean table/list with param names + examples
 *  - Optional inputs (collapsed by default)
 *  - Detected parameters from the user's prompt (non-internal)
 *  - An "Upload files" CTA that wires into the parent's upload handler
 */

import React, { useState } from 'react';

// ---------------------------------------------------------------------------
// Types (mirror _build_needs_inputs_response backend output)
// ---------------------------------------------------------------------------

export interface InputSpec {
  name: string;
  description?: string;
  example?: string;
}

export interface NeedsInputsCardProps {
  /** Full output object from /execute */
  output: Record<string, unknown>;
  /** Called when the user clicks the upload CTA */
  onUpload?: () => void;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function extractNeedsInputsData(output: Record<string, unknown>): {
  displayName: string;
  description: string;
  required: InputSpec[];
  optional: InputSpec[];
  detected: Record<string, string>;
} {
  // Structured data lives in data.results (set by _build_needs_inputs_response)
  const results =
    (output?.data as any)?.results ??
    (output?.result as any)?.result ??
    (output?.result as any) ??
    {};

  const required: InputSpec[] = Array.isArray(results.required_inputs)
    ? results.required_inputs
    : [];
  const optional: InputSpec[] = Array.isArray(results.optional_inputs)
    ? results.optional_inputs
    : [];
  const rawDetected = results.detected_parameters ?? {};
  // Exclude any remaining internal keys that may have slipped through
  const SKIP = new Set(['router_reasoning', 'session_id', 'needs_inputs']);
  const detected: Record<string, string> = {};
  for (const [k, v] of Object.entries(rawDetected)) {
    if (SKIP.has(k) || v === null || v === undefined || v === '') continue;
    detected[k] = String(v);
  }

  // Derive display name: prefer tool_name from results, fall back to output.tool
  const rawTool = (results.tool_name ?? output?.tool ?? '') as string;
  const displayName = rawTool.replace(/_/g, ' ').replace(/\b\w/g, (c) => c.toUpperCase());

  // Description: first non-empty sentence of the markdown text, up to period
  const rawText = ((output?.text ?? results.text ?? '') as string).trim();
  const afterBold = rawText.replace(/^\*\*[^*]+\*\*\n*/, '').trim();
  const description = afterBold.split('\n')[0].trim();

  return { displayName, description, required, optional, detected };
}

// ---------------------------------------------------------------------------
// Sub-components
// ---------------------------------------------------------------------------

const InputRow: React.FC<{ spec: InputSpec; required?: boolean }> = ({ spec, required }) => (
  <div
    style={{
      display: 'flex',
      alignItems: 'flex-start',
      gap: 10,
      padding: '7px 0',
      borderBottom: '1px dashed #f1f5f9',
    }}
  >
    <code
      style={{
        fontSize: '0.78rem',
        background: required ? '#fef3c7' : '#f1f5f9',
        color: required ? '#92400e' : '#475569',
        border: `1px solid ${required ? '#fde68a' : '#e2e8f0'}`,
        borderRadius: 4,
        padding: '2px 6px',
        flexShrink: 0,
        whiteSpace: 'nowrap',
      }}
    >
      {spec.name}
    </code>
    <div style={{ flex: 1 }}>
      {spec.description && (
        <div style={{ fontSize: '0.82rem', color: '#334155', lineHeight: 1.4 }}>
          {spec.description}
        </div>
      )}
      {spec.example && (
        <div style={{ fontSize: '0.75rem', color: '#94a3b8', marginTop: 2 }}>
          e.g. <code style={{ fontSize: '0.73rem' }}>{spec.example}</code>
        </div>
      )}
    </div>
    {required && (
      <span
        style={{
          fontSize: '0.65rem',
          fontWeight: 700,
          color: '#b45309',
          background: '#fef3c7',
          border: '1px solid #fde68a',
          borderRadius: 3,
          padding: '1px 5px',
          flexShrink: 0,
          textTransform: 'uppercase',
          letterSpacing: '0.04em',
        }}
      >
        required
      </span>
    )}
  </div>
);

// ---------------------------------------------------------------------------
// Main component
// ---------------------------------------------------------------------------

export const NeedsInputsCard: React.FC<NeedsInputsCardProps> = ({ output, onUpload }) => {
  const [showOptional, setShowOptional] = useState(false);
  const { displayName, description, required, optional, detected } =
    extractNeedsInputsData(output);

  const detectedEntries = Object.entries(detected);

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
      {/* ── Header ────────────────────────────────────────────────────────── */}
      <div style={{ display: 'flex', alignItems: 'flex-start', gap: 10, marginBottom: 14 }}>
        <span style={{ fontSize: '1.3rem', lineHeight: 1, flexShrink: 0 }}>📂</span>
        <div style={{ flex: 1 }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: 8, flexWrap: 'wrap' }}>
            <span style={{ fontWeight: 700, fontSize: '0.97rem', color: '#1e293b' }}>
              {displayName}
            </span>
            <span
              style={{
                fontSize: '0.68rem',
                fontWeight: 600,
                color: '#0369a1',
                background: '#e0f2fe',
                border: '1px solid #bae6fd',
                borderRadius: 4,
                padding: '1px 6px',
                textTransform: 'uppercase',
                letterSpacing: '0.04em',
              }}
            >
              Needs inputs
            </span>
          </div>
          {description && (
            <div style={{ fontSize: '0.83rem', color: '#475569', marginTop: 3, lineHeight: 1.45 }}>
              {description}
            </div>
          )}
        </div>
      </div>

      {/* ── Required inputs ───────────────────────────────────────────────── */}
      {required.length > 0 && (
        <div style={{ marginBottom: 12 }}>
          <div
            style={{
              fontSize: '0.72rem',
              fontWeight: 700,
              color: '#64748b',
              textTransform: 'uppercase',
              letterSpacing: '0.06em',
              marginBottom: 6,
            }}
          >
            Required inputs
          </div>
          {required.map((spec, i) => (
            <InputRow key={i} spec={spec} required />
          ))}
        </div>
      )}

      {/* ── Optional inputs (collapsed) ───────────────────────────────────── */}
      {optional.length > 0 && (
        <div style={{ marginBottom: 12 }}>
          <button
            onClick={() => setShowOptional((p) => !p)}
            style={{
              background: 'none',
              border: 'none',
              padding: 0,
              cursor: 'pointer',
              fontSize: '0.72rem',
              fontWeight: 700,
              color: '#64748b',
              textTransform: 'uppercase',
              letterSpacing: '0.06em',
              display: 'flex',
              alignItems: 'center',
              gap: 4,
              marginBottom: showOptional ? 6 : 0,
            }}
          >
            {showOptional ? '▾' : '▸'} Optional inputs ({optional.length})
          </button>
          {showOptional &&
            optional.map((spec, i) => <InputRow key={i} spec={spec} />)}
        </div>
      )}

      {/* ── Detected parameters ───────────────────────────────────────────── */}
      {detectedEntries.length > 0 && (
        <div
          style={{
            background: '#f0fdf4',
            border: '1px solid #bbf7d0',
            borderRadius: 7,
            padding: '8px 12px',
            marginBottom: 14,
            fontSize: '0.8rem',
          }}
        >
          <div
            style={{
              fontWeight: 700,
              color: '#166534',
              fontSize: '0.72rem',
              textTransform: 'uppercase',
              letterSpacing: '0.06em',
              marginBottom: 4,
            }}
          >
            ✓ Detected from your prompt
          </div>
          {detectedEntries.map(([k, v]) => (
            <div key={k} style={{ color: '#15803d', lineHeight: 1.6 }}>
              <code style={{ fontSize: '0.75rem' }}>{k}</code>{' '}
              <span style={{ color: '#475569' }}>→</span>{' '}
              <code style={{ fontSize: '0.75rem', color: '#166534' }}>{v}</code>
            </div>
          ))}
        </div>
      )}

      {/* ── CTA ───────────────────────────────────────────────────────────── */}
      <div style={{ display: 'flex', alignItems: 'center', gap: 12, flexWrap: 'wrap' }}>
        {onUpload && (
          <button
            onClick={onUpload}
            style={{
              background: '#1d4ed8',
              color: '#fff',
              border: 'none',
              borderRadius: 7,
              padding: '8px 18px',
              fontWeight: 700,
              fontSize: '0.88rem',
              cursor: 'pointer',
              display: 'flex',
              alignItems: 'center',
              gap: 6,
            }}
          >
            📁 Upload files
          </button>
        )}
        <span style={{ fontSize: '0.73rem', color: '#94a3b8' }}>
          Provide {required.length > 0 ? `${required.length} required input${required.length !== 1 ? 's' : ''}` : 'inputs'}{' '}
          — analysis runs immediately after upload
        </span>
      </div>
    </div>
  );
};

export default NeedsInputsCard;
