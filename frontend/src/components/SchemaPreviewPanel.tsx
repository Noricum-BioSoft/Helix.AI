/**
 * SchemaPreviewPanel
 *
 * Shows the upload-time profile returned by the backend's file_intelligence
 * profiler.  Handles all supported format families:
 *   tabular      → column list with dtype badges + numeric stats
 *   sequence     → read/record count + length distribution
 *   variant      → variant count + chromosome breakdown
 *   genomic_interval → interval count + feature types
 *   single_cell  → n_cells × n_genes + embeddings
 *   alignment    → read count + mapped %
 *   unknown      → graceful fallback
 *
 * Also renders clickable suggested-question chips that pre-fill the parent's
 * command textarea via the onSuggest callback.
 */

import React, { useState } from 'react';
import type { SchemaPreview, ColumnStat } from '../services/helixApi';

interface SchemaPreviewPanelProps {
  filename: string;
  preview: SchemaPreview;
  /** Called when the user clicks a suggested question chip. */
  onSuggest?: (question: string) => void;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

const FAMILY_LABELS: Record<string, string> = {
  tabular: '📊 Tabular',
  sequence: '🧬 Sequence',
  variant: '🔬 Variant',
  genomic_interval: '📍 Genomic Intervals',
  single_cell: '🔵 Single Cell',
  alignment: '🗺️ Alignment',
  unknown: '📄 File',
};

const DTYPE_COLORS: Record<string, string> = {
  object: '#6c757d',
  int64: '#0d6efd',
  int32: '#0d6efd',
  float64: '#198754',
  float32: '#198754',
  bool: '#fd7e14',
  datetime64: '#6f42c1',
};

function dtypeColor(dtype: string): string {
  const base = dtype.replace(/\[.*\]/, '').trim();
  return DTYPE_COLORS[base] ?? '#6c757d';
}

function formatCount(n: number | null | undefined): string {
  if (n == null) return '—';
  return n.toLocaleString();
}

function summaryStat(label: string, value: unknown): React.ReactNode {
  if (value == null) return null;
  return (
    <span className="me-3 text-muted" style={{ fontSize: '0.8rem' }}>
      <strong>{label}:</strong> {String(value)}
    </span>
  );
}

// ---------------------------------------------------------------------------
// Sub-renderers per family
// ---------------------------------------------------------------------------

function TabularBody({ preview }: { preview: SchemaPreview }) {
  const cols: ColumnStat[] = preview.schema?.columns ?? [];
  const summary = preview.summary as Record<string, unknown> | undefined;
  const sheets = preview.available_sheets;
  const [expanded, setExpanded] = useState(false);
  const visible = expanded ? cols : cols.slice(0, 6);

  return (
    <>
      <div className="mb-2">
        {summaryStat('Rows', summary?.n_rows)}
        {summaryStat('Cols', summary?.n_cols)}
        {summary?.source_sheet && summaryStat('Sheet', summary.source_sheet as string)}
        {sheets && sheets.length > 1 && (
          <span className="badge bg-secondary me-2" style={{ fontSize: '0.75rem' }}>
            {sheets.length} sheets
          </span>
        )}
      </div>

      {cols.length > 0 && (
        <div style={{ overflowX: 'auto' }}>
          <table
            className="table table-sm table-borderless mb-1"
            style={{ fontSize: '0.78rem', minWidth: 280 }}
          >
            <thead>
              <tr>
                <th className="text-muted fw-normal py-0">Column</th>
                <th className="text-muted fw-normal py-0">Type</th>
                <th className="text-muted fw-normal py-0">Range / Top values</th>
              </tr>
            </thead>
            <tbody>
              {visible.map((col) => (
                <tr key={col.name}>
                  <td className="py-0 text-truncate" style={{ maxWidth: 140 }}>
                    <code style={{ fontSize: '0.78rem' }}>{col.name}</code>
                  </td>
                  <td className="py-0">
                    <span
                      className="badge"
                      style={{
                        background: dtypeColor(col.dtype),
                        color: '#fff',
                        fontSize: '0.7rem',
                      }}
                    >
                      {col.dtype}
                    </span>
                  </td>
                  <td className="py-0 text-muted" style={{ fontSize: '0.75rem' }}>
                    {col.min != null && col.max != null
                      ? `${col.min} – ${col.max}`
                      : col.top_values
                      ? Object.keys(col.top_values).slice(0, 3).join(', ')
                      : '—'}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
          {cols.length > 6 && (
            <button
              className="btn btn-link btn-sm p-0"
              style={{ fontSize: '0.75rem' }}
              onClick={() => setExpanded((e) => !e)}
            >
              {expanded ? '▲ Show less' : `▼ +${cols.length - 6} more columns`}
            </button>
          )}
        </div>
      )}
    </>
  );
}

function SequenceBody({ preview }: { preview: SchemaPreview }) {
  const summary = preview.summary as Record<string, unknown> | undefined;
  return (
    <div>
      {summaryStat('Records', preview.n_records)}
      {summaryStat('Min length', summary?.length_min as number)}
      {summaryStat('Max length', summary?.length_max as number)}
      {summaryStat('Mean length', summary?.length_mean as number)}
    </div>
  );
}

function VariantBody({ preview }: { preview: SchemaPreview }) {
  const summary = preview.summary as Record<string, unknown> | undefined;
  const chroms = summary?.chromosomes as string[] | undefined;
  return (
    <div>
      {summaryStat('Variants', preview.n_records)}
      {chroms && chroms.length > 0 && (
        <span className="text-muted" style={{ fontSize: '0.8rem' }}>
          <strong>Chromosomes:</strong> {chroms.slice(0, 6).join(', ')}
          {chroms.length > 6 ? ` +${chroms.length - 6} more` : ''}
        </span>
      )}
    </div>
  );
}

function GenomicIntervalBody({ preview }: { preview: SchemaPreview }) {
  const summary = preview.summary as Record<string, unknown> | undefined;
  const features = summary?.feature_types as string[] | undefined;
  return (
    <div>
      {summaryStat('Intervals', preview.n_records)}
      {features && features.length > 0 && (
        <span className="text-muted" style={{ fontSize: '0.8rem' }}>
          <strong>Features:</strong> {features.slice(0, 5).join(', ')}
        </span>
      )}
    </div>
  );
}

function SingleCellBody({ preview }: { preview: SchemaPreview }) {
  const summary = preview.summary as Record<string, unknown> | undefined;
  const embeddings = summary?.embeddings as string[] | undefined;
  return (
    <div>
      {summaryStat('Cells', summary?.n_cells as number)}
      {summaryStat('Genes', summary?.n_genes as number)}
      {embeddings && embeddings.length > 0 && (
        <span className="me-3 text-muted" style={{ fontSize: '0.8rem' }}>
          <strong>Embeddings:</strong> {embeddings.join(', ')}
        </span>
      )}
    </div>
  );
}

function AlignmentBody({ preview }: { preview: SchemaPreview }) {
  const summary = preview.summary as Record<string, unknown> | undefined;
  return (
    <div>
      {summaryStat('Reads sampled', summary?.n_reads_sampled as number)}
      {summary?.pct_mapped != null &&
        summaryStat('Mapped', `${summary.pct_mapped as number}%`)}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Suggested questions per family
// ---------------------------------------------------------------------------

function suggestionsFor(filename: string, preview: SchemaPreview): string[] {
  const f = preview.family ?? 'unknown';
  const cols: ColumnStat[] = preview.schema?.columns ?? [];
  const numCols = cols.filter((c) => ['float64', 'float32', 'int64', 'int32'].includes(c.dtype));
  const firstNum = numCols[0]?.name;
  const name = filename.replace(/\.[^.]+$/, '');

  switch (f) {
    case 'tabular': {
      const base = [
        `Show me the first 10 rows of ${name}`,
        `How many rows and columns does ${name} have?`,
      ];
      if (firstNum) base.push(`What is the distribution of ${firstNum}?`);
      if (numCols.length >= 2)
        base.push(`Correlate ${numCols[0].name} and ${numCols[1].name}`);
      base.push(`Are there any missing values in ${name}?`);
      return base;
    }
    case 'sequence':
      return [
        `How many sequences are in ${name}?`,
        `What is the average sequence length in ${name}?`,
        `Perform multiple sequence alignment on ${name}`,
      ];
    case 'variant':
      return [
        `How many variants are in ${name}?`,
        `What is the variant distribution across chromosomes in ${name}?`,
        `Filter high-quality variants from ${name}`,
      ];
    case 'genomic_interval':
      return [
        `How many intervals are in ${name}?`,
        `What feature types are present in ${name}?`,
      ];
    case 'single_cell':
      return [
        `How many cells and genes are in ${name}?`,
        `What cell types are present in ${name}?`,
        `Run dimensionality reduction on ${name}`,
      ];
    case 'alignment':
      return [
        `What is the mapping rate for ${name}?`,
        `Summarize alignment statistics for ${name}`,
      ];
    default:
      return [`Analyze ${name}`];
  }
}

// ---------------------------------------------------------------------------
// Main component
// ---------------------------------------------------------------------------

export function SchemaPreviewPanel({
  filename,
  preview,
  onSuggest,
}: SchemaPreviewPanelProps) {
  const [open, setOpen] = useState(true);
  const family = preview.family ?? 'unknown';
  const label = FAMILY_LABELS[family] ?? FAMILY_LABELS.unknown;
  const hasError = Boolean(preview.profiler_error);
  const suggestions = suggestionsFor(filename, preview);

  return (
    <div
      style={{
        background: '#f8f9fb',
        border: '1px solid #dee2e6',
        borderRadius: 6,
        marginTop: 6,
        fontSize: '0.82rem',
      }}
    >
      {/* Header row */}
      <div
        className="d-flex justify-content-between align-items-center px-2 py-1"
        style={{ cursor: 'pointer', userSelect: 'none' }}
        onClick={() => setOpen((o) => !o)}
      >
        <span>
          <span className="me-2">{label}</span>
          {!hasError && preview.n_records != null && (
            <span className="text-muted" style={{ fontSize: '0.78rem' }}>
              {formatCount(preview.n_records)} records
            </span>
          )}
          {hasError && (
            <span className="text-warning" style={{ fontSize: '0.78rem' }}>
              ⚠ profiling limited
            </span>
          )}
        </span>
        <span className="text-muted" style={{ fontSize: '0.75rem' }}>
          {open ? '▲' : '▼'}
        </span>
      </div>

      {/* Collapsible body */}
      {open && (
        <div className="px-2 pb-2" style={{ borderTop: '1px solid #dee2e6' }}>
          {hasError ? (
            <div className="text-muted mt-1" style={{ fontSize: '0.78rem' }}>
              {preview.profiler_error}
            </div>
          ) : (
            <div className="mt-1">
              {family === 'tabular' && <TabularBody preview={preview} />}
              {family === 'sequence' && <SequenceBody preview={preview} />}
              {family === 'variant' && <VariantBody preview={preview} />}
              {family === 'genomic_interval' && <GenomicIntervalBody preview={preview} />}
              {family === 'single_cell' && <SingleCellBody preview={preview} />}
              {family === 'alignment' && <AlignmentBody preview={preview} />}
            </div>
          )}

          {/* Suggested questions */}
          {onSuggest && suggestions.length > 0 && (
            <div className="mt-2">
              <div className="text-muted mb-1" style={{ fontSize: '0.75rem' }}>
                Ask about this file:
              </div>
              <div className="d-flex flex-wrap gap-1">
                {suggestions.map((q) => (
                  <button
                    key={q}
                    className="btn btn-outline-primary btn-sm"
                    style={{ fontSize: '0.72rem', padding: '1px 7px', borderRadius: 12 }}
                    onClick={() => onSuggest(q)}
                  >
                    {q}
                  </button>
                ))}
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
