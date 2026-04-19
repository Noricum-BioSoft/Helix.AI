import { describe, it, expect, vi } from 'vitest';
import { render, screen, fireEvent } from '@testing-library/react';
import { SchemaPreviewPanel } from './SchemaPreviewPanel';
import type { SchemaPreview } from '../services/helixApi';

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

const tabularPreview: SchemaPreview = {
  format: 'csv',
  family: 'tabular',
  n_records: 100,
  summary: { n_rows: 100, n_cols: 3, size_bytes: 2048 },
  schema: {
    columns: [
      { name: 'gene', dtype: 'object', n_unique: 100 },
      { name: 'log2fc', dtype: 'float64', min: -3.2, max: 2.56, mean: 0.1 },
      { name: 'padj', dtype: 'float64', min: 0.0001, max: 0.99, mean: 0.3 },
    ],
  },
  sample: [{ gene: 'BRCA1', log2fc: 1.41, padj: 0.001 }],
  profiler_error: null,
};

const sequencePreview: SchemaPreview = {
  format: 'fasta',
  family: 'sequence',
  n_records: 42,
  summary: { length_min: 100, length_max: 500, length_mean: 250 },
  profiler_error: null,
};

const singleCellPreview: SchemaPreview = {
  format: 'h5ad',
  family: 'single_cell',
  n_records: 3000,
  summary: { n_cells: 3000, n_genes: 20000, embeddings: ['X_umap', 'X_pca'] },
  schema: { obs_columns: ['cell_type', 'batch'], var_columns: ['gene_name'] },
  profiler_error: null,
};

const errorPreview: SchemaPreview = {
  format: 'csv.gz',
  family: 'tabular',
  n_records: null,
  profiler_error: 'gzip-compressed csv files are not yet supported. Please decompress first.',
};

// ---------------------------------------------------------------------------
// Tests: rendering
// ---------------------------------------------------------------------------

describe('SchemaPreviewPanel', () => {
  it('renders family label for tabular', () => {
    render(<SchemaPreviewPanel filename="data.csv" preview={tabularPreview} />);
    expect(screen.getByText(/Tabular/)).toBeTruthy();
  });

  it('renders record count', () => {
    render(<SchemaPreviewPanel filename="data.csv" preview={tabularPreview} />);
    expect(screen.getByText(/100 records/)).toBeTruthy();
  });

  it('renders all column names', () => {
    render(<SchemaPreviewPanel filename="data.csv" preview={tabularPreview} />);
    expect(screen.getByText('gene')).toBeTruthy();
    expect(screen.getByText('log2fc')).toBeTruthy();
    expect(screen.getByText('padj')).toBeTruthy();
  });

  it('renders dtype badges', () => {
    render(<SchemaPreviewPanel filename="data.csv" preview={tabularPreview} />);
    const badges = screen.getAllByText('float64');
    expect(badges.length).toBeGreaterThanOrEqual(2);
  });

  it('renders numeric min-max range', () => {
    render(<SchemaPreviewPanel filename="data.csv" preview={tabularPreview} />);
    expect(screen.getByText(/-3.2.*2.56/)).toBeTruthy();
  });

  it('renders sequence family label and stats', () => {
    render(<SchemaPreviewPanel filename="seqs.fasta" preview={sequencePreview} />);
    expect(screen.getByText(/Sequence/)).toBeTruthy();
    expect(screen.getByText(/42 records/)).toBeTruthy();
  });

  it('renders single-cell cells and genes', () => {
    render(<SchemaPreviewPanel filename="cells.h5ad" preview={singleCellPreview} />);
    expect(screen.getByText(/Single Cell/)).toBeTruthy();
    expect(screen.getByText(/3,000/)).toBeTruthy();
  });

  it('shows profiler error instead of schema when present', () => {
    render(<SchemaPreviewPanel filename="data.csv.gz" preview={errorPreview} />);
    expect(screen.getByText(/please decompress/i)).toBeTruthy();
  });

  it('shows "profiling limited" badge when error present', () => {
    render(<SchemaPreviewPanel filename="data.csv.gz" preview={errorPreview} />);
    expect(screen.getByText(/profiling limited/)).toBeTruthy();
  });

  // ---------------------------------------------------------------------------
  // Tests: collapse/expand
  // ---------------------------------------------------------------------------

  it('collapses when header is clicked', () => {
    render(<SchemaPreviewPanel filename="data.csv" preview={tabularPreview} />);
    const header = screen.getByText(/Tabular/).closest('div')!;
    // Initially open — column names visible
    expect(screen.getByText('gene')).toBeTruthy();
    fireEvent.click(header);
    // After collapse — column names gone
    expect(screen.queryByText('gene')).toBeNull();
  });

  it('re-expands when header clicked again', () => {
    render(<SchemaPreviewPanel filename="data.csv" preview={tabularPreview} />);
    const header = screen.getByText(/Tabular/).closest('div')!;
    fireEvent.click(header);
    fireEvent.click(header);
    expect(screen.getByText('gene')).toBeTruthy();
  });

  // ---------------------------------------------------------------------------
  // Tests: suggested questions
  // ---------------------------------------------------------------------------

  it('renders suggested question chips when onSuggest is provided', () => {
    render(<SchemaPreviewPanel filename="data.csv" preview={tabularPreview} onSuggest={vi.fn()} />);
    expect(screen.getByText(/Ask about this file/)).toBeTruthy();
  });

  it('calls onSuggest with the question text when chip is clicked', () => {
    const onSuggest = vi.fn();
    render(
      <SchemaPreviewPanel
        filename="expr.csv"
        preview={tabularPreview}
        onSuggest={onSuggest}
      />
    );
    const chip = screen.getByText(/Show me the first 10 rows/);
    fireEvent.click(chip);
    expect(onSuggest).toHaveBeenCalledWith(expect.stringContaining('first 10 rows'));
  });

  it('does not render suggestions when onSuggest is not provided', () => {
    render(<SchemaPreviewPanel filename="data.csv" preview={tabularPreview} />);
    // "Ask about this file" should not appear without the callback
    expect(screen.queryByText(/Ask about this file/)).toBeNull();
  });

  // ---------------------------------------------------------------------------
  // Tests: column expand/collapse for wide tables
  // ---------------------------------------------------------------------------

  it('shows "more columns" button when table has more than 6 columns', () => {
    const widePreview: SchemaPreview = {
      ...tabularPreview,
      schema: {
        columns: Array.from({ length: 10 }, (_, i) => ({
          name: `col_${i}`,
          dtype: 'float64',
          min: 0,
          max: 1,
        })),
      },
    };
    render(<SchemaPreviewPanel filename="wide.csv" preview={widePreview} onSuggest={vi.fn()} />);
    expect(screen.getByText(/\+4 more columns/)).toBeTruthy();
  });

  it('expands all columns when "more columns" is clicked', () => {
    const widePreview: SchemaPreview = {
      ...tabularPreview,
      schema: {
        columns: Array.from({ length: 10 }, (_, i) => ({
          name: `col_${i}`,
          dtype: 'float64',
          min: 0,
          max: 1,
        })),
      },
    };
    render(<SchemaPreviewPanel filename="wide.csv" preview={widePreview} onSuggest={vi.fn()} />);
    fireEvent.click(screen.getByText(/\+4 more columns/));
    // All 10 columns should now be visible
    for (let i = 0; i < 10; i++) {
      expect(screen.getByText(`col_${i}`)).toBeTruthy();
    }
  });
});
