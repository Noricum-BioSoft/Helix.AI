import React, { useState } from 'react';

const defaultParams = {
  target_property: 'thermal_stability',
  library_size: 30,
  strategy: 'rational',
  num_cycles: 1,
  assay_type: 'thermal_stability',
};

export default function DirectedEvolutionDemo() {
  const [params, setParams] = useState(defaultParams);
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<any>(null);
  const [error, setError] = useState<string | null>(null);

  const handleChange = (e: React.ChangeEvent<HTMLInputElement | HTMLSelectElement>) => {
    const { name, value } = e.target;
    setParams((prev) => ({ ...prev, [name]: name === 'library_size' || name === 'num_cycles' ? Number(value) : value }));
  };

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setLoading(true);
    setError(null);
    setResult(null);
    try {
      const command = `Run a directed evolution workflow to optimize ${params.target_property} with a library size of ${params.library_size} using ${params.strategy} strategy for ${params.num_cycles} cycle${params.num_cycles > 1 ? 's' : ''}.`;
      const res = await fetch('/execute', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ command }),
      });
      const data = await res.json();
      if (data.output && data.output.status === 'success') {
        setResult(data.output);
      } else {
        setError(data.output?.message || 'Unknown error');
      }
    } catch (err: any) {
      setError(err.message || 'Request failed');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div style={{ maxWidth: 600, margin: '2rem auto', padding: 24, background: '#fff', borderRadius: 8, boxShadow: '0 2px 8px #0001' }}>
      <h2>Directed Evolution Demo</h2>
      <form onSubmit={handleSubmit} style={{ marginBottom: 24 }}>
        <div style={{ marginBottom: 12 }}>
          <label>Property to Optimize:&nbsp;
            <select name="target_property" value={params.target_property} onChange={handleChange}>
              <option value="thermal_stability">Thermal Stability</option>
              <option value="activity">Activity</option>
              <option value="expression">Expression</option>
            </select>
          </label>
        </div>
        <div style={{ marginBottom: 12 }}>
          <label>Library Size:&nbsp;
            <input type="number" name="library_size" min={5} max={100} value={params.library_size} onChange={handleChange} />
          </label>
        </div>
        <div style={{ marginBottom: 12 }}>
          <label>Strategy:&nbsp;
            <select name="strategy" value={params.strategy} onChange={handleChange}>
              <option value="rational">Rational</option>
              <option value="random">Random</option>
            </select>
          </label>
        </div>
        <div style={{ marginBottom: 12 }}>
          <label>Number of Cycles:&nbsp;
            <input type="number" name="num_cycles" min={1} max={5} value={params.num_cycles} onChange={handleChange} />
          </label>
        </div>
        <button type="submit" disabled={loading} style={{ padding: '8px 24px', fontWeight: 600 }}>
          {loading ? 'Running...' : 'Run Workflow'}
        </button>
      </form>
      {error && <div style={{ color: 'red', marginBottom: 16 }}>{error}</div>}
      {result && (
        <div style={{ background: '#f6f8fa', padding: 16, borderRadius: 6 }}>
          <h3>Results</h3>
          <div><b>Status:</b> {result.status}</div>
          <div><b>Message:</b> {result.message}</div>
          {result.results && (
            <>
              <div><b>Best Mutant:</b> {result.results.best_mutant}</div>
              <div><b>Improvement:</b> {result.results.improvement?.toFixed(2)}</div>
              <div><b>Total Mutants:</b> {result.results.total_mutants}</div>
              <div><b>Target Property:</b> {result.results.target_property}</div>
              <div><b>Assay Type:</b> {result.results.assay_type}</div>
            </>
          )}
          {result.workflow_info && (
            <div style={{ marginTop: 12 }}>
              <b>Workflow ID:</b> {result.workflow_info.workflow_id}<br />
              <b>Started:</b> {result.workflow_info.start_time}<br />
              <b>Status:</b> {result.workflow_info.status}<br />
              <b>Total Steps:</b> {result.workflow_info.total_steps}
            </div>
          )}
        </div>
      )}
    </div>
  );
} 