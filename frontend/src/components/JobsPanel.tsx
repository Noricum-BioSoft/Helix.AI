import React, { useEffect, useMemo, useState } from 'react';
import Offcanvas from 'react-bootstrap/Offcanvas';
import Button from 'react-bootstrap/Button';
import Badge from 'react-bootstrap/Badge';
import Spinner from 'react-bootstrap/Spinner';
import { mcpApi } from '../services/mcpApi';

type JobStatus = {
  job_id: string;
  status?: string;
  tool_name?: string;
  job_type?: string;
  submitted_at?: string;
  updated_at?: string;
  error?: string | null;
};

interface JobsPanelProps {
  show: boolean;
  onHide: () => void;
  jobIds: string[];
}

export const JobsPanel: React.FC<JobsPanelProps> = ({ show, onHide, jobIds }) => {
  const [loading, setLoading] = useState(false);
  const [jobs, setJobs] = useState<Record<string, JobStatus>>({});

  const uniqueJobIds = useMemo(() => Array.from(new Set(jobIds)).slice(0, 50), [jobIds]);

  useEffect(() => {
    if (!show) return;
    let cancelled = false;

    const pollOnce = async () => {
      if (uniqueJobIds.length === 0) return;
      setLoading(true);
      try {
        const results = await Promise.all(
          uniqueJobIds.map(async (jobId) => {
            try {
              const job = await mcpApi.getJob(jobId);
              return [jobId, job] as const;
            } catch (e) {
              return [jobId, { job_id: jobId, status: 'error', error: String(e) }] as const;
            }
          })
        );
        if (cancelled) return;
        const next: Record<string, JobStatus> = {};
        for (const [jobId, job] of results) next[jobId] = job as JobStatus;
        setJobs(next);
      } finally {
        if (!cancelled) setLoading(false);
      }
    };

    pollOnce();
    const t = window.setInterval(pollOnce, 5000);
    return () => {
      cancelled = true;
      window.clearInterval(t);
    };
  }, [show, uniqueJobIds]);

  const badgeVariant = (status?: string) => {
    const s = (status || '').toLowerCase();
    if (s.includes('completed') || s.includes('success')) return 'success';
    if (s.includes('failed') || s.includes('error')) return 'danger';
    if (s.includes('running')) return 'primary';
    if (s.includes('submitted') || s.includes('pending') || s.includes('queued')) return 'warning';
    return 'secondary';
  };

  return (
    <Offcanvas show={show} onHide={onHide} placement="end">
      <Offcanvas.Header closeButton>
        <Offcanvas.Title>
          Jobs{' '}
          <Badge bg="secondary" className="ms-2">
            {uniqueJobIds.length}
          </Badge>
        </Offcanvas.Title>
      </Offcanvas.Header>
      <Offcanvas.Body>
        {uniqueJobIds.length === 0 ? (
          <div className="text-muted">No async jobs yet. Run a long operation (e.g. FastQC) to see it here.</div>
        ) : (
          <>
            <div className="d-flex justify-content-between align-items-center mb-3">
              <div className="text-muted small">Auto-refreshes every 5s.</div>
              <Button variant="outline-secondary" size="sm" onClick={() => setJobs({})}>
                Clear cache
              </Button>
            </div>

            {loading && (
              <div className="d-flex align-items-center gap-2 mb-3">
                <Spinner animation="border" size="sm" />
                <span className="text-muted">Refreshing…</span>
              </div>
            )}

            <div className="d-flex flex-column gap-2">
              {uniqueJobIds.map((jobId) => {
                const j = jobs[jobId];
                return (
                  <div key={jobId} className="border rounded p-2">
                    <div className="d-flex justify-content-between align-items-start gap-2">
                      <div style={{ minWidth: 0 }}>
                        <div className="fw-semibold text-truncate">{jobId}</div>
                        <div className="text-muted small">
                          {(j?.tool_name || j?.job_type) ? `${j?.tool_name || j?.job_type}` : '—'}
                        </div>
                      </div>
                      <Badge bg={badgeVariant(j?.status)}>{j?.status || 'unknown'}</Badge>
                    </div>
                    {j?.error && <div className="text-danger small mt-2">{j.error}</div>}
                  </div>
                );
              })}
            </div>
          </>
        )}
      </Offcanvas.Body>
    </Offcanvas>
  );
};

export default JobsPanel;





