import React, { useEffect, useMemo, useState } from 'react';
import Offcanvas from 'react-bootstrap/Offcanvas';
import Button from 'react-bootstrap/Button';
import Badge from 'react-bootstrap/Badge';
import Spinner from 'react-bootstrap/Spinner';
import Modal from 'react-bootstrap/Modal';
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
  const [selectedJobId, setSelectedJobId] = useState<string | null>(null);
  const [jobResults, setJobResults] = useState<any>(null);
  const [loadingResults, setLoadingResults] = useState(false);
  const [resultsError, setResultsError] = useState<string | null>(null);

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

  const isCompleted = (status?: string) => {
    const s = (status || '').toLowerCase();
    return s.includes('completed') || s.includes('success');
  };

  const handleJobClick = async (jobId: string, status?: string) => {
    if (!isCompleted(status)) {
      return; // Only show results for completed jobs
    }

    setSelectedJobId(jobId);
    setJobResults(null);
    setResultsError(null);
    setLoadingResults(true);

    try {
      const response = await mcpApi.getJobResults(jobId);
      setJobResults(response.results || response);
    } catch (error: any) {
      setResultsError(error?.response?.data?.detail || error?.message || 'Failed to load job results');
      console.error('Error loading job results:', error);
    } finally {
      setLoadingResults(false);
    }
  };

  const handleCloseResults = () => {
    setSelectedJobId(null);
    setJobResults(null);
    setResultsError(null);
  };

  const formatResults = (results: any): string => {
    if (!results) return '';
    if (typeof results === 'string') return results;
    if (typeof results === 'object') {
      try {
        return JSON.stringify(results, null, 2);
      } catch {
        return String(results);
      }
    }
    return String(results);
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
                const completed = isCompleted(j?.status);
                const clickable = completed;
                return (
                  <div 
                    key={jobId} 
                    className={`border rounded p-2 ${clickable ? 'cursor-pointer' : ''}`}
                    style={{ 
                      cursor: clickable ? 'pointer' : 'default',
                      transition: 'background-color 0.2s'
                    }}
                    onClick={() => clickable && handleJobClick(jobId, j?.status)}
                    onMouseEnter={(e) => {
                      if (clickable) {
                        e.currentTarget.style.backgroundColor = 'rgba(0, 0, 0, 0.05)';
                      }
                    }}
                    onMouseLeave={(e) => {
                      if (clickable) {
                        e.currentTarget.style.backgroundColor = '';
                      }
                    }}
                  >
                    <div className="d-flex justify-content-between align-items-start gap-2">
                      <div style={{ minWidth: 0 }}>
                        <div className="fw-semibold text-truncate">{jobId}</div>
                        <div className="text-muted small">
                          {(j?.tool_name || j?.job_type) ? `${j?.tool_name || j?.job_type}` : '—'}
                        </div>
                      </div>
                      <Badge 
                        bg={badgeVariant(j?.status)}
                        style={clickable ? { cursor: 'pointer' } : {}}
                      >
                        {j?.status || 'unknown'}
                      </Badge>
                    </div>
                    {j?.error && <div className="text-danger small mt-2">{j.error}</div>}
                    {completed && (
                      <div className="text-muted small mt-1" style={{ fontStyle: 'italic' }}>
                        Click to view results
                      </div>
                    )}
                  </div>
                );
              })}
            </div>
          </>
        )}
      </Offcanvas.Body>

      {/* Results Modal */}
      <Modal show={selectedJobId !== null} onHide={handleCloseResults} size="lg">
        <Modal.Header closeButton>
          <Modal.Title>
            Job Results: {selectedJobId}
            {jobs[selectedJobId || ''] && (
              <Badge bg={badgeVariant(jobs[selectedJobId || ''].status)} className="ms-2">
                {jobs[selectedJobId || ''].status}
              </Badge>
            )}
          </Modal.Title>
        </Modal.Header>
        <Modal.Body>
          {loadingResults ? (
            <div className="d-flex align-items-center justify-content-center py-4">
              <Spinner animation="border" className="me-2" />
              <span>Loading results...</span>
            </div>
          ) : resultsError ? (
            <div className="alert alert-danger">
              <strong>Error:</strong> {resultsError}
            </div>
          ) : jobResults ? (
            <div>
              {/* Show validation warnings if present */}
              {jobResults.output_validation && (
                <div className="mb-3">
                  {!jobResults.output_validation.all_files_exist && (
                    <div className="alert alert-warning">
                      <strong>⚠️ Output File Validation Warning:</strong>
                      <ul className="mb-0 mt-2">
                        {jobResults.output_validation.validation_warnings?.map((warning: string, idx: number) => (
                          <li key={idx}>{warning}</li>
                        ))}
                      </ul>
                      {jobResults.output_validation.output_files_found && jobResults.output_validation.output_files_found.length > 0 && (
                        <div className="mt-3">
                          <strong>Output Files Checked:</strong>
                          <ul className="mb-0 mt-2">
                            {jobResults.output_validation.output_files_found.map((file: any, idx: number) => (
                              <li key={idx}>
                                <span className={file.exists ? 'text-success' : 'text-danger'}>
                                  {file.exists ? '✅' : '❌'} {file.location}: {file.path}
                                </span>
                                {file.error && <small className="text-muted d-block ms-3">({file.error})</small>}
                              </li>
                            ))}
                          </ul>
                        </div>
                      )}
                    </div>
                  )}
                  {jobResults.output_validation.all_files_exist && jobResults.output_validation.output_files_found && jobResults.output_validation.output_files_found.length > 0 && (
                    <div className="alert alert-success">
                      <strong>✅ All Output Files Verified:</strong>
                      <ul className="mb-0 mt-2">
                        {jobResults.output_validation.output_files_found.map((file: any, idx: number) => (
                          <li key={idx}>{file.location}: {file.path}</li>
                        ))}
                      </ul>
                    </div>
                  )}
                </div>
              )}
              
              {/* Show top-level warnings */}
              {jobResults.warning && (
                <div className="alert alert-warning mb-3">
                  <strong>Warning:</strong> {jobResults.warning}
                </div>
              )}
              
              <h6 className="mb-3">Results:</h6>
              <pre 
                className="bg-light p-3 rounded border"
                style={{ 
                  maxHeight: '60vh', 
                  overflow: 'auto',
                  whiteSpace: 'pre-wrap',
                  wordBreak: 'break-word',
                  fontSize: '0.9em'
                }}
              >
                {formatResults(jobResults)}
              </pre>
            </div>
          ) : null}
        </Modal.Body>
        <Modal.Footer>
          <Button variant="secondary" onClick={handleCloseResults}>
            Close
          </Button>
        </Modal.Footer>
      </Modal>
    </Offcanvas>
  );
};

export default JobsPanel;






