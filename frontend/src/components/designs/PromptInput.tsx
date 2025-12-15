import React, { useMemo } from 'react';
import OverlayTrigger from 'react-bootstrap/OverlayTrigger';
import Tooltip from 'react-bootstrap/Tooltip';
import Button from 'react-bootstrap/Button';
import { getRecommendedButton } from '../../utils/buttonRecommendation';

interface PromptInputProps {
  command: string;
  placeholder?: string;
  loading: boolean;
  agentLoading: boolean;
  onAgentSubmit: () => void;
  dragActive: boolean;
  onCommandChange: (value: string) => void;
  onSubmit: () => void;
  examplesOpen?: boolean;
  onToggleExamples?: () => void;
  jobsOpen?: boolean;
  onToggleJobs?: () => void;
  jobsCount?: number;
  onDropZoneDragOver: (event: React.DragEvent) => void;
  onDropZoneDragLeave: (event: React.DragEvent) => void;
  onDropZoneDrop: (event: React.DragEvent) => void;
  onBrowseClick: () => void;
}

export const PromptInput: React.FC<PromptInputProps> = ({
  command,
  placeholder,
  loading,
  agentLoading,
  onAgentSubmit,
  dragActive,
  onCommandChange,
  onSubmit,
  examplesOpen,
  onToggleExamples,
  jobsOpen,
  onToggleJobs,
  jobsCount,
  onDropZoneDragOver,
  onDropZoneDragLeave,
  onDropZoneDrop,
  onBrowseClick,
}) => {
  const handlePrimarySubmit = () => {
    // Single-button UX: route to agent or direct path based on heuristic recommendation.
    if (recommendedButton === 'agent') {
      onAgentSubmit();
    } else {
      onSubmit();
    }
  };

  const handleKeyDown = (event: React.KeyboardEvent<HTMLTextAreaElement>) => {
    if ((event.ctrlKey || event.metaKey) && event.key === 'Enter') {
      event.preventDefault();
      handlePrimarySubmit();
    }
  };

  // Determine which button to recommend based on command
  const recommendedButton = useMemo(() => getRecommendedButton(command), [command]);

  return (
    <div className="prompt-shell">
      <div className="prompt-hero text-center mb-4">
        <h2 className="prompt-hero-title mb-2">What data do you want to process today?</h2>
        <div className="prompt-hero-subtitle text-muted">
          Ask anything, drop a FASTA file, or paste your sequence instructions.
        </div>
      </div>

      <div
        className={`prompt-input-wrapper ${dragActive ? 'drag-active' : ''}`}
        onDragOver={onDropZoneDragOver}
        onDragLeave={onDropZoneDragLeave}
        onDrop={onDropZoneDrop}
      >
        <button
          type="button"
          className="prompt-icon-button"
          onClick={onBrowseClick}
          aria-label="Add attachments"
        >
          +
        </button>

        <textarea
          className="prompt-textarea"
          value={command}
          onChange={(e) => onCommandChange(e.target.value)}
          rows={Math.min(15, Math.max(1, command.split('\n').length))}
          placeholder={placeholder || 'Ask anything'}
          onKeyDown={handleKeyDown}
          wrap="soft"
        />

        <div className="prompt-actions d-flex align-items-center gap-2">
          {(onToggleExamples || onToggleJobs) && (
            <div className="prompt-top-icons d-inline-flex align-items-center gap-2">
              {onToggleExamples && (
                <OverlayTrigger
                  placement="top"
                  overlay={
                    <Tooltip id="prompt-examples-tooltip">
                      {examplesOpen ? 'Hide examples' : 'Show examples'}
                    </Tooltip>
                  }
                >
                  <button
                    type="button"
                    className={`prompt-mini-icon-button ${examplesOpen ? 'is-active' : ''}`}
                    onClick={onToggleExamples}
                    aria-label={examplesOpen ? 'Hide examples' : 'Show examples'}
                  >
                    ✨
                  </button>
                </OverlayTrigger>
              )}

              {onToggleJobs && (
                <OverlayTrigger
                  placement="top"
                  overlay={
                    <Tooltip id="prompt-jobs-tooltip">
                      {jobsOpen ? 'Hide jobs' : 'Show jobs'}
                    </Tooltip>
                  }
                >
                  <button
                    type="button"
                    className={`prompt-mini-icon-button ${jobsOpen ? 'is-active' : ''}`}
                    onClick={onToggleJobs}
                    aria-label={jobsOpen ? 'Hide jobs' : 'Show jobs'}
                  >
                    <span className="prompt-mini-icon">🧾</span>
                    {typeof jobsCount === 'number' && jobsCount > 0 && (
                      <span className="prompt-mini-badge" aria-hidden="true">
                        {jobsCount}
                      </span>
                    )}
                  </button>
                </OverlayTrigger>
              )}
            </div>
          )}

          <OverlayTrigger
            placement="top"
            overlay={
              <Tooltip id="prompt-primary-tooltip">
                {recommendedButton === 'agent'
                  ? 'Run with Helix agent (recommended)'
                  : 'Run directly (recommended)'}
              </Tooltip>
            }
          >
            <span className="d-inline-flex">
              <Button
                variant="primary"
                className={`prompt-primary-button ${recommendedButton === 'agent' ? 'recommended-button' : ''}`}
                onClick={handlePrimarySubmit}
                disabled={loading || agentLoading || !command.trim()}
                aria-label={loading || agentLoading ? 'Processing' : 'Run'}
              >
                {(loading || agentLoading) ? (
                  <span className="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span>
                ) : (
                  <span aria-hidden="true">🤖</span>
                )}
              </Button>
            </span>
          </OverlayTrigger>
        </div>

        {dragActive && (
          <div className="prompt-drop-overlay">Drop your file to add it</div>
        )}
      </div>

      <div className="prompt-hint text-muted mt-3">
        Tip: Press <kbd>Ctrl</kbd> + <kbd>Enter</kbd> (or <kbd>⌘</kbd> + <kbd>Enter</kbd>) to submit instantly.
      </div>
    </div>
  );
};

export default PromptInput;
