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
  onDropZoneDragOver,
  onDropZoneDragLeave,
  onDropZoneDrop,
  onBrowseClick,
}) => {
  const handleKeyDown = (event: React.KeyboardEvent<HTMLTextAreaElement>) => {
    if ((event.ctrlKey || event.metaKey) && event.key === 'Enter') {
      event.preventDefault();
      onSubmit();
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
          <OverlayTrigger
            placement="top"
            overlay={
              <Tooltip id="prompt-primary-tooltip">
                Run with Helix agent (single-button UX)
              </Tooltip>
            }
          >
            <span className="d-inline-flex">
              <Button
                variant="primary"
                className="prompt-primary-button"
                onClick={onAgentSubmit}
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
