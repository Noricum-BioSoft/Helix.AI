import React from 'react';
import Card from 'react-bootstrap/Card';
import Button from 'react-bootstrap/Button';
import { PromptDesignProps } from './types';
import { PromptInput } from './PromptInput';

export const DesignOptionOne: React.FC<PromptDesignProps> = ({
  command,
  onCommandChange,
  onSubmit,
  loading,
  agentLoading,
  placeholder,
  onAgentSubmit,
  dragActive,
  uploadedFiles,
  onFileRemove,
  onDropZoneDragOver,
  onDropZoneDragLeave,
  onDropZoneDrop,
  onBrowseClick,
  onToggleExamples,
  jobsOpen,
  onToggleJobs,
  jobsCount,
  workflowContextContent,
  historyContent,
}) => {
  return (
    <div className="design-option design-option-one">
      <div className="d-flex flex-column gap-3">
        <div className="prompt-toolbar">
          <Button
            variant="outline-secondary"
            className="prompt-toolbar-button"
            onClick={onToggleExamples}
            aria-label="Toggle examples"
          >
            📚 Examples
          </Button>
          <Button
            variant="outline-secondary"
            className="prompt-toolbar-button"
            onClick={onToggleJobs}
            aria-label="Toggle jobs"
          >
            🧾 Jobs{typeof jobsCount === 'number' && jobsCount > 0 ? ` (${jobsCount})` : ''}
          </Button>
        </div>

        <Card className="shadow-sm border-0">
          <Card.Body>
            <PromptInput
              command={command}
              onCommandChange={onCommandChange}
              onSubmit={onSubmit}
              loading={loading}
              agentLoading={agentLoading}
              onAgentSubmit={onAgentSubmit}
              placeholder={placeholder}
              dragActive={dragActive}
              onDropZoneDragOver={onDropZoneDragOver}
              onDropZoneDragLeave={onDropZoneDragLeave}
              onDropZoneDrop={onDropZoneDrop}
              onBrowseClick={onBrowseClick}
            />

            {uploadedFiles.length > 0 && (
              <Card className="mt-3 border-0 bg-blue-subtle">
                <Card.Body>
                  <div className="d-flex justify-content-between align-items-center mb-2">
                    <div className="fw-semibold">Uploaded Files ({uploadedFiles.length})</div>
                    <Button variant="outline-secondary" size="sm" onClick={() => onFileRemove()}>
                      Clear All
                    </Button>
                  </div>
                  {uploadedFiles.map((file, index) => (
                    <div key={index} className="d-flex justify-content-between align-items-center mb-2 pb-2 border-bottom">
                      <div>
                        <div className="fw-semibold small">{file.name}</div>
                        <div className="text-muted small">{file.content.length.toLocaleString()} characters</div>
                      </div>
                      <Button variant="outline-secondary" size="sm" onClick={() => onFileRemove(index)}>
                        Remove
                      </Button>
                    </div>
                  ))}
                </Card.Body>
              </Card>
            )}
          </Card.Body>
        </Card>

        {workflowContextContent}
        {historyContent}
      </div>
    </div>
  );
};

export default DesignOptionOne;
