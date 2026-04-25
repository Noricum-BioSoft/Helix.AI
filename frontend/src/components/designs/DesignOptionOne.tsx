import React from 'react';
import Card from 'react-bootstrap/Card';
import Button from 'react-bootstrap/Button';
import { PromptDesignProps } from './types';
import { PromptInput } from './PromptInput';
import { SchemaPreviewPanel } from '../SchemaPreviewPanel';

export const DesignOptionOne: React.FC<PromptDesignProps> = ({
  command,
  onCommandChange,
  onSubmit,
  onSuggestSubmit,
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
  workflowContextContent,
  historyContent,
}) => {
  return (
    <div className="design-option design-option-one">
      <div className="d-flex flex-column gap-3">
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
                    <div key={index} className="mb-2 pb-2 border-bottom">
                      <div className="d-flex justify-content-between align-items-center">
                        <div>
                          <div className="fw-semibold small">{file.name}</div>
                          <div className="text-muted small">
                            {(file.size / (1024 * 1024)).toFixed(2)} MB
                            {file.status === 'uploading' && ' • uploading…'}
                            {file.status === 'uploaded' && ' • ready'}
                            {file.status === 'failed' && ` • ${file.error ?? 'failed'}`}
                          </div>
                        </div>
                        <Button variant="outline-secondary" size="sm" onClick={() => onFileRemove(index)}>
                          Remove
                        </Button>
                      </div>
                      {file.schema_preview && file.status === 'uploaded' && (
                        <SchemaPreviewPanel
                          filename={file.name}
                          preview={file.schema_preview}
                          onSuggest={onCommandChange}
                          onSuggestSubmit={onSuggestSubmit}
                        />
                      )}
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
