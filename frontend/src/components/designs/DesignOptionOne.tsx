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
  onRemoveAllFiles,
  onDropZoneDragOver,
  onDropZoneDragLeave,
  onDropZoneDrop,
  onBrowseClick,
  examplesPanel,
  quickExamples,
  onExampleClick,
  examplesOpen,
  onToggleExamples,
  workflowContextContent,
  historyContent,
}) => {
  const renderToggleButton = (isOpen: boolean, onToggle: () => void) => (
    <Button
      variant="outline-primary"
      size="sm"
      onClick={onToggle}
      aria-label={isOpen ? 'Collapse examples' : 'Expand examples'}
      className={`collapsed-sidebar-button${isOpen ? ' is-open' : ''}`}
    >
      ⟩
    </Button>
  );

  const renderCollapsedHandle = (label: string, onToggle: () => void) => (
    <div className="collapsed-sidebar-panel">
      {renderToggleButton(false, onToggle)}
      <div className="collapsed-sidebar-label">{label}</div>
    </div>
  );

  const renderSection = (
    label: string,
    isOpen: boolean,
    onToggle: () => void,
    content: React.ReactNode,
  ) => {
    if (!isOpen) {
      return (
        <div className="sidebar-section sidebar-section--collapsed">
          {renderCollapsedHandle(label, onToggle)}
        </div>
      );
    }

    return (
      <div className="sidebar-section">
        <div className="expanded-sidebar-header">
          {renderToggleButton(true, onToggle)}
          <div className="expanded-sidebar-label">{label}</div>
        </div>
        <div className="expanded-sidebar-content">{content}</div>
      </div>
    );
  };

  const renderExamplesContent = () => {
    if (examplesPanel) {
      return examplesPanel;
    }

    return (
      <div className="p-3 d-flex flex-column gap-2">
        {quickExamples.map((example, idx) => (
          <Button
            key={idx}
            variant="outline-secondary"
            size="sm"
            className="text-start"
            onClick={() => onExampleClick(example.command)}
            style={{ whiteSpace: 'normal' }}
          >
            <div className="fw-semibold">{example.title}</div>
            {example.description && (
              <div className="text-muted small">{example.description}</div>
            )}
          </Button>
        ))}
      </div>
    );
  };

  const examplesCard = (
    <Card className="shadow-sm border-0 mt-3">
      <Card.Body className="p-0" style={{ maxHeight: '65vh', overflowY: 'auto' }}>
        {renderExamplesContent()}
      </Card.Body>
    </Card>
  );

  return (
    <div className="design-option design-option-one">
      <div className="d-flex flex-column flex-lg-row gap-4">
        <aside className="flex-shrink-0" style={{ width: examplesOpen ? '320px' : 'auto', flex: '0 0 auto' }}>
          <div className="d-flex flex-column align-items-center gap-3">
            {renderSection('Examples', examplesOpen, onToggleExamples, examplesCard)}
          </div>
        </aside>

        <div className="flex-grow-1 d-flex flex-column gap-4" style={{ minWidth: 0 }}>
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
                <div className="mt-3">
                  <div className="d-flex justify-content-between align-items-center mb-2">
                    <div className="fw-semibold">Uploaded Files ({uploadedFiles.length})</div>
                    <Button variant="outline-secondary" size="sm" onClick={onRemoveAllFiles}>
                      Remove All
                    </Button>
                  </div>
                  {uploadedFiles.map((file, index) => (
                    <Card key={index} className="mb-2 border-0 bg-blue-subtle">
                      <Card.Body className="d-flex justify-content-between align-items-center">
                        <div>
                          <div className="fw-semibold">
                            {file.name}
                            {file.name.match(/[._-]R?1[._-]/i) && <span className="badge bg-primary ms-2">R1</span>}
                            {file.name.match(/[._-]R?2[._-]/i) && <span className="badge bg-success ms-2">R2</span>}
                          </div>
                          <div className="text-muted small">{file.content.length.toLocaleString()} characters</div>
                        </div>
                        <Button variant="outline-secondary" size="sm" onClick={() => onFileRemove(index)}>
                          Remove
                        </Button>
                      </Card.Body>
                    </Card>
                  ))}
                </div>
              )}
            </Card.Body>
          </Card>

          {workflowContextContent}
          {historyContent}
        </div>
      </div>
    </div>
  );
};

export default DesignOptionOne;
