import React, { useState } from 'react';
import Card from 'react-bootstrap/Card';
import Button from 'react-bootstrap/Button';
import Nav from 'react-bootstrap/Nav';
import Tab from 'react-bootstrap/Tab';
import { PromptDesignProps } from './types';
import { PromptInput } from './PromptInput';

export const DesignOptionThree: React.FC<PromptDesignProps> = ({
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
  examplesPanel,
  quickExamples,
  onExampleClick,
  examplesOpen,
  onToggleExamples,
  workflowContextContent,
  historyContent,
}) => {
  const [activeTab, setActiveTab] = useState<'prompt' | 'upload'>('prompt');
  const examplesToggleLabel = examplesOpen ? 'Collapse Examples' : 'Expand Examples';

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

  const examplesCard = (
    <Card className="shadow-sm border-0 mt-3">
      <Card.Body className="p-0" style={{ maxHeight: '60vh', overflowY: 'auto' }}>
        {examplesPanel ?? (
          <div className="p-3 d-flex flex-column gap-2">
            {quickExamples.map((example, idx) => (
              <Button
                key={idx}
                variant="outline-secondary"
                size="sm"
                className="text-start"
                style={{ whiteSpace: 'normal' }}
                onClick={() => onExampleClick(example.command)}
              >
                <div className="fw-semibold">{example.title}</div>
                {example.description && (
                  <div className="text-muted small">{example.description}</div>
                )}
              </Button>
            ))}
          </div>
        )}
      </Card.Body>
    </Card>
  );

  return (
    <div className="design-option design-option-three d-flex flex-column gap-4">
      <div className="d-flex justify-content-between align-items-center">
        <h2 className="h4 mb-0">Workspace Tabs</h2>
        <Button variant="outline-primary" size="sm" onClick={onToggleExamples}>
          {examplesToggleLabel}
        </Button>
      </div>

      <Card className="shadow-sm border-0">
        <Card.Body>
          <Tab.Container activeKey={activeTab} onSelect={(key) => key && setActiveTab(key as 'prompt' | 'upload')}>
            <div className="d-flex justify-content-between align-items-center flex-wrap gap-2 mb-3">
              <Nav variant="pills">
                <Nav.Item>
                  <Nav.Link eventKey="prompt">Prompt</Nav.Link>
                </Nav.Item>
                <Nav.Item>
                  <Nav.Link eventKey="upload">Upload</Nav.Link>
                </Nav.Item>
              </Nav>
              <Button
                variant="primary"
                onClick={onSubmit}
                disabled={loading || agentLoading || !command.trim()}
              >
                {loading ? (
                  <>
                    <span className="spinner-border spinner-border-sm me-2" role="status" aria-hidden="true"></span>
                    Running…
                  </>
                ) : (
                  'Run Command'
                )}
              </Button>
            </div>

            <Tab.Content>
              <Tab.Pane eventKey="prompt">
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
              </Tab.Pane>

              <Tab.Pane eventKey="upload">
                <h5>Upload Sequence Files</h5>
                <p className="text-muted small">
                  Drag and drop FASTA/CSV files or pick one from your computer. Uploaded files stay in context for subsequent commands.
                </p>
                <div
                  className={`p-5 text-center rounded-4 border ${dragActive ? 'border-brand-blue bg-blue-subtle' : 'border-brand-gold bg-gold-subtle'}`}
                  style={{ borderStyle: 'dashed', transition: 'all 0.2s ease-in-out' }}
                  onDragOver={onDropZoneDragOver}
                  onDragLeave={onDropZoneDragLeave}
                  onDrop={onDropZoneDrop}
                  onClick={onBrowseClick}
                >
                  <div className="mb-3" style={{ fontSize: '3rem' }}>🗂️</div>
                  <h6 className="fw-semibold">Drop FASTA or CSV files here</h6>
                  <p className="text-muted small mb-3">{dragActive ? 'Release to upload your file.' : 'Supports FASTA (.fasta, .fa, .fas), CSV (.csv), TXT (.txt) formats.'}</p>
                  <Button variant="outline-primary" onClick={onBrowseClick}>Browse Files</Button>
                </div>

                {uploadedFiles.length > 0 && (
                  <Card className="mt-3 border-0 bg-light">
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
              </Tab.Pane>
            </Tab.Content>
          </Tab.Container>
        </Card.Body>
      </Card>

      <div className="d-flex flex-column flex-lg-row gap-4">
        <div style={{ width: examplesOpen ? '320px' : 'auto' }}>
          <div className="d-flex flex-column align-items-center gap-3">
            {renderSection('Examples', examplesOpen, onToggleExamples, examplesCard)}
          </div>
        </div>
      </div>

      {workflowContextContent}
      {historyContent}
    </div>
  );
};

export default DesignOptionThree;
