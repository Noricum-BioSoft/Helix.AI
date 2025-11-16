import React from 'react';
import Card from 'react-bootstrap/Card';
import Button from 'react-bootstrap/Button';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';
import { PromptDesignProps } from './types';
import { PromptInput } from './PromptInput';

export const DesignOptionTwo: React.FC<PromptDesignProps> = ({
  command,
  onCommandChange,
  onSubmit,
  loading,
  agentLoading,
  placeholder,
  onAgentSubmit,
  dragActive,
  uploadedFile,
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
  const examplesToggleLabel = examplesOpen ? 'Collapse Examples' : 'Expand Examples';
  const sidebarExpanded = examplesOpen;

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

  const quickExamplesCard = (
    <Card className="shadow-sm border-0 mt-3">
      <Card.Body>
        <div className="d-flex flex-column gap-2">
          {quickExamples.map((example, idx) => (
            <Button
              key={idx}
              variant="outline-secondary"
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
      </Card.Body>
    </Card>
  );

  const examplesLibraryCard = (
    <Card className="shadow-sm border-0 mt-3">
      <Card.Body className="p-0" style={{ maxHeight: '65vh', overflowY: 'auto' }}>
        {examplesPanel ?? (
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
        )}
      </Card.Body>
    </Card>
  );

  return (
    <div className="design-option design-option-two">
      <div className="d-flex justify-content-between align-items-center mb-3">
        <h2 className="h4 mb-0">Split Workspace</h2>
        <Button variant="outline-primary" size="sm" onClick={onToggleExamples}>
          {examplesToggleLabel}
        </Button>
      </div>

      <Row className="g-4">
        <Col
          lg={sidebarExpanded ? 8 : 11}
          className="d-flex flex-column gap-4"
          style={{ minWidth: 0 }}
        >
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

              {uploadedFile && (
                <Card className="mt-3 border-0 bg-light">
                  <Card.Body className="d-flex justify-content-between align-items-center">
                    <div>
                      <div className="fw-semibold">{uploadedFile.name}</div>
                      <div className="text-muted small">{uploadedFile.content.length.toLocaleString()} characters</div>
                    </div>
                    <Button variant="outline-secondary" size="sm" onClick={onFileRemove}>
                      Remove
                    </Button>
                  </Card.Body>
                </Card>
              )}
            </Card.Body>
          </Card>

          {workflowContextContent}
          {historyContent}
        </Col>

        <Col xs="auto" lg={sidebarExpanded ? 4 : 'auto'} style={{ maxWidth: sidebarExpanded ? '320px' : 'auto', flex: sidebarExpanded ? '0 0 320px' : '0 0 auto' }}>
          <div className="d-flex flex-column align-items-center gap-4">
            {renderSection('Examples', examplesOpen, onToggleExamples, (
              <>
                {quickExamplesCard}
                {examplesLibraryCard}
              </>
            ))}
          </div>
        </Col>
      </Row>
    </div>
  );
};

export default DesignOptionTwo;
