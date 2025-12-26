import React from 'react';
import Offcanvas from 'react-bootstrap/Offcanvas';
import Badge from 'react-bootstrap/Badge';
import { ExampleCommandsPanel } from './ExampleCommandsPanel';

interface ExamplesPanelProps {
  show: boolean;
  onHide: () => void;
  onSelect: (command: string) => void;
}

export const ExamplesPanel: React.FC<ExamplesPanelProps> = ({ show, onHide, onSelect }) => {
  return (
    <Offcanvas show={show} onHide={onHide} placement="end">
      <Offcanvas.Header closeButton>
        <Offcanvas.Title>
          Examples <Badge bg="secondary" className="ms-2">new</Badge>
        </Offcanvas.Title>
      </Offcanvas.Header>
      <Offcanvas.Body>
        <ExampleCommandsPanel
          onCommandSelect={(cmd) => {
            onSelect(cmd);
            onHide();
          }}
        />
      </Offcanvas.Body>
    </Offcanvas>
  );
};

export default ExamplesPanel;





