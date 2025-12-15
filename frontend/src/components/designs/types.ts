import React from 'react';

export interface QuickExample {
  title: string;
  description?: string;
  command: string;
}

export interface PromptDesignProps {
  command: string;
  onCommandChange: (value: string) => void;
  onSubmit: () => void;
  onAgentSubmit: () => void;
  loading: boolean;
  agentLoading: boolean;
  placeholder?: string;
  dragActive: boolean;
  uploadedFiles: Array<{ name: string; content: string }>;
  onFileRemove: (index?: number) => void;
  onDropZoneDragOver: (event: React.DragEvent) => void;
  onDropZoneDragLeave: (event: React.DragEvent) => void;
  onDropZoneDrop: (event: React.DragEvent) => void;
  onBrowseClick: () => void;
  examplesPanel: React.ReactNode;
  quickExamples: QuickExample[];
  onExampleClick: (command: string) => void;
  examplesOpen: boolean;
  onToggleExamples: () => void;
  jobsOpen: boolean;
  onToggleJobs: () => void;
  jobsCount?: number;
  workflowContextContent: React.ReactNode;
  historyContent: React.ReactNode;
}
