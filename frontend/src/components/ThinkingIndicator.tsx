import React from 'react';
import '../styles/ThinkingIndicator.css';

interface ThinkingIndicatorProps {
  message?: string;
  variant?: 'primary' | 'agent' | 'processing';
}

export const ThinkingIndicator: React.FC<ThinkingIndicatorProps> = ({ 
  message = 'Processing...', 
  variant = 'primary' 
}) => {
  const getVariantClass = () => {
    switch (variant) {
      case 'agent':
        return 'thinking-indicator-agent';
      case 'processing':
        return 'thinking-indicator-processing';
      default:
        return 'thinking-indicator-primary';
    }
  };

  const getVariantMessage = () => {
    if (message !== 'Processing...') return message;
    
    switch (variant) {
      case 'agent':
        return 'Agent is thinking...';
      case 'processing':
        return 'Analyzing data...';
      default:
        return 'Processing your request...';
    }
  };

  return (
    <div className={`thinking-indicator ${getVariantClass()}`}>
      <div className="thinking-indicator-content">
        <div className="thinking-dots">
          <div className="thinking-dot"></div>
          <div className="thinking-dot"></div>
          <div className="thinking-dot"></div>
        </div>
        <span className="thinking-message">{getVariantMessage()}</span>
      </div>
      <div className="thinking-pulse"></div>
    </div>
  );
};

interface LoadingOverlayProps {
  message?: string;
  variant?: 'primary' | 'agent' | 'processing';
  fullScreen?: boolean;
}

export const LoadingOverlay: React.FC<LoadingOverlayProps> = ({ 
  message, 
  variant = 'primary',
  fullScreen = false 
}) => {
  return (
    <div className={`loading-overlay ${fullScreen ? 'loading-overlay-fullscreen' : ''}`}>
      <ThinkingIndicator message={message} variant={variant} />
    </div>
  );
};

interface InlineLoadingProps {
  message?: string;
  size?: 'small' | 'medium' | 'large';
}

export const InlineLoading: React.FC<InlineLoadingProps> = ({ 
  message = 'Loading...', 
  size = 'medium' 
}) => {
  return (
    <div className={`inline-loading inline-loading-${size}`}>
      <div className="spinner-border" role="status">
        <span className="visually-hidden">Loading...</span>
      </div>
      <span className="inline-loading-message">{message}</span>
    </div>
  );
};

interface ActivityIndicatorProps {
  activities: Array<{
    id: string;
    message: string;
    status: 'active' | 'completed' | 'error';
    timestamp?: Date;
  }>;
  position?: 'top-right' | 'bottom-right' | 'top-left' | 'bottom-left';
}

export const ActivityIndicator: React.FC<ActivityIndicatorProps> = ({ 
  activities, 
  position = 'bottom-right' 
}) => {
  if (activities.length === 0) return null;

  return (
    <div className={`activity-indicator activity-indicator-${position}`}>
      <div className="activity-indicator-header">
        <div className="activity-pulse-icon"></div>
        <span>Backend Activity</span>
      </div>
      <div className="activity-list">
        {activities.map((activity) => (
          <div 
            key={activity.id} 
            className={`activity-item activity-item-${activity.status}`}
          >
            <div className="activity-status-dot"></div>
            <div className="activity-content">
              <span className="activity-message">{activity.message}</span>
              {activity.timestamp && (
                <span className="activity-timestamp">
                  {activity.timestamp.toLocaleTimeString()}
                </span>
              )}
            </div>
          </div>
        ))}
      </div>
    </div>
  );
};




